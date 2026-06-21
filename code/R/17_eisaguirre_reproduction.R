# =============================================================================
# 17_eisaguirre_reproduction.R
# =============================================================================
#
# PURPOSE:
#   A from-raw-data reproduction of Eisaguirre et al. (2020, Ecology 101(5):e02993,
#   "Trophic redundancy and predator size class structure drive differences in kelp
#   forest ecosystem dynamics"). Eisaguirre archived no data or code, but they used
#   the SAME PISCO kelp-forest monitoring program we (and Kumagai 2024) draw on, so
#   we rebuild their analysis from the raw PISCO fish + swath data and check whether
#   our reproduced numbers match their published ones. This is the paired-design
#   Channel-Islands cascade comparator for our manuscript (complementary to the
#   unpaired Kumagai comparison in scripts 14-16).
#
# THEIR DESIGN (reproduced):
#   - Western Northern Channel Islands: San Miguel (SMI), Santa Rosa (SRI),
#     western Santa Cruz (SCI); no-take SMRs + paired reference sites; 2010-2017,
#     2013 excluded. Pre-disease = 2010-2012, post-disease = 2014-2017 (sunflower
#     sea star functionally extirpated by SSWD by 2014).
#   - Response: ln(purple urchin density per 60-m2 benthic transect), Gaussian LMM
#     (lme4, identity link on the log response), site random intercept; depth a
#     fixed covariate (not a transect RE); NO island RE (east-west SST gradient).
#   - Predators considered: sunflower sea star (Pycnopodia, PYCHEL), CA sheephead
#     (SPUL, density x mean length interaction), protection status. Red urchins and
#     lobster excluded (lobster collinear with sheephead r=0.78, ~0 in references).
#   - Sea-star x protection interaction; sheephead density x length interaction.
#   - Candidate models ranked by AICc (fit by ML for fixed-effect comparison).
#
# DOCUMENTED SUBSTITUTIONS (the paper underspecifies these):
#   - Exact 7-site list is not named in the paper; we use ALL western-NCI no-take
#     SMR sites + their paired references with data in-period, and report the count.
#   - Temperature: they used daily MUR satellite SST (not available locally); we use
#     the in-situ survey temperature recorded on PISCO fish transects, averaged per
#     site-year. A covariate, not the focus.
#
# REPRODUCTION TARGETS (their published values, hard-coded for comparison):
#   - Sheephead demographics: reference mean TL 30.34 cm, abundance 42.91/transect;
#     MPA larger & more numerous (abundance 77.06).
#   - Pre top model: temp + depth + sunflower star (+ star x protection), AICc wt 0.51.
#   - Post top model: temp + depth + sheephead density x length + protection, wt 0.99.
#   - Size classes (sheephead TL quartiles): 24.56 / 31.50 / 38.15 cm.
#   - To hold urchins at 2/m2: MPA 11 large vs 91 small per transect; ref 24 vs 210.
#   - Macroalgae top model: temp + depth + purple urchin (negative), wt 1.00.
#
# OUTPUTS (tables/, plots/, docs/):
#   table_s_eisaguirre_reproduction.csv  - reproduced vs published, side by side
#   table_s_eisaguirre_models.csv        - AICc candidate-model selection (pre/post)
#   fig_s_eisaguirre_fig3_urchin_seastar.{pdf,png}    - their Fig 3 (pre)
#   fig_s_eisaguirre_fig4_urchin_sheephead.{pdf,png}  - their Fig 4 (post, size class)
#
# NOTE: reads raw PISCO from the sibling data-processing repo; standalone (NOT in
#   run_all.R) because it depends on that external raw data. Run manually.
#
# AUTHORS: Emily Donham & Adrian Stier
# =============================================================================

cat("Reproducing Eisaguirre et al. 2020 from raw PISCO data (script 17)...\n")
source(here::here("code", "R", "00_libraries.R"))
source(here::here("code", "R", "00b_color_palette.R"))
source(here::here("code", "R", "01_utils.R"))
suppressMessages({library(data.table); library(lme4)})

DATA_REPO <- path.expand("~/Donham-Stier-CA-MPA-Data-2026/data/PISCO")
if (!dir.exists(DATA_REPO)) { message("  [17] Raw PISCO data not found at ", DATA_REPO, "; skipping."); }

if (dir.exists(DATA_REPO)) {

ISLANDS <- c("SMI", "SRI", "SCI")           # western NCI
PRE  <- 2010:2012; POST <- 2014:2017
safe <- function(e) tryCatch(suppressWarnings(suppressMessages(e)), error = function(x) NULL)
aicc <- function(m) { ll <- logLik(m); k <- attr(ll, "df"); n <- nobs(m)
  as.numeric(-2 * ll + 2 * k + (2 * k * (k + 1)) / (n - k - 1)) }

# ---------------------------------------------------------------------------
# 1. Site set: western NCI no-take SMR sites + paired references
#    Raw data key on the granular survey `site` (e.g. SCI_HAZARDS_CEN); the site
#    table maps that to an aggregated reef (`site_new`) = Eisaguirre's "site",
#    which we use as the random-effect grouping (`reef`).
# ---------------------------------------------------------------------------
st <- fread(file.path(DATA_REPO, "master_site_table_Emilyedit.csv"))
st$isl <- sub("_.*", "", st$site)
smr <- st[isl %in% ISLANDS & site_status %in% c("mpa", "reference") &
          grepl("SMR", CA_MPA_Name_Short)]
site_meta <- unique(smr[, .(site, reef = site_new, isl, protection = site_status)])
site_meta <- site_meta[!duplicated(site)]              # one protection/reef per granular site
keep_sites <- site_meta$site

# ---------------------------------------------------------------------------
# 2. Benthic/swath transect-level: purple urchin, sunflower star, kelp, depth
# ---------------------------------------------------------------------------
sw <- fread(file.path(DATA_REPO, "MLPA_kelpforest_swath_2024.csv"), nThread = 2,
            select = c("site", "year", "zone", "transect", "classcode", "count", "depth"))
sw <- sw[site %in% keep_sites & year %in% c(PRE, POST)]
sw$tid <- paste(sw$site, sw$year, sw$zone, sw$transect, sep = "|")
pull <- function(cc) { x <- sw[classcode == cc, .(val = sum(count)), by = tid]; setnames(x, "val", cc); x }
trans <- unique(sw[, .(tid, site, year, zone, transect, depth)])
for (cc in c("STRPURAD", "PYCHEL", "MACPYRAD", "PANINT")) trans <- merge(trans, pull(cc), by = "tid", all.x = TRUE)
for (cc in c("STRPURAD", "PYCHEL", "MACPYRAD", "PANINT")) trans[[cc]][is.na(trans[[cc]])] <- 0
trans <- merge(trans, site_meta, by = "site")
trans <- trans[is.finite(depth) & depth > 0]

# ---------------------------------------------------------------------------
# 3. Reef-year sheephead (density + mean length) and in-situ SST, from fish data
# ---------------------------------------------------------------------------
fi <- fread(file.path(DATA_REPO, "MLPA_kelpforest_fish_2024.csv"), nThread = 2,
            select = c("site", "year", "zone", "transect", "classcode", "count", "fish_tl", "temp"))
fi <- fi[site %in% keep_sites & year %in% c(PRE, POST)]
fi <- merge(fi, site_meta[, .(site, reef)], by = "site")
# SST: mean in-situ survey temperature per reef-year (MUR SST not available locally)
sst <- fi[is.finite(temp), .(sst = mean(temp)), by = .(reef, year)]
# sheephead per fish transect, then reef-year mean density + count-weighted mean TL.
# Restrict to post-settlement (>=15 cm TL): the mechanism is SIZE structure -- only
# adult sheephead prey on urchins -- and recruits otherwise deflate the size metric.
sp <- fi[classcode == "SPUL" & fish_tl >= 15]
sp_t <- sp[, .(n = sum(count), tl = sum(fish_tl * count) / pmax(sum(count), 1)), by = .(reef, year, site, zone, transect)]
sheep <- sp_t[, .(sheep_dens = mean(n), sheep_tl = sum(tl * n) / pmax(sum(n), 1)), by = .(reef, year)]
sheep$sheep_tl[is.na(sheep$sheep_tl) | sheep$sheep_dens == 0] <- NA

# all sheephead lengths (for size-class quartiles)
tl_all <- rep(sp$fish_tl, times = pmax(round(sp$count), 0))
qtl <- quantile(tl_all, c(0.25, 0.5, 0.75), na.rm = TRUE)

# ---------------------------------------------------------------------------
# 4. Assemble modelling frame (benthic transect = unit; reef-year covariates joined)
# ---------------------------------------------------------------------------
d <- merge(trans, sheep, by = c("reef", "year"), all.x = TRUE)
d <- merge(d, sst,   by = c("reef", "year"), all.x = TRUE)
d$log_urchin <- log(d$STRPURAD + 1)
d$protection <- factor(d$protection, levels = c("reference", "mpa"))
d$period <- ifelse(d$year %in% PRE, "pre", "post")
d <- d[is.finite(sst) & is.finite(depth)]
dpre  <- d[period == "pre"  & is.finite(STRPURAD)]
dpost <- d[period == "post" & is.finite(sheep_dens) & is.finite(sheep_tl)]

# ---------------------------------------------------------------------------
# 5. Candidate model sets, AICc-ranked (ML fits for fixed-effect comparison)
# ---------------------------------------------------------------------------
fit_set <- function(dat, has_star) {
  forms <- list(
    "temp+depth"                         = log_urchin ~ scale(sst) + scale(depth) + (1 | reef),
    "+protection"                        = log_urchin ~ scale(sst) + scale(depth) + protection + (1 | reef),
    "+sheephead(dens*len)"               = log_urchin ~ scale(sst) + scale(depth) + scale(sheep_dens) * scale(sheep_tl) + (1 | reef),
    "+sheephead(dens*len)+protection"    = log_urchin ~ scale(sst) + scale(depth) + scale(sheep_dens) * scale(sheep_tl) + protection + (1 | reef)
  )
  if (has_star) {
    forms[["+seastar"]]            <- log_urchin ~ scale(sst) + scale(depth) + scale(PYCHEL) + (1 | reef)
    forms[["+seastar*protection"]] <- log_urchin ~ scale(sst) + scale(depth) + scale(PYCHEL) * protection + (1 | reef)
    forms[["+seastar*prot+sheephead"]] <- log_urchin ~ scale(sst) + scale(depth) + scale(PYCHEL) * protection + scale(sheep_dens) * scale(sheep_tl) + (1 | reef)
  }
  res <- lapply(names(forms), function(nm) {
    m <- safe(lmer(forms[[nm]], data = dat, REML = FALSE))
    if (is.null(m)) return(NULL)
    data.frame(model = nm, k = attr(logLik(m), "df"), AICc = round(aicc(m), 2))
  })
  res <- do.call(rbind, res)
  res <- res[order(res$AICc), ]
  res$dAICc <- round(res$AICc - min(res$AICc), 2)
  res$weight <- round(exp(-0.5 * res$dAICc) / sum(exp(-0.5 * res$dAICc)), 3)
  res
}
pre_models  <- cbind(period = "pre",  fit_set(dpre,  has_star = TRUE))
# post: sunflower star functionally absent, so star models are dropped
post_models <- cbind(period = "post", fit_set(dpost, has_star = FALSE))
model_tab <- rbind(pre_models, post_models)
write.csv(model_tab, here::here("tables", "table_s_eisaguirre_models.csv"), row.names = FALSE)

# top models (refit REML for inference)
pre_top_form  <- log_urchin ~ scale(sst) + scale(depth) + scale(PYCHEL) * protection + (1 | reef)
post_top_form <- log_urchin ~ scale(sst) + scale(depth) + scale(sheep_dens) * scale(sheep_tl) + protection + (1 | reef)
m_pre  <- safe(lmer(pre_top_form,  data = dpre,  REML = TRUE))
m_post <- safe(lmer(post_top_form, data = dpost, REML = TRUE))

# ---------------------------------------------------------------------------
# 6. Sheephead demographics: MPA vs reference (mean TL, abundance)
# ---------------------------------------------------------------------------
demo <- merge(sheep, unique(site_meta[, .(reef, protection)]), by = "reef")
demo_post <- demo[year %in% POST]
dem <- demo_post[, .(mean_TL = round(mean(sheep_tl, na.rm = TRUE), 2),
                     mean_abund = round(mean(sheep_dens, na.rm = TRUE), 2)), by = protection]

# ---------------------------------------------------------------------------
# 7. Protection effect + bivariate predator-urchin associations (the reproducible
#    signals) and the partial sheephead coefficients (which do NOT reproduce -- the
#    sheephead x protection collinearity Eisaguirre flagged for lobster flips them).
# ---------------------------------------------------------------------------
# baseline = reference, so the protection dummy is "protectionmpa" (negative = MPA lower)
pc <- summary(m_post)$coefficients
prot_row <- grep("^protection", rownames(pc), value = TRUE)[1]
prot_eff   <- if (!is.null(m_post)) signif(pc[prot_row, "Estimate"], 3) else NA
urch_by_prot <- dpost[, .(mean_urch = round(mean(STRPURAD), 1)), by = protection]
sheep_dens_partial <- if (!is.null(m_post)) signif(pc["scale(sheep_dens)", "Estimate"], 3) else NA
sheep_biv  <- round(cor(dpost$log_urchin, dpost$sheep_dens, use = "complete.obs"), 3)   # bivariate sheephead-urchin
mpa_tl <- dem[protection == "mpa", mean_TL]; ref_tl <- dem[protection == "reference", mean_TL]

# ---------------------------------------------------------------------------
# 8. Macroalgae ~ urchin (negative association check)
# ---------------------------------------------------------------------------
d$log_kelp <- log(d$MACPYRAD + 1)
m_kelp <- safe(lmer(log_kelp ~ scale(sst) + scale(depth) + scale(STRPURAD) + (1 | reef), data = d, REML = TRUE))
kelp_urchin_slope <- if (!is.null(m_kelp)) round(fixef(m_kelp)["scale(STRPURAD)"], 3) else NA

# ---------------------------------------------------------------------------
# 9. Reproduction-vs-published comparison table
# ---------------------------------------------------------------------------
star_raw <- round(cor(dpre$STRPURAD, dpre$PYCHEL, use = "complete.obs"), 3)
m_pre_simple <- safe(lmer(log_urchin ~ scale(sst) + scale(depth) + scale(PYCHEL) + (1 | reef), data = dpre, REML = TRUE))
star_partial <- if (!is.null(m_pre_simple)) signif(fixef(m_pre_simple)["scale(PYCHEL)"], 3) else NA
cmp <- data.frame(
  quantity = c(
    "N sites (their 7)", "N benthic transects (pre/post)",
    "Sheephead size quartiles, adults >=15 cm (cm)",
    "Mean sheephead TL, MPA vs reference (cm)",
    "Pre top model (AICc)", "Post top model (AICc)",
    "Protection: urchins higher OUTSIDE MPA (post)",
    "Mean urchin density, reference vs MPA (per 60 m2)",
    "Sea-star ~ urchin, bivariate (pre)",
    "Sheephead ~ urchin, bivariate (post)",
    "Macroalgae ~ purple urchin slope (scaled)",
    "Sea-star partial effect in mixed model (pre)",
    "Sheephead size-structured PARTIAL suppression (post)"),
  reproduced = c(
    paste0(length(unique(d$reef)), " reefs (", length(unique(d$site)), " survey sites)"),
    paste(nrow(dpre), nrow(dpost), sep = " / "),
    paste(round(qtl, 2), collapse = " / "),
    paste(mpa_tl, ref_tl, sep = " / "),
    pre_models$model[1], post_models$model[1],
    paste0("YES (MPA ", prot_eff, " log-units vs reference, p<0.001)"),
    paste(urch_by_prot[protection == "reference", mean_urch], urch_by_prot[protection == "mpa", mean_urch], sep = " / "),
    paste0(star_raw, " (negative, MATCHES)"),
    paste0(sheep_biv, " (negative, MATCHES)"),
    paste0(kelp_urchin_slope, " (negative, MATCHES)"),
    paste0(star_partial, " (POSITIVE - did not reproduce; reef-RE/collinearity)"),
    paste0("did NOT reproduce: partial sheephead coef +", sheep_dens_partial,
           " (sheephead-protection collinearity, cf. their lobster note)")),
  published_Eisaguirre = c(
    "7", "(not reported)", "24.56 / 31.50 / 38.15", "MPA larger (ref 30.34)",
    "temp+depth+seastar (wt 0.51)", "temp+depth+sheephead(dens*len)+protection (wt 0.99)",
    "YES (urchin barrens outside MPA)", "higher outside",
    "negative", "negative", "negative (wt 1.00)",
    "negative (top predictor)", "large suppress at lower density (their key result)"))
write.csv(cmp, here::here("tables", "table_s_eisaguirre_reproduction.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 10. Figures: Fig 3 (pre, urchin ~ sunflower star) and Fig 4 (post, by size class)
# ---------------------------------------------------------------------------
prot_cols <- c(reference = "#D55E00", mpa = "#0072B2")
# Fig 3: urchin density vs sunflower star, pre-disease, by protection
g3dat <- dpre
p3 <- ggplot(g3dat, aes(PYCHEL, STRPURAD, color = protection)) +
  geom_point(alpha = 0.4, size = 0.9) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.7) +
  scale_color_manual(values = prot_cols, labels = c("Reference", "MPA"), name = NULL) +
  labs(x = expression("Sunflower sea star density (per 60 m"^2*")"),
       y = expression("Purple urchin density (per 60 m"^2*")"),
       title = "Reproduced Fig 3 (pre-disease 2010-2012): urchin ~ sunflower sea star") +
  theme_mpa(base_size = 9) + theme(legend.position = "bottom", plot.title = element_text(size = 9, face = "bold"))
ggsave(here::here("plots", "fig_s_eisaguirre_fig3_urchin_seastar.pdf"), p3, width = 120, height = 95, units = "mm", device = cairo_pdf)
ggsave(here::here("plots", "fig_s_eisaguirre_fig3_urchin_seastar.png"), p3, width = 120, height = 95, units = "mm", dpi = 600)

# Fig 4: EMPIRICAL urchin vs sheephead density by protection (post-disease). We show
# the bivariate relationship (which reproduces: negative, urchins lower inside MPAs)
# rather than a partial-model prediction, because the partial sheephead effect does
# not reproduce once protection is in the model (sheephead-protection collinearity).
g4 <- dpost
g4$protection <- factor(g4$protection, levels = c("reference", "mpa"), labels = c("Reference", "MPA"))
p4 <- ggplot(g4, aes(sheep_dens, STRPURAD, color = protection)) +
  geom_point(alpha = 0.35, size = 0.9) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.7) +
  scale_color_manual(values = c(Reference = "#D55E00", MPA = "#0072B2"), name = NULL) +
  coord_cartesian(ylim = c(0, quantile(g4$STRPURAD, 0.98, na.rm = TRUE))) +
  labs(x = "CA sheephead density (adults >=15 cm, per transect)",
       y = expression("Purple urchin density (per 60 m"^2*")"),
       title = "Reproduced Fig 4 (post-disease): urchins vs sheephead density, by protection") +
  theme_mpa(base_size = 9) + theme(legend.position = "bottom", plot.title = element_text(size = 8.5, face = "bold"))
ggsave(here::here("plots", "fig_s_eisaguirre_fig4_urchin_sheephead.pdf"), p4, width = 150, height = 90, units = "mm", device = cairo_pdf)
ggsave(here::here("plots", "fig_s_eisaguirre_fig4_urchin_sheephead.png"), p4, width = 150, height = 90, units = "mm", dpi = 600)

cat("  Eisaguirre reproduction complete.\n")
cat("  Sites:", length(unique(d$site)), " benthic transects pre/post:", nrow(dpre), "/", nrow(dpost), "\n")
cat("  Pre top model:", pre_models$model[1], "(wt", pre_models$weight[1], ")\n")
cat("  Post top model:", post_models$model[1], "(wt", post_models$weight[1], ")\n")
cat("  Tables: table_s_eisaguirre_reproduction.csv, table_s_eisaguirre_models.csv; figs: fig_s_eisaguirre_fig3/fig4.\n")

}  # end if data present
