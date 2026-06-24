# =============================================================================
# 20_compound_disturbance.R
# =============================================================================
#
# PURPOSE:
#   Evaluate MPA resilience to the 2013-14 sea-star wasting disease (SSWD) as a
#   complement to the heat-only framing of Kumagai (2024) / scripts 14 & 19 and the
#   disease-focused Eisaguirre (2020). SSWD removed the SUNFLOWER STAR (Pycnopodia
#   helianthoides), the one sea star that is a major urchin predator.
#
# CRITICAL CAVEAT (this version):
#   The sunflower star is a COLD-WATER species. In the Southern California Bight it
#   was present mainly at the colder, northern / island reefs and essentially ABSENT
#   at most warmer (mainland, Catalina) reefs. Pre-SSWD (2008-2012), only ~1/3 of
#   our MPAs had meaningful Pycnopodia; ~half had ~0. An earlier draft lumped
#   Pycnopodia with Pisaster giganteus -- but Pisaster (which dominated total
#   sea-star density) is NOT a primary urchin predator, so that conflated the wrong
#   species. Here we use Pycnopodia SPECIFICALLY, document its patchy distribution,
#   and -- because pre-SSWD Pycnopodia VARIES across reserves while the MHW hit them
#   all -- use that variation to partly DISENTANGLE the keystone-loss (SSWD) signal
#   from the heatwave: if the cascade reorganization is stronger where more sunflower
#   star was lost, that is an SSWD effect beyond the MHW.
#
# ANALYSES:
#   (A) Pycnopodia distribution + crash where present (raw PISCO swath; Pisaster
#       reported separately for transparency, not as the urchin predator).
#   (B) Keystone moderation: is the MPA cascade response (kelp lnRR up, urchin lnRR
#       down) across before/compound/after stronger at HIGH- vs LOW-Pycnopodia
#       reserves? Stratified period effects + a period x Pycnopodia interaction.
#   (C) Overall resilience (lnRR, AR1) for context.
#   (D) Cross-study synthesis (reframed: the keystone mechanism is geographically
#       restricted to the cold-water reefs that had sunflower stars).
#
#   SCOPE: Southern California Bight (lat <= 34.45 N). period = before(<=2013) /
#   compound(2014-16, SSWD+MHW) / after(>=2017).
#
# OUTPUTS (tables/, plots/):
#   table_compound_pycnopodia_distribution.csv  - per-MPA pre-SSWD Pycnopodia/Pisaster, group
#   table_compound_pycnopodia_crash.csv         - Pycnopodia density by year at high-Pyc reefs
#   table_compound_keystone_moderation.csv      - cascade response, high vs low Pycnopodia
#   table_compound_resilience.csv               - lnRR MPA effect by period (AR1)
#   table_compound_crossstudy.csv               - 3-study synthesis
#   fig_compound_disturbance.{pdf,png}          - Pycnopodia distribution + stratified response
#
# AUTHORS: Emily Donham & Adrian Stier
# =============================================================================

cat("Running compound-disturbance analysis (SSWD/Pycnopodia + MHW; Southern California)...\n")
source(here::here("code", "R", "00_libraries.R"))
source(here::here("code", "R", "00b_color_palette.R"))
source(here::here("code", "R", "01_utils.R"))
suppressMessages({library(glmmTMB); library(emmeans); library(lme4); library(lmerTest); library(data.table)})
safe <- function(e) tryCatch(suppressWarnings(suppressMessages(e)), error = function(x) NULL)
period3 <- function(y) factor(ifelse(y <= 2013, "before", ifelse(y <= 2016, "compound", "after")),
                              levels = c("before", "compound", "after"))
PYC_HI <- 0.5   # pre-SSWD sunflower-star density (per 60 m2) threshold for "had sunflower stars"
pick_resp <- function(tx) if (tx == "Macrocystis pyrifera") "Bio" else "Den"  # one response per taxon (avoid Bio+Den pooling)

# ---------------------------------------------------------------------------
# (A) Pycnopodia distribution + crash (raw PISCO swath, SoCal)
# ---------------------------------------------------------------------------
DR <- path.expand("~/Donham-Stier-CA-MPA-Data-2026/data/PISCO")
have_swath <- dir.exists(DR)
pyc_grp <- NULL
if (have_swath) {
  st <- fread(file.path(DR, "master_site_table_Emilyedit.csv"))
  soc <- unique(st[latitude <= 34.45 & site_status %in% c("mpa", "reference"),
                   .(site, lat = latitude, mpa = CA_MPA_Name_Short)])[!duplicated(site)]
  sw <- fread(file.path(DR, "MLPA_kelpforest_swath_2024.csv"), nThread = 2,
              select = c("site", "year", "zone", "transect", "classcode", "count"))
  sw <- merge(sw, soc, by = "site")[year %in% 2008:2022]
  sw$tid <- paste(sw$site, sw$year, sw$zone, sw$transect)
  tr <- unique(sw[, .(tid, site, mpa, lat, year)])
  for (cc in c("PYCHEL", "PISGIG")) { x <- sw[classcode == cc, .(v = sum(count)), by = tid]
    setnames(x, "v", cc); tr <- merge(tr, x, by = "tid", all.x = TRUE); tr[[cc]][is.na(tr[[cc]])] <- 0 }
  # per-MPA pre-SSWD baseline (2008-2012)
  pre <- tr[year <= 2012, .(pyc_pre = round(mean(PYCHEL), 3), pis_pre = round(mean(PISGIG), 2),
                            lat = round(mean(lat), 3)), by = mpa]
  pre$pyc_group <- ifelse(pre$pyc_pre >= PYC_HI, "high", "low/zero")
  pre <- pre[order(-pyc_pre)]
  write.csv(pre, here::here("tables", "table_compound_pycnopodia_distribution.csv"), row.names = FALSE)
  pyc_grp <- pre[, .(mpa, pyc_pre, pyc_group)]
  # Pycnopodia crash through time at HIGH-Pyc reefs (where there was something to lose)
  hi_mpa <- pre[pyc_group == "high", mpa]
  crash <- tr[mpa %in% hi_mpa, .(pyc = round(mean(PYCHEL), 3), n = .N), by = year][order(year)]
  write.csv(crash, here::here("tables", "table_compound_pycnopodia_crash.csv"), row.names = FALSE)
  pyc_share <- round(sum(pre$pyc_pre) / sum(pre$pyc_pre + pre$pis_pre), 3)
}

# ---------------------------------------------------------------------------
# (B) Keystone moderation: cascade response high vs low Pycnopodia reserves
# ---------------------------------------------------------------------------
rr <- read.csv(here::here("data", "harmonized", "harmonized_response_ratios.csv"), stringsAsFactors = FALSE)
rr <- subset(rr, !is.na(lnDiff) & year >= 2002); rr$period <- period3(rr$year); rr$yrf <- factor(rr$year)
fit_ar1 <- function(form_fixed, data) {
  nsrc <- length(unique(data$source))
  re <- if (nsrc > 1) "+ (1|source) + ar1(yrf + 0 | CA_MPA_Name_Short)" else "+ ar1(yrf + 0 | CA_MPA_Name_Short)"
  m <- safe(glmmTMB(as.formula(paste(form_fixed, "+ (1|CA_MPA_Name_Short)", re)), data = data))
  if (!is.null(m) && isTRUE(m$sdr$pdHess)) m else NULL
}
ab_contrasts <- function(m) {
  em <- emmeans(m, ~ period)
  cc <- as.data.frame(contrast(em, method = list("compound-before" = c(-1, 1, 0), "after-before" = c(-1, 0, 1))))
  list(cb = cc$estimate[1], pcb = cc[[ncol(cc)]][1], ab = cc$estimate[2], pab = cc[[ncol(cc)]][2])
}
mod_rows <- list()
if (!is.null(pyc_grp)) {
  rrm <- merge(rr, pyc_grp, by.x = "CA_MPA_Name_Short", by.y = "mpa")
  for (tx in c("Macrocystis pyrifera", "Strongylocentrotus purpuratus", "Mesocentrotus franciscanus")) {
    for (grp in c("high", "low/zero")) {
      s <- subset(rrm, y == tx & resp == pick_resp(tx) & pyc_group == grp)
      if (length(unique(s$CA_MPA_Name_Short)) < 3) next
      m <- fit_ar1("lnDiff ~ period", s)
      if (is.null(m)) m <- safe(lmer(lnDiff ~ period + (1 | CA_MPA_Name_Short), data = s, REML = TRUE))
      if (is.null(m)) next
      ct <- ab_contrasts(m)
      mod_rows[[paste(tx, grp)]] <- data.frame(taxon = tx, pyc_group = grp,
        n_mpa = length(unique(s$CA_MPA_Name_Short)),
        compound_before = round(ct$cb, 3), p_compound = signif(ct$pcb, 3),
        after_before = round(ct$ab, 3), p_after = signif(ct$pab, 3))
    }
    # period x Pycnopodia-group interaction (does the response DIFFER by group?)
    s <- subset(rrm, y == tx & resp == pick_resp(tx)); s$pyc_group <- factor(s$pyc_group, levels = c("low/zero", "high"))
    mi <- fit_ar1("lnDiff ~ period * pyc_group", s)
    if (!is.null(mi)) { cc <- summary(mi)$coefficients$cond
      ix <- grep("periodafter:pyc_grouphigh", rownames(cc))
      if (length(ix)) mod_rows[[paste(tx, "INTERACTION")]] <- data.frame(taxon = tx,
        pyc_group = "after x high (interaction)", n_mpa = length(unique(s$CA_MPA_Name_Short)),
        compound_before = NA, p_compound = NA,
        after_before = round(cc[ix, 1], 3), p_after = signif(cc[ix, 4], 3)) }
  }
}
mod_tab <- do.call(rbind, mod_rows)
write.csv(mod_tab, here::here("tables", "table_compound_keystone_moderation.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# (C) Overall resilience (lnRR, AR1) for context + figure emmeans
# ---------------------------------------------------------------------------
taxa <- c("Panulirus interruptus", "Semicossyphus pulcher", "Strongylocentrotus purpuratus",
          "Mesocentrotus franciscanus", "Macrocystis pyrifera")
short <- c("Panulirus interruptus" = "P. interruptus", "Semicossyphus pulcher" = "S. pulcher",
           "Strongylocentrotus purpuratus" = "S. purpuratus",
           "Mesocentrotus franciscanus" = "M. franciscanus", "Macrocystis pyrifera" = "M. pyrifera")
res_rows <- list()
for (tx in taxa) {
  sub <- subset(rr, y == tx & resp == pick_resp(tx)); m <- fit_ar1("lnDiff ~ period", sub)
  if (is.null(m)) m <- safe(lmer(lnDiff ~ period + (1 | CA_MPA_Name_Short), data = sub, REML = TRUE))
  ct <- ab_contrasts(m)
  res_rows[[tx]] <- data.frame(taxon = tx, compound_before = round(ct$cb, 3), p_compound = signif(ct$pcb, 3),
                               after_before = round(ct$ab, 3), p_after = signif(ct$pab, 3))
}
write.csv(do.call(rbind, res_rows), here::here("tables", "table_compound_resilience.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# (D) Cross-study synthesis (reframed)
# ---------------------------------------------------------------------------
cross <- data.frame(
  study = c("Eisaguirre et al. 2020", "Kumagai et al. 2024", "This study"),
  region = c("W. N. Channel Islands (cold)", "S. + Central California", "S. California Bight"),
  disturbance = c("SSWD (sunflower-star loss)", "2014-16 marine heatwave", "2014-16 MHW (SSWD keystone target mostly absent in SoCal)"),
  sunflower_star = c("abundant (cold reefs) -> strong keystone loss",
                     "context only (not modeled)",
                     "PATCHY: abundant at ~1/3 (cold/island) reefs, ~0 at the rest"),
  inference = c("protected predators compensate for sunflower-star loss",
                "MPAs preserve cascade -> kelp resilience to heat",
                "kelp resilient to compound stressor; the keystone-loss mechanism applies only where sunflower stars occurred -- elsewhere the response is heat + fishing protection"))
write.csv(cross, here::here("tables", "table_compound_crossstudy.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# Figure: (a) per-MPA pre-SSWD Pycnopodia distribution; (b) cascade response by group
# ---------------------------------------------------------------------------
plots <- list()
if (have_swath) {
  pd <- pre; pd$Group <- factor(pd$pyc_group, levels = c("high", "low/zero"), labels = c("had sunflower stars", "~none"))
  pd$mpa_f <- factor(pd$mpa, levels = pd$mpa[order(pd$pyc_pre)])
  plots$a <- ggplot(pd, aes(pyc_pre + 0.02, mpa_f, color = Group)) +
    geom_vline(xintercept = PYC_HI, linetype = "dotted", color = "grey50") +
    geom_point(size = 1.6) +
    scale_color_manual(values = c("had sunflower stars" = "#0072B2", "~none" = "#D55E00"), name = NULL) +
    labs(x = expression("Pre-SSWD sunflower-star density (per 60 m"^2*", log)"), y = NULL,
         title = "a  Sunflower star was patchy: present at ~1/3 of SoCal reserves, ~0 at the rest") +
    scale_x_log10() +
    theme_mpa(base_size = 8) + theme(legend.position = "bottom", axis.text.y = element_text(size = 5.5),
                                     plot.title = element_text(size = 8, face = "bold"))
}
if (!is.null(mod_tab)) {
  mt <- mod_tab[!grepl("INTERACTION|interaction", mod_tab$pyc_group), ]
  mt$lab <- short[mt$taxon]; mt$Group <- factor(mt$pyc_group, levels = c("high", "low/zero"),
                                                labels = c("had sunflower stars", "~none"))
  plots$b <- ggplot(mt, aes(after_before, lab, color = Group)) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
    geom_point(position = position_dodge(0.5), size = 2.2) +
    scale_color_manual(values = c("had sunflower stars" = "#0072B2", "~none" = "#D55E00"), name = NULL) +
    labs(x = "After-before change in MPA effect (lnRR)", y = NULL,
         title = "b  Is the cascade response stronger where sunflower stars were lost?") +
    theme_mpa(base_size = 8.5) + theme(legend.position = "bottom", axis.text.y = element_text(face = "italic"),
                                       plot.title = element_text(size = 8.5, face = "bold"))
}
fig <- if (length(plots) == 2) patchwork::wrap_plots(plots$a, plots$b, ncol = 1, heights = c(1.4, 1)) else plots[[1]]
ggsave(here::here("plots", "fig_compound_disturbance.pdf"), fig, width = 180, height = 190, units = "mm", device = cairo_pdf)
ggsave(here::here("plots", "fig_compound_disturbance.png"), fig, width = 180, height = 190, units = "mm", dpi = 600)

cat("\n=== Pre-SSWD sunflower-star (Pycnopodia): patchy distribution ===\n")
if (have_swath) { cat("Pycnopodia = only", round(100*pyc_share,1), "% of large-star density (rest is Pisaster, not an urchin predator)\n")
  cat("MPAs with sunflower stars (>=", PYC_HI, "):", sum(pre$pyc_group=="high"), "of", nrow(pre), "\n") }
cat("\n=== Cascade response: high vs low Pycnopodia reserves ===\n")
print(mod_tab, row.names = FALSE)
cat("\n  Tables + figure written.\n")
