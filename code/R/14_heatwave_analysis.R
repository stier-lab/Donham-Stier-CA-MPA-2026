# =============================================================================
# 14_heatwave_analysis.R
# =============================================================================
#
# PURPOSE:
#   Test whether, and how, the MPA-driven trophic cascade changed with the
#   2014-2016 marine heatwave (MHW) in Southern California, with a like-for-like
#   comparison to Kumagai et al. (2024, Global Change Biology).
#
# APPROACH (rigorous specification):
#   Protection is encoded in our proportion-based response ratio
#   (lnRR = ln[MPA / Reference]). Annual lnRR within an MPA is strongly
#   temporally autocorrelated (lag-1 ACF up to ~0.6 for kelp) and pools three
#   monitoring programs, so the PRIMARY models are glmmTMB mixed models with a
#   site random intercept, a Source random intercept (PISCO/KFM/LTER), and a
#   first-order autoregressive AR1(year | MPA) residual structure. Naive
#   lme4 models (no AR1) are retained only as a sensitivity to show the
#   autocorrelation-driven inflation. Diagnostics (DHARMa uniformity/dispersion/
#   outliers; within-MPA residual ACF) are written to a table.
#
#   Heatwave periods: before (<=2013) / during (2014-2016) / after (>=2017),
#   matching Kumagai et al.; MHW catalog = Hobday (2016)/heatwaveR on daily OISST
#   for the SBC region (sbc-oceanography source-of-truth; data/heatwave_exposure_
#   SBC_annual.csv). Per-MPA exposure = Kumagai's 1-km MHW intensity sampled at
#   MPA centroids (data/per_mpa_mhw_exposure.csv).
#
#   SCOPE: Southern California Bight only (all MPAs < 34.45 N; asserted).
#
# KEY RESULT: only giant kelp shows a heatwave-specific effect beyond the MPA
#   recovery trajectory; predator/urchin components are recovery-driven and, for
#   urchins, much weaker once autocorrelation is accounted for.
#
# OUTPUTS (tables/, plots/):
#   table_heatwave_period_effects.csv, table_heatwave_contrasts.csv (AR1 primary,
#   naive sensitivity), table_heatwave_diagnostics.csv, table_heatwave_cascade_
#   regression.csv, table_heatwave_per_mpa_exposure.csv, table_heatwave_pbacips_
#   integrated.csv, table_heatwave_sensitivity.csv; fig_heatwave_cascade.{pdf,png},
#   fig_heatwave_cascade_regression.{pdf,png}
#
# AUTHORS: Emily Donham & Adrian Stier
# =============================================================================

cat("Running heatwave analysis (Southern California; AR1 + Source RE)...\n")
source(here::here("code", "R", "00_libraries.R"))
source(here::here("code", "R", "00b_color_palette.R"))
source(here::here("code", "R", "01_utils.R"))
suppressMessages({library(glmmTMB); library(emmeans); library(lme4); library(lmerTest); library(DHARMa)})

taxa <- c("Panulirus interruptus", "Semicossyphus pulcher",
          "Strongylocentrotus purpuratus", "Mesocentrotus franciscanus",
          "Macrocystis pyrifera")
role <- c("Panulirus interruptus" = "Predator", "Semicossyphus pulcher" = "Predator",
          "Strongylocentrotus purpuratus" = "Herbivore",
          "Mesocentrotus franciscanus" = "Herbivore", "Macrocystis pyrifera" = "Producer")
# Use ONE response per taxon (giant kelp = biomass; animals = density) so each MPA-year
# contributes a single lnRR row -- avoids pooling the non-independent Bio/Den/RD rows and
# feeding duplicate years into ar1(year|MPA). (The cascade regression below stays density.)
pick_resp <- function(tx) if (tx == "Macrocystis pyrifera") "Bio" else "Den"

# ---------------------------------------------------------------------------
# 1. Load + join
# ---------------------------------------------------------------------------
rr   <- read.csv(here::here("data", "harmonized", "harmonized_response_ratios.csv"), stringsAsFactors = FALSE)
hw   <- read.csv(here::here("data", "heatwave_exposure_SBC_annual.csv"), stringsAsFactors = FALSE)
meta <- read.csv(here::here("data", "harmonized", "harmonized_site_metadata.csv"), stringsAsFactors = FALSE)
stopifnot("Non-Southern-California sites (lat > 34.45 N)" = all(meta$Lat <= 34.45, na.rm = TRUE))

rr <- merge(rr, hw[, c("year", "mhw_days", "mhw_icum", "period")], by = "year", all.x = TRUE)
rr <- subset(rr, !is.na(period) & !is.na(lnDiff))
rr$period <- factor(rr$period, levels = c("before", "during", "after"))
rr$yrf <- factor(rr$year)

# AR1 + Source-RE fitter (drops Source RE when a taxon has a single program)
fit_ar1 <- function(form_fixed, data) {
  nsrc <- length(unique(data$source))
  re <- if (nsrc > 1) "+ (1|source) + ar1(yrf + 0 | CA_MPA_Name_Short)" else
                      "+ ar1(yrf + 0 | CA_MPA_Name_Short)"
  f <- as.formula(paste(form_fixed, "+ (1|CA_MPA_Name_Short)", re))
  m <- tryCatch(suppressWarnings(glmmTMB(f, data = data)), error = function(e) NULL)
  if (!is.null(m) && isTRUE(m$sdr$pdHess)) m else NULL
}
pval <- function(m, term) { cc <- summary(m)$coefficients$cond
  if (term %in% rownames(cc)) cc[term, ncol(cc)] else NA }

# ---------------------------------------------------------------------------
# 2. Period model (PRIMARY: AR1 + Source RE) + emmeans contrasts; naive sensitivity
# ---------------------------------------------------------------------------
em_list <- list(); ct_list <- list()
for (tx in taxa) {
  sub <- subset(rr, y == tx & resp == pick_resp(tx))
  m  <- fit_ar1("lnDiff ~ period", sub)
  mn <- suppressWarnings(lmer(lnDiff ~ period + (1 | CA_MPA_Name_Short), data = sub, REML = TRUE))
  if (is.null(m)) { message("  [14] AR1 period model failed for ", tx, "; using naive."); m <- mn }
  em <- emmeans(m, ~ period)
  emdf <- as.data.frame(em)
  names(emdf)[names(emdf) %in% c("asymp.LCL", "lower.CL")] <- "lower"
  names(emdf)[names(emdf) %in% c("asymp.UCL", "upper.CL")] <- "upper"
  em_list[[tx]] <- transform(emdf[, c("period", "emmean", "SE", "lower", "upper")],
                             taxon = tx, n = nrow(sub),
                             n_mpa = length(unique(sub$CA_MPA_Name_Short)))
  cc  <- as.data.frame(contrast(em, method = list("during-before" = c(-1, 1, 0), "after-before" = c(-1, 0, 1))))
  ccn <- as.data.frame(contrast(emmeans(mn, ~ period),
                                method = list("during-before" = c(-1, 1, 0), "after-before" = c(-1, 0, 1))))
  ct_list[[tx]] <- data.frame(taxon = tx, contrast = cc$contrast, dlnRR = round(cc$estimate, 3),
                              p_AR1 = signif(cc[[ncol(cc)]], 3), p_naive = signif(ccn[[ncol(ccn)]], 3))
}
emm <- do.call(rbind, em_list)
ctr <- do.call(rbind, ct_list)
write.csv(ctr, here::here("tables", "table_heatwave_contrasts.csv"), row.names = FALSE)

tab <- reshape(data.frame(taxon = emm$taxon, period = emm$period, RR = round(exp(emm$emmean), 2)),
               idvar = "taxon", timevar = "period", direction = "wide")
names(tab) <- gsub("RR\\.", "RR_", names(tab)); tab$Role <- role[tab$taxon]
write.csv(tab[match(taxa, tab$taxon), c("taxon", "Role", "RR_before", "RR_during", "RR_after")],
          here::here("tables", "table_heatwave_period_effects.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 3. Diagnostics (DHARMa on the naive LMM + within-MPA residual ACF)
# ---------------------------------------------------------------------------
diag_rows <- list()
for (tx in taxa) {
  sub <- subset(rr, y == tx & resp == pick_resp(tx))
  mn <- suppressWarnings(lmer(lnDiff ~ period + (1 | CA_MPA_Name_Short), data = sub, REML = TRUE))
  sim <- simulateResiduals(mn, n = 400, plot = FALSE)
  sub$r <- resid(mn)
  acf1 <- tapply(seq_len(nrow(sub)), sub$CA_MPA_Name_Short, function(i) {
    z <- sub[i, ]; z <- z[order(z$year), ]
    if (nrow(z) > 4) acf(z$r, plot = FALSE, lag.max = 1)$acf[2] else NA_real_ })
  diag_rows[[tx]] <- data.frame(taxon = tx,
    KS_uniformity_p = signif(testUniformity(sim, plot = FALSE)$p.value, 3),
    dispersion_p = signif(testDispersion(sim, plot = FALSE)$p.value, 3),
    outlier_p = signif(testOutliers(sim, plot = FALSE)$p.value, 3),
    mean_lag1_resid_ACF = round(mean(acf1, na.rm = TRUE), 2))
}
write.csv(do.call(rbind, diag_rows), here::here("tables", "table_heatwave_diagnostics.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 4. Figure (AR1 emmeans; project theme + RR-scaled axis)
# ---------------------------------------------------------------------------
short <- c("Panulirus interruptus" = "P. interruptus", "Semicossyphus pulcher" = "S. pulcher",
           "Strongylocentrotus purpuratus" = "S. purpuratus",
           "Mesocentrotus franciscanus" = "M. franciscanus", "Macrocystis pyrifera" = "M. pyrifera")
emm$lab <- factor(short[emm$taxon], levels = short[taxa])
emm$period <- factor(emm$period, levels = c("before", "during", "after"),
                     labels = c("Before\n(≤2013)", "During\n(2014-16)", "After\n(≥2017)"))
yr <- range(c(emm$lower, emm$upper))
p <- ggplot(emm, aes(period, emmean, color = lab, group = lab)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  geom_line(linewidth = 0.5, alpha = 0.6) +
  geom_pointrange(aes(ymin = lower, ymax = upper), linewidth = 0.6, size = 0.35) +
  facet_wrap(~ lab, nrow = 1) + scale_color_taxa() +
  scale_y_rr(yr[1], yr[2], name = "MPA effect on abundance  (RR = MPA / Reference)") +
  labs(x = NULL) + theme_mpa(base_size = 9) +
  theme(legend.position = "none", axis.text.x = element_text(size = 7))
ggsave(here::here("plots", "fig_heatwave_cascade.pdf"), p, width = 180, height = 80, units = "mm", device = cairo_pdf)
ggsave(here::here("plots", "fig_heatwave_cascade.png"), p, width = 180, height = 80, units = "mm", dpi = 600)

# ---------------------------------------------------------------------------
# 5. Cascade regressions across MPA-years (density lnRR; AR1 + Source RE; cf.
#    Kumagai 2024 Fig 7)
# ---------------------------------------------------------------------------
den <- subset(rr, resp == "Den")
agg <- aggregate(lnDiff ~ CA_MPA_Name_Short + year + y + source, data = den, FUN = mean)
wide <- reshape(aggregate(lnDiff ~ CA_MPA_Name_Short + year + y, data = den, FUN = mean),
                idvar = c("CA_MPA_Name_Short", "year"), timevar = "y", direction = "wide")
names(wide) <- gsub("lnDiff\\.", "", names(wide)); wide$yrf <- factor(wide$year)
casc_fit <- function(yv, xv, ylab, xlab) {
  s <- wide[!is.na(wide[[yv]]) & !is.na(wide[[xv]]), ]; s$Y <- s[[yv]]; s$X <- s[[xv]]
  m <- tryCatch(suppressWarnings(glmmTMB(Y ~ X + (1 | CA_MPA_Name_Short) + ar1(yrf + 0 | CA_MPA_Name_Short), data = s)),
                error = function(e) NULL)
  if (is.null(m)) return(NULL)
  co <- summary(m)$coefficients$cond["X", ]
  data.frame(response = ylab, predictor = xlab, n = nrow(s), slope = round(co[1], 3),
             SE = round(co[2], 3), p = signif(co[4], 3))
}
casc <- do.call(rbind, list(
  casc_fit("Macrocystis pyrifera", "Strongylocentrotus purpuratus", "kelp", "purple urchin"),
  casc_fit("Macrocystis pyrifera", "Mesocentrotus franciscanus", "kelp", "red urchin"),
  casc_fit("Strongylocentrotus purpuratus", "Panulirus interruptus", "purple urchin", "lobster"),
  casc_fit("Strongylocentrotus purpuratus", "Semicossyphus pulcher", "purple urchin", "sheephead"),
  casc_fit("Macrocystis pyrifera", "Panulirus interruptus", "kelp", "lobster"),
  casc_fit("Macrocystis pyrifera", "Semicossyphus pulcher", "kelp", "sheephead")))
write.csv(casc, here::here("tables", "table_heatwave_cascade_regression.csv"), row.names = FALSE)

# cascade scatter figure
cs <- data.frame(x = c(wide[["Strongylocentrotus purpuratus"]], wide[["Panulirus interruptus"]]),
                 y = c(wide[["Macrocystis pyrifera"]], wide[["Strongylocentrotus purpuratus"]]),
                 panel = rep(c("kelp ~ purple urchin", "purple urchin ~ lobster"), each = nrow(wide)))
cs <- cs[is.finite(cs$x) & is.finite(cs$y), ]
pc <- ggplot(cs, aes(x, y)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_point(alpha = 0.35, size = 0.9, colour = "#0072B2") +
  geom_smooth(method = "lm", se = TRUE, colour = "#D55E00", linewidth = 0.7) +
  facet_wrap(~ panel, scales = "free") +
  labs(x = "Predictor lnRR (MPA / Reference)", y = "Response lnRR (MPA / Reference)") +
  theme_mpa(base_size = 9)
ggsave(here::here("plots", "fig_heatwave_cascade_regression.pdf"), pc, width = 140, height = 75, units = "mm", device = cairo_pdf)
ggsave(here::here("plots", "fig_heatwave_cascade_regression.png"), pc, width = 140, height = 75, units = "mm", dpi = 600)

# ---------------------------------------------------------------------------
# 6. Per-MPA continuous exposure + integrated recovery-vs-heatwave model (AR1)
# ---------------------------------------------------------------------------
expo_path <- here::here("data", "per_mpa_mhw_exposure.csv")
if (file.exists(expo_path)) {
  expo <- read.csv(expo_path, stringsAsFactors = FALSE)
  je <- merge(rr, expo, by = c("CA_MPA_Name_Short", "year"))
  je <- subset(je, !is.na(mhw_icum_mpa) & !is.na(lnDiff))

  ex_rows <- list(); int_rows <- list()
  for (tx in taxa) {
    s <- subset(je, y == tx & resp == pick_resp(tx)); s$z <- as.numeric(scale(s$mhw_icum_mpa))
    m <- fit_ar1("lnDiff ~ z", s)
    if (!is.null(m)) { co <- summary(m)$coefficients$cond["z", ]
      ex_rows[[tx]] <- data.frame(Taxon = tx, Role = role[tx], n = nrow(s),
        slope_per_SD = round(co[1], 3), SE = round(co[2], 3), p = signif(co[4], 3)) }
    si <- subset(s, BA == "After" & !is.na(time))
    mi <- fit_ar1("lnDiff ~ time + z", si)
    if (!is.null(mi)) { co <- summary(mi)$coefficients$cond
      int_rows[[tx]] <- data.frame(Taxon = tx, Role = role[tx], n = nrow(si),
        recovery_slope = round(co["time", 1], 3), recovery_p = signif(co["time", 4], 3),
        heatwave_slope_perSD = round(co["z", 1], 3), heatwave_p = signif(co["z", 4], 3)) }
  }
  write.csv(do.call(rbind, ex_rows), here::here("tables", "table_heatwave_per_mpa_exposure.csv"), row.names = FALSE)
  write.csv(do.call(rbind, int_rows), here::here("tables", "table_heatwave_pbacips_integrated.csv"), row.names = FALSE)
}

# ---------------------------------------------------------------------------
# 7. Sensitivity of the after-vs-before contrast (AR1) to MPA set / panel balance
# ---------------------------------------------------------------------------
sheephead_only <- c("Blue Cavern Onshore SMCA", "Farnsworth Onshore SMCA", "Swamis SMCA",
                    "Long Point SMR", "Point Dume SMR", "Point Dume SMCA", "Painted Cave SMCA",
                    "Cat Harbor SMCA", "Dana Point SMCA", "Anacapa Island SMCA", "Santa Barbara Island SMR")
ab_p <- function(sub) { m <- fit_ar1("lnDiff ~ period", sub); if (is.null(m)) return(c(NA, NA))
  cc <- as.data.frame(contrast(emmeans(m, ~ period), method = list("after-before" = c(-1, 0, 1))))
  c(cc$estimate, cc[[ncol(cc)]]) }
scen <- list(all_MPAs = function(x) x,
  balanced_3periods = function(x) { keep <- names(which(tapply(x$period, x$CA_MPA_Name_Short, function(z) length(unique(z))) == 3)); x[x$CA_MPA_Name_Short %in% keep, ] },
  SMR_only = function(x) x[grepl("SMR", x$CA_MPA_Name_Short), ],
  drop_sheephead_only = function(x) x[!x$CA_MPA_Name_Short %in% sheephead_only, ]
)
sens_rows <- list()
for (tx in c("Panulirus interruptus", "Strongylocentrotus purpuratus", "Macrocystis pyrifera"))
  for (nm in names(scen)) { sub <- scen[[nm]](subset(rr, y == tx & resp == pick_resp(tx))); r <- ab_p(sub)
    sens_rows[[paste(tx, nm)]] <- data.frame(Taxon = tx, scenario = nm, n = nrow(sub),
      n_mpa = length(unique(sub$CA_MPA_Name_Short)), after_before_dlnRR = round(r[1], 3), p = signif(r[2], 3)) }
write.csv(do.call(rbind, sens_rows), here::here("tables", "table_heatwave_sensitivity.csv"), row.names = FALSE)

cat("  Heatwave analysis complete (AR1 + Source RE primary). Tables + figures written to tables/ and plots/.\n")
