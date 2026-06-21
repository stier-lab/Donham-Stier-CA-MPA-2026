# =============================================================================
# 15_methods_comparison.R
# =============================================================================
#
# PURPOSE:
#   Formal, quantitative comparison of HOW the method of analysis drives
#   variation in MPA x marine-heatwave effect sizes, contrasting our
#   paired-reference / proportion-based log-response-ratio + multilevel
#   inverse-variance meta-analysis against the Kumagai et al. (2024, Global
#   Change Biology) pooled-protection / raw-density GLMM approach.
#
#   This is the supporting-information "analytical multiverse". It separates
#   variation attributable to METHOD from variation attributable to DATA by:
#
#   (A) METHOD BRIDGE  -- on a single common dataset (OURS), walk from the
#       Kumagai-style endpoint to our endpoint, flipping exactly ONE analytical
#       choice per step. The change in effect size at each step = the
#       contribution of that one methodological axis.
#
#         WP0  pooled  | raw density | Tweedie GLMM (log link)   = Kumagai-style
#         WP1  pooled  | proportion  | Tweedie GLMM (log link)   (flip: response scale)
#         WP2  pooled  | proportion  | Gaussian LMM on log       (flip: distribution/link)
#         WP3  paired  | proportion  | LMM on lnRR (no AR1)      (flip: reference/pairing)
#         WP4  paired  | proportion  | LMM on lnRR + AR1(year|MPA) + Source RE
#                                                                (flip: temporal autocorrelation)
#                                                                = our heatwave endpoint (script 14)
#         WP5  paired  | proportion  | inverse-variance meta-analysis (rma.mv)
#                                                                (flip: weighting) = our main-paper synthesis
#
#   The WP3->WP4 autocorrelation flip is the single largest lever in our own
#   analysis (script 14): annual lnRR is strongly autocorrelated within MPAs, and
#   ignoring it inflates significance by many orders of magnitude.
#
#   (B) DATA EFFECT  -- run the pooled waypoints (WP0-WP2, the only ones the
#       Kumagai dataset structurally supports) on BOTH substrates:
#         OURS   = our harmonized PISCO+KFM+LTER data (MPA-group level)
#         THEIRS = Kumagai MLPA_data_summarized_wo_siteblocks.csv (site level,
#                  Southern Coast, Full vs Reference)
#       The OURS-vs-THEIRS gap at a matched waypoint = the DATA contribution.
#
#   Common currency for every cell: the during-before and after-before CHANGE
#   in ln(MPA / Reference) abundance, per taxon, with 95% CI. Temporal contrast
#   is fixed to the heatwave windows (before <=2013 / during 2014-16 / after
#   >=2017, matching Kumagai) so all cells answer the same question. Density
#   only (THEIRS is all density), Southern California only.
#
# OUTPUTS (tables/, plots/, docs/):
#   table_s_methods_crosswalk.csv      - qualitative axis-by-axis crosswalk
#   table_s_methods_multiverse.csv     - every waypoint x substrate x taxon estimate
#   table_s_methods_decomposition.csv  - per-axis contribution (method) + data effect
#   fig_s_methods_multiverse.{pdf,png} - bridge waterfall + cross-substrate panel
#   docs/methods_comparison_supplement.md is maintained separately (prose/legend)
#
# AUTHORS: Emily Donham & Adrian Stier
# =============================================================================

cat("Running methods comparison / analytical multiverse (script 15)...\n")
source(here::here("code", "R", "00_libraries.R"))
source(here::here("code", "R", "00b_color_palette.R"))
source(here::here("code", "R", "01_utils.R"))
suppressMessages({library(glmmTMB); library(lme4); library(lmerTest); library(metafor)})

taxa <- c("Panulirus interruptus", "Semicossyphus pulcher",
          "Strongylocentrotus purpuratus", "Mesocentrotus franciscanus",
          "Macrocystis pyrifera")
role <- c("Panulirus interruptus" = "Predator", "Semicossyphus pulcher" = "Predator",
          "Strongylocentrotus purpuratus" = "Herbivore",
          "Mesocentrotus franciscanus" = "Herbivore", "Macrocystis pyrifera" = "Producer")
short <- c("Panulirus interruptus" = "P. interruptus", "Semicossyphus pulcher" = "S. pulcher",
           "Strongylocentrotus purpuratus" = "S. purpuratus",
           "Mesocentrotus franciscanus" = "M. franciscanus", "Macrocystis pyrifera" = "M. pyrifera")

# Heatwave-period coding identical to Kumagai et al. (calendar windows)
hw_period <- function(year) factor(ifelse(year <= 2013, "before",
                                    ifelse(year <= 2016, "during", "after")),
                                   levels = c("before", "during", "after"))
# small constant for log() of a vector that may contain zeros (half min nonzero)
log_eps <- function(x) { nz <- x[x > 0 & is.finite(x)]; e <- if (length(nz)) min(nz)/2 else 1e-6; log(x + e) }

# ---------------------------------------------------------------------------
# Coefficient extractors -> tidy (during-before, after-before) change in lnRR
#   POOLED models are coded ~ status*period with status = {reference, mpa} and
#   period = {before, during, after}; the interaction coefficients ARE the
#   during/after-vs-before changes in lnRR, and `statusmpa` is the before lnRR.
#   PAIRED models are coded ~ period on the per-MPA lnRR; the period
#   coefficients are the changes and the intercept is the before lnRR.
# ---------------------------------------------------------------------------
tidy_pooled <- function(cmat, est = 1, se = 2, p = ncol(cmat)) {
  g <- function(r, j) if (r %in% rownames(cmat)) cmat[r, j] else NA_real_
  b  <- g("statusmpa", est)
  data.frame(
    lnRR_before = b,
    lnRR_during = b + g("statusmpa:periodduring", est),
    lnRR_after  = b + g("statusmpa:periodafter",  est),
    dDB = g("statusmpa:periodduring", est), dDB_se = g("statusmpa:periodduring", se), dDB_p = g("statusmpa:periodduring", p),
    dAB = g("statusmpa:periodafter",  est), dAB_se = g("statusmpa:periodafter",  se), dAB_p = g("statusmpa:periodafter",  p))
}
tidy_paired <- function(cmat, est = 1, se = 2, p = ncol(cmat)) {
  g <- function(r, j) if (r %in% rownames(cmat)) cmat[r, j] else NA_real_
  b <- g("(Intercept)", est)
  data.frame(
    lnRR_before = b,
    lnRR_during = b + g("periodduring", est),
    lnRR_after  = b + g("periodafter",  est),
    dDB = g("periodduring", est), dDB_se = g("periodduring", se), dDB_p = g("periodduring", p),
    dAB = g("periodafter",  est), dAB_se = g("periodafter",  se), dAB_p = g("periodafter",  p))
}
safe <- function(expr) tryCatch(suppressWarnings(suppressMessages(expr)), error = function(e) NULL)

# ===========================================================================
# Build the two substrates in a common long format:
#   columns: taxon, MPA, site, year, period, status (reference/mpa), raw, prop
#   (one row per observational unit x status). OURS is at MPA-group level
#   (mpa vs reference), THEIRS is at site level (Full vs Reference).
# ===========================================================================

# --- OURS ------------------------------------------------------------------
# RAW: harmonized_raw_responses (status = mpa/reference, value = raw density)
oraw <- read.csv(here::here("data", "harmonized", "harmonized_raw_responses.csv"), stringsAsFactors = FALSE)
oraw <- subset(oraw, resp == "Den" & year >= 2002 & taxon_name %in% taxa)
oraw <- aggregate(value ~ CA_MPA_Name_Short + year + taxon_name + status, data = oraw, FUN = mean)
# PROP + paired lnRR: harmonized_response_ratios (mpa, reference = proportions; lnDiff = paired lnRR)
orr <- read.csv(here::here("data", "harmonized", "harmonized_response_ratios.csv"), stringsAsFactors = FALSE)
orr <- subset(orr, resp == "Den" & year >= 2002 & y %in% taxa)
orr$period <- hw_period(orr$year)

ours_pooled_raw <- transform(oraw, period = hw_period(year),
                             status = factor(status, levels = c("reference", "mpa")),
                             taxon = taxon_name, val = value)
# pooled proportion: melt the mpa/reference proportion columns to long status
ours_pooled_prop <- rbind(
  data.frame(taxon = orr$y, CA_MPA_Name_Short = orr$CA_MPA_Name_Short, year = orr$year,
             period = orr$period, status = "mpa",       val = orr$mpa),
  data.frame(taxon = orr$y, CA_MPA_Name_Short = orr$CA_MPA_Name_Short, year = orr$year,
             period = orr$period, status = "reference", val = orr$reference))
ours_pooled_prop$status <- factor(ours_pooled_prop$status, levels = c("reference", "mpa"))
# paired proportion lnRR (one value per MPA-year-taxon); keep source for AR1/Source-RE
ours_paired <- data.frame(taxon = orr$y, MPA = orr$CA_MPA_Name_Short, year = orr$year,
                          period = orr$period, source = orr$source, lnRR = orr$lnDiff)
ours_paired <- subset(ours_paired, is.finite(lnRR))

# --- THEIRS (Kumagai South Coast; Full vs Reference) -----------------------
theirs_path <- "~/kumagai2024-comparison/repo/Processed_data/MLPA_data_summarized_wo_siteblocks.csv"
theirs_cols <- c("Panulirus interruptus" = "PANINT_d", "Semicossyphus pulcher" = "SPUL_d",
                 "Strongylocentrotus purpuratus" = "STRPURAD_d",
                 "Mesocentrotus franciscanus" = "MESFRAAD_d", "Macrocystis pyrifera" = "MACPYRAD_d")
have_theirs <- file.exists(path.expand(theirs_path))
if (have_theirs) {
  kd <- read.csv(path.expand(theirs_path), stringsAsFactors = FALSE)
  kd <- subset(kd, region == "South_Coast" & mpa_status %in% c("Full", "Reference") & year >= 2002)
  kd$period <- hw_period(kd$year)
  kd$status <- factor(ifelse(kd$mpa_status == "Full", "mpa", "reference"),
                      levels = c("reference", "mpa"))
  kd$site_total <- rowSums(kd[, theirs_cols], na.rm = TRUE)   # for proportion scale
  theirs_long <- do.call(rbind, lapply(taxa, function(tx) {
    data.frame(taxon = tx, site = kd$site, year = kd$year, period = kd$period,
               status = kd$status, raw = kd[[theirs_cols[tx]]],
               prop = ifelse(kd$site_total > 0, kd[[theirs_cols[tx]]] / kd$site_total, NA_real_))
  }))
} else message("  [15] Kumagai data not found at ", theirs_path, "; cross-substrate panel skipped.")

# ===========================================================================
# (A) METHOD BRIDGE on OURS  -- one analytical flip per waypoint
# ===========================================================================
fit_tweedie_pooled <- function(d) {  # WP0/WP1: pooled Tweedie GLMM on `val` (raw or proportion)
  m <- safe(glmmTMB(val ~ status * period + (1 | CA_MPA_Name_Short), data = d, family = tweedie(link = "log")))
  if (is.null(m) || !isTRUE(m$sdr$pdHess)) return(NULL)
  tidy_pooled(summary(m)$coefficients$cond)
}
fit_lmm_pooled <- function(d) {  # WP2: pooled Gaussian LMM on log(val)
  d$lv <- log_eps(d$val)
  m <- safe(lmer(lv ~ status * period + (1 | CA_MPA_Name_Short), data = d, REML = TRUE))
  if (is.null(m)) return(NULL)
  tidy_pooled(coef(summary(m)))
}
fit_lmm_paired <- function(d) {  # WP3: paired LMM on lnRR (no AR1)
  m <- safe(lmer(lnRR ~ period + (1 | MPA), data = d, REML = TRUE))
  if (is.null(m)) return(NULL)
  tidy_paired(coef(summary(m)))
}
fit_ar1_paired <- function(d) {  # WP4: paired LMM on lnRR + AR1(year|MPA) + Source RE (script-14 spec)
  d$yrf <- factor(d$year)
  nsrc <- length(unique(d$source))
  re <- if (nsrc > 1) "+ (1|source) + ar1(yrf + 0 | MPA)" else "+ ar1(yrf + 0 | MPA)"
  f <- as.formula(paste("lnRR ~ period + (1|MPA)", re))
  m <- safe(glmmTMB(f, data = d))
  if (is.null(m) || !isTRUE(m$sdr$pdHess)) return(NULL)
  tidy_paired(summary(m)$coefficients$cond)
}
fit_meta_paired <- function(d) {  # WP5: paired, inverse-variance multilevel meta-analysis
  ag <- do.call(rbind, lapply(split(d, list(d$MPA, d$period), drop = TRUE), function(s) {
    if (nrow(s) < 2) return(NULL)
    data.frame(MPA = s$MPA[1], period = s$period[1], yi = mean(s$lnRR), vi = var(s$lnRR) / nrow(s))
  }))
  ag <- subset(ag, is.finite(yi) & is.finite(vi) & vi > 0)
  if (is.null(ag) || nrow(ag) < 4 || length(unique(ag$period)) < 2) return(NULL)
  ag$period <- factor(ag$period, levels = c("before", "during", "after"))
  m <- safe(rma.mv(yi = yi, V = vi, mods = ~ period, random = ~ 1 | MPA, data = ag, method = "REML"))
  if (is.null(m)) return(NULL)
  cmat <- cbind(est = as.numeric(m$b), se = m$se, p = m$pval)
  rownames(cmat) <- rownames(m$b)
  rownames(cmat)[rownames(cmat) == "intrcpt"] <- "(Intercept)"
  tidy_paired(cmat, est = 1, se = 2, p = 3)
}

waypoints <- c("WP0 pooled|raw|Tweedie-GLMM", "WP1 pooled|prop|Tweedie-GLMM",
               "WP2 pooled|prop|LMM-log", "WP3 paired|prop|LMM-lnRR",
               "WP4 paired|prop|AR1-LMM", "WP5 paired|prop|IV-meta")
bridge <- list()
for (tx in taxa) {
  d0 <- subset(ours_pooled_raw,  taxon == tx)              # val = raw density
  d1 <- subset(ours_pooled_prop, taxon == tx)              # val = proportion
  d3 <- subset(ours_paired,      taxon == tx)              # lnRR
  res <- list(WP0 = fit_tweedie_pooled(d0), WP1 = fit_tweedie_pooled(d1),
              WP2 = fit_lmm_pooled(d1), WP3 = fit_lmm_paired(d3),
              WP4 = fit_ar1_paired(d3), WP5 = fit_meta_paired(d3))
  for (k in seq_along(waypoints)) {
    r <- res[[k]]; if (is.null(r)) next
    bridge[[length(bridge) + 1]] <- cbind(substrate = "OURS", waypoint = waypoints[k],
                                          taxon = tx, role = role[tx], r)
  }
}
bridge <- do.call(rbind, bridge)

# ===========================================================================
# (B) DATA EFFECT  -- pooled waypoints on THEIRS (the cells their data supports)
# ===========================================================================
fit_theirs_tweedie <- function(d, col) {
  d$v <- d[[col]]
  m <- safe(glmmTMB(v ~ status * period + (1 | site), data = d, family = tweedie(link = "log")))
  if (is.null(m) || !isTRUE(m$sdr$pdHess)) return(NULL)
  tidy_pooled(summary(m)$coefficients$cond)
}
fit_theirs_lmm <- function(d, col) {
  d$lv <- log_eps(d[[col]])
  m <- safe(lmer(lv ~ status * period + (1 | site), data = d, REML = TRUE))
  if (is.null(m)) return(NULL)
  tidy_pooled(coef(summary(m)))
}
theirs_rows <- list()
if (have_theirs) {
  twp <- c("WP0 pooled|raw|Tweedie-GLMM", "WP1 pooled|prop|Tweedie-GLMM", "WP2 pooled|prop|LMM-log")
  for (tx in taxa) {
    d <- subset(theirs_long, taxon == tx & is.finite(raw))
    dp <- subset(d, is.finite(prop))
    res <- list(fit_theirs_tweedie(d, "raw"), fit_theirs_tweedie(dp, "prop"),
                fit_theirs_lmm(dp, "prop"))
    for (k in seq_along(twp)) {
      r <- res[[k]]; if (is.null(r)) next
      theirs_rows[[length(theirs_rows) + 1]] <- cbind(substrate = "THEIRS", waypoint = twp[k],
                                                      taxon = tx, role = role[tx], r)
    }
  }
}
multiverse <- if (length(theirs_rows)) rbind(bridge, do.call(rbind, theirs_rows)) else bridge
rownames(multiverse) <- NULL

# rounded export
num <- sapply(multiverse, is.numeric)
mv_out <- multiverse; mv_out[num] <- lapply(mv_out[num], function(x) signif(x, 4))
write.csv(mv_out, here::here("tables", "table_s_methods_multiverse.csv"), row.names = FALSE)

# ===========================================================================
# (C) DECOMPOSITION -- per-axis contribution to the after-before change in lnRR
# ===========================================================================
ab <- multiverse[, c("substrate", "waypoint", "taxon", "dAB", "dAB_p")]
ab$wp <- substr(ab$waypoint, 1, 3)
oab <- subset(ab, substrate == "OURS")
wide_ab <- reshape(oab[, c("taxon", "wp", "dAB")],   idvar = "taxon", timevar = "wp", direction = "wide")
wide_p  <- reshape(oab[, c("taxon", "wp", "dAB_p")], idvar = "taxon", timevar = "wp", direction = "wide")
names(wide_ab) <- gsub("dAB\\.", "", names(wide_ab))
names(wide_p)  <- gsub("dAB_p\\.", "", names(wide_p))
# per axis: mean/max |change in point estimate|, and # taxa whose significance (a=0.05) flips
method_axis <- function(a, b) {
  d <- abs(wide_ab[[b]] - wide_ab[[a]]); d <- d[is.finite(d)]
  sa <- wide_p[[a]] < 0.05; sb <- wide_p[[b]] < 0.05
  flips <- sum(xor(sa, sb), na.rm = TRUE)
  c(mean = mean(d), max = max(d), flips = flips)
}
axes <- rbind(
  data.frame(contribution = "METHOD: response scale (raw density -> proportion)",        t(method_axis("WP0", "WP1"))),
  data.frame(contribution = "METHOD: distribution/link (Tweedie GLMM -> Gaussian log-LMM)", t(method_axis("WP1", "WP2"))),
  data.frame(contribution = "METHOD: reference (pooled -> paired)",                       t(method_axis("WP2", "WP3"))),
  data.frame(contribution = "METHOD: temporal autocorrelation (no AR1 -> AR1+SourceRE)",  t(method_axis("WP3", "WP4"))),
  data.frame(contribution = "METHOD: weighting (AR1 LMM -> IV meta-analysis)",            t(method_axis("WP4", "WP5")))
)
# DATA effect: OURS vs THEIRS gap at matched pooled waypoints (WP0-WP2)
if (length(theirs_rows)) {
  g <- merge(subset(ab, substrate == "OURS"), subset(ab, substrate == "THEIRS"),
             by = c("taxon", "wp"), suffixes = c("_ours", "_theirs"))
  gd <- abs(g$dAB_ours - g$dAB_theirs); gd <- gd[is.finite(gd)]
  gflip <- sum(xor(g$dAB_p_ours < 0.05, g$dAB_p_theirs < 0.05), na.rm = TRUE)
  axes <- rbind(axes, data.frame(contribution = "DATA: substrate (ours vs Kumagai) at matched pooled waypoints",
                                 mean = mean(gd), max = max(gd), flips = gflip))
}
axes$mean <- round(axes$mean, 3); axes$max <- round(axes$max, 3)
names(axes) <- c("contribution", "mean_abs_delta_lnRR", "max_abs_delta_lnRR", "n_taxa_significance_flips")
write.csv(axes, here::here("tables", "table_s_methods_decomposition.csv"), row.names = FALSE)

# ===========================================================================
# (D) QUALITATIVE CROSSWALK TABLE
# ===========================================================================
crosswalk <- data.frame(
  axis = c("Protection coding", "Reference structure", "Response scale", "Effect metric",
           "Error model / link", "Temporal contrast", "Synthesis / weighting",
           "Heterogeneity", "Kelp resilience", "Multiple testing"),
  kumagai_2024 = c(
    "3-level Reference / Partial / Full",
    "Pooled protected-vs-unprotected within region (no pairing)",
    "Raw density (count per 60 m2) with offset(log n_transects)",
    "Model marginal-mean contrast on the link scale",
    "glmmTMB nbinom1/nbinom2/Tweedie, log link; (1 + year_std | site)",
    "Heatwave windows before/during/after (also TWFE on years-since-establishment)",
    "Single one-stage GLMM (population-level, conditional on site membership)",
    "Absorbed by site random intercept/slope; MPAs exchangeable within region",
    "Separate pixel-level permutation test on satellite canopy (10,000 perms)",
    "AIC model selection; Wald / permutation p-values"),
  donham_stier = c(
    "Binary MPA vs designated reference",
    "Each MPA paired to its own reference (within-pair difference)",
    "Proportion of site total (composition-standardized)",
    "Log response ratio ln(MPA / Reference)",
    "Gaussian on lnRR; glmmTMB with AR1(year|MPA) + Source RE (script 14)",
    "pBACIPS time-since-establishment (main) and heatwave windows (script 14)",
    "Two-stage inverse-variance multilevel meta-analysis (rma.mv)",
    "Explicit (1|MPA),(1|Source) variance components; tau2, I2 reported",
    "Same lnRR framework as all other taxa (no separate test)",
    "Benjamini-Hochberg FDR across meta-analytic p-values"),
  expected_consequence = c(
    "Collapsing Partial changes the contrast set but not direction",
    "Pairing removes among-site composition; lowers pre-heatwave baseline RR for lobster/kelp",
    "Proportion nets out shared site-level abundance trends; shrinks extreme RRs",
    "lnRR is symmetric and unit-free; comparable across taxa/programs",
    "Count vs Gaussian-on-log mainly affects SE/zeros, less the point estimate",
    "Establishment vs heatwave timing separates recovery from acute MHW effect",
    "IV weighting downweights noisy small-n MPAs; widens CIs vs one-stage GLMM",
    "Explicit components reveal high between-MPA heterogeneity (large I2)",
    "Unified framework trades permutation robustness for cross-taxa comparability",
    "FDR is more conservative than uncorrected Wald p-values"))
write.csv(crosswalk, here::here("tables", "table_s_methods_crosswalk.csv"), row.names = FALSE)

# ===========================================================================
# (E) FIGURE -- bridge waterfall (after-before change) + cross-substrate panel
# ===========================================================================
mv <- multiverse
mv$wp <- substr(mv$waypoint, 1, 3)
mv$lab <- factor(short[mv$taxon], levels = short[taxa])
mv$lo <- mv$dAB - 1.96 * mv$dAB_se
mv$hi <- mv$dAB + 1.96 * mv$dAB_se

# Panel a: the method bridge on OURS (WP0 -> WP4)
br <- subset(mv, substrate == "OURS")
br$wp <- factor(br$wp, levels = paste0("WP", 0:5))
endpts <- ifelse(br$wp == "WP0", "Kumagai-style",
                 ifelse(br$wp %in% c("WP4", "WP5"), "Our endpoint", "intermediate"))
br$endpoint <- factor(endpts, levels = c("Kumagai-style", "intermediate", "Our endpoint"))
yr1 <- range(c(br$lo, br$hi), na.rm = TRUE)
pa <- ggplot(br, aes(wp, dAB, group = lab)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  geom_line(linewidth = 0.4, color = "grey60") +
  geom_pointrange(aes(ymin = lo, ymax = hi, color = endpoint), linewidth = 0.5, size = 0.3) +
  facet_wrap(~ lab, nrow = 1) +
  scale_color_manual(values = c("Kumagai-style" = "#D55E00", "intermediate" = "grey50",
                                "Our endpoint" = "#0072B2"), name = NULL) +
  scale_y_rr(yr1[1], yr1[2], name = "After-before change in MPA effect  (x-fold change in RR)") +
  labs(x = NULL, title = "a  Method bridge (common dataset): one analytical flip per step") +
  theme_mpa(base_size = 9) +
  theme(legend.position = "bottom", axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
        plot.title = element_text(size = 9, face = "bold"))

# Panel b: cross-substrate at the shared pooled waypoints (data effect)
plist <- list(pa)
if (length(theirs_rows)) {
  cs <- subset(mv, wp %in% c("WP0", "WP1", "WP2"))
  cs$wp <- factor(cs$wp, levels = c("WP0", "WP1", "WP2"))
  pb <- ggplot(cs, aes(wp, dAB, color = substrate, group = substrate)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
    geom_line(position = position_dodge(0.3), linewidth = 0.4) +
    geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(0.3),
                    linewidth = 0.5, size = 0.3) +
    facet_wrap(~ lab, nrow = 1) +
    scale_color_manual(values = c("OURS" = "#0072B2", "THEIRS" = "#009E73"), name = NULL) +
    scale_y_rr(min(c(cs$lo, 0), na.rm = TRUE), max(c(cs$hi, 0), na.rm = TRUE),
               name = "After-before change in MPA effect  (x-fold change in RR)") +
    labs(x = NULL, title = "b  Data effect: same pooled methods on each dataset") +
    theme_mpa(base_size = 9) +
    theme(legend.position = "bottom", axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
          plot.title = element_text(size = 9, face = "bold"))
  plist <- list(pa, pb)
}
fig <- if (length(plist) == 2) patchwork::wrap_plots(plist, ncol = 1) else plist[[1]]
ggsave(here::here("plots", "fig_s_methods_multiverse.pdf"), fig,
       width = 180, height = if (length(plist) == 2) 170 else 95, units = "mm", device = cairo_pdf)
ggsave(here::here("plots", "fig_s_methods_multiverse.png"), fig,
       width = 180, height = if (length(plist) == 2) 170 else 95, units = "mm", dpi = 600)

cat("  Methods comparison complete. Multiverse:", nrow(multiverse), "cells across",
    length(unique(multiverse$waypoint)), "waypoints.\n")
cat("  Tables: table_s_methods_{crosswalk,multiverse,decomposition}.csv; figure: fig_s_methods_multiverse.\n")
