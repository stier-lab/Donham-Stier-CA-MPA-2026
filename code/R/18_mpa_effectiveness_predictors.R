# =============================================================================
# 18_mpa_effectiveness_predictors.R
# =============================================================================
#
# PURPOSE:
#   "Dig in hard on relevant predictors and predict the variation in MPA
#   effectiveness." Script 16 screened environmental moderators one at a time
#   (null after FDR). Here we assemble the FULLEST per-MPA predictor set
#   (environmental + reserve design + trophic context) and ask, with a multi-
#   predictor model AND an honest leave-one-out cross-validated R^2, how much of
#   the among-MPA variation in effectiveness we can actually explain and predict.
#
# EFFECTIVENESS TARGETS (per-MPA effect sizes, lnRR = ln[MPA / Reference]):
#   - kelp_eff   : giant kelp (M. pyrifera) BIOMASS lnRR  -- the foundation-species
#                  outcome, the cleanest "did this reserve restore kelp?" metric (k~31)
#   - purp_supp  : purple urchin (S. purpuratus) DENSITY lnRR (negative = suppression)
#
# PREDICTORS (per MPA):
#   Environmental : MHW intensity during 2014-16, cold-spell intensity, latitude
#   Reserve design: log area (Hectares), SMR vs SMCA, island vs mainland, establishment year
#   Trophic       : per-MPA predator/urchin density lnRR (lobster, sheephead, urchins)
#                   -- tests whether the cascade itself predicts kelp effectiveness
#
# MODELS / HONESTY:
#   (1) univariate inverse-variance meta-regression screen (Knapp-Hartung) per predictor
#   (2) multi-predictor meta-regression (env + design) with a pseudo-R^2
#   (3) random forest variable importance + out-of-bag R^2 (nonlinear, all predictors)
#   (4) LEAVE-ONE-OUT CV R^2 of a parsimonious model -- the real test of whether we
#       can *predict* (not just fit) MPA effectiveness. With k~20-31 MPAs this is
#       the arbiter; in-sample fits will look better than they predict.
#
# OUTPUTS:
#   tables/table_s_effectiveness_predictors.csv  - univariate screen (both targets)
#   tables/table_s_effectiveness_models.csv      - multi-predictor + RF + LOO-CV summary
#   plots/fig_s_effectiveness_predictors.{pdf,png} - variable importance + obs vs LOO-pred
#
# NOTE: standalone (reads the Kumagai cold-spell grid like script 16; skips cs if absent).
#
# AUTHORS: Emily Donham & Adrian Stier
# =============================================================================

cat("Modeling predictors of MPA effectiveness (script 18)...\n")
source(here::here("code", "R", "00_libraries.R"))
source(here::here("code", "R", "00b_color_palette.R"))
source(here::here("code", "R", "00c_analysis_constants.R"))
source(here::here("code", "R", "01_utils.R"))
suppressMessages({library(metafor); library(data.table)})
have_rf <- requireNamespace("randomForest", quietly = TRUE)
safe <- function(e) tryCatch(suppressWarnings(suppressMessages(e)), error = function(x) NULL)
external_file <- function(envvar, repo_rel, fallback) {
  override <- Sys.getenv(envvar, unset = "")
  if (nzchar(override)) return(path.expand(override))
  local <- here::here(repo_rel)
  if (file.exists(local)) return(local)
  path.expand(fallback)
}

# ---------------------------------------------------------------------------
# 1. Assemble the per-MPA predictor + target table
# ---------------------------------------------------------------------------
ss <- as.data.table(read.csv(here::here("data", "sumstats_final.csv"), stringsAsFactors = FALSE))
ss <- ss[is.finite(Mean) & is.finite(SE) & SE > 0]
wide_eff <- function(taxon, resp) {
  x <- ss[Taxa == taxon & Resp == resp, .(MPA, Mean, SE)]
  x[!duplicated(MPA)]
}
kelp  <- wide_eff("M. pyrifera", "Bio");  setnames(kelp,  c("Mean","SE"), c("kelp_eff","kelp_se"))
purp  <- wide_eff("S. purpuratus", "Den"); setnames(purp,  c("Mean","SE"), c("purp_eff","purp_se"))
lob   <- wide_eff("P. interruptus", "Den")[, .(MPA, lob_eff = Mean)]
sheep <- wide_eff("S. pulcher", "Den")[, .(MPA, sheep_eff = Mean)]
red   <- wide_eff("M. franciscanus", "Den")[, .(MPA, red_eff = Mean)]

# design + location from sumstats (one row per MPA)
meta <- ss[!duplicated(MPA), .(MPA, lat = Lat, type, MPA_Start, island = as.integer(!is.na(ChannelIsland)))]
meta[, `:=`(smr = as.integer(type == "SMR"), est_year = MPA_Start)]

# reserve size
rr <- as.data.table(load_harmonized_rr())
size <- unique(rr[, .(MPA = CA_MPA_Name_Short, log_size = log(Hectares))])[!duplicated(MPA)]

# environmental: MHW intensity during MHW1 (2014-16)
mhw <- as.data.table(read.csv(here::here("data", "per_mpa_mhw_exposure.csv"), stringsAsFactors = FALSE))
mhw_d <- mhw[year %in% MHW1_YEARS, .(mhw_during = mean(mhw_icum_mpa)), by = .(MPA = CA_MPA_Name_Short)]
# cold-spell climatology (nearest 1-km cell), if Kumagai grid present
cs_path <- external_file(
  "KUMAGAI_CS_GRID",
  file.path("data", "external", "kumagai2024", "CS_cummulative_intensity_1km.rds"),
  "~/kumagai2024-comparison/repo/Processed_data/SST/CS_cummulative_intensity_1km.rds"
)
meta[, lon := ss[!duplicated(MPA)]$Lon]   # lon for cs nearest-cell match
if (file.exists(cs_path)) {
  cs <- as.data.table(readRDS(cs_path)); csagg <- cs[, .(cs = mean(CS_cummulative)), by = .(long, lat)]
  meta[, cs_mean := vapply(seq_len(.N), function(i) {
    d <- (csagg$lat - lat[i])^2 + (csagg$long - lon[i])^2; csagg$cs[which.min(d)] }, numeric(1))]
} else meta[, cs_mean := NA_real_]

# merge all
dt <- Reduce(function(a, b) merge(a, b, by = "MPA", all.x = TRUE),
             list(kelp, purp, lob, sheep, red, meta, size, mhw_d))
# environmental + design predictors are available for most MPAs (the main model);
# trophic effect-size predictors are only co-measured at a few MPAs, so they enter
# the univariate screen only (multi-predictor + CV models would collapse to k~9).
env_design <- c("mhw_during", "cs_mean", "lat", "log_size", "smr", "island", "est_year")
trophic    <- c("lob_eff", "sheep_eff", "purp_eff", "red_eff")
predictors <- c(env_design, trophic)

# ---------------------------------------------------------------------------
# 2. Univariate inverse-variance meta-regression screen (per target)
# ---------------------------------------------------------------------------
screen <- function(target, target_se, preds, data) {
  rows <- lapply(preds, function(p) {
    d <- data[is.finite(get(target)) & is.finite(get(target_se)) & is.finite(get(p))]
    if (nrow(d) < 7 || length(unique(d[[p]])) < 2) return(NULL)
    d$z <- as.numeric(scale(d[[p]])); d$.yi <- d[[target]]; d$.sei <- d[[target_se]]
    m <- safe(rma(yi = .yi, sei = .sei, mods = ~ z, data = d, method = "REML", test = "knha"))
    if (is.null(m)) return(NULL)
    data.frame(target = target, predictor = p, k = nrow(d),
               slope = as.numeric(m$b[2]), ci_lo = m$ci.lb[2], ci_hi = m$ci.ub[2], p = m$pval[2])
  })
  do.call(rbind, rows)
}
scr <- rbind(screen("kelp_eff", "kelp_se", predictors, dt),
             screen("purp_eff", "purp_se", setdiff(predictors, "purp_eff"), dt))
scr$FDR <- p.adjust(scr$p, "BH")
for (cl in c("slope","ci_lo","ci_hi","p","FDR")) scr[[cl]] <- signif(as.numeric(scr[[cl]]), 3)
write.csv(scr, here::here("tables", "table_s_effectiveness_predictors.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 3-4. Multi-predictor meta-regression, random forest, and LEAVE-ONE-OUT CV R^2
# ---------------------------------------------------------------------------
analyze_target <- function(target, target_se, preds, data, label) {
  out <- list(target = label, n_pred = length(preds))
  # rank predictors by the univariate screen; keep the usable ones (>=3 unique, n>=7)
  sc <- screen(target, target_se, preds, data)
  if (is.null(sc) || nrow(sc) == 0) { out$note <- "no usable predictors"; out$k <- 0; return(out) }
  usable <- sc[order(sc$p), "predictor"]
  top <- head(usable, 3)                          # parsimonious set, given small k
  d <- data[is.finite(get(target)) & is.finite(get(target_se))]
  d <- d[, c(target, target_se, usable), with = FALSE]; d <- d[complete.cases(d)]
  out$k <- nrow(d); out$usable <- paste(usable, collapse = ", ")
  if (nrow(d) < 10) { out$note <- "too few complete MPAs for multivariate"; }

  # (3a) multi-predictor inverse-variance meta-regression (usable set) + pseudo-R^2
  if (nrow(d) >= 10) {
    zd <- copy(d); for (p in usable) zd[[p]] <- as.numeric(scale(zd[[p]]))
    zd$.yi <- zd[[target]]; zd$.sei <- zd[[target_se]]
    f <- as.formula(paste("~", paste(usable, collapse = " + ")))
    m_full <- safe(rma(yi = .yi, sei = .sei, mods = f, data = zd, method = "REML"))
    m_null <- safe(rma(yi = .yi, sei = .sei, data = zd, method = "REML"))
    out$metareg_R2 <- if (!is.null(m_full) && !is.null(m_null) && m_null$tau2 > 0)
      round(max(0, (m_null$tau2 - m_full$tau2) / m_null$tau2), 3) else NA

    # (3b) random forest (unweighted) variable importance + OOB R^2
    if (have_rf) {
      rfd <- as.data.frame(d[, c(target, usable), with = FALSE]); names(rfd)[1] <- "y"
      set.seed(1)
      rf <- safe(randomForest::randomForest(y ~ ., data = rfd, importance = TRUE, ntree = 2000))
      if (!is.null(rf)) {
        out$rf_oob_R2 <- round(rf$rsq[length(rf$rsq)], 3)
        imp <- randomForest::importance(rf, type = 1)
        out$rf_importance <- imp[order(-imp[, 1]), , drop = FALSE]
      }
    }
  }

  # (4) leave-one-out CV R^2 of a PARSIMONIOUS model (top-3 screened predictors)
  cvd <- as.data.frame(d[, c(target, target_se, top), with = FALSE]); cvd <- cvd[complete.cases(cvd), ]
  if (nrow(cvd) >= 10) {
    y <- cvd[[target]]; wt <- 1 / cvd[[target_se]]^2
    fcv <- as.formula(paste(target, "~", paste(top, collapse = " + ")))
    pred <- vapply(seq_len(nrow(cvd)), function(i) {
      tr <- cvd[-i, , drop = FALSE]; tr$.w <- wt[-i]
      mi <- safe(lm(fcv, data = tr, weights = .w))
      if (is.null(mi)) NA_real_ else as.numeric(predict(mi, newdata = cvd[i, , drop = FALSE]))
    }, numeric(1))
    ok <- is.finite(pred)
    ss_res <- sum((y[ok] - pred[ok])^2); ss_tot <- sum((y[ok] - mean(y[ok]))^2)
    out$loo_predictors <- paste(top, collapse = " + ")
    out$loo_cv_R2 <- round(1 - ss_res / ss_tot, 3)
    out$loo_pred <- data.frame(MPA_obs = y[ok], MPA_pred = pred[ok])
  }
  out
}
res_kelp <- analyze_target("kelp_eff", "kelp_se", env_design, dt, "Giant kelp effectiveness (biomass lnRR)")
res_purp <- analyze_target("purp_eff", "purp_se", env_design, dt, "Purple urchin suppression (density lnRR)")

# model summary table
mk_row <- function(r) data.frame(
  target = r$target, k = r$k, n_predictors = r$n_pred,
  metareg_pseudoR2 = ifelse(is.null(r$metareg_R2), NA, r$metareg_R2),
  rf_oob_R2 = ifelse(is.null(r$rf_oob_R2), NA, r$rf_oob_R2),
  loo_cv_R2 = ifelse(is.null(r$loo_cv_R2), NA, r$loo_cv_R2),
  loo_predictors = ifelse(is.null(r$loo_predictors), NA, r$loo_predictors))
model_tab <- rbind(mk_row(res_kelp), mk_row(res_purp))
write.csv(model_tab, here::here("tables", "table_s_effectiveness_models.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 5. Figure: RF variable importance (kelp) + observed-vs-LOO-predicted (kelp)
# ---------------------------------------------------------------------------
plots <- list()
if (!is.null(res_kelp$rf_importance)) {
  imp <- res_kelp$rf_importance
  impdf <- data.frame(predictor = rownames(imp), importance = imp[, 1])
  impdf$predictor <- factor(impdf$predictor, levels = impdf$predictor[order(impdf$importance)])
  plots$imp <- ggplot(impdf, aes(importance, predictor)) +
    geom_col(fill = "#0072B2", width = 0.7) +
    labs(x = "Variable importance (% increase in MSE)", y = NULL,
         title = "a  What predicts kelp effectiveness? (random forest)") +
    theme_mpa(base_size = 9) + theme(plot.title = element_text(size = 9, face = "bold"))
}
if (!is.null(res_kelp$loo_pred) && nrow(res_kelp$loo_pred) > 0) {
  lp <- res_kelp$loo_pred; rng <- range(c(lp$MPA_obs, lp$MPA_pred))
  plots$cv <- ggplot(lp, aes(MPA_pred, MPA_obs)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey50") +
    geom_point(alpha = 0.6, size = 1.4, color = "#D55E00") +
    coord_equal(xlim = rng, ylim = rng) +
    labs(x = "Leave-one-out predicted kelp lnRR", y = "Observed kelp lnRR",
         title = paste0("b  Can we predict it? LOO-CV R2 = ",
                        ifelse(is.null(res_kelp$loo_cv_R2), "NA", res_kelp$loo_cv_R2))) +
    theme_mpa(base_size = 9) + theme(plot.title = element_text(size = 9, face = "bold"))
}
if (length(plots) == 2) {
  fig <- patchwork::wrap_plots(plots$imp, plots$cv, nrow = 1)
  ggsave(here::here("plots", "fig_s_effectiveness_predictors.pdf"), fig, width = 180, height = 85, units = "mm", device = cairo_pdf)
  ggsave(here::here("plots", "fig_s_effectiveness_predictors.png"), fig, width = 180, height = 85, units = "mm", dpi = 600)
}

cat("\n=== MPA effectiveness predictability ===\n"); print(model_tab)
cat("\nKelp RF importance (top):\n"); if (!is.null(res_kelp$rf_importance)) print(round(head(res_kelp$rf_importance, 6), 2))
cat("\nUnivariate predictors surviving FDR:", sum(scr$FDR < 0.05, na.rm = TRUE), "of", nrow(scr), "\n")
cat("  Script 18 complete. Tables + figure written.\n")
