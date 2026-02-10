# =============================================================================
# 08_effect_sizes.R
# =============================================================================
#
# PURPOSE:
#   Run pBACIPS analyses and calculate effect sizes for all combinations of
#   data source, taxa, response type, and MPA.
#
# WHAT THIS SCRIPT DOES:
#   1. For each data subset (taxa x source x response x MPA):
#      a. Tests for linear trends in the "Before" period (violates BACI assumptions)
#      b. Runs pBACIPS model selection (step, linear, asymptotic, sigmoid)
#      c. Calculates effect size (mean difference on log-ratio scale)
#      d. Records the best-fitting model and effect size with SE
#   2. Handles three types of analyses:
#      - pBACIPS: Full before-after comparison with model selection
#      - BACI: Simple step-change model (for sites with linear before period)
#      - CI: Control-Impact only (for sites without before period data)
#   3. Compiles all effect sizes into SumStats dataframe
#   4. Filters to final publication-ready subset (SumStats.Final)
#
# EFFECT SIZE INTERPRETATION:
#   - Effect size = change in log response ratio after MPA implementation
#   - Positive values: MPA site increased relative to reference
#   - Negative values: MPA site decreased relative to reference
#   - SE allows calculation of confidence intervals and meta-analysis weights
#
# STATISTICAL NOTES:
#   - Linear before period: If there's a significant trend before MPA,
#     we cannot attribute post-MPA changes to protection (confounded)
#   - pBACIPS vs BACI: pBACIPS allows progressive change; BACI assumes step
#   - CI analysis: Used when no before-period data exists
#
# INPUTS:
#   - All.RR.sub.trans: Combined response ratio data
#   - Site, Sites2: Site metadata with MPA start years
#   - ProgressiveChangeBACIPS function from 02_pBACIPS_function.R
#
# OUTPUTS:
#   - SumStats: All effect size estimates with metadata
#   - SumStats.Final: Filtered subset for publication
#     (excludes MPAs with data quality issues)
#
# DEPENDENCIES:
#   Requires 00-07 scripts to be sourced first
#
# AUTHORS: Emily Donham & Adrian Stier
# PROJECT: CA MPA Kelp Forest pBACIPS Analysis
# =============================================================================

####################################################################################################
## Source pBACIPS function and set up ##############################################################
####################################################################################################

source(here::here("code", "R", "02_pBACIPS_function.R"))

# Color palette for diagnostic plots
# Note: `cols` is already defined in 00b_color_palette.R as col_taxa_short

# Rename Site columns for MPA_implement reference
MPA_implement <- Site
names(MPA_implement)[names(MPA_implement) == "Site"] <- "Site_ID"
names(MPA_implement)[names(MPA_implement) == "CA_MPA_Name_Short"] <- "Site"

####################################################################################################
## DHARMa Model Diagnostics Setup ##################################################################
####################################################################################################

# Check if DHARMa is available
DHARMA_AVAILABLE <- requireNamespace("DHARMa", quietly = TRUE)
if (DHARMA_AVAILABLE) {
  library(DHARMa)
}

# Initialize diagnostic tracking dataframe
ModelDiagnostics <- data.frame(
 Taxa = character(),
  MPA = character(),
  Source = character(),
  Model_Type = character(),
  N_Obs = integer(),
  Uniformity_p = numeric(),
  Dispersion_p = numeric(),
  Outlier_p = numeric(),
  Shapiro_p = numeric(),
  Hetero_Cor = numeric(),
  R_Squared = numeric(),
  Pass_All = logical(),
  Notes = character(),
  stringsAsFactors = FALSE
)

#' Run DHARMa diagnostics on a linear model
#'
#' @param model Fitted lm object
#' @param taxa Character taxa name
#' @param mpa Character MPA name
#' @param source Character data source
#' @param model_type Character model type (Step, Linear, etc.)
#' @return One-row dataframe with diagnostic results
run_dharma_diagnostics <- function(model, taxa, mpa, source, model_type) {
  if (!DHARMA_AVAILABLE || is.null(model)) {
    return(data.frame(
      Taxa = taxa, MPA = mpa, Source = source, Model_Type = model_type,
      N_Obs = NA, Uniformity_p = NA, Dispersion_p = NA, Outlier_p = NA,
      Shapiro_p = NA, Hetero_Cor = NA, R_Squared = NA, Pass_All = NA,
      Notes = ifelse(DHARMA_AVAILABLE, "NULL model", "DHARMa not installed"),
      stringsAsFactors = FALSE
    ))
  }

  tryCatch({
    # Get basic model info
    n_obs <- length(residuals(model))
    r_sq <- summary(model)$r.squared

    # Shapiro-Wilk test on residuals
    resid <- residuals(model)
    shapiro_p <- if (n_obs >= 3 && n_obs <= 5000) {
      shapiro.test(resid)$p.value
    } else {
      NA
    }

    # Heteroscedasticity check
    hetero_cor <- cor(abs(resid), fitted(model), use = "complete.obs")

    # DHARMa simulation-based tests
    sim_resid <- simulateResiduals(model, n = 250, plot = FALSE)
    uniformity <- testUniformity(sim_resid, plot = FALSE)
    dispersion <- testDispersion(sim_resid, plot = FALSE)
    outliers <- testOutliers(sim_resid, plot = FALSE)

    # Determine pass/fail
    pass_all <- (uniformity$p.value > 0.05) &&
                (dispersion$p.value > 0.05) &&
                (outliers$p.value > 0.05) &&
                (is.na(shapiro_p) || shapiro_p > 0.05) &&
                (abs(hetero_cor) < 0.5)

    notes <- c()
    if (uniformity$p.value <= 0.05) notes <- c(notes, "non-uniform")
    if (dispersion$p.value <= 0.05) notes <- c(notes, "overdispersion")
    if (outliers$p.value <= 0.05) notes <- c(notes, "outliers")
    if (!is.na(shapiro_p) && shapiro_p <= 0.05) notes <- c(notes, "non-normal")
    if (abs(hetero_cor) >= 0.5) notes <- c(notes, "heteroscedastic")

    data.frame(
      Taxa = taxa, MPA = mpa, Source = source, Model_Type = model_type,
      N_Obs = n_obs,
      Uniformity_p = round(uniformity$p.value, 4),
      Dispersion_p = round(dispersion$p.value, 4),
      Outlier_p = round(outliers$p.value, 4),
      Shapiro_p = round(shapiro_p, 4),
      Hetero_Cor = round(hetero_cor, 3),
      R_Squared = round(r_sq, 3),
      Pass_All = pass_all,
      Notes = if (length(notes) > 0) paste(notes, collapse = "; ") else "OK",
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      Taxa = taxa, MPA = mpa, Source = source, Model_Type = model_type,
      N_Obs = NA, Uniformity_p = NA, Dispersion_p = NA, Outlier_p = NA,
      Shapiro_p = NA, Hetero_Cor = NA, R_Squared = NA, Pass_All = NA,
      Notes = paste("Error:", e$message),
      stringsAsFactors = FALSE
    )
  })
}

#' Run diagnostics on NLS (nonlinear) models
#'
#' DHARMa doesn't work directly with NLS, so use custom diagnostics
#'
#' @param model Fitted nls or similar object
#' @param taxa Character taxa name
#' @param mpa Character MPA name
#' @param source Character data source
#' @param model_type Character model type
#' @return One-row dataframe with diagnostic results
run_nls_diagnostics <- function(model, taxa, mpa, source, model_type) {
  if (is.null(model)) {
    return(data.frame(
      Taxa = taxa, MPA = mpa, Source = source, Model_Type = model_type,
      N_Obs = NA, Uniformity_p = NA, Dispersion_p = NA, Outlier_p = NA,
      Shapiro_p = NA, Hetero_Cor = NA, R_Squared = NA, Pass_All = NA,
      Notes = "NULL model",
      stringsAsFactors = FALSE
    ))
  }

  tryCatch({
    resid <- residuals(model)
    fitted_vals <- fitted(model)
    n_obs <- length(resid)

    # Shapiro-Wilk test
    shapiro_p <- if (n_obs >= 3 && n_obs <= 5000) {
      shapiro.test(resid)$p.value
    } else {
      NA
    }

    # Heteroscedasticity
    hetero_cor <- cor(abs(resid), fitted_vals, use = "complete.obs")

    # R-squared approximation
    ss_res <- sum(resid^2)
    ss_tot <- sum((fitted_vals + resid - mean(fitted_vals + resid))^2)
    r_sq <- 1 - ss_res / ss_tot

    # Outlier count
    sigma <- tryCatch(summary(model)$sigma, error = function(e) sd(resid))
    n_outliers <- sum(abs(resid) > 3 * sigma)

    # Pass criteria for NLS
    pass_all <- (is.na(shapiro_p) || shapiro_p > 0.01) &&  # More lenient for NLS
                (abs(hetero_cor) < 0.6) &&
                (n_outliers <= 1)

    notes <- c()
    if (!is.na(shapiro_p) && shapiro_p <= 0.01) notes <- c(notes, "non-normal")
    if (abs(hetero_cor) >= 0.6) notes <- c(notes, "heteroscedastic")
    if (n_outliers > 1) notes <- c(notes, paste0(n_outliers, " outliers"))

    data.frame(
      Taxa = taxa, MPA = mpa, Source = source, Model_Type = model_type,
      N_Obs = n_obs,
      Uniformity_p = NA,  # Not applicable for NLS
      Dispersion_p = NA,
      Outlier_p = NA,
      Shapiro_p = round(shapiro_p, 4),
      Hetero_Cor = round(hetero_cor, 3),
      R_Squared = round(r_sq, 3),
      Pass_All = pass_all,
      Notes = if (length(notes) > 0) paste(notes, collapse = "; ") else "OK",
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      Taxa = taxa, MPA = mpa, Source = source, Model_Type = model_type,
      N_Obs = NA, Uniformity_p = NA, Dispersion_p = NA, Outlier_p = NA,
      Shapiro_p = NA, Hetero_Cor = NA, R_Squared = NA, Pass_All = NA,
      Notes = paste("Error:", e$message),
      stringsAsFactors = FALSE
    )
  })
}

#' Run Durbin-Watson test for residual autocorrelation
#'
#' Tests whether residuals from a linear model show significant temporal
#' autocorrelation. Important for ecological time series spanning 10-20+ years,
#' where correlated residuals can inflate Type I error by underestimating SEs.
#'
#' @param model Fitted lm object
#' @return Named list with DW_stat (Durbin-Watson statistic) and DW_pval (p-value).
#'   DW_stat near 2 indicates no autocorrelation; <2 positive, >2 negative.
#'   Returns list(DW_stat = NA, DW_pval = NA) if test cannot be run.
run_dw_test <- function(model) {
  if (is.null(model) || !inherits(model, "lm")) {
    return(list(DW_stat = NA, DW_pval = NA))
  }
  tryCatch({
    if (requireNamespace("lmtest", quietly = TRUE)) {
      dw <- lmtest::dwtest(model)
      list(DW_stat = round(dw$statistic, 4), DW_pval = round(dw$p.value, 4))
    } else {
      list(DW_stat = NA, DW_pval = NA)
    }
  }, error = function(e) {
    warning("Durbin-Watson test failed: ", e$message)
    list(DW_stat = NA, DW_pval = NA)
  })
}

####################################################################################################
## Initialize SumStats dataframe ###################################################################
####################################################################################################

SumStats <- data.frame(
  Taxa = NA, MPA = NA, Mean = NA, SE = NA, SD = NA, CI = NA,
  Model = NA, Source = NA, Resp = NA, BA = NA, Primary = NA,
  Type = NA, LinearBefore = NA, N = NA, DW_stat = NA, DW_pval = NA
)

####################################################################################################
## Helper functions for repeated analysis patterns #################################################
####################################################################################################

#' Add time.model and time.true columns to a subset of RR data
#'
#' time.model: 0 for Before period, sequential starting at 0 for After
#' time.true: sequential index across all time points per MPA
#'
#' @param dat Dataframe with CA_MPA_Name_Short, BA, year columns
#' @return Modified dataframe with time.model and time.true columns
add_time_columns <- function(dat) {
  dat$time.model <- NA
  dat$time.true <- NA
  dat <- dat[order(as.numeric(dat$year), decreasing = FALSE), ]
  mpas <- unique(dat$CA_MPA_Name_Short)
  for (i in seq_along(mpas)) {
    idx <- which(dat$CA_MPA_Name_Short == mpas[i])
    n_before <- length(which(dat$BA[idx] == "Before"))
    n_after <- length(which(dat$BA[idx] == "After"))
    dat$time.model[idx] <- c(rep(0, n_before), seq(0, n_after - 1))
    dat$time.true[idx] <- seq(1, length(idx))
  }
  dat
}

#' Test for linear trend in the Before period for each MPA
#'
#' METHODOLOGICAL NOTE:
#' The BACI/pBACIPS design assumes that MPA and reference sites have parallel
#' trajectories before intervention. A significant pre-existing trend in the
#' log response ratio violates this assumption and can confound estimates.
#'
#' DECISION RATIONALE:
#' Sites with significant (p < 0.05) linear trends in the before period are
#' excluded from pBACIPS analysis and analyzed separately as CI (Control-Impact)
#' sites, which don't require the parallel trend assumption.
#'
#' LIMITATIONS AND CAVEATS:
#' 1. Low statistical power: With only 3-5 before-period points, power to detect
#'    trends is low. Sites may pass this test simply due to noise.
#' 2. Selection bias: This approach may favor inclusion of noisy sites over sites
#'    with detectable patterns, potentially biasing effect estimates.
#' 3. Arbitrary threshold: The p < 0.05 cutoff is conventional but arbitrary.
#'
#' POTENTIAL SENSITIVITY ANALYSES:
#' - Compare results using different p-value thresholds (0.01, 0.10)
#' - Include all sites and use detrending approaches
#' - Report results both with and without exclusions
#'
#' @param dat Dataframe with CA_MPA_Name_Short, lnDiff, time.true, BA columns
#' @return Dataframe with pInt (intercept p-value), pSlope (slope p-value),
#'         and site columns. Sites with pSlope < 0.05 have significant before trends.
test_linear_before <- function(dat) {
  mpas <- unique(dat$CA_MPA_Name_Short)
  LinearBefore <- data.frame(pInt = NA, pSlope = NA)
  valid_mpas <- c()

  for (j in mpas) {
    temp <- subset(dat, CA_MPA_Name_Short == j & BA == "Before",
                   select = c(CA_MPA_Name_Short, mpa, reference, lnDiff, time.true, time.model, BA))

    # Skip if insufficient before data
    if (nrow(temp) < 3) {
      warning("test_linear_before: skipping ", j, " (insufficient Before data: n=", nrow(temp), ")")
      next
    }

    result <- tryCatch({
      lmBefore <- lm(lnDiff ~ time.true, data = temp)
      coef_summary <- summary(lmBefore)$coefficients
      if (nrow(coef_summary) < 2) {
        warning("test_linear_before: no slope for ", j)
        return(NULL)
      }
      coef_summary[, 4]
    }, error = function(e) {
      warning("test_linear_before failed for ", j, ": ", e$message)
      NULL
    })

    if (!is.null(result)) {
      LinearBefore <- rbind(LinearBefore, result)
      valid_mpas <- c(valid_mpas, j)
    }
  }

  LinearBefore <- na.omit(LinearBefore)
  if (nrow(LinearBefore) > 0) {
    LinearBefore$site <- valid_mpas
  }
  LinearBefore
}

#' Run ProgressiveChangeBACIPS across MPAs and collect model weights/p-values
#'
#' @param dat Dataframe with lnDiff, lnmpa, lnreference, time.true, time.model
#' @return List with likelihood.poly, pvals.poly, Norm.poly dataframes
run_pbacips_loop <- function(dat) {
  mpas <- unique(dat$CA_MPA_Name_Short)
  likelihood.poly <- data.frame(matrix(ncol = 4, nrow = 0))
  pvals.poly <- data.frame(matrix(ncol = 9, nrow = 0))
  Norm.poly <- data.frame(matrix(ncol = 9, nrow = 0))
  valid_mpas <- c()

  dat$lnmpa <- log(dat$mpa)
  dat$lnreference <- log(dat$reference)

  for (j in mpas) {
    temp <- subset(dat, CA_MPA_Name_Short == j,
                   select = c(CA_MPA_Name_Short, lnmpa, lnreference, lnDiff, time.true, time.model))

    # Skip if insufficient data
    if (nrow(temp) < 5) {
      warning("run_pbacips_loop: skipping ", j, " (insufficient data: n=", nrow(temp), ")")
      next
    }

    result <- tryCatch({
      ProgressiveChangeBACIPS(
        control = temp$lnreference,
        impact = temp$lnmpa,
        time.true = temp$time.true,
        time.model = temp$time.model
      )
    }, error = function(e) {
      warning("run_pbacips_loop: pBACIPS failed for ", j, ": ", e$message)
      NULL
    })

    if (is.null(result)) next

    # Safely extract p-values with error handling
    pvals <- tryCatch({
      ps <- if (!is.null(result$step)) summary(result$step)[[1]][["Pr(>F)"]][1] else NA
      pl <- if (!is.null(result$linear)) summary(result$linear)$coefficients[, 4] else c(NA, NA)
      psi <- if (!is.null(result$sigmoid)) summary(result$sigmoid)$coefficients[, 4] else rep(NA, 4)
      pa <- if (!is.null(result$asymptotic)) summary(result$asymptotic)$coefficients[, 4] else rep(NA, 3)
      c(ps, pl, psi, pa)
    }, error = function(e) {
      warning("run_pbacips_loop: p-value extraction failed for ", j, ": ", e$message)
      rep(NA, 9)
    })

    likelihood.poly <- rbind(likelihood.poly, result$weights)
    pvals.poly <- rbind(pvals.poly, pvals)

    # Normality test
    norm_p <- tryCatch({
      shapiro.test(temp$lnDiff)$p.value
    }, error = function(e) NA)

    Norm.poly <- rbind(Norm.poly, norm_p)
    valid_mpas <- c(valid_mpas, j)
  }

  colnames(likelihood.poly) <- c("Step", "Linear", "Asymptotic", "Sigmoid")
  colnames(pvals.poly) <- c("Step", "Linear", "Asymptotic M", "Asymptotic B", "Asymptotic L",
                             "Sigmoid M", "Sigmoid B", "Sigmoid L", "Sigmoid K")
  if (nrow(likelihood.poly) > 0) {
    likelihood.poly <- cbind(likelihood.poly, mpas = valid_mpas)
  }

  list(likelihood = likelihood.poly, pvals = pvals.poly, normality = Norm.poly)
}

#' Run CI analysis: calculate both linear and mean effect sizes
#'
#' For each MPA in the data, fits lnDiff ~ time, calculates effect size at time 0 vs t=11,
#' and also calculates the simple mean effect size. Returns two rows per MPA.
#'
#' The effect size represents the predicted change in lnRR from time 0 (MPA establishment)
#' to 11 years post-implementation. We use a STANDARDIZED time point (t=11) rather than
#' the maximum observed time to ensure effect sizes are comparable across MPAs with
#' different establishment dates. The value of 11 years corresponds to the age of
#' our youngest MPA in 2023 (MLPA South Coast MPAs implemented in 2012). This controls
#' for differences in effect size that would arise from longer protection duration.
#'
#' @param dat Dataframe subset with lnDiff, time columns
#' @param taxa_name Character, species name for SumStats
#' @param source_name Character, data source for SumStats
#' @param resp_name Character, response type ("Den" or "Bio")
#' @param time_var Character, name of the time variable to use (default "time")
#' @param time_after Numeric, time point for "after" prediction (default = EFFECT_SIZE_TIME_YEARS = 11)
#' @param run_diagnostics Logical, whether to run DHARMa diagnostics (default TRUE)
#' @return Dataframe with two rows per MPA (linear and mean effect sizes)
run_ci_analysis <- function(dat, taxa_name, source_name, resp_name, time_var = "time",
                             time_after = NULL, run_diagnostics = TRUE) {
  # Validate input
  if (is.null(dat) || nrow(dat) == 0) {
    warning("run_ci_analysis: empty dataset for ", taxa_name)
    return(NULL)
  }

  mpas <- unique(dat$CA_MPA_Name_Short)
  rows <- list()

  for (i in seq_along(mpas)) {
    idx <- which(dat$CA_MPA_Name_Short == mpas[i])
    sub_dat <- dat[idx, ]

    # Skip if insufficient data
    if (nrow(sub_dat) < 3) {
      warning("run_ci_analysis: skipping ", mpas[i], " (n=", nrow(sub_dat), ")")
      next
    }

    # Skip if no variance in time
    if (length(unique(sub_dat[[time_var]])) < 2) {
      warning("run_ci_analysis: skipping ", mpas[i], " (no time variance)")
      next
    }

    # Get sample size before tryCatch
    n_obs <- nrow(sub_dat)

    # Wrap in tryCatch for robust error handling
    result <- tryCatch({
      # Fit linear model
      formula_str <- as.formula(paste("lnDiff ~", time_var))
      Lm.Ab <- lm(formula_str, data = sub_dat)

      # Durbin-Watson test for residual autocorrelation
      dw_result <- run_dw_test(Lm.Ab)

      # Run DHARMa diagnostics and add to global tracking
      if (run_diagnostics && exists("ModelDiagnostics", envir = .GlobalEnv)) {
        diag_result <- run_dharma_diagnostics(Lm.Ab, taxa_name, mpas[i], source_name, "CI_Linear")
        ModelDiagnostics <<- rbind(ModelDiagnostics, diag_result)
      }

      # Check model has slope coefficient
      coef_summary <- summary(Lm.Ab)$coefficients
      if (nrow(coef_summary) < 2) {
        warning("run_ci_analysis: no slope for ", mpas[i])
        return(NULL)
      }

      p <- coef_summary[2, 4]

      # Calculate linear effect size using contrast (properly handles covariance)
      # This uses emmeans::pairs() which correctly calculates Var(A-B) = Var(A) + Var(B) - 2*Cov(A,B)
      # Use standardized t=11 years to ensure effect sizes are comparable across MPAs
      # RATIONALE: Using a fixed time point controls for MPA age effects - older MPAs would
      # naturally show larger effect sizes simply due to longer protection duration.
      # t=11 corresponds to the youngest MPA age in 2023 (MLPA South Coast, implemented 2012).
      actual_time_after <- if (is.null(time_after)) EFFECT_SIZE_TIME_YEARS else time_after
      es <- calculate_effect_size_from_contrast(Lm.Ab, time_var, time_before = 0, time_after = actual_time_after)
      mean_val <- es$mean

      # --- EDGE CASE: Validate effect size is reasonable ---
      if (!is.finite(mean_val) || !is.finite(es$SE) || es$SE <= 0) {
        warning("run_ci_analysis: invalid effect size for ", mpas[i], " (mean=", mean_val, ", SE=", es$SE, ")")
        return(NULL)
      }

      # Calculate simple mean CI
      # EDGE CASE: summarySE requires n >= 2 for SD calculation
      if (n_obs < 2) {
        warning("run_ci_analysis: single observation for ", mpas[i], ", using linear model SE only")
        CP.mean <- data.frame(lnDiff = mean(sub_dat$lnDiff), sd = NA, se = NA, ci = NA)
      } else {
        CP.mean <- summarySE(sub_dat, measurevar = "lnDiff")
      }

      # Build result rows
      # Note: calculate_effect_size_from_contrast returns SE with proper covariance handling
      # SD is calculated from SE and df for backward compatibility
      es_sd <- es$SE * sqrt(es$df + 1)  # Approximate SD from SE and df
      linear_row <- data.frame(
        Taxa = taxa_name, MPA = mpas[i], Mean = mean_val, SE = es$SE, SD = es_sd, CI = es$CI,
        Model = "Linear", Source = source_name, Resp = resp_name, BA = "N",
        Primary = if (p <= 0.05) "Y" else "N",
        Type = "CI", LinearBefore = "NA", N = n_obs,
        DW_stat = dw_result$DW_stat, DW_pval = dw_result$DW_pval,
        stringsAsFactors = FALSE
      )

      mean_row <- data.frame(
        Taxa = taxa_name, MPA = mpas[i], Mean = CP.mean$lnDiff, SE = CP.mean$se,
        SD = CP.mean$sd, CI = CP.mean$ci, Model = "Mean", Source = source_name, Resp = resp_name,
        BA = "N", Primary = if (p <= 0.05) "N" else "Y",
        Type = "CI", LinearBefore = "NA", N = n_obs,
        DW_stat = dw_result$DW_stat, DW_pval = dw_result$DW_pval,
        stringsAsFactors = FALSE
      )

      rbind(linear_row, mean_row)
    }, error = function(e) {
      warning("run_ci_analysis failed for ", mpas[i], ": ", e$message)
      NULL
    })

    if (!is.null(result)) {
      rows[[length(rows) + 1]] <- result
    }
  }

  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

# =============================================================================
# NOTE ON SE ESTIMATION METHODOLOGY:
# =============================================================================
# - Step and linear models: SE from emmeans::pairs() contrast (uses t-distribution,
#   accounts for model covariance structure). This is the preferred approach as it
#   correctly handles Var(A-B) = Var(A) + Var(B) - 2*Cov(A,B).
# - NLS models (asymptotic/sigmoid): SE back-calculated from investr::predFit()
#   confidence intervals as CI_width / (2 * 1.96). This assumes normality (z-distribution)
#   rather than t-distribution, making it slightly anti-conservative for small samples.
#   The independence assumption between t=0 and t=11 predictions makes the SE
#   slightly conservative (overestimated).
# - Net effect: These biases approximately cancel, but readers should note the
#   methodological difference. All effect sizes are downstream weighted by 1/SE^2
#   in the meta-analysis (09_meta_analysis.R), so any SE inflation for NLS models
#   results in slightly lower weight, which is conservative.
# =============================================================================

#' Add a single step-model (BACI) effect size row
#'
#' @param dat Dataframe for a single MPA with lnDiff and BA columns
#' @param taxa_name Character species name
#' @param mpa_name Character MPA name
#' @param source_name Character data source
#' @param resp_name Character response type
#' @param run_diagnostics Logical, whether to run DHARMa diagnostics (default TRUE)
#' @return One-row dataframe with the same columns as SumStats
add_step_effect_size <- function(dat, taxa_name, mpa_name, source_name, resp_name,
                                  run_diagnostics = TRUE) {
  # Validate input data
  if (is.null(dat) || nrow(dat) == 0) {
    warning("add_step_effect_size: empty dataset for ", mpa_name)
    return(NULL)
  }

  # Check BA levels - need both Before and After with at least 2 obs each
  ba_counts <- table(dat$BA)
  n_before <- ifelse("Before" %in% names(ba_counts), ba_counts["Before"], 0)
  n_after <- ifelse("After" %in% names(ba_counts), ba_counts["After"], 0)

  if (length(ba_counts) < 2 || n_before < 2 || n_after < 2) {
    warning("add_step_effect_size: insufficient BA levels for ", mpa_name,
            " (Before=", n_before, ", After=", n_after, ")")
    return(NULL)
  }

  # Get sample size
  n_obs <- nrow(dat)

  # Wrap model fitting in tryCatch
  tryCatch({
    mod1 <- lm(data = dat, lnDiff ~ BA)

    # Durbin-Watson test for residual autocorrelation
    dw_result <- run_dw_test(mod1)

    # Run DHARMa diagnostics and add to global tracking
    if (run_diagnostics && exists("ModelDiagnostics", envir = .GlobalEnv)) {
      diag_result <- run_dharma_diagnostics(mod1, taxa_name, mpa_name, source_name, "Step")
      ModelDiagnostics <<- rbind(ModelDiagnostics, diag_result)
    }

    # Use pairs() for proper covariance-aware contrast estimation
    # This correctly calculates Var(After - Before) accounting for shared model variance
    # Note: pairs() works on emmeans objects (method dispatch), but can't be called as emmeans::pairs()
    em <- emmeans::emmeans(mod1, ~ BA)
    contrast_result <- pairs(em, reverse = TRUE)  # After - Before
    contrast_summary <- summary(contrast_result)

    # Extract effect size with proper SE from contrast
    mean_es <- contrast_summary$estimate[1]
    pSE <- contrast_summary$SE[1]
    df_es <- contrast_summary$df[1]
    pCI <- pSE * qt(0.975, df_es)  # Use t-distribution CI

    # Approximate SD from SE and df for backward compatibility
    # SD ≈ SE * sqrt(n), where n ≈ df + 1 for simple models
    pSD <- pSE * sqrt(df_es + 1)

    data.frame(
      Taxa = taxa_name, MPA = mpa_name, Mean = mean_es, SE = pSE, SD = pSD, CI = pCI,
      Model = "Step", Source = source_name, Resp = resp_name, BA = "Y", Primary = "Y",
      Type = "BACI", LinearBefore = "N", N = n_obs,
      DW_stat = dw_result$DW_stat, DW_pval = dw_result$DW_pval,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    warning("add_step_effect_size failed for ", mpa_name, ": ", e$message)
    NULL
  })
}

#' Add a single linear-model (pBACIPS) effect size row
#'
#' Calculates effect size as the predicted change from t=0 (before MPA) to t=11 years
#' post-implementation. We use a STANDARDIZED time point (t=11) rather than the maximum
#' observed time to ensure effect sizes are comparable across MPAs with different
#' establishment dates. This controls for the confounding effect of protection duration.
#'
#' @param dat Dataframe for a single MPA with lnDiff and time.model columns
#' @param taxa_name Character species name
#' @param mpa_name Character MPA name
#' @param source_name Character data source
#' @param resp_name Character response type
#' @param time_var Character, time variable name (default "time.model")
#' @param time_after Numeric, the time point for "after" prediction (default = EFFECT_SIZE_TIME_YEARS = 11)
#' @param run_diagnostics Logical, whether to run DHARMa diagnostics (default TRUE)
#' @return One-row dataframe with the same columns as SumStats
add_linear_effect_size <- function(dat, taxa_name, mpa_name, source_name, resp_name,
                                    time_var = "time.model", time_after = NULL,
                                    run_diagnostics = TRUE) {
  # Validate input data
  if (is.null(dat) || nrow(dat) < 3) {
    warning("add_linear_effect_size: insufficient data for ", mpa_name, " (n=", nrow(dat), ")")
    return(NULL)
  }

  # Check time variable has variance
  time_values <- dat[[time_var]]
  if (length(unique(time_values)) < 2) {
    warning("add_linear_effect_size: no variance in time for ", mpa_name)
    return(NULL)
  }

  # Get sample size
  n_obs <- nrow(dat)

  # Wrap model fitting in tryCatch
  tryCatch({
    # Use standardized t=11 years if not specified
    # RATIONALE: Using a fixed time point ensures effect sizes are comparable across MPAs
    # with different establishment dates. t=11 corresponds to the youngest MPA age in 2023
    # (MLPA South Coast MPAs, implemented 2012). This controls for the confounding effect
    # of protection duration - older MPAs would naturally show larger effect sizes.
    if (is.null(time_after)) {
      time_after <- EFFECT_SIZE_TIME_YEARS
    }
    formula_str <- as.formula(paste("lnDiff ~", time_var))
    mod1 <- lm(formula_str, data = dat)

    # Durbin-Watson test for residual autocorrelation
    dw_result <- run_dw_test(mod1)

    # Run DHARMa diagnostics and add to global tracking
    if (run_diagnostics && exists("ModelDiagnostics", envir = .GlobalEnv)) {
      diag_result <- run_dharma_diagnostics(mod1, taxa_name, mpa_name, source_name, "Linear")
      ModelDiagnostics <<- rbind(ModelDiagnostics, diag_result)
    }

    # Use covariance-aware effect size calculation
    # This properly accounts for correlation between predictions at different time points
    es <- calculate_effect_size_from_contrast(mod1, time_var, time_before = 0, time_after = time_after)

    # Approximate SD from SE and df for backward compatibility
    es_sd <- es$SE * sqrt(es$df + 1)

    data.frame(
      Taxa = taxa_name, MPA = mpa_name, Mean = es$mean, SE = es$SE, SD = es_sd, CI = es$CI,
      Model = "Linear", Source = source_name, Resp = resp_name, BA = "Y", Primary = "Y",
      Type = "pBACIPS", LinearBefore = "N", N = n_obs,
      DW_stat = dw_result$DW_stat, DW_pval = dw_result$DW_pval,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    warning("add_linear_effect_size failed for ", mpa_name, ": ", e$message)
    NULL
  })
}

#' Add a mean-only effect size row (for cases without enough data for pBACIPS)
#'
#' @param dat Dataframe for a single MPA
#' @param taxa_name Character species name
#' @param mpa_name Character MPA name
#' @param source_name Character data source
#' @param resp_name Character response type
#' @return One-row dataframe with the same columns as SumStats
add_mean_effect_size <- function(dat, taxa_name, mpa_name, source_name, resp_name) {
  # Validate input data
  if (is.null(dat) || nrow(dat) < 2) {
    warning("add_mean_effect_size: insufficient data for ", mpa_name, " (n=", nrow(dat), ")")
    return(NULL)
  }

  # Check for valid lnDiff values
  if (all(is.na(dat$lnDiff))) {
    warning("add_mean_effect_size: all lnDiff values are NA for ", mpa_name)
    return(NULL)
  }

  # Get sample size
  n_obs <- nrow(dat)

  # Wrap in tryCatch
  tryCatch({
    CP.mean <- summarySE(dat, measurevar = "lnDiff")
    data.frame(
      Taxa = taxa_name, MPA = mpa_name, Mean = CP.mean$lnDiff, SE = CP.mean$se,
      SD = CP.mean$sd, CI = CP.mean$ci, Model = "Mean", Source = source_name, Resp = resp_name,
      BA = "N", Primary = "Y", Type = "CI", LinearBefore = "NA", N = n_obs,
      DW_stat = NA, DW_pval = NA,  # No linear model fit for mean-only effect sizes
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    warning("add_mean_effect_size failed for ", mpa_name, ": ", e$message)
    NULL
  })
}

####################################################################################################
## LTER ANALYSES ###################################################################################
####################################################################################################

#----- S. purpuratus Density (LTER) -----#
LTER.purps <- subset(All.RR.sub.trans, source == "LTER" &
                       y == "Strongylocentrotus purpuratus" & resp == "Den")
LTER.purps <- add_time_columns(LTER.purps)
LinearBefore.purps <- test_linear_before(LTER.purps)
# Note: Linear in before period -- excluded from pBACIPS (conservative: exclude den & bio)

#----- M. franciscanus Density (LTER) -----#
LTER.reds <- subset(All.RR.sub.trans, source == "LTER" &
                      y == "Mesocentrotus franciscanus" & resp == "Den")
LTER.reds <- add_time_columns(LTER.reds)
LinearBefore.reds <- test_linear_before(LTER.reds)
pbacips.reds.den <- run_pbacips_loop(LTER.reds)

# Only Naples (Campus Point is linear in before); step was best fit model
LTER.redsNap <- subset(LTER.reds, CA_MPA_Name_Short == "Naples SMCA")
SumStats <- rbind(SumStats, add_step_effect_size(LTER.redsNap, "M. franciscanus", "Naples SMCA", "LTER", "Den"))

#----- M. franciscanus Biomass (LTER) -----#
LTER.reds.bio <- subset(All.RR.sub.trans, source == "LTER" &
                           y == "Mesocentrotus franciscanus" & resp == "Bio")
LTER.reds.bio <- add_time_columns(LTER.reds.bio)

# Naples only since Campus Point was linear in density in the before
LTER.nap.red <- dplyr::filter(LTER.reds.bio, CA_MPA_Name_Short == "Naples SMCA" & time > 0)
SumStats <- rbind(SumStats, add_mean_effect_size(LTER.nap.red, "M. franciscanus", "Naples SMCA", "LTER", "Bio"))

#----- M. pyrifera Biomass (LTER) -----#
# Note: only biomass since density is the exact same proportional data (bio is a multiplier of density)
LTER.macro <- subset(All.RR.sub.trans, source == "LTER" &
                       y == "Macrocystis pyrifera" & resp == "Bio")
LTER.macro <- add_time_columns(LTER.macro)
LinearBefore.macro <- test_linear_before(LTER.macro)
pbacips.macro <- run_pbacips_loop(LTER.macro)

# Naples
LTER.nap.macro <- dplyr::filter(LTER.macro, CA_MPA_Name_Short == "Naples SMCA")
SumStats <- rbind(SumStats, add_linear_effect_size(LTER.nap.macro, "M. pyrifera", "Naples SMCA", "LTER", "Bio"))

# Campus Point
LTER.cp.macro <- dplyr::filter(LTER.macro, CA_MPA_Name_Short == "Campus Point SMCA")
SumStats <- rbind(SumStats, add_linear_effect_size(LTER.cp.macro, "M. pyrifera", "Campus Point SMCA", "LTER", "Bio"))

#----- S. pulcher Density (LTER) -----#
LTER.SPUL.den <- subset(All.RR.sub.trans, source == "LTER" &
                           y == "Semicossyphus pulcher" & resp == "Den")
LTER.SPUL.den <- add_time_columns(LTER.SPUL.den)
LinearBefore.spul.den <- test_linear_before(LTER.SPUL.den)
pbacips.spul.den <- run_pbacips_loop(LTER.SPUL.den)

# Naples was linear in the before; Campus Point only
LTER.cp.SPUL.den <- dplyr::filter(LTER.SPUL.den, CA_MPA_Name_Short == "Campus Point SMCA")
SumStats <- rbind(SumStats, add_linear_effect_size(LTER.cp.SPUL.den, "S. pulcher", "Campus Point SMCA", "LTER", "Den"))

#----- S. pulcher Biomass (LTER) -----#
LTER.SPUL.bio <- subset(All.RR.sub.trans, source == "LTER" &
                           y == "Semicossyphus pulcher" & resp == "Bio")
LTER.SPUL.bio <- add_time_columns(LTER.SPUL.bio)
LinearBefore.spul.bio <- test_linear_before(LTER.SPUL.bio)
pbacips.spul.bio <- run_pbacips_loop(LTER.SPUL.bio)

# Campus Point (Naples biomass was linear in the before)
LTER.cp.SPUL.bio <- dplyr::filter(LTER.SPUL.bio, CA_MPA_Name_Short == "Campus Point SMCA")
SumStats <- rbind(SumStats, add_linear_effect_size(LTER.cp.SPUL.bio, "S. pulcher", "Campus Point SMCA", "LTER", "Bio"))

#----- P. interruptus Density (LTER) -- CI data only -----#
LTER.PANINT.den <- subset(All.RR.sub.trans, source == "LTER" &
                             y == "Panulirus interruptus" & resp == "Den")

# Campus Point
LTER.cp.LOB <- dplyr::filter(LTER.PANINT.den, CA_MPA_Name_Short == "Campus Point SMCA")
SumStats <- rbind(SumStats, add_mean_effect_size(LTER.cp.LOB, "P. interruptus", "Campus Point SMCA", "LTER", "Den"))

# Naples
LTER.nap.LOB <- dplyr::filter(LTER.PANINT.den, CA_MPA_Name_Short == "Naples SMCA")
SumStats <- rbind(SumStats, add_mean_effect_size(LTER.nap.LOB, "P. interruptus", "Naples SMCA", "LTER", "Den"))

#----- P. interruptus Biomass (LTER) -----#
LTER.PANINT.bio <- subset(All.RR.sub.trans, source == "LTER" &
                             y == "Panulirus interruptus" & resp == "Bio")

# Campus Point (linear)
LTER.cp.LOB.bio <- dplyr::filter(LTER.PANINT.bio, CA_MPA_Name_Short == "Campus Point SMCA")
SumStats <- rbind(SumStats, add_linear_effect_size(LTER.cp.LOB.bio, "P. interruptus", "Campus Point SMCA", "LTER", "Bio",
                       time_var = "time"))

# Naples (mean only)
LTER.nap.LOB.bio <- dplyr::filter(LTER.PANINT.bio, CA_MPA_Name_Short == "Naples SMCA")
SumStats <- rbind(SumStats, add_mean_effect_size(LTER.nap.LOB.bio, "P. interruptus", "Naples SMCA", "LTER", "Bio"))


####################################################################################################
## PISCO ANALYSES ##################################################################################
####################################################################################################

# The PISCO loop pattern: for each taxa x response, loop over MPAs and calculate both
# linear and mean effect sizes using run_ci_analysis()

#----- M. franciscanus Density (PISCO) -----#
Mes.den <- subset(All.RR.sub.trans, y == "Mesocentrotus franciscanus" &
                    source == "PISCO" & resp == "Den")
SumStats <- rbind(SumStats, run_ci_analysis(Mes.den, "M. franciscanus", "PISCO", "Den"))

#----- M. franciscanus Biomass (PISCO) -----#
Mes.bio <- subset(All.RR.sub.trans, y == "Mesocentrotus franciscanus" &
                    source == "PISCO" & resp == "Bio")
Mes.bio <- subset(Mes.bio, CA_MPA_Name_Short != "Matlahuayl SMR")
SumStats <- rbind(SumStats, run_ci_analysis(Mes.bio, "M. franciscanus", "PISCO", "Bio"))

#----- S. purpuratus Density (PISCO) -----#
Str.den <- subset(All.RR.sub.trans, y == "Strongylocentrotus purpuratus" &
                    source == "PISCO" & resp == "Den")
SumStats <- rbind(SumStats, run_ci_analysis(Str.den, "S. purpuratus", "PISCO", "Den"))

#----- S. purpuratus Biomass (PISCO) -----#
Str.bio <- subset(All.RR.sub.trans, y == "Strongylocentrotus purpuratus" &
                    source == "PISCO" & resp == "Bio")
Str.bio <- subset(Str.bio, CA_MPA_Name_Short != "Matlahuayl SMR")
SumStats <- rbind(SumStats, run_ci_analysis(Str.bio, "S. purpuratus", "PISCO", "Bio"))

#----- M. pyrifera Biomass (PISCO) -----#
Mac <- subset(All.RR.sub.trans, y == "Macrocystis pyrifera" &
                source == "PISCO" & resp == "Bio")
SumStats <- rbind(SumStats, run_ci_analysis(Mac, "M. pyrifera", "PISCO", "Bio"))

#----- P. interruptus Density (PISCO) -----#
Pan.den <- subset(All.RR.sub.trans, y == "Panulirus interruptus" &
                    source == "PISCO" & resp == "Den")
SumStats <- rbind(SumStats, run_ci_analysis(Pan.den, "P. interruptus", "PISCO", "Den"))

#----- P. interruptus Biomass (PISCO) -----#
Pan.bio <- subset(All.RR.sub.trans, y == "Panulirus interruptus" &
                    resp == "Bio" & source == "PISCO")
# Exclude MPAs with insufficient lobster biomass data
Pan.bio <- subset(Pan.bio,
  CA_MPA_Name_Short != "Gull Island SMR" &
  CA_MPA_Name_Short != "Harris Point SMR" &
  CA_MPA_Name_Short != "South Point SMR" &
  CA_MPA_Name_Short != "Scorpion SMR" &
  CA_MPA_Name_Short != "Anacapa Island SMR 2003" &
  CA_MPA_Name_Short != "Santa Barbara Island SMR" &
  CA_MPA_Name_Short != "Point Dume SMCA"
)
SumStats <- rbind(SumStats, run_ci_analysis(Pan.bio, "P. interruptus", "PISCO", "Bio"))

#----- S. pulcher Density (PISCO) -----#
Sheep.den <- subset(All.RR.sub.trans, y == "Semicossyphus pulcher" &
                      resp == "Den" & source == "PISCO")
SumStats <- rbind(SumStats, run_ci_analysis(Sheep.den, "S. pulcher", "PISCO", "Den"))

#----- S. pulcher Biomass (PISCO) -----#
Sheep.bio <- subset(All.RR.sub.trans, y == "Semicossyphus pulcher" &
                      resp == "Bio" & source == "PISCO")
SumStats <- rbind(SumStats, run_ci_analysis(Sheep.bio, "S. pulcher", "PISCO", "Bio"))


####################################################################################################
## KFM ANALYSES ####################################################################################
####################################################################################################

# KFM sites are split into:
# BA sites (Scorpion, SBI, Harris Point, Gull Island) -> pBACIPS with before/after data
# CI sites (Anacapa 2003, South Point) -> CI-only linear regression

# KFM BA site lists
KFM_BA_SITES <- c("Scorpion SMR", "Santa Barbara Island SMR",
                   "Harris Point SMR", "Gull Island SMR")
KFM_CI_SITES <- c("Anacapa Island SMR 2003", "South Point SMR")

#----- S. purpuratus Density (KFM) -----#
KFM.purps.den <- subset(All.RR.sub, source == "KFM" &
                           y == "Strongylocentrotus purpuratus" & resp == "Den")
KFM.purps.den <- add_time_columns(KFM.purps.den)

# BA sites: test for linear before
KFM.purps.BA <- subset(KFM.purps.den, CA_MPA_Name_Short %in% KFM_BA_SITES)
KFM.purps.CI <- subset(KFM.purps.den, CA_MPA_Name_Short %in% KFM_CI_SITES)
LinearBefore.kfm.purps <- test_linear_before(KFM.purps.BA)

# Santa Barbara Island SMR is only site not linear in before period
KFM.SB.purp <- dplyr::filter(KFM.purps.den, CA_MPA_Name_Short == "Santa Barbara Island SMR")
SumStats <- rbind(SumStats, add_step_effect_size(KFM.SB.purp, "S. purpuratus", "Santa Barbara Island SMR", "KFM", "Den"))

# CI sites
SumStats <- rbind(SumStats, run_ci_analysis(KFM.purps.CI, "S. purpuratus", "KFM", "Den"))

#----- S. purpuratus Biomass (KFM) -----#
KFM.purps.bio <- subset(All.RR.sub, source == "KFM" &
                           y == "Strongylocentrotus purpuratus" & resp == "Bio")
KFM.purps.bio <- add_time_columns(KFM.purps.bio)

# All BA sites with before data were linear in before for density, so removed
# CI sites only
KFM.purps.bio.CI <- subset(KFM.purps.bio, CA_MPA_Name_Short %in% KFM_CI_SITES)
SumStats <- rbind(SumStats, run_ci_analysis(KFM.purps.bio.CI, "S. purpuratus", "KFM", "Bio"))

#----- M. franciscanus Density (KFM) -----#
KFM.reds.den <- subset(All.RR.sub, source == "KFM" &
                          y == "Mesocentrotus franciscanus" & resp == "Den")
KFM.reds.den <- add_time_columns(KFM.reds.den)

KFM.reds.den.BA <- subset(KFM.reds.den, CA_MPA_Name_Short %in% KFM_BA_SITES)
KFM.reds.den.CI <- subset(KFM.reds.den, CA_MPA_Name_Short %in% KFM_CI_SITES)
LinearBefore.kfm.reds <- test_linear_before(KFM.reds.den.BA)

# Scorpion SMR - step model
SumStats <- rbind(SumStats, add_step_effect_size(
  dplyr::filter(KFM.reds.den, CA_MPA_Name_Short == "Scorpion SMR"),
  "M. franciscanus", "Scorpion SMR", "KFM", "Den"
))

# Gull Island SMR - step model
SumStats <- rbind(SumStats, add_step_effect_size(
  dplyr::filter(KFM.reds.den, CA_MPA_Name_Short == "Gull Island SMR"),
  "M. franciscanus", "Gull Island SMR", "KFM", "Den"
))

# Harris Point SMR - linear model
SumStats <- rbind(SumStats, add_linear_effect_size(
  subset(KFM.reds.den, CA_MPA_Name_Short == "Harris Point SMR"),
  "M. franciscanus", "Harris Point SMR", "KFM", "Den"
))

# CI sites
SumStats <- rbind(SumStats, run_ci_analysis(KFM.reds.den.CI, "M. franciscanus", "KFM", "Den"))

#----- M. franciscanus Biomass (KFM) -----#
KFM.reds.bio <- subset(All.RR.sub, source == "KFM" &
                          y == "Mesocentrotus franciscanus" & resp == "Bio")
KFM.reds.bio <- add_time_columns(KFM.reds.bio)

KFM.reds.bio.BA <- subset(KFM.reds.bio, CA_MPA_Name_Short %in% KFM_BA_SITES)
KFM.reds.bio.CI <- subset(KFM.reds.bio, CA_MPA_Name_Short %in% KFM_CI_SITES)

# Run pBACIPS on BA sites
pbacips.kfm.reds.bio <- run_pbacips_loop(KFM.reds.bio.BA)

# Scorpion SMR - step model
SumStats <- rbind(SumStats, add_step_effect_size(
  dplyr::filter(KFM.reds.bio.BA, CA_MPA_Name_Short == "Scorpion SMR"),
  "M. franciscanus", "Scorpion SMR", "KFM", "Bio"
))

# Gull Island SMR - linear model
SumStats <- rbind(SumStats, add_linear_effect_size(
  dplyr::filter(KFM.reds.bio.BA, CA_MPA_Name_Short == "Gull Island SMR"),
  "M. franciscanus", "Gull Island SMR", "KFM", "Bio"
))

# Harris Point SMR - sigmoid model
KFM.mfran.harris <- subset(KFM.reds.bio.BA, CA_MPA_Name_Short == "Harris Point SMR")
if (nrow(KFM.mfran.harris) >= 5) {
  time.model.of.impact <- max(which(KFM.mfran.harris$time.model == 0))
  time.true <- KFM.mfran.harris$time.true
  sigmoid.Model <- tryCatch(
    mySIGfun_standalone(
      delta = KFM.mfran.harris$lnDiff,
      time.model = KFM.mfran.harris$time.model,
      time.model.of.impact = time.model.of.impact,
      time.true = time.true
    ),
    error = function(e) {
      warning("Harris Point M. franciscanus sigmoid failed: ", e$message)
      NULL
    }
  )
  if (!is.null(sigmoid.Model)) {
    # STATISTICAL FIX (2026-02-06): Extract effect size at standardized time point (t=11)
    # rather than at maximum observed time to ensure comparability across MPAs.
    # Create prediction data at time.model = 0 (before) and time.model = 11 (after)
    n_obs_original <- nrow(KFM.mfran.harris)
    time.model.of.impact_original <- max(which(KFM.mfran.harris$time.model == 0))

    # Create standardized prediction points: time 0 (before) and time 11 (after)
    pred_data <- data.frame(
      time.model = c(0, EFFECT_SIZE_TIME_YEARS),
      time.true = c(time.model.of.impact_original, time.model.of.impact_original + EFFECT_SIZE_TIME_YEARS)
    )

    # Get predictions at standardized times with confidence intervals
    interval <- tryCatch({
      data.frame(predFit(sigmoid.Model, newdata = pred_data, interval = "confidence"))
    }, error = function(e) NULL)

    if (!is.null(interval) && nrow(interval) == 2) {
      # Calculate effect size as difference between t=11 and t=0
      mean_es <- interval$fit[2] - interval$fit[1]

      # Calculate SE from confidence interval width (95% CI)
      se_before <- abs((interval$lwr[1] - interval$fit[1]) / 1.96)
      se_after <- abs((interval$lwr[2] - interval$fit[2]) / 1.96)

      # Pooled SE for the difference (conservative: assumes independence)
      # Note: This may slightly overestimate SE but is conservative
      pSE <- sqrt(se_before^2 + se_after^2)
      pCI <- pSE * 1.96

      # Approximate SD from SE (for backward compatibility with meta-analysis)
      # Use original sample size for df approximation
      n_params <- length(coef(sigmoid.Model))
      df_approx <- n_obs_original - n_params
      pSD <- pSE * sqrt(df_approx + 1)

      SumStats[nrow(SumStats) + 1, ] <- c("M. franciscanus", "Harris Point SMR", mean_es, pSE, pSD, pCI,
                                            "Sigmoid", "KFM", "Bio", "Y", "Y", "pBACIPS", "N", n_obs_original,
                                            NA, NA)  # DW_stat, DW_pval: not applicable for NLS models
    }
  }
}

# CI sites
SumStats <- rbind(SumStats, run_ci_analysis(KFM.reds.bio.CI, "M. franciscanus", "KFM", "Bio"))

#----- M. pyrifera Biomass (KFM) -----#
KFM.macro <- subset(All.RR.sub.trans, source == "KFM" &
                      y == "Macrocystis pyrifera" & resp == "Bio")
KFM.macro <- add_time_columns(KFM.macro)

KFM.macro.BA <- subset(KFM.macro, CA_MPA_Name_Short == "Scorpion SMR" |
                          CA_MPA_Name_Short == "Harris Point SMR" |
                          CA_MPA_Name_Short == "Gull Island SMR")
KFM.macro.CI <- subset(KFM.macro, CA_MPA_Name_Short %in%
                          c("Anacapa Island SMR 2003", "South Point SMR", "Santa Barbara Island SMR"))

LinearBefore.kfm.macro <- test_linear_before(KFM.macro.BA)
# Harris Point is linear in before
pbacips.kfm.macro <- run_pbacips_loop(KFM.macro.BA)

# Gull Island SMR - linear
SumStats <- rbind(SumStats, add_linear_effect_size(
  dplyr::filter(KFM.macro.BA, CA_MPA_Name_Short == "Gull Island SMR"),
  "M. pyrifera", "Gull Island SMR", "KFM", "Bio"
))

# Scorpion SMR - sigmoid
KFM.macro.stipe <- dplyr::filter(KFM.macro.BA, CA_MPA_Name_Short == "Scorpion SMR")
if (nrow(KFM.macro.stipe) >= 5) {
  time.model.of.impact <- max(which(KFM.macro.stipe$time.model == 0))
  time.true <- KFM.macro.stipe$time.true
  sigmoid.Model <- tryCatch(
    mySIGfun_standalone(
      delta = KFM.macro.stipe$lnDiff,
      time.model = KFM.macro.stipe$time.model,
      time.model.of.impact = time.model.of.impact,
      time.true = time.true
    ),
    error = function(e) {
      warning("Scorpion M. pyrifera sigmoid failed: ", e$message)
      NULL
    }
  )
  if (!is.null(sigmoid.Model)) {
    # STATISTICAL FIX (2026-02-06): Extract effect size at standardized time point (t=11)
    # rather than at maximum observed time to ensure comparability across MPAs.
    n_obs_original <- nrow(KFM.macro.stipe)
    time.model.of.impact_original <- max(which(KFM.macro.stipe$time.model == 0))

    # Create standardized prediction points: time 0 (before) and time 11 (after)
    pred_data <- data.frame(
      time.model = c(0, EFFECT_SIZE_TIME_YEARS),
      time.true = c(time.model.of.impact_original, time.model.of.impact_original + EFFECT_SIZE_TIME_YEARS)
    )

    # Get predictions at standardized times with confidence intervals
    interval <- tryCatch({
      data.frame(predFit(sigmoid.Model, newdata = pred_data, interval = "confidence"))
    }, error = function(e) NULL)

    if (!is.null(interval) && nrow(interval) == 2) {
      # Calculate effect size as difference between t=11 and t=0
      mean_es <- interval$fit[2] - interval$fit[1]

      # Calculate SE from confidence interval width
      se_before <- abs((interval$lwr[1] - interval$fit[1]) / 1.96)
      se_after <- abs((interval$lwr[2] - interval$fit[2]) / 1.96)

      # Pooled SE for the difference (conservative: assumes independence)
      pSE <- sqrt(se_before^2 + se_after^2)
      pCI <- pSE * 1.96

      # Approximate SD from SE (for backward compatibility)
      n_params <- length(coef(sigmoid.Model))
      df_approx <- n_obs_original - n_params
      pSD <- pSE * sqrt(df_approx + 1)

      SumStats[nrow(SumStats) + 1, ] <- c("M. pyrifera", "Scorpion SMR", mean_es, pSE, pSD, pCI,
                                            "Sigmoid", "KFM", "Bio", "Y", "Y", "pBACIPS", "N", n_obs_original,
                                            NA, NA)  # DW_stat, DW_pval: not applicable for NLS models
    }
  }
}

# CI sites
SumStats <- rbind(SumStats, run_ci_analysis(KFM.macro.CI, "M. pyrifera", "KFM", "Bio"))

#----- P. interruptus Density (KFM) -----#
KFM.lob <- subset(All.RR.sub.trans, source == "KFM" &
                     y == "Panulirus interruptus")
KFM.lob$CA_MPA_Name_Short <- as.character(KFM.lob$CA_MPA_Name_Short)
KFM.lob <- add_time_columns(KFM.lob)

# Dropped Harris, no lobsters in before suggests not good sampling
KFM.lob.BA <- subset(KFM.lob, CA_MPA_Name_Short %in%
                        c("Scorpion SMR", "Santa Barbara Island SMR", "Gull Island SMR"))
KFM.lob.CI <- subset(KFM.lob, CA_MPA_Name_Short %in% KFM_CI_SITES)

LinearBefore.kfm.lob <- test_linear_before(KFM.lob.BA)

# Run pBACIPS on subset (exclude Harris and SBI from pBACIPS loop)
KFM.lob.sub <- subset(KFM.lob.BA, CA_MPA_Name_Short != "Harris Point SMR" &
                         CA_MPA_Name_Short != "Santa Barbara Island SMR")
pbacips.kfm.lob <- run_pbacips_loop(KFM.lob.sub)

# Gull Island SMR - sigmoid
KFM.lob.gull <- dplyr::filter(KFM.lob.sub, CA_MPA_Name_Short == "Gull Island SMR")
if (nrow(KFM.lob.gull) >= 5) {
  time.model.of.impact <- max(which(KFM.lob.gull$time == 0))
  time.true <- KFM.lob.gull$time
  sigmoid.Model <- tryCatch(
    mySIGfun_standalone(
      delta = KFM.lob.gull$lnDiff,
      time.model = KFM.lob.gull$time.model,
      time.model.of.impact = time.model.of.impact,
      time.true = time.true
    ),
    error = function(e) {
      warning("Gull Island P. interruptus sigmoid failed: ", e$message)
      NULL
    }
  )
  if (!is.null(sigmoid.Model)) {
    # STATISTICAL FIX (2026-02-06): Extract effect size at standardized time point (t=11)
    # rather than at maximum observed time to ensure comparability across MPAs.
    n_obs_original <- nrow(KFM.lob.gull)
    time.model.of.impact_original <- max(which(KFM.lob.gull$time == 0))

    # Create standardized prediction points: time 0 (before) and time 11 (after)
    pred_data <- data.frame(
      time.model = c(0, EFFECT_SIZE_TIME_YEARS),
      time.true = c(time.model.of.impact_original, time.model.of.impact_original + EFFECT_SIZE_TIME_YEARS)
    )

    # Get predictions at standardized times with confidence intervals
    interval <- tryCatch({
      data.frame(predFit(sigmoid.Model, newdata = pred_data, interval = "confidence"))
    }, error = function(e) NULL)

    if (!is.null(interval) && nrow(interval) == 2) {
      # Calculate effect size as difference between t=11 and t=0
      mean_es <- interval$fit[2] - interval$fit[1]

      # Calculate SE from confidence interval width
      se_before <- abs((interval$lwr[1] - interval$fit[1]) / 1.96)
      se_after <- abs((interval$lwr[2] - interval$fit[2]) / 1.96)

      # Pooled SE for the difference (conservative: assumes independence)
      pSE <- sqrt(se_before^2 + se_after^2)
      pCI <- pSE * 1.96

      # Approximate SD from SE (for backward compatibility)
      n_params <- length(coef(sigmoid.Model))
      df_approx <- n_obs_original - n_params
      pSD <- pSE * sqrt(df_approx + 1)

      SumStats[nrow(SumStats) + 1, ] <- c("P. interruptus", "Gull Island SMR", mean_es, pSE, pSD, pCI,
                                            "Sigmoid", "KFM", "Den", "Y", "Y", "pBACIPS", "N", n_obs_original,
                                            NA, NA)  # DW_stat, DW_pval: not applicable for NLS models
    }
  }
}

# Scorpion SMR - linear
SumStats <- rbind(SumStats, add_linear_effect_size(
  dplyr::filter(KFM.lob.sub, CA_MPA_Name_Short == "Scorpion SMR"),
  "P. interruptus", "Scorpion SMR", "KFM", "Den"
))

# Santa Barbara Island SMR - linear
KFM.lobs.sbi <- subset(KFM.lob.BA, CA_MPA_Name_Short == "Santa Barbara Island SMR")
SumStats <- rbind(SumStats, add_linear_effect_size(KFM.lobs.sbi, "P. interruptus", "Santa Barbara Island SMR", "KFM", "Den"))

# CI sites
SumStats <- rbind(SumStats, run_ci_analysis(KFM.lob.CI, "P. interruptus", "KFM", "Den"))

#----- S. pulcher Density (KFM) -----#
Sheep.den.kfm <- subset(All.RR.sub.trans, y == "Semicossyphus pulcher" &
                           resp == "Den" & source == "KFM")
Sheep.den.kfm <- add_time_columns(Sheep.den.kfm)

KFM.sheep.BA <- subset(Sheep.den.kfm, CA_MPA_Name_Short %in% KFM_BA_SITES)
KFM.sheep.CI <- subset(Sheep.den.kfm, CA_MPA_Name_Short %in% KFM_CI_SITES)

LinearBefore.kfm.sheep <- test_linear_before(KFM.sheep.BA)

# Harris Point SMR - linear
SumStats <- rbind(SumStats, add_linear_effect_size(
  dplyr::filter(KFM.sheep.BA, CA_MPA_Name_Short == "Harris Point SMR"),
  "S. pulcher", "Harris Point SMR", "KFM", "Den"
))

# Scorpion SMR - step
SumStats <- rbind(SumStats, add_step_effect_size(
  dplyr::filter(KFM.sheep.BA, CA_MPA_Name_Short == "Scorpion SMR"),
  "S. pulcher", "Scorpion SMR", "KFM", "Den"
))

# CI sites
SumStats <- rbind(SumStats, run_ci_analysis(KFM.sheep.CI, "S. pulcher", "KFM", "Den"))


####################################################################################################
## LANDSAT Kelp Biomass Effect Sizes ###############################################################
####################################################################################################

# Landsat data processing is now handled by 06b_landsat_processing.R which produces Landsat.RR.
# Here we only run the effect size analysis loop over the pre-processed response ratios.

if (exists("Landsat.RR") && nrow(Landsat.RR) > 0) {
  mpas_landsat <- unique(Landsat.RR$CA_MPA_Name_Short)
  for (mpa in mpas_landsat) {
    mpa_dat <- dplyr::filter(Landsat.RR, CA_MPA_Name_Short == mpa)
    if (nrow(mpa_dat) > 3) {
      SumStats <- rbind(SumStats, add_linear_effect_size(mpa_dat, "M. pyrifera", mpa, "Landsat", "Bio",
                             time_var = "time"))
    }
  }
}

####################################################################################################
## Finalize SumStats ###############################################################################
####################################################################################################

# Convert numeric columns
SumStats$Mean <- as.numeric(SumStats$Mean)
SumStats$SE <- as.numeric(SumStats$SE)
SumStats$CI <- as.numeric(SumStats$CI)
SumStats$SD <- as.numeric(SumStats$SD)
SumStats$DW_stat <- as.numeric(SumStats$DW_stat)
SumStats$DW_pval <- as.numeric(SumStats$DW_pval)

# Remove incomplete rows (exclude DW columns from completeness check since
# DW is NA for mean-only and NLS models by design)
core_cols <- setdiff(names(SumStats), c("DW_stat", "DW_pval"))
SumStats.sub <- SumStats[complete.cases(SumStats[, core_cols]), ]

# Set factor levels for plotting
SumStats.sub$Taxa <- factor(SumStats.sub$Taxa,
  levels = c("S. purpuratus", "M. franciscanus", "M. pyrifera", "P. interruptus", "S. pulcher"),
  ordered = TRUE
)
SumStats.sub$Resp <- factor(SumStats.sub$Resp, levels = c("Bio", "Den"), ordered = TRUE)
SumStats.sub$Source <- factor(SumStats.sub$Source,
  levels = c("KFM", "LTER", "PISCO", "Landsat"), ordered = TRUE
)
SumStats.sub$Type <- factor(SumStats.sub$Type,
  levels = c("pBACIPS", "BACI", "CI"), ordered = TRUE
)

# Merge with site metadata
# Rename Type column to AnalysisType before merge to avoid collision with Site$type
names(SumStats.sub)[names(SumStats.sub) == "Type"] <- "AnalysisType"
SumStats.sub <- merge(SumStats.sub, Site, by.x = "MPA", by.y = "CA_MPA_Name_Short")

####################################################################################################
## Filter SumStats.Final for publication ###########################################################
####################################################################################################

# Include pBACIPS results, plus CI results where before period was not linear and is primary
SumStats.Final <- subset(SumStats.sub,
  AnalysisType == "pBACIPS" | (AnalysisType == "CI" & LinearBefore == "NA" & Primary == "Y")
)

# Remove problematic MPAs:
# Painted Cave SMCA (lobster take allowed), San Miguel Island SC (weird location),
# Arrow Point (finfish fishing allowed), Judith Rk SMR (overlaps San Miguel Island SC),
# Point Conception SMR
SumStats.Final <- subset(SumStats.Final,
  MPA != "Painted Cave SMCA" &
  MPA != "San Miguel Island SC" &
  MPA != "Arrow Point to Lion Head Point SMCA" &
  MPA != "Judith Rk SMR" &
  MPA != "Point Conception SMR"
)

cat("Effect size calculation complete.\n")
cat("SumStats rows:", nrow(SumStats), "\n")
cat("SumStats.Final rows:", nrow(SumStats.Final), "\n")

####################################################################################################
## FILTERING AUDIT: Track exactly what is excluded and why ########################################
####################################################################################################

cat("\n")
cat("====================================\n")
cat("FILTERING AUDIT: Effect Size Data\n")
cat("====================================\n")

# Start with all complete cases from SumStats (exclude DW columns from completeness check)
SumStats_complete <- SumStats[complete.cases(SumStats[, core_cols]), ]
cat("Step 0: Complete cases in SumStats:", nrow(SumStats_complete), "\n")

# Create audit dataframe tracking each observation
FilterAudit <- SumStats_complete[, c("Taxa", "MPA", "Source", "Resp", "Mean", "SE", "Model", "Primary", "Type", "LinearBefore")]
FilterAudit$Step0_Complete <- TRUE

# Step 1: Check if MPA exists in Site metadata (merge success)
FilterAudit$Step1_SiteMerge <- FilterAudit$MPA %in% Site$CA_MPA_Name_Short
FilterAudit$Step1_Reason <- ifelse(FilterAudit$Step1_SiteMerge, "OK", "MPA not in Site metadata")

# Step 2: Analysis type filter (pBACIPS OR (CI + LinearBefore==NA + Primary==Y))
FilterAudit$Step2_TypeFilter <- (
  FilterAudit$Type == "pBACIPS" |
  (FilterAudit$Type == "CI" & FilterAudit$LinearBefore == "NA" & FilterAudit$Primary == "Y")
)
FilterAudit$Step2_Reason <- ifelse(FilterAudit$Step2_TypeFilter, "OK",
  ifelse(FilterAudit$Type == "CI" & FilterAudit$LinearBefore != "NA",
         paste0("CI but LinearBefore=", FilterAudit$LinearBefore),
  ifelse(FilterAudit$Type == "CI" & FilterAudit$Primary != "Y",
         "CI but Primary=N (not the selected model)",
  ifelse(FilterAudit$Type == "BACI", "BACI type not included",
         paste0("Unknown type: ", FilterAudit$Type)))))

# Step 3: Excluded MPA filter (uses EXCLUDED_MPAS from 01_utils.R)
FilterAudit$Step3_MPAFilter <- !(FilterAudit$MPA %in% EXCLUDED_MPAS)
FilterAudit$Step3_Reason <- ifelse(FilterAudit$Step3_MPAFilter, "OK",
  ifelse(FilterAudit$MPA == "Painted Cave SMCA", "Excluded: lobster take allowed",
  ifelse(FilterAudit$MPA == "San Miguel Island SC", "Excluded: unusual location",
  ifelse(FilterAudit$MPA == "Arrow Point to Lion Head Point SMCA", "Excluded: finfish fishing allowed",
  ifelse(FilterAudit$MPA == "Judith Rk SMR", "Excluded: overlaps San Miguel Island SC",
  ifelse(FilterAudit$MPA == "Point Conception SMR", "Excluded: data quality",
         "Excluded: unknown reason"))))))

# Overall: passes all filters
FilterAudit$Passes_All <- FilterAudit$Step1_SiteMerge & FilterAudit$Step2_TypeFilter & FilterAudit$Step3_MPAFilter

# Determine exclusion reason (first failing step)
FilterAudit$Exclusion_Reason <- ifelse(FilterAudit$Passes_All, "INCLUDED",
  ifelse(!FilterAudit$Step1_SiteMerge, FilterAudit$Step1_Reason,
  ifelse(!FilterAudit$Step2_TypeFilter, FilterAudit$Step2_Reason,
         FilterAudit$Step3_Reason)))

# Summary statistics
cat("\n--- Filtering Summary ---\n")
cat("Total complete effect sizes:", nrow(FilterAudit), "\n")
cat("  Pass Site merge:", sum(FilterAudit$Step1_SiteMerge), "\n")
cat("  Pass Type filter:", sum(FilterAudit$Step2_TypeFilter), "\n")
cat("  Pass MPA filter:", sum(FilterAudit$Step3_MPAFilter), "\n")
cat("  Pass ALL filters:", sum(FilterAudit$Passes_All), "\n")
cat("  Excluded:", sum(!FilterAudit$Passes_All), "\n")

# Breakdown by exclusion reason
cat("\n--- Exclusion Reasons ---\n")
reason_table <- table(FilterAudit$Exclusion_Reason)
for (reason in names(sort(reason_table, decreasing = TRUE))) {
  cat(sprintf("  %s: %d\n", reason, reason_table[reason]))
}

# Breakdown by taxa
cat("\n--- By Taxa (Passes All) ---\n")
taxa_summary <- aggregate(Passes_All ~ Taxa + Resp, data = FilterAudit,
                          FUN = function(x) c(total = length(x), included = sum(x)))
for (i in seq_len(nrow(taxa_summary))) {
  cat(sprintf("  %s (%s): %d/%d included\n",
              taxa_summary$Taxa[i], taxa_summary$Resp[i],
              taxa_summary$Passes_All[i, "included"],
              taxa_summary$Passes_All[i, "total"]))
}

# Write detailed audit to CSV
audit_file <- here::here("outputs", "filter_audit_effect_sizes.csv")
if (!dir.exists(here::here("outputs"))) {
  dir.create(here::here("outputs"), recursive = TRUE)
}
write.csv(FilterAudit, audit_file, row.names = FALSE)
cat("\nDetailed filter audit saved to:", audit_file, "\n")

# Create taxa-specific summary
TaxaSummary <- FilterAudit %>%
  dplyr::group_by(Taxa, Resp) %>%
  dplyr::summarise(
    Total_Generated = dplyr::n(),
    Passed_SiteMerge = sum(Step1_SiteMerge),
    Passed_TypeFilter = sum(Step2_TypeFilter),
    Passed_MPAFilter = sum(Step3_MPAFilter),
    Final_Included = sum(Passes_All),
    Excluded = sum(!Passes_All),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Taxa, Resp)

taxa_summary_file <- here::here("outputs", "filter_summary_by_taxa.csv")
write.csv(TaxaSummary, taxa_summary_file, row.names = FALSE)
cat("Taxa summary saved to:", taxa_summary_file, "\n")

# Focus on lobster for the specific question
cat("\n--- LOBSTER (P. interruptus) DETAIL ---\n")
lobster_audit <- subset(FilterAudit, Taxa == "P. interruptus")
if (nrow(lobster_audit) > 0) {
  cat("Total lobster effect sizes generated:", nrow(lobster_audit), "\n")
  cat("By response type:\n")
  for (resp in unique(lobster_audit$Resp)) {
    lob_resp <- subset(lobster_audit, Resp == resp)
    cat(sprintf("  %s: %d total, %d included\n", resp, nrow(lob_resp), sum(lob_resp$Passes_All)))
  }
  cat("\nLobster exclusions:\n")
  lob_excluded <- subset(lobster_audit, !Passes_All)
  if (nrow(lob_excluded) > 0) {
    for (i in seq_len(nrow(lob_excluded))) {
      cat(sprintf("  - %s %s (%s): %s\n",
                  lob_excluded$MPA[i], lob_excluded$Resp[i],
                  lob_excluded$Source[i], lob_excluded$Exclusion_Reason[i]))
    }
  } else {
    cat("  (none excluded)\n")
  }
  cat("\nLobster inclusions:\n")
  lob_included <- subset(lobster_audit, Passes_All)
  if (nrow(lob_included) > 0) {
    for (i in seq_len(nrow(lob_included))) {
      cat(sprintf("  + %s %s (%s): Effect=%.3f\n",
                  lob_included$MPA[i], lob_included$Resp[i],
                  lob_included$Source[i], as.numeric(lob_included$Mean[i])))
    }
  }
}

####################################################################################################
## Model Diagnostics Summary #######################################################################
####################################################################################################

if (exists("ModelDiagnostics") && nrow(ModelDiagnostics) > 0) {
  cat("\n")
  cat("============================\n")
  cat("MODEL DIAGNOSTICS SUMMARY\n")
  cat("============================\n")

  # Remove any rows with all NA
  ModelDiagnostics <- ModelDiagnostics[!is.na(ModelDiagnostics$Model_Type), ]

  if (nrow(ModelDiagnostics) > 0) {
    # Summary statistics
    n_total <- nrow(ModelDiagnostics)
    n_pass <- sum(ModelDiagnostics$Pass_All == TRUE, na.rm = TRUE)
    n_fail <- sum(ModelDiagnostics$Pass_All == FALSE, na.rm = TRUE)
    n_na <- sum(is.na(ModelDiagnostics$Pass_All))

    cat(sprintf("Total models diagnosed: %d\n", n_total))
    cat(sprintf("  Pass all tests: %d (%.1f%%)\n", n_pass, 100 * n_pass / n_total))
    cat(sprintf("  Fail one or more: %d (%.1f%%)\n", n_fail, 100 * n_fail / n_total))
    if (n_na > 0) cat(sprintf("  Could not diagnose: %d\n", n_na))

    # Breakdown by model type
    cat("\nBy Model Type:\n")
    type_summary <- aggregate(Pass_All ~ Model_Type, data = ModelDiagnostics,
                               FUN = function(x) c(n = length(x), pass = sum(x, na.rm = TRUE)))
    for (i in seq_len(nrow(type_summary))) {
      mt <- type_summary$Model_Type[i]
      n <- type_summary$Pass_All[i, "n"]
      pass <- type_summary$Pass_All[i, "pass"]
      cat(sprintf("  %s: %d/%d pass (%.1f%%)\n", mt, pass, n, 100 * pass / n))
    }

    # List models with issues
    failed_models <- ModelDiagnostics[ModelDiagnostics$Pass_All == FALSE & !is.na(ModelDiagnostics$Pass_All), ]
    if (nrow(failed_models) > 0) {
      cat("\nModels with diagnostic issues:\n")
      # Show up to 10 failures
      show_n <- min(10, nrow(failed_models))
      for (i in seq_len(show_n)) {
        row <- failed_models[i, ]
        cat(sprintf("  - %s %s (%s): %s\n", row$Taxa, row$MPA, row$Model_Type, row$Notes))
      }
      if (nrow(failed_models) > 10) {
        cat(sprintf("  ... and %d more\n", nrow(failed_models) - 10))
      }
    }

    # Save diagnostics to file
    diag_file <- here::here("data", "model_diagnostics.csv")
    write.csv(ModelDiagnostics, diag_file, row.names = FALSE)
    cat(sprintf("\nDiagnostics saved to: %s\n", diag_file))

    # Create diagnostic plots directory if needed
    diag_dir <- here::here("plots", "diagnostics")
    if (!dir.exists(diag_dir)) {
      dir.create(diag_dir, recursive = TRUE)
    }

    # Generate summary plot if ggplot2 is available
    if (requireNamespace("ggplot2", quietly = TRUE) && n_total > 0) {
      # Summary bar plot
      diag_summary <- data.frame(
        Status = c("Pass", "Fail", "NA"),
        Count = c(n_pass, n_fail, n_na)
      )
      diag_summary <- diag_summary[diag_summary$Count > 0, ]

      p <- ggplot2::ggplot(diag_summary, ggplot2::aes(x = Status, y = Count, fill = Status)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::scale_fill_manual(values = c("Pass" = "#4CAF50", "Fail" = "#F44336", "NA" = "#9E9E9E")) +
        ggplot2::labs(title = "Model Diagnostic Results",
                      subtitle = sprintf("%d models tested with DHARMa", n_total),
                      x = "", y = "Number of Models") +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none")

      ggplot2::ggsave(file.path(diag_dir, "diagnostic_summary.pdf"), p, width = 6, height = 4)
      ggplot2::ggsave(file.path(diag_dir, "diagnostic_summary.png"), p, width = 6, height = 4, dpi = 150)
      cat(sprintf("Summary plot saved to: %s\n", file.path(diag_dir, "diagnostic_summary.pdf")))
    }
  }
} else {
  cat("\nNo model diagnostics recorded (DHARMa may not be installed).\n")
}

####################################################################################################
## Durbin-Watson Autocorrelation Summary ###########################################################
####################################################################################################

cat("\n")
cat("==========================================\n")
cat("DURBIN-WATSON AUTOCORRELATION DIAGNOSTICS\n")
cat("==========================================\n")

# Summarize DW results across all effect size models
dw_available <- SumStats.sub[!is.na(SumStats.sub$DW_stat), ]
dw_na <- SumStats.sub[is.na(SumStats.sub$DW_stat), ]

cat(sprintf("Models with DW test: %d\n", nrow(dw_available)))
cat(sprintf("Models without DW test (Mean/NLS): %d\n", nrow(dw_na)))

if (nrow(dw_available) > 0) {
  n_sig <- sum(dw_available$DW_pval < 0.05, na.rm = TRUE)
  n_nonsig <- sum(dw_available$DW_pval >= 0.05, na.rm = TRUE)
  cat(sprintf("\nSignificant autocorrelation (p < 0.05): %d/%d (%.1f%%)\n",
              n_sig, nrow(dw_available), 100 * n_sig / nrow(dw_available)))
  cat(sprintf("No significant autocorrelation: %d/%d (%.1f%%)\n",
              n_nonsig, nrow(dw_available), 100 * n_nonsig / nrow(dw_available)))
  cat(sprintf("Mean DW statistic: %.3f (values near 2 indicate no autocorrelation)\n",
              mean(dw_available$DW_stat, na.rm = TRUE)))

  # List models with significant autocorrelation
  if (n_sig > 0) {
    sig_models <- dw_available[dw_available$DW_pval < 0.05, ]
    cat("\nModels with significant autocorrelation:\n")
    for (i in seq_len(nrow(sig_models))) {
      row <- sig_models[i, ]
      cat(sprintf("  - %s %s (%s, %s): DW=%.3f, p=%.4f\n",
                  row$Taxa, row$MPA, row$Source, row$Model,
                  row$DW_stat, row$DW_pval))
    }
    cat("\nNote: Significant autocorrelation may inflate Type I error rates by\n")
    cat("underestimating SEs. However, the meta-analytic framework (rma.mv in\n")
    cat("09_meta_analysis.R) partially mitigates this by modeling heterogeneity\n")
    cat("across studies. Consider GLS or Newey-West corrections for a sensitivity analysis.\n")
  } else {
    cat("\nNo models show significant temporal autocorrelation. The independence\n")
    cat("assumption for residuals appears reasonable across all effect size models.\n")
  }
}

####################################################################################################
## Memory cleanup
####################################################################################################
# Remove intermediate analysis objects to free memory
# Keep only: SumStats.Final, ModelDiagnostics (if needed)

# LTER analysis intermediates
rm(list = intersect(ls(), c(
  "LTER.purps", "LTER.reds", "LTER.reds.bio", "LTER.macro", "LTER.SPUL.den", "LTER.SPUL.bio",
  "LTER.PANINT.den", "LTER.PANINT.bio", "LTER.nap.red", "LTER.cp.macro", "LTER.cp.SPUL.den",
  "LTER.cp.SPUL.bio", "LTER.cp.LOB", "LTER.nap.LOB", "LTER.cp.LOB.bio", "LTER.nap.LOB.bio",
  "LTER.redsNap", "LinearBefore.purps", "LinearBefore.reds", "LinearBefore.macro",
  "LinearBefore.spul.den", "LinearBefore.spul.bio", "pbacips.reds.den", "pbacips.macro",
  "pbacips.spul.den", "pbacips.spul.bio"
)))

# PISCO analysis intermediates
rm(list = intersect(ls(), c(
  "Mes.den", "Mes.bio", "Str.den", "Str.bio", "Mac", "Pan.den", "Pan.bio",
  "Sheep.den", "Sheep.bio"
)))

# KFM analysis intermediates
rm(list = intersect(ls(), c(
  "KFM.purps.den", "KFM.purps.BA", "KFM.purps.CI", "KFM.purps.bio", "KFM.purps.bio.CI",
  "KFM.reds.den", "KFM.reds.den.BA", "KFM.reds.den.CI", "KFM.reds.bio",
  "KFM.macro", "KFM.macro.BA", "KFM.macro.CI", "KFM.macro.stipe",
  "KFM.lob", "KFM.lob.BA", "KFM.lob.CI", "KFM.lob.sub", "KFM.lob.gull", "KFM.lobs.sbi",
  "KFM.mfran.harris", "Sheep.den.kfm", "KFM.sheep.BA", "KFM.sheep.CI",
  "LinearBefore.kfm.purps", "LinearBefore.kfm.reds", "LinearBefore.kfm.macro",
  "LinearBefore.kfm.lob", "LinearBefore.kfm.sheep",
  "pbacips.kfm.reds.bio", "pbacips.kfm.lob", "KFM.SB.purp"
)))

# Model objects and intermediates
rm(list = intersect(ls(), c(
  "sigmoid.Model", "interval", "mpas_landsat", "mpa_dat"
)))

# Force garbage collection
gc(verbose = FALSE)

cat("\nEffect size calculation complete.\n")
cat("  Output: SumStats.Final (", nrow(SumStats.Final), " effect sizes)\n", sep = "")
