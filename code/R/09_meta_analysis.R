# =============================================================================
# 09_meta_analysis.R
# =============================================================================
#
# PURPOSE:
#   Perform multilevel meta-analysis of MPA effects across all taxa, data sources,
#   and MPAs. This aggregates the individual effect sizes from 08_effect_sizes.R
#   to estimate overall MPA effects by taxa.
#
# WHAT THIS SCRIPT DOES:
#   1. Prepares effect size data for meta-analysis (calculates sampling variance)
#   2. Fits joint multilevel models (for sensitivity comparison & variance components)
#   3. Fits PER-TAXA meta-analyses (PRIMARY) with within-taxon outlier detection
#   4. Produces Table 2: Meta-analytic estimates by taxa (from per-taxa models)
#   5. Produces outlier sensitivity table (full data vs per-taxa vs joint Cook's D)
#   6. Produces leave-one-out analysis for borderline kelp result
#   7. Produces Table 3: Cross-taxa relationships (trophic cascade tests)
#
# PRIMARY META-ANALYSIS STRUCTURE (Per-Taxa):
#   For each taxon-response separately:
#   yi ~ 1, random = ~1|MPA   (for k >= 5 and >= 3 MPAs)
#   yi ~ 1                     (simple RE for k < 5)
#   Cook's D applied within each taxon at threshold 4/k
#
# SENSITIVITY: JOINT MODEL STRUCTURE:
#   yi ~ Taxa - 1, random = list(~1|MPA, ~1|Source)
#   Cook's D at 4/n across all taxa (retained for comparison)
#
# WHY PER-TAXA INSTEAD OF JOINT MODEL:
#   The joint model applies Cook's D globally, flagging observations that differ
#   from OTHER taxa, not just unusual within their own taxon. This resulted in
#   62% outlier removal in the joint model. Per-taxa Cook's D correctly asks
#   "Is this observation influential for THIS taxon's estimate?"
#
# WHY SEPARATE MODELS FOR BIOMASS AND DENSITY:
#   Biomass and density estimates for the same taxa/MPA/source are non-independent
#   (they come from the same surveys). Running them together would violate
#   independence assumptions. By splitting, we can properly estimate variance.
#
# METHODS REFERENCE:
#   "We fit multilevel meta-analysis models with restricted maximum-likelihood
#    estimation (REML) using effect size estimates (ES) from pBACIPS and/or CI
#    analyses for each taxa using the rma.mv function in the metafor package.
#    Taxa was a fixed effect (moderator) in meta-analysis models with MPA as
#    a random effect."
#
# INPUTS:
#   - SumStats.Final: Effect size estimates from 08_effect_sizes.R
#
# OUTPUTS:
#   - meta_biomass: rma.mv model object for biomass
#   - meta_density: rma.mv model object for density
#   - Table2: Summary dataframe (also exported as CSV)
#   - table_02_meta_analysis.csv: Written to project root
#
# DEPENDENCIES:
#   Requires 00-08 scripts to be sourced first
#
# AUTHORS: Emily Donham & Adrian Stier
# PROJECT: CA MPA Kelp Forest pBACIPS Analysis
# =============================================================================
#
# This script implements the multilevel meta-analysis described in the Methods:
#'
#'   "To test for an effect of MPA implementation on species biomass/density, we fit
#'    multilevel meta-analysis models with restricted maximum-likelihood estimation
#'    (REML) using effect size estimates (ES) from pBACIPS and/or CI analyses for each
#'    taxa using the rma.mv function in the metafor package. Taxa was a fixed effect
#'    (moderator) in meta-analysis models with MPA as a random effect. We ran models
#'    separately for the responses of density and biomass since estimates were
#'    non-independent. Initially, models were fit using all data. We then assessed the
#'    influence of outliers using Cook's distance and removed influential outliers
#'    before rerunning analyses."
#'
#' Inputs:
#'   SumStats.Final  - Filtered effect size estimates from 08_effect_sizes.R
#'                     Columns used: Taxa, MPA, Mean, SE, Resp
#'
#' Outputs:
#'   meta_biomass    - rma.mv model object for biomass (after outlier removal)
#'   meta_density    - rma.mv model object for density (after outlier removal)
#'   Table2          - Summary dataframe matching Table 2 in the paper
#'   Table2 CSV      - Exported to project root as "table_02_meta_analysis.csv"

####################################################################################################
## Prepare data for meta-analysis ##################################################################
####################################################################################################

# Compute sampling variance (vi = SE^2)
# SE is the standard error of the effect size estimate
SumStats.Final$vi <- as.numeric(SumStats.Final$SE)^2

# Ensure Mean and vi are numeric (they may have been stored as character from 08)
SumStats.Final$Mean <- as.numeric(SumStats.Final$Mean)
SumStats.Final$vi   <- as.numeric(SumStats.Final$vi)

# --- OUTPUT VALIDATION: Check variance calculation ---
n_invalid_se <- sum(SumStats.Final$SE <= 0, na.rm = TRUE)
n_invalid_vi <- sum(!is.finite(SumStats.Final$vi))
if (n_invalid_se > 0) {
  warning("Found ", n_invalid_se, " effect sizes with SE <= 0")
}
if (n_invalid_vi > 0) {
  warning("Found ", n_invalid_vi, " invalid variance values (NA/NaN/Inf)")
  cat("  Removing", n_invalid_vi, "rows with invalid variance\n")
  SumStats.Final <- SumStats.Final[is.finite(SumStats.Final$vi), ]
}
stopifnot("No valid effect sizes remaining" = nrow(SumStats.Final) > 0)

# --- OUTPUT VALIDATION: Check effect size range ---
extreme_es <- sum(abs(SumStats.Final$Mean) > 5, na.rm = TRUE)
if (extreme_es > 0) {
  warning("Found ", extreme_es, " effect sizes with |Mean| > 5 (may be outliers)")
}

# Ensure Taxa is a factor for the moderator specification
SumStats.Final$Taxa <- factor(SumStats.Final$Taxa)

# Ensure Source is a factor for random effects nesting
SumStats.Final$Source <- factor(SumStats.Final$Source)

# Split into biomass and density subsets (non-independent, so modeled separately)
biomass_data <- subset(SumStats.Final, Resp == "Bio")
density_data <- subset(SumStats.Final, Resp == "Den")

# --- OUTPUT VALIDATION: Ensure both subsets have sufficient data ---
stopifnot("No biomass data for meta-analysis" = nrow(biomass_data) > 0)
stopifnot("No density data for meta-analysis" = nrow(density_data) > 0)

cat("Biomass observations:", nrow(biomass_data), "\n")
cat("Density observations:", nrow(density_data), "\n")

####################################################################################################
## JOINT MODEL: Biomass meta-analysis (Sensitivity / Variance Components) ##########################
####################################################################################################

# Step 1: Fit initial model with all biomass data
# - yi = Mean: the effect size (mean difference on log-ratio scale)
# - V = vi: sampling variance (SE^2)
# - mods = ~ Taxa - 1: taxa as a fixed-effect moderator (cell-means parameterization,
#   removing the intercept so each taxon gets its own estimate)
# - random = ~ 1 | MPA: random intercept for MPA to account for non-independence
#   of effect sizes within the same MPA
# - method = "REML": restricted maximum-likelihood estimation
meta_biomass_full <- tryCatch(
  rma.mv(
    yi     = Mean,
    V      = vi,
    mods   = ~ Taxa - 1,
    random = list(~ 1 | MPA, ~ 1 | Source),
    data   = biomass_data,
    method = "REML",
    test   = "t"
  ),
  error = function(e) {
    stop("Biomass meta-analysis model failed to converge: ", e$message)
  }
)

# --- OUTPUT VALIDATION: Check model convergence ---
if (!is.null(meta_biomass_full$not.converged) && meta_biomass_full$not.converged) {
  warning("Biomass meta-analysis model did not converge")
}
stopifnot("Model coefficients contain NA" = !any(is.na(coef(meta_biomass_full))))
stopifnot("Model coefficients not finite" = all(is.finite(coef(meta_biomass_full))))

cat("\n--- Biomass model (full data) ---\n")
print(summary(meta_biomass_full))

# Step 2: Cook's distance outlier detection
# Cook's distance measures the influence of each observation on the fitted model.
# Observations with Cook's distance exceeding the conventional threshold (cutoff based
# on the number of observations) are considered influential outliers.
cooks_bio <- cooks.distance(meta_biomass_full)
n_bio <- nrow(biomass_data)
p_bio <- length(coef(meta_biomass_full))

# Standard threshold: values exceeding 4/n or the median + 3*IQR
# We use the 4/n rule which is standard in meta-analytic practice
cooks_threshold_bio <- 4 / n_bio
outliers_bio <- which(cooks_bio > cooks_threshold_bio)

cat("\nBiomass Cook's distance threshold (4/n):", round(cooks_threshold_bio, 4), "\n")
cat("Biomass outliers detected:", length(outliers_bio), "\n")
if (length(outliers_bio) > 0) {
  cat("Outlier rows (Taxa, MPA, Mean):\n")
  print(biomass_data[outliers_bio, c("Taxa", "MPA", "Mean", "SE")])
}

# Step 3: Remove outliers and refit
if (length(outliers_bio) > 0) {
  # --- OUTPUT VALIDATION: Check outlier removal is not excessive ---
  pct_removed <- 100 * length(outliers_bio) / n_bio
  cat("Removing", length(outliers_bio), "outliers (", round(pct_removed, 1), "% of data)\n")
  if (pct_removed > 20) {
    warning("Removing >20% of biomass data as outliers - consider reviewing threshold")
  }
  biomass_clean <- biomass_data[-outliers_bio, ]
} else {
  biomass_clean <- biomass_data
}

# --- OUTPUT VALIDATION: Ensure sufficient data remains ---
stopifnot("Insufficient biomass data after outlier removal" = nrow(biomass_clean) >= 10)

meta_biomass <- tryCatch(
  rma.mv(
    yi     = Mean,
    V      = vi,
    mods   = ~ Taxa - 1,
    random = list(~ 1 | MPA, ~ 1 | Source),
    data   = biomass_clean,
    method = "REML",
    test   = "t"
  ),
  error = function(e) {
    stop("Biomass meta-analysis (clean) failed: ", e$message)
  }
)

cat("\n--- Biomass model (outliers removed) ---\n")
print(summary(meta_biomass))

# Report heterogeneity statistics for biomass model
# tau² (tau-squared): Between-study variance components
# I²: Percentage of total variability due to heterogeneity (vs sampling error)
cat("\n--- Biomass Heterogeneity Statistics ---\n")
cat("Between-MPA variance (tau²_MPA):", round(meta_biomass$sigma2[1], 4), "\n")
cat("Between-Source variance (tau²_Source):", round(meta_biomass$sigma2[2], 4), "\n")
# Calculate pseudo-I² for multilevel models
# Total heterogeneity variance = sum of random effect variances
# Typical within-study variance approximated by mean sampling variance
total_hetero_bio <- sum(meta_biomass$sigma2)
# Use harmonic mean of sampling variances (more robust to outlier variances;
# see Nakagawa & Santos 2012 for multilevel I² approximation)
typical_v_bio <- 1 / mean(1 / biomass_clean$vi)
pseudo_I2_bio <- 100 * total_hetero_bio / (total_hetero_bio + typical_v_bio)
cat("Pseudo-I² (total):", round(pseudo_I2_bio, 1), "%\n")
cat("Interpretation: ",
    ifelse(pseudo_I2_bio < 25, "Low heterogeneity",
           ifelse(pseudo_I2_bio < 75, "Moderate heterogeneity", "High heterogeneity")), "\n")

####################################################################################################
## JOINT MODEL: Density meta-analysis (Sensitivity / Variance Components) ##########################
####################################################################################################

# Step 1: Fit initial model with all density data
meta_density_full <- rma.mv(
  yi     = Mean,
  V      = vi,
  mods   = ~ Taxa - 1,
  random = list(~ 1 | MPA, ~ 1 | Source),
  data   = density_data,
  method = "REML",
  test   = "t"
)

cat("\n--- Density model (full data) ---\n")
print(summary(meta_density_full))

# Step 2: Cook's distance outlier detection
cooks_den <- cooks.distance(meta_density_full)
n_den <- nrow(density_data)
cooks_threshold_den <- 4 / n_den
outliers_den <- which(cooks_den > cooks_threshold_den)

cat("\nDensity Cook's distance threshold (4/n):", round(cooks_threshold_den, 4), "\n")
cat("Density outliers detected:", length(outliers_den), "\n")
if (length(outliers_den) > 0) {
  cat("Outlier rows (Taxa, MPA, Mean):\n")
  print(density_data[outliers_den, c("Taxa", "MPA", "Mean", "SE")])
}

# Step 3: Remove outliers and refit
if (length(outliers_den) > 0) {
  density_clean <- density_data[-outliers_den, ]
} else {
  density_clean <- density_data
}

meta_density <- rma.mv(
  yi     = Mean,
  V      = vi,
  mods   = ~ Taxa - 1,
  random = list(~ 1 | MPA, ~ 1 | Source),
  data   = density_clean,
  method = "REML",
  test   = "t"
)

cat("\n--- Density model (outliers removed) ---\n")
print(summary(meta_density))

# Report heterogeneity statistics for density model
cat("\n--- Density Heterogeneity Statistics ---\n")
cat("Between-MPA variance (tau²_MPA):", round(meta_density$sigma2[1], 4), "\n")
cat("Between-Source variance (tau²_Source):", round(meta_density$sigma2[2], 4), "\n")
# Calculate pseudo-I² for multilevel models
total_hetero_den <- sum(meta_density$sigma2)
# Use harmonic mean of sampling variances (consistent with biomass model)
typical_v_den <- 1 / mean(1 / density_clean$vi)
pseudo_I2_den <- 100 * total_hetero_den / (total_hetero_den + typical_v_den)
cat("Pseudo-I² (total):", round(pseudo_I2_den, 1), "%\n")
cat("Interpretation: ",
    ifelse(pseudo_I2_den < 25, "Low heterogeneity",
           ifelse(pseudo_I2_den < 75, "Moderate heterogeneity", "High heterogeneity")), "\n")

####################################################################################################
## SENSITIVITY ANALYSIS: Models with/without Source random effect #################################
####################################################################################################
# NOTE: The Source random effect accounts for program-level variation (PISCO, KFM, LTER).
# However, with only 3-4 Source levels (below the recommended 5-6 minimum), variance
# estimates may be uncertain. This sensitivity analysis compares results with and without
# the Source random effect to assess its impact on conclusions.

cat("\n")
cat("====================================\n")
cat("SENSITIVITY ANALYSIS: Source Effect\n")
cat("====================================\n")

# Fit alternative models WITHOUT Source random effect
meta_biomass_no_source <- tryCatch({
  rma.mv(
    yi     = Mean,
    V      = vi,
    mods   = ~ Taxa - 1,
    random = list(~ 1 | MPA),  # Only MPA, no Source
    data   = biomass_clean,
    method = "REML",
    test   = "t"
  )
}, error = function(e) {
  warning("Could not fit biomass model without Source: ", e$message)
  NULL
})

meta_density_no_source <- tryCatch({
  rma.mv(
    yi     = Mean,
    V      = vi,
    mods   = ~ Taxa - 1,
    random = list(~ 1 | MPA),  # Only MPA, no Source
    data   = density_clean,
    method = "REML",
    test   = "t"
  )
}, error = function(e) {
  warning("Could not fit density model without Source: ", e$message)
  NULL
})

# Model comparison (AIC/BIC)
cat("\n--- Model Comparison: With vs Without Source Random Effect ---\n")

if (!is.null(meta_biomass_no_source)) {
  comparison_biomass <- data.frame(
    Model = c("With Source", "Without Source"),
    AIC = c(AIC(meta_biomass), AIC(meta_biomass_no_source)),
    BIC = c(BIC(meta_biomass), BIC(meta_biomass_no_source)),
    logLik = c(as.numeric(logLik(meta_biomass)), as.numeric(logLik(meta_biomass_no_source))),
    tau2_MPA = c(meta_biomass$sigma2[1], meta_biomass_no_source$sigma2[1]),
    tau2_Source = c(meta_biomass$sigma2[2], NA)
  )
  comparison_biomass$deltaAIC <- comparison_biomass$AIC - min(comparison_biomass$AIC)
  cat("\nBIOMASS Model Comparison:\n")
  print(comparison_biomass)
  cat("  -> Preferred model (lower AIC):", comparison_biomass$Model[which.min(comparison_biomass$AIC)], "\n")
}

if (!is.null(meta_density_no_source)) {
  comparison_density <- data.frame(
    Model = c("With Source", "Without Source"),
    AIC = c(AIC(meta_density), AIC(meta_density_no_source)),
    BIC = c(BIC(meta_density), BIC(meta_density_no_source)),
    logLik = c(as.numeric(logLik(meta_density)), as.numeric(logLik(meta_density_no_source))),
    tau2_MPA = c(meta_density$sigma2[1], meta_density_no_source$sigma2[1]),
    tau2_Source = c(meta_density$sigma2[2], NA)
  )
  comparison_density$deltaAIC <- comparison_density$AIC - min(comparison_density$AIC)
  cat("\nDENSITY Model Comparison:\n")
  print(comparison_density)
  cat("  -> Preferred model (lower AIC):", comparison_density$Model[which.min(comparison_density$AIC)], "\n")
}

# Compare coefficient estimates
cat("\n--- Effect of Source on Coefficient Estimates ---\n")
if (!is.null(meta_biomass_no_source)) {
  coef_with <- coef(summary(meta_biomass))[, c("estimate", "se", "pval")]
  coef_without <- coef(summary(meta_biomass_no_source))[, c("estimate", "se", "pval")]
  coef_diff_bio <- data.frame(
    Taxa = gsub("^Taxa", "", rownames(coef_with)),
    Est_with = round(coef_with$estimate, 3),
    Est_without = round(coef_without$estimate, 3),
    SE_with = round(coef_with$se, 3),
    SE_without = round(coef_without$se, 3),
    pval_with = round(coef_with$pval, 4),
    pval_without = round(coef_without$pval, 4)
  )
  cat("\nBIOMASS - Coefficient comparison:\n")
  print(coef_diff_bio)
}

if (!is.null(meta_density_no_source)) {
  coef_with <- coef(summary(meta_density))[, c("estimate", "se", "pval")]
  coef_without <- coef(summary(meta_density_no_source))[, c("estimate", "se", "pval")]
  coef_diff_den <- data.frame(
    Taxa = gsub("^Taxa", "", rownames(coef_with)),
    Est_with = round(coef_with$estimate, 3),
    Est_without = round(coef_without$estimate, 3),
    SE_with = round(coef_with$se, 3),
    SE_without = round(coef_without$se, 3),
    pval_with = round(coef_with$pval, 4),
    pval_without = round(coef_without$pval, 4)
  )
  cat("\nDENSITY - Coefficient comparison:\n")
  print(coef_diff_den)
}

####################################################################################################
## Variance Component Confidence Intervals ########################################################
####################################################################################################
# NOTE: With only 3-4 Source levels, variance component estimates have high uncertainty.
# Reporting confidence intervals provides transparency about this limitation.

cat("\n")
cat("==========================================\n")
cat("VARIANCE COMPONENT CONFIDENCE INTERVALS\n")
cat("==========================================\n")
cat("(Using profile likelihood method from metafor::confint)\n")

# Extract confidence intervals for variance components
ci_biomass <- tryCatch({
  confint(meta_biomass)
}, error = function(e) {
  warning("Could not compute variance component CIs for biomass: ", e$message)
  NULL
})

ci_density <- tryCatch({
  confint(meta_density)
}, error = function(e) {
  warning("Could not compute variance component CIs for density: ", e$message)
  NULL
})

if (!is.null(ci_biomass)) {
  cat("\n--- BIOMASS Variance Component CIs ---\n")
  print(ci_biomass)
}

if (!is.null(ci_density)) {
  cat("\n--- DENSITY Variance Component CIs ---\n")
  print(ci_density)
}

# Create summary table for supplementary materials
tau2_summary <- data.frame(
  Response = character(),
  Component = character(),
  Estimate = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  stringsAsFactors = FALSE
)

if (!is.null(ci_biomass)) {
  # Extract from confint output (structure depends on metafor version)
  if ("random" %in% names(ci_biomass)) {
    for (i in seq_along(ci_biomass$random)) {
      comp_name <- names(ci_biomass$random)[i]
      tau2_summary <- rbind(tau2_summary, data.frame(
        Response = "Biomass",
        Component = comp_name,
        Estimate = ci_biomass$random[[i]]["estimate"],
        CI_Lower = ci_biomass$random[[i]]["ci.lb"],
        CI_Upper = ci_biomass$random[[i]]["ci.ub"],
        stringsAsFactors = FALSE
      ))
    }
  }
}

if (!is.null(ci_density)) {
  if ("random" %in% names(ci_density)) {
    for (i in seq_along(ci_density$random)) {
      comp_name <- names(ci_density$random)[i]
      tau2_summary <- rbind(tau2_summary, data.frame(
        Response = "Density",
        Component = comp_name,
        Estimate = ci_density$random[[i]]["estimate"],
        CI_Lower = ci_density$random[[i]]["ci.lb"],
        CI_Upper = ci_density$random[[i]]["ci.ub"],
        stringsAsFactors = FALSE
      ))
    }
  }
}

if (nrow(tau2_summary) > 0) {
  cat("\n--- Summary Table: Variance Components with CIs ---\n")
  print(tau2_summary)

  # Export to supplementary materials
  write.csv(tau2_summary, here::here("data", "table_s_variance_components.csv"), row.names = FALSE)
  cat("\nVariance component summary exported to: data/table_s_variance_components.csv\n")
}

cat("\n")
cat("NOTE: The Source random effect has only 3-4 levels (PISCO, KFM, LTER, Landsat),\n")
cat("which is below the recommended minimum of 5-6 for reliable variance estimation.\n")
cat("Interpret tau²_Source estimates with caution; wide CIs are expected.\n")

####################################################################################################
## PER-TAXA META-ANALYSIS (Primary Analysis) ######################################################
####################################################################################################
#
# Why per-taxa instead of joint model?
# The joint model (rma.mv with Taxa moderator) applies Cook's distance globally across
# all taxa. An observation can be flagged as an "outlier" because it differs from OTHER
# taxa's effect sizes, not because it's unusual within its OWN taxon. With the joint
# model, the Cook's D threshold of 4/n (where n = all observations) is extremely aggressive
# when taxa have naturally different effect magnitudes.
#
# Per-taxa models apply Cook's distance within each taxon separately, which is the
# correct question: "Is this observation influential for THIS taxon's estimate?"
#
# For taxa with k >= 5 and >= 3 MPAs: rma.mv() with MPA random effect
# For taxa with 2 <= k < 5: rma() simple random-effects model
# For taxa with k < 2: report descriptively (no model)
#
# The joint model results are retained above as a sensitivity comparison.
#

#' Fit per-taxa meta-analysis with within-taxon outlier detection
#'
#' @param data Data frame with columns: Taxa, MPA, Source, Mean, vi
#' @param response_label Character: "Biomass" or "Density"
#' @return List with: table (summary df), clean_data (filtered data), audit (all data with flags)
run_per_taxa_meta <- function(data, response_label) {
  taxa_list <- sort(unique(as.character(data$Taxa)))
  results_list <- list()
  clean_data_list <- list()
  audit_list <- list()

  for (taxon in taxa_list) {
    sub <- data[data$Taxa == taxon, ]
    k_input <- nrow(sub)

    cat(sprintf("\n  %s (%s): k_input = %d\n", taxon, response_label, k_input))

    # --- Handle k < 2: report descriptively ---
    if (k_input < 2) {
      cat(sprintf("    -> k < 2: reporting descriptively\n"))
      mean_val <- as.numeric(sub$Mean[1])
      se_val <- as.numeric(sub$SE[1])
      z <- mean_val / se_val
      results_list[[taxon]] <- data.frame(
        Taxa = taxon, k = k_input,
        Estimate = round(mean_val, 4), SE = round(se_val, 4),
        tval = round(z, 4),
        pval = 2 * pnorm(-abs(z)),
        CI_lower = round(mean_val - 1.96 * se_val, 4),
        CI_upper = round(mean_val + 1.96 * se_val, 4),
        Response = response_label, stringsAsFactors = FALSE
      )
      sub$Is_Outlier <- FALSE
      sub$Cooks_Distance <- NA_real_
      sub$Cooks_Threshold <- NA_real_
      sub$In_Final_Analysis <- TRUE
      sub$Exclusion_Reason <- "INCLUDED"
      audit_list[[taxon]] <- sub
      clean_data_list[[taxon]] <- sub
      next
    }

    # --- Fit initial model ---
    model_full <- tryCatch({
      if (k_input >= 5 && length(unique(sub$MPA)) >= 3) {
        rma.mv(yi = Mean, V = vi, random = list(~ 1 | MPA),
               data = sub, method = "REML", test = "t")
      } else {
        rma(yi = Mean, vi = vi, data = sub, method = "REML", test = "t")
      }
    }, error = function(e) {
      warning(sprintf("Model failed for %s %s: %s. Falling back to FE.",
                      taxon, response_label, e$message))
      rma(yi = Mean, vi = vi, data = sub, method = "FE")
    })

    # --- Per-taxa Cook's distance ---
    cooks_vals <- tryCatch(cooks.distance(model_full), error = function(e) rep(0, k_input))
    threshold <- 4 / k_input
    outlier_idx <- which(cooks_vals > threshold)

    cat(sprintf("    -> Cook's D threshold (4/%d): %.4f, outliers: %d\n",
                k_input, threshold, length(outlier_idx)))

    # Track audit info
    sub$Cooks_Distance <- cooks_vals
    sub$Cooks_Threshold <- threshold

    # --- Remove outliers and refit (if any, and if enough data remains) ---
    if (length(outlier_idx) > 0 && (k_input - length(outlier_idx)) >= 2) {
      sub_clean <- sub[-outlier_idx, ]
      k_clean <- nrow(sub_clean)

      model <- tryCatch({
        if (k_clean >= 5 && length(unique(sub_clean$MPA)) >= 3) {
          rma.mv(yi = Mean, V = vi, random = list(~ 1 | MPA),
                 data = sub_clean, method = "REML", test = "t")
        } else {
          rma(yi = Mean, vi = vi, data = sub_clean, method = "REML", test = "t")
        }
      }, error = function(e) {
        warning(sprintf("Clean model failed for %s %s: %s. Using full model.",
                        taxon, response_label, e$message))
        model_full
      })

      sub$Is_Outlier <- seq_len(nrow(sub)) %in% outlier_idx
    } else {
      sub_clean <- sub
      model <- model_full
      if (length(outlier_idx) > 0) {
        cat(sprintf("    -> Cannot remove outliers: would leave k < 2\n"))
      }
      sub$Is_Outlier <- FALSE
    }

    sub$In_Final_Analysis <- !sub$Is_Outlier
    sub$Exclusion_Reason <- ifelse(
      sub$In_Final_Analysis, "INCLUDED",
      paste0("Outlier (Cook's D=", round(sub$Cooks_Distance, 4),
             " > ", round(sub$Cooks_Threshold, 4), ")")
    )

    audit_list[[taxon]] <- sub
    clean_data_list[[taxon]] <- sub[sub$In_Final_Analysis, ]

    # --- Extract results ---
    coef_tbl <- coef(summary(model))
    stat_col <- if ("tval" %in% colnames(coef_tbl)) "tval" else "zval"

    results_list[[taxon]] <- data.frame(
      Taxa = taxon, k = model$k,
      Estimate = round(coef_tbl[1, "estimate"], 4),
      SE = round(coef_tbl[1, "se"], 4),
      tval = round(coef_tbl[1, stat_col], 4),
      pval = coef_tbl[1, "pval"],
      CI_lower = round(coef_tbl[1, "ci.lb"], 4),
      CI_upper = round(coef_tbl[1, "ci.ub"], 4),
      Response = response_label, stringsAsFactors = FALSE
    )

    cat(sprintf("    -> Final k=%d, estimate=%.3f, SE=%.3f, p=%.4f\n",
                model$k, coef_tbl[1, "estimate"], coef_tbl[1, "se"], coef_tbl[1, "pval"]))
  }

  table_out <- do.call(rbind, results_list)
  rownames(table_out) <- NULL
  clean_out <- do.call(rbind, clean_data_list)
  audit_out <- do.call(rbind, audit_list)

  list(table = table_out, clean_data = clean_out, audit = audit_out)
}

# --- Run per-taxa meta-analysis ---
cat("\n")
cat("============================================\n")
cat("PER-TAXA META-ANALYSIS (Primary Analysis)\n")
cat("============================================\n")

pertaxa_biomass <- run_per_taxa_meta(biomass_data, "Biomass")
pertaxa_density <- run_per_taxa_meta(density_data, "Density")

# Summary comparison: per-taxa vs joint outlier removal
cat("\n--- Outlier Removal Comparison ---\n")
pertaxa_bio_removed <- sum(!pertaxa_biomass$audit$In_Final_Analysis)
pertaxa_den_removed <- sum(!pertaxa_density$audit$In_Final_Analysis)
joint_bio_removed <- length(outliers_bio)
joint_den_removed <- length(outliers_den)
cat(sprintf("  Biomass: Joint model removed %d/%d (%.0f%%), Per-taxa removed %d/%d (%.0f%%)\n",
            joint_bio_removed, n_bio, 100 * joint_bio_removed / n_bio,
            pertaxa_bio_removed, n_bio, 100 * pertaxa_bio_removed / n_bio))
cat(sprintf("  Density: Joint model removed %d/%d (%.0f%%), Per-taxa removed %d/%d (%.0f%%)\n",
            joint_den_removed, n_den, 100 * joint_den_removed / n_den,
            pertaxa_den_removed, n_den, 100 * pertaxa_den_removed / n_den))
cat(sprintf("  TOTAL:   Joint removed %d/%d (%.0f%%), Per-taxa removed %d/%d (%.0f%%)\n",
            joint_bio_removed + joint_den_removed, n_bio + n_den,
            100 * (joint_bio_removed + joint_den_removed) / (n_bio + n_den),
            pertaxa_bio_removed + pertaxa_den_removed, n_bio + n_den,
            100 * (pertaxa_bio_removed + pertaxa_den_removed) / (n_bio + n_den)))

# Export per-taxa audit
pertaxa_audit <- rbind(pertaxa_biomass$audit, pertaxa_density$audit)
write.csv(pertaxa_audit, here::here("outputs", "filter_audit_pertaxa_meta.csv"), row.names = FALSE)
cat("Per-taxa audit saved to: outputs/filter_audit_pertaxa_meta.csv\n")

# --- Outlier sensitivity table: 3 approaches side by side ---
cat("\n")
cat("============================================\n")
cat("OUTLIER SENSITIVITY TABLE\n")
cat("============================================\n")

# Approach 1: No outlier removal (full data)
extract_meta_table_helper <- function(model, data, response_label) {
  coef_tbl <- coef(summary(model))
  stat_col <- if ("tval" %in% colnames(coef_tbl)) "tval" else "zval"
  taxa_names <- gsub("^Taxa", "", rownames(coef_tbl))
  k_per_taxa <- vapply(taxa_names, function(t) sum(data$Taxa == t), integer(1))
  data.frame(
    Taxa = taxa_names, k = k_per_taxa,
    Estimate = round(coef_tbl[, "estimate"], 4),
    SE = round(coef_tbl[, "se"], 4),
    pval = coef_tbl[, "pval"],
    Response = response_label, row.names = NULL, stringsAsFactors = FALSE
  )
}

full_bio_tbl <- extract_meta_table_helper(meta_biomass_full, biomass_data, "Biomass")
full_den_tbl <- extract_meta_table_helper(meta_density_full, density_data, "Density")
full_tbl <- rbind(full_bio_tbl, full_den_tbl)
full_tbl$Method <- "No outlier removal"

# Approach 2: Per-taxa Cook's D (primary)
pertaxa_tbl <- rbind(pertaxa_biomass$table, pertaxa_density$table)
pertaxa_tbl_sens <- pertaxa_tbl[, c("Taxa", "k", "Estimate", "SE", "pval", "Response")]
pertaxa_tbl_sens$Method <- "Per-taxa Cook's D (4/k)"

# Approach 3: Joint Cook's D (current/legacy)
joint_bio_tbl <- extract_meta_table_helper(meta_biomass, biomass_clean, "Biomass")
joint_den_tbl <- extract_meta_table_helper(meta_density, density_clean, "Density")
joint_tbl <- rbind(joint_bio_tbl, joint_den_tbl)
joint_tbl$Method <- "Joint Cook's D (4/n)"

sensitivity_table <- rbind(full_tbl, pertaxa_tbl_sens, joint_tbl)
sensitivity_table <- sensitivity_table[order(sensitivity_table$Response,
                                              sensitivity_table$Taxa,
                                              sensitivity_table$Method), ]

write.csv(sensitivity_table, here::here("data", "table_s_outlier_sensitivity.csv"), row.names = FALSE)
cat("Outlier sensitivity table saved to: data/table_s_outlier_sensitivity.csv\n")
print(sensitivity_table)

# --- Leave-one-out analysis for kelp biomass ---
cat("\n")
cat("============================================\n")
cat("LEAVE-ONE-OUT: M. pyrifera biomass\n")
cat("============================================\n")

kelp_clean <- pertaxa_biomass$clean_data[pertaxa_biomass$clean_data$Taxa == "M. pyrifera", ]
if (nrow(kelp_clean) >= 3) {
  kelp_model <- rma(yi = Mean, vi = vi, data = kelp_clean, method = "REML")
  loo_kelp <- leave1out(kelp_model)

  cat(sprintf("Full model: estimate = %.4f, p = %.4f, k = %d\n",
              kelp_model$beta[1], kelp_model$pval[1], kelp_model$k))
  cat("\nLeave-one-out results:\n")
  loo_df <- data.frame(
    Removed_MPA = kelp_clean$MPA,
    Removed_Source = kelp_clean$Source,
    Removed_ES = round(as.numeric(kelp_clean$Mean), 3),
    Estimate = round(loo_kelp$estimate, 4),
    SE = round(loo_kelp$se, 4),
    pval = round(loo_kelp$pval, 4)
  )
  loo_df$Significant <- ifelse(loo_df$pval < 0.05, "YES", "no")
  print(loo_df)

  n_sig <- sum(loo_df$pval < 0.05)
  cat(sprintf("\nKelp is significant (p < 0.05) in %d of %d leave-one-out permutations (%.0f%%)\n",
              n_sig, nrow(loo_df), 100 * n_sig / nrow(loo_df)))

  write.csv(loo_df, here::here("data", "table_s_kelp_leave1out.csv"), row.names = FALSE)
  cat("Leave-one-out results saved to: data/table_s_kelp_leave1out.csv\n")
} else {
  cat("Insufficient kelp data for leave-one-out analysis\n")
}

####################################################################################################
## Build Table 2: Summary of meta-analysis results ################################################
####################################################################################################

#' Extract a tidy summary table from an rma.mv model
#'
#' Returns a dataframe with one row per taxon containing the estimate, SE,
#' t-value, p-value, 95% confidence interval bounds, and sample size (k).
#'
#' @param model An rma.mv model object fit with mods = ~ Taxa - 1
#' @param data The data used to fit the model (for calculating sample sizes)
#' @param response Character label for the response type ("Biomass" or "Density")
#' @return A dataframe with columns: Taxa, k, Estimate, SE, tval, pval, CI_lower, CI_upper, Response
extract_meta_table <- function(model, data, response) {
  coef_table <- coef(summary(model))
  # Column name is "zval" in rma.mv (not "tval" unless test="t" is used)
  stat_col <- if ("tval" %in% colnames(coef_table)) "tval" else "zval"

  # Calculate sample size (k = number of effect sizes) per taxon
  taxa_names <- gsub("^Taxa", "", rownames(coef_table))
  k_per_taxa <- vapply(taxa_names, function(taxon) {
    sum(data$Taxa == taxon, na.rm = TRUE)
  }, integer(1))

  data.frame(
    Taxa      = taxa_names,
    k         = k_per_taxa,  # Number of effect sizes
    Estimate  = round(coef_table[, "estimate"], 4),
    SE        = round(coef_table[, "se"], 4),
    tval      = round(coef_table[, stat_col], 4),
    pval      = coef_table[, "pval"],  # Full precision for FDR correction
    CI_lower  = round(coef_table[, "ci.lb"], 4),
    CI_upper  = round(coef_table[, "ci.ub"], 4),
    Response  = response,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

# Build Table 2 from PER-TAXA models (primary analysis)
# Joint model results are retained as Table2_joint for sensitivity comparison
Table2_joint_biomass <- extract_meta_table(meta_biomass, biomass_clean, "Biomass")
Table2_joint_density <- extract_meta_table(meta_density, density_clean, "Density")
Table2_joint <- rbind(Table2_joint_biomass, Table2_joint_density)

# Primary Table 2: per-taxa meta-analysis results
Table2 <- rbind(pertaxa_biomass$table, pertaxa_density$table)

# Use per-taxa clean data for sample size reporting
biomass_clean_pertaxa <- pertaxa_biomass$clean_data
density_clean_pertaxa <- pertaxa_density$clean_data

# Report sample size summary
cat("\n--- Sample Size Summary (Per-Taxa Analysis) ---\n")
cat("Total biomass effect sizes (after per-taxa outlier removal):", nrow(biomass_clean_pertaxa), "\n")
cat("  - Number of unique MPAs:", length(unique(biomass_clean_pertaxa$MPA)), "\n")
cat("  - Number of unique data sources:", length(unique(biomass_clean_pertaxa$Source)), "\n")
if ("N" %in% names(biomass_clean_pertaxa)) {
  cat("  - Total underlying observations (N):", sum(as.numeric(biomass_clean_pertaxa$N), na.rm = TRUE), "\n")
}
cat("Total density effect sizes (after per-taxa outlier removal):", nrow(density_clean_pertaxa), "\n")
cat("  - Number of unique MPAs:", length(unique(density_clean_pertaxa$MPA)), "\n")
cat("  - Number of unique data sources:", length(unique(density_clean_pertaxa$Source)), "\n")
if ("N" %in% names(density_clean_pertaxa)) {
  cat("  - Total underlying observations (N):", sum(as.numeric(density_clean_pertaxa$N), na.rm = TRUE), "\n")
}

# Order taxa to match manuscript Table 2 presentation
taxa_order <- c("S. purpuratus", "M. franciscanus", "M. pyrifera",
                "P. interruptus", "S. pulcher")
Table2$Taxa <- factor(Table2$Taxa, levels = taxa_order)
Table2 <- Table2[order(Table2$Response, Table2$Taxa), ]
rownames(Table2) <- NULL

# Flag taxa with very few effect sizes (k < 5)
Table2$Flag <- ifelse(Table2$k < 5, "preliminary (k<5)", "")

low_k <- Table2[Table2$k < 5, ]
if (nrow(low_k) > 0) {
  cat("\nWARNING: The following taxa have k < 5 effect sizes (interpret with caution):\n")
  print(low_k[, c("Taxa", "Response", "k", "Estimate", "pval")])
}
# Additional warning for k < 3: cannot meaningfully estimate heterogeneity
very_low_k <- Table2[Table2$k < 3, ]
if (nrow(very_low_k) > 0) {
  cat("\nCAUTION: Taxa with k < 3 have insufficient replication for reliable meta-analysis.\n")
  cat("  These should be reported descriptively, not as formal meta-analytic estimates.\n")
  cat("  Affected: ", paste(paste0(very_low_k$Taxa, " (", very_low_k$Response, ", k=", very_low_k$k, ")"),
                            collapse = "; "), "\n")
}

# Back-transform to Response Ratio scale for interpretability
# RR = exp(lnRR); RR = 1 means no difference, RR > 1 = higher inside MPA
Table2$RR        <- round(exp(Table2$Estimate), 3)
Table2$RR_lower  <- round(exp(Table2$CI_lower), 3)
Table2$RR_upper  <- round(exp(Table2$CI_upper), 3)
Table2$Pct_Change <- round((exp(Table2$Estimate) - 1) * 100, 1)

# FDR-corrected p-values for multiple testing across all taxa-response tests
Table2$pval_fdr <- p.adjust(Table2$pval, method = "fdr")

cat("\n--- Multiple Testing Correction (FDR) ---\n")
cat("Tests conducted:", nrow(Table2), "\n")
sig_uncorrected <- sum(Table2$pval < 0.05)
sig_fdr <- sum(Table2$pval_fdr < 0.05)
cat("Significant at p<0.05 (uncorrected):", sig_uncorrected, "\n")
cat("Significant at p<0.05 (FDR-corrected):", sig_fdr, "\n")

cat("\n============================\n")
cat("TABLE 2: Meta-analysis results\n")
cat("============================\n")
print(Table2)

# --- OUTPUT VALIDATION: Ensure RR columns are present before export ---
required_cols <- c("Taxa", "k", "Estimate", "SE", "tval", "pval",
                   "CI_lower", "CI_upper", "Response", "Flag",
                   "RR", "RR_lower", "RR_upper", "Pct_Change", "pval_fdr")
missing_cols <- setdiff(required_cols, names(Table2))
if (length(missing_cols) > 0) {
  stop("Table2 is missing expected columns before export: ",
       paste(missing_cols, collapse = ", "))
}

# Export Table 2 as CSV to data directory (with explicit column order)
write.csv(Table2[, required_cols], here::here("data", "table_02_meta_analysis.csv"), row.names = FALSE)
cat("\nTable 2 exported to:", here::here("data", "table_02_meta_analysis.csv"), "\n")
cat("  Columns:", paste(required_cols, collapse = ", "), "\n")

####################################################################################################
## META-ANALYSIS FILTERING AUDIT: Track exactly what enters meta-analysis #########################
####################################################################################################

cat("\n")
cat("==========================================\n")
cat("META-ANALYSIS FILTERING AUDIT\n")
cat("==========================================\n")

# Create comprehensive audit of what enters meta-analysis
# Start with all data entering this script (SumStats.Final)

# Track biomass data through filtering
cat("\n--- BIOMASS Meta-Analysis Filtering ---\n")
cat("Input to meta-analysis (SumStats.Final, Resp='Bio'):", nrow(biomass_data), "\n")

# Create biomass audit
BiomassAudit <- biomass_data[, c("Taxa", "MPA", "Source", "Mean", "SE", "vi")]
BiomassAudit$Input_Order <- seq_len(nrow(BiomassAudit))

# Track outlier status
if (length(outliers_bio) > 0) {
  BiomassAudit$Is_Outlier <- seq_len(nrow(biomass_data)) %in% outliers_bio
  BiomassAudit$Cooks_Distance <- cooks_bio
} else {
  BiomassAudit$Is_Outlier <- FALSE
  BiomassAudit$Cooks_Distance <- cooks_bio
}
BiomassAudit$Cooks_Threshold <- cooks_threshold_bio
BiomassAudit$In_Final_Analysis <- !BiomassAudit$Is_Outlier

cat("Outliers removed:", sum(BiomassAudit$Is_Outlier), "\n")
cat("Final biomass observations:", sum(BiomassAudit$In_Final_Analysis), "\n")

# k-values by taxa
cat("\nBiomass k-values (observations per taxa):\n")
bio_k <- aggregate(In_Final_Analysis ~ Taxa, data = BiomassAudit,
                   FUN = function(x) c(input = length(x), final = sum(x)))
for (i in seq_len(nrow(bio_k))) {
  cat(sprintf("  %s: %d input -> %d final (k=%d)\n",
              bio_k$Taxa[i], bio_k$In_Final_Analysis[i, "input"],
              bio_k$In_Final_Analysis[i, "final"],
              bio_k$In_Final_Analysis[i, "final"]))
}

# Track density data through filtering
cat("\n--- DENSITY Meta-Analysis Filtering ---\n")
cat("Input to meta-analysis (SumStats.Final, Resp='Den'):", nrow(density_data), "\n")

# Create density audit
DensityAudit <- density_data[, c("Taxa", "MPA", "Source", "Mean", "SE", "vi")]
DensityAudit$Input_Order <- seq_len(nrow(DensityAudit))

# Track outlier status
if (length(outliers_den) > 0) {
  DensityAudit$Is_Outlier <- seq_len(nrow(density_data)) %in% outliers_den
  DensityAudit$Cooks_Distance <- cooks_den
} else {
  DensityAudit$Is_Outlier <- FALSE
  DensityAudit$Cooks_Distance <- cooks_den
}
DensityAudit$Cooks_Threshold <- cooks_threshold_den
DensityAudit$In_Final_Analysis <- !DensityAudit$Is_Outlier

cat("Outliers removed:", sum(DensityAudit$Is_Outlier), "\n")
cat("Final density observations:", sum(DensityAudit$In_Final_Analysis), "\n")

# k-values by taxa
cat("\nDensity k-values (observations per taxa):\n")
den_k <- aggregate(In_Final_Analysis ~ Taxa, data = DensityAudit,
                   FUN = function(x) c(input = length(x), final = sum(x)))
for (i in seq_len(nrow(den_k))) {
  cat(sprintf("  %s: %d input -> %d final (k=%d)\n",
              den_k$Taxa[i], den_k$In_Final_Analysis[i, "input"],
              den_k$In_Final_Analysis[i, "final"],
              den_k$In_Final_Analysis[i, "final"]))
}

# Combine and write audit
BiomassAudit$Response <- "Biomass"
DensityAudit$Response <- "Density"
MetaAnalysisAudit <- rbind(BiomassAudit, DensityAudit)

# Add exclusion reason
MetaAnalysisAudit$Exclusion_Reason <- ifelse(MetaAnalysisAudit$In_Final_Analysis,
                                              "INCLUDED",
                                              paste0("Outlier (Cook's D=",
                                                     round(MetaAnalysisAudit$Cooks_Distance, 4),
                                                     " > ", round(MetaAnalysisAudit$Cooks_Threshold, 4), ")"))

# Write to CSV
meta_audit_file <- here::here("outputs", "filter_audit_meta_analysis.csv")
write.csv(MetaAnalysisAudit, meta_audit_file, row.names = FALSE)
cat("\nMeta-analysis filter audit saved to:", meta_audit_file, "\n")

# LOBSTER SPECIFIC DETAIL
cat("\n--- LOBSTER (P. interruptus) META-ANALYSIS DETAIL ---\n")
lob_bio <- subset(BiomassAudit, Taxa == "P. interruptus")
lob_den <- subset(DensityAudit, Taxa == "P. interruptus")

cat("\nLobster Biomass:\n")
cat(sprintf("  Total in meta-analysis input: %d\n", nrow(lob_bio)))
cat(sprintf("  Outliers removed: %d\n", sum(lob_bio$Is_Outlier)))
cat(sprintf("  Final k: %d\n", sum(lob_bio$In_Final_Analysis)))
if (nrow(lob_bio) > 0) {
  cat("  Detail:\n")
  for (i in seq_len(nrow(lob_bio))) {
    status <- ifelse(lob_bio$In_Final_Analysis[i], "+", "-")
    cat(sprintf("    %s %s (%s): Effect=%.3f, SE=%.3f, Cook's=%.4f %s\n",
                status, lob_bio$MPA[i], lob_bio$Source[i],
                as.numeric(lob_bio$Mean[i]), as.numeric(lob_bio$SE[i]),
                lob_bio$Cooks_Distance[i],
                ifelse(lob_bio$Is_Outlier[i], "[OUTLIER]", "")))
  }
}

cat("\nLobster Density:\n")
cat(sprintf("  Total in meta-analysis input: %d\n", nrow(lob_den)))
cat(sprintf("  Outliers removed: %d\n", sum(lob_den$Is_Outlier)))
cat(sprintf("  Final k: %d\n", sum(lob_den$In_Final_Analysis)))
if (nrow(lob_den) > 0) {
  cat("  Detail:\n")
  for (i in seq_len(nrow(lob_den))) {
    status <- ifelse(lob_den$In_Final_Analysis[i], "+", "-")
    cat(sprintf("    %s %s (%s): Effect=%.3f, SE=%.3f, Cook's=%.4f %s\n",
                status, lob_den$MPA[i], lob_den$Source[i],
                as.numeric(lob_den$Mean[i]), as.numeric(lob_den$SE[i]),
                lob_den$Cooks_Distance[i],
                ifelse(lob_den$Is_Outlier[i], "[OUTLIER]", "")))
  }
}

# Create combined summary across both filtering stages
cat("\n--- COMPLETE DATA FLOW SUMMARY ---\n")
cat("This traces data from effect size calculation through meta-analysis.\n")
cat("For full detail, see:\n")
cat("  1. outputs/filter_audit_effect_sizes.csv (08_effect_sizes.R filtering)\n")
cat("  2. outputs/filter_audit_meta_analysis.csv (09_meta_analysis.R outlier removal)\n")
cat("  3. outputs/filter_summary_by_taxa.csv (taxa-level summary)\n")

####################################################################################################
## Table 3: Meta-regressions between taxa (trophic cascade relationships) ##########################
####################################################################################################

# Table 3 tests the trophic cascade relationship: whether urchin density effect sizes
# predict kelp biomass effect sizes across MPAs.
# Uses metafor::rma() meta-regression instead of OLS because:
#   1. Both X and Y are estimated effect sizes with measurement error
#   2. OLS biases slopes toward zero when the predictor has error (attenuation bias)
#   3. rma() properly weights by the response's sampling variance (vi = SE^2)

# Merge effect sizes by MPA to test cross-taxa relationships
# We need MPA-level estimates for different taxa in the same row

# Get primary effect sizes per MPA for each taxa x response combination
# Pivot both Mean and SE so we can use sampling variance in meta-regression
es_wide_mean <- SumStats.Final %>%
  dplyr::select(Taxa, MPA, Mean, Resp) %>%
  dplyr::mutate(Mean = as.numeric(Mean)) %>%
  tidyr::unite("Taxa_Resp", Taxa, Resp, sep = "_") %>%
  tidyr::pivot_wider(names_from = Taxa_Resp, values_from = Mean, values_fn = mean)

es_wide_se <- SumStats.Final %>%
  dplyr::select(Taxa, MPA, SE, Resp) %>%
  dplyr::mutate(SE = as.numeric(SE)) %>%
  tidyr::unite("Taxa_Resp", Taxa, Resp, sep = "_SE_") %>%
  tidyr::pivot_wider(names_from = Taxa_Resp, values_from = SE, values_fn = mean)

es_wide <- dplyr::left_join(es_wide_mean, es_wide_se, by = "MPA")

# Clean column names for formula use (replace spaces and dots)
names(es_wide) <- gsub(" ", "_", names(es_wide))
names(es_wide) <- gsub("\\.", "_", names(es_wide))

cat("\n============================\n")
cat("TABLE 3: Cross-taxa meta-regressions (metafor::rma)\n")
cat("============================\n")

#' Helper to print rma meta-regression results in a compact, readable format
#' @param model An rma model object
#' @param label Character string describing the model
print_meta_reg <- function(model, label) {
  cat("\n---", label, "---\n")
  cat(sprintf("k = %d effect sizes\n", model$k))
  # Moderator (slope) coefficient
  coef_tbl <- coef(summary(model))
  # Row 1 = intercept, Row 2 = moderator slope
  cat(sprintf("Intercept:  estimate = %.4f, SE = %.4f, p = %.4f\n",
              coef_tbl[1, "estimate"], coef_tbl[1, "se"], coef_tbl[1, "pval"]))
  cat(sprintf("Slope:      estimate = %.4f, SE = %.4f, p = %.4f\n",
              coef_tbl[2, "estimate"], coef_tbl[2, "se"], coef_tbl[2, "pval"]))
  cat(sprintf("  95%% CI for slope: [%.4f, %.4f]\n",
              coef_tbl[2, "ci.lb"], coef_tbl[2, "ci.ub"]))
  # Model fit statistics
  cat(sprintf("Test of moderator: QM(df=%d) = %.4f, p = %.4f\n",
              model$m, model$QM, model$QMp))
  cat(sprintf("Residual heterogeneity: QE(df=%d) = %.4f, p = %.4f\n",
              model$k - model$p, model$QE, model$QEp))
  cat(sprintf("tau² = %.4f, I² = %.1f%%\n", model$tau2, model$I2))
  cat(sprintf("R² (variance explained by moderator) = %.1f%%\n",
              ifelse(is.null(model$R2), NA, model$R2)))
  cat("\n")
}

# Model 1: S. purpuratus density predicts M. pyrifera biomass
# (urchin grazing hypothesis: more urchins -> less kelp)
if (all(c("S_purpuratus_Den", "M_pyrifera_Bio") %in% names(es_wide))) {
  w_col <- "M_pyrifera_SE_Bio"
  req_cols <- c("S_purpuratus_Den", "M_pyrifera_Bio", w_col)
  if (w_col %in% names(es_wide)) {
    es_purp_macro <- es_wide[complete.cases(es_wide[, req_cols]), ]
    if (nrow(es_purp_macro) >= 3) {
      meta_purp_macro <- rma(yi = es_purp_macro$M_pyrifera_Bio,
                             vi = as.numeric(es_purp_macro[[w_col]])^2,
                             mods = ~ S_purpuratus_Den,
                             data = es_purp_macro, method = "REML")
      print_meta_reg(meta_purp_macro, "S. purpuratus density -> M. pyrifera biomass")
    }
  }
}

# Model 2: M. franciscanus density predicts M. pyrifera biomass
if (all(c("M_franciscanus_Den", "M_pyrifera_Bio") %in% names(es_wide))) {
  w_col <- "M_pyrifera_SE_Bio"
  req_cols <- c("M_franciscanus_Den", "M_pyrifera_Bio", w_col)
  if (w_col %in% names(es_wide)) {
    es_reds_macro <- es_wide[complete.cases(es_wide[, req_cols]), ]
    if (nrow(es_reds_macro) >= 3) {
      meta_reds_macro <- rma(yi = es_reds_macro$M_pyrifera_Bio,
                             vi = as.numeric(es_reds_macro[[w_col]])^2,
                             mods = ~ M_franciscanus_Den,
                             data = es_reds_macro, method = "REML")
      print_meta_reg(meta_reds_macro, "M. franciscanus density -> M. pyrifera biomass")
    }
  }
}

# Model 3: P. interruptus density predicts S. purpuratus density
# (predator-prey: more lobsters -> fewer urchins)
if (all(c("P_interruptus_Den", "S_purpuratus_Den") %in% names(es_wide))) {
  w_col <- "S_purpuratus_SE_Den"
  req_cols <- c("P_interruptus_Den", "S_purpuratus_Den", w_col)
  if (w_col %in% names(es_wide)) {
    es_lob_purp <- es_wide[complete.cases(es_wide[, req_cols]), ]
    if (nrow(es_lob_purp) >= 3) {
      meta_lob_purp <- rma(yi = es_lob_purp$S_purpuratus_Den,
                           vi = as.numeric(es_lob_purp[[w_col]])^2,
                           mods = ~ P_interruptus_Den,
                           data = es_lob_purp, method = "REML")
      print_meta_reg(meta_lob_purp, "P. interruptus density -> S. purpuratus density")
    }
  }
}

# Model 4: S. pulcher density predicts S. purpuratus density
# (predator-prey: more sheephead -> fewer urchins)
if (all(c("S_pulcher_Den", "S_purpuratus_Den") %in% names(es_wide))) {
  w_col <- "S_purpuratus_SE_Den"
  req_cols <- c("S_pulcher_Den", "S_purpuratus_Den", w_col)
  if (w_col %in% names(es_wide)) {
    es_sheep_purp <- es_wide[complete.cases(es_wide[, req_cols]), ]
    if (nrow(es_sheep_purp) >= 3) {
      meta_sheep_purp <- rma(yi = es_sheep_purp$S_purpuratus_Den,
                             vi = as.numeric(es_sheep_purp[[w_col]])^2,
                             mods = ~ S_pulcher_Den,
                             data = es_sheep_purp, method = "REML")
      print_meta_reg(meta_sheep_purp, "S. pulcher density -> S. purpuratus density")
    }
  }
}

# Model 5: P. interruptus biomass predicts S. purpuratus biomass
if (all(c("P_interruptus_Bio", "S_purpuratus_Bio") %in% names(es_wide))) {
  w_col <- "S_purpuratus_SE_Bio"
  req_cols <- c("P_interruptus_Bio", "S_purpuratus_Bio", w_col)
  if (w_col %in% names(es_wide)) {
    es_lob_purp_bio <- es_wide[complete.cases(es_wide[, req_cols]), ]
    if (nrow(es_lob_purp_bio) >= 3) {
      meta_lob_purp_bio <- rma(yi = es_lob_purp_bio$S_purpuratus_Bio,
                               vi = as.numeric(es_lob_purp_bio[[w_col]])^2,
                               mods = ~ P_interruptus_Bio,
                               data = es_lob_purp_bio, method = "REML")
      print_meta_reg(meta_lob_purp_bio, "P. interruptus biomass -> S. purpuratus biomass")
    }
  }
}

# Model 6: S. pulcher biomass predicts S. purpuratus biomass
if (all(c("S_pulcher_Bio", "S_purpuratus_Bio") %in% names(es_wide))) {
  w_col <- "S_purpuratus_SE_Bio"
  req_cols <- c("S_pulcher_Bio", "S_purpuratus_Bio", w_col)
  if (w_col %in% names(es_wide)) {
    es_sheep_purp_bio <- es_wide[complete.cases(es_wide[, req_cols]), ]
    if (nrow(es_sheep_purp_bio) >= 3) {
      meta_sheep_purp_bio <- rma(yi = es_sheep_purp_bio$S_purpuratus_Bio,
                                 vi = as.numeric(es_sheep_purp_bio[[w_col]])^2,
                                 mods = ~ S_pulcher_Bio,
                                 data = es_sheep_purp_bio, method = "REML")
      print_meta_reg(meta_sheep_purp_bio, "S. pulcher biomass -> S. purpuratus biomass")
    }
  }
}

####################################################################################################
## Done ############################################################################################
####################################################################################################

cat("\nMeta-analysis complete.\n")
cat("Objects available: meta_biomass, meta_density, Table2\n")
cat("table_02_meta_analysis.csv written to:", here::here("data", "table_02_meta_analysis.csv"), "\n")
