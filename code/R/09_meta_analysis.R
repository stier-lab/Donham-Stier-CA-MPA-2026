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
#   2. Fits multilevel random-effects models separately for biomass and density
#   3. Detects and removes influential outliers using Cook's distance
#   4. Produces Table 2: Meta-analytic estimates by taxa
#   5. Produces Table 3: Cross-taxa relationships (trophic cascade tests)
#
# META-ANALYSIS MODEL STRUCTURE:
#   yi ~ Taxa - 1, random = list(~1|MPA, ~1|Source)
#   - yi: Effect size (mean difference on log-ratio scale)
#   - V: Sampling variance (SE^2)
#   - Taxa: Fixed effect moderator (estimates each taxa's mean effect)
#   - MPA: Random intercept (accounts for non-independence within MPAs)
#   - Source: Random intercept (accounts for program-level variation)
#   - Method: REML (Restricted Maximum Likelihood)
#
# WHY SEPARATE MODELS FOR BIOMASS AND DENSITY:
#   Biomass and density estimates for the same taxa/MPA/source are non-independent
#   (they come from the same surveys). Running them together would violate
#   independence assumptions. By splitting, we can properly estimate variance.
#
# OUTLIER DETECTION:
#   Cook's distance identifies observations that disproportionately influence
#   model estimates. Threshold is 4/n (standard in meta-analysis).
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
## Biomass meta-analysis ###########################################################################
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
typical_v_bio <- mean(biomass_clean$vi)
pseudo_I2_bio <- 100 * total_hetero_bio / (total_hetero_bio + typical_v_bio)
cat("Pseudo-I² (total):", round(pseudo_I2_bio, 1), "%\n")
cat("Interpretation: ",
    ifelse(pseudo_I2_bio < 25, "Low heterogeneity",
           ifelse(pseudo_I2_bio < 75, "Moderate heterogeneity", "High heterogeneity")), "\n")

####################################################################################################
## Density meta-analysis ###########################################################################
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
typical_v_den <- mean(density_clean$vi)
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
    pval      = round(coef_table[, "pval"], 4),
    CI_lower  = round(coef_table[, "ci.lb"], 4),
    CI_upper  = round(coef_table[, "ci.ub"], 4),
    Response  = response,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

# Build Table 2 from both models
Table2_biomass <- extract_meta_table(meta_biomass, biomass_clean, "Biomass")
Table2_density <- extract_meta_table(meta_density, density_clean, "Density")
Table2 <- rbind(Table2_biomass, Table2_density)

# Report sample size summary
cat("\n--- Sample Size Summary ---\n")
cat("Total biomass effect sizes (after outlier removal):", nrow(biomass_clean), "\n")
cat("  - Number of unique MPAs:", length(unique(biomass_clean$MPA)), "\n")
cat("  - Number of unique data sources:", length(unique(biomass_clean$Source)), "\n")
if ("N" %in% names(biomass_clean)) {
  cat("  - Total underlying observations (N):", sum(as.numeric(biomass_clean$N), na.rm = TRUE), "\n")
}
cat("Total density effect sizes (after outlier removal):", nrow(density_clean), "\n")
cat("  - Number of unique MPAs:", length(unique(density_clean$MPA)), "\n")
cat("  - Number of unique data sources:", length(unique(density_clean$Source)), "\n")
if ("N" %in% names(density_clean)) {
  cat("  - Total underlying observations (N):", sum(as.numeric(density_clean$N), na.rm = TRUE), "\n")
}

# Order taxa to match manuscript Table 2 presentation
taxa_order <- c("S. purpuratus", "M. franciscanus", "M. pyrifera",
                "P. interruptus", "S. pulcher")
Table2$Taxa <- factor(Table2$Taxa, levels = taxa_order)
Table2 <- Table2[order(Table2$Response, Table2$Taxa), ]
rownames(Table2) <- NULL

cat("\n============================\n")
cat("TABLE 2: Meta-analysis results\n")
cat("============================\n")
print(Table2)

# Export Table 2 as CSV to data directory
write.csv(Table2, here::here("data", "table_02_meta_analysis.csv"), row.names = FALSE)
cat("\nTable 2 exported to:", here::here("data", "table_02_meta_analysis.csv"), "\n")

####################################################################################################
## Table 3: Linear relationships between taxa (urchin density ~ kelp biomass) ######################
####################################################################################################

# Table 3 tests the trophic cascade relationship: whether urchin density effect sizes
# predict kelp biomass effect sizes across MPAs.
# This uses simple linear regression on the MPA-level effect sizes.

# Merge effect sizes by MPA to test cross-taxa relationships
# We need MPA-level estimates for different taxa in the same row

# Get primary effect sizes per MPA for each taxa x response combination
# Pivot both Mean and SE so we can use inverse-variance weights in Table 3
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
cat("TABLE 3: Cross-taxa linear models\n")
cat("============================\n")

# Model 1: S. purpuratus density predicts M. pyrifera biomass
# (urchin grazing hypothesis: more urchins -> less kelp)
if (all(c("S_purpuratus_Den", "M_pyrifera_Bio") %in% names(es_wide))) {
  es_purp_macro <- es_wide[complete.cases(es_wide[, c("S_purpuratus_Den", "M_pyrifera_Bio")]), ]
  if (nrow(es_purp_macro) >= 3) {
    # Inverse-variance weights from SE of the response variable
    w_col <- "M_pyrifera_SE_Bio"
    w <- if (w_col %in% names(es_purp_macro)) 1 / as.numeric(es_purp_macro[[w_col]])^2 else NULL
    lm_purp_macro <- lm(M_pyrifera_Bio ~ S_purpuratus_Den, data = es_purp_macro, weights = w)
    cat("\n--- S. purpuratus density -> M. pyrifera biomass ---\n")
    print(summary(lm_purp_macro))
  }
}

# Model 2: M. franciscanus density predicts M. pyrifera biomass
if (all(c("M_franciscanus_Den", "M_pyrifera_Bio") %in% names(es_wide))) {
  es_reds_macro <- es_wide[complete.cases(es_wide[, c("M_franciscanus_Den", "M_pyrifera_Bio")]), ]
  if (nrow(es_reds_macro) >= 3) {
    w_col <- "M_pyrifera_SE_Bio"
    w <- if (w_col %in% names(es_reds_macro)) 1 / as.numeric(es_reds_macro[[w_col]])^2 else NULL
    lm_reds_macro <- lm(M_pyrifera_Bio ~ M_franciscanus_Den, data = es_reds_macro, weights = w)
    cat("\n--- M. franciscanus density -> M. pyrifera biomass ---\n")
    print(summary(lm_reds_macro))
  }
}

# Model 3: P. interruptus density predicts S. purpuratus density
# (predator-prey: more lobsters -> fewer urchins)
if (all(c("P_interruptus_Den", "S_purpuratus_Den") %in% names(es_wide))) {
  es_lob_purp <- es_wide[complete.cases(es_wide[, c("P_interruptus_Den", "S_purpuratus_Den")]), ]
  if (nrow(es_lob_purp) >= 3) {
    w_col <- "S_purpuratus_SE_Den"
    w <- if (w_col %in% names(es_lob_purp)) 1 / as.numeric(es_lob_purp[[w_col]])^2 else NULL
    lm_lob_purp <- lm(S_purpuratus_Den ~ P_interruptus_Den, data = es_lob_purp, weights = w)
    cat("\n--- P. interruptus density -> S. purpuratus density ---\n")
    print(summary(lm_lob_purp))
  }
}

# Model 4: S. pulcher density predicts S. purpuratus density
# (predator-prey: more sheephead -> fewer urchins)
if (all(c("S_pulcher_Den", "S_purpuratus_Den") %in% names(es_wide))) {
  es_sheep_purp <- es_wide[complete.cases(es_wide[, c("S_pulcher_Den", "S_purpuratus_Den")]), ]
  if (nrow(es_sheep_purp) >= 3) {
    w_col <- "S_purpuratus_SE_Den"
    w <- if (w_col %in% names(es_sheep_purp)) 1 / as.numeric(es_sheep_purp[[w_col]])^2 else NULL
    lm_sheep_purp <- lm(S_purpuratus_Den ~ S_pulcher_Den, data = es_sheep_purp, weights = w)
    cat("\n--- S. pulcher density -> S. purpuratus density ---\n")
    print(summary(lm_sheep_purp))
  }
}

# Model 5: P. interruptus biomass predicts S. purpuratus biomass
if (all(c("P_interruptus_Bio", "S_purpuratus_Bio") %in% names(es_wide))) {
  es_lob_purp_bio <- es_wide[complete.cases(es_wide[, c("P_interruptus_Bio", "S_purpuratus_Bio")]), ]
  if (nrow(es_lob_purp_bio) >= 3) {
    w_col <- "S_purpuratus_SE_Bio"
    w <- if (w_col %in% names(es_lob_purp_bio)) 1 / as.numeric(es_lob_purp_bio[[w_col]])^2 else NULL
    lm_lob_purp_bio <- lm(S_purpuratus_Bio ~ P_interruptus_Bio, data = es_lob_purp_bio, weights = w)
    cat("\n--- P. interruptus biomass -> S. purpuratus biomass ---\n")
    print(summary(lm_lob_purp_bio))
  }
}

# Model 6: S. pulcher biomass predicts S. purpuratus biomass
if (all(c("S_pulcher_Bio", "S_purpuratus_Bio") %in% names(es_wide))) {
  es_sheep_purp_bio <- es_wide[complete.cases(es_wide[, c("S_pulcher_Bio", "S_purpuratus_Bio")]), ]
  if (nrow(es_sheep_purp_bio) >= 3) {
    w_col <- "S_purpuratus_SE_Bio"
    w <- if (w_col %in% names(es_sheep_purp_bio)) 1 / as.numeric(es_sheep_purp_bio[[w_col]])^2 else NULL
    lm_sheep_purp_bio <- lm(S_purpuratus_Bio ~ S_pulcher_Bio, data = es_sheep_purp_bio, weights = w)
    cat("\n--- S. pulcher biomass -> S. purpuratus biomass ---\n")
    print(summary(lm_sheep_purp_bio))
  }
}

####################################################################################################
## Done ############################################################################################
####################################################################################################

cat("\nMeta-analysis complete.\n")
cat("Objects available: meta_biomass, meta_density, Table2\n")
cat("table_02_meta_analysis.csv written to:", here::here("data", "table_02_meta_analysis.csv"), "\n")
