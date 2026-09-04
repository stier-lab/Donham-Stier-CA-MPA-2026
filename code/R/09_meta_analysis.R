# =============================================================================
# 09_meta_analysis.R
# =============================================================================
#
# PURPOSE:
#   This script answers the central question: "On average, how much do MPAs
#   change the biomass and density of each kelp forest species?"
#
#   It takes the individual MPA-level effect sizes computed in 08_effect_sizes.R
#   and synthesizes them into overall (network-wide) estimates for each of the
#   five focal taxa (lobster, sheephead, purple urchin, red urchin, giant kelp).
#   It also tests whether effect sizes on different trophic levels are related
#   to each other, which is the statistical test for the trophic cascade hypothesis.
#
# WHAT THIS SCRIPT DOES (in order):
#   1. Data prep: converts SE to sampling variance (vi = SE^2) for meta-analysis
#   2. PRIMARY ANALYSIS -- Joint multilevel models (one for biomass, one for density):
#      Estimates the mean MPA effect for each taxon while sharing variance
#      components across taxa, which stabilizes estimates for taxa with few MPAs.
#   3. SENSITIVITY -- Source random effect: checks whether including the monitoring
#      program (PISCO/KFM/LTER) as a random effect changes conclusions.
#   4. SENSITIVITY -- Variance components: reports how much of the variation in
#      effect sizes comes from differences among MPAs vs. monitoring programs.
#   5. SENSITIVITY -- Per-taxa models: fits separate models for each taxon
#      (instead of a joint model) to show results are robust to modeling choice.
#   6. SENSITIVITY -- Outlier diagnostics: compares 4 outlier-handling approaches
#      to show results are not driven by any single extreme effect size.
#   7. SENSITIVITY -- Leave-one-out for kelp: tests whether the kelp biomass
#      result depends on any single MPA (it's a borderline result).
#   8. TABLE 2 (manuscript): The main results table -- meta-analytic mean effect
#      sizes with 95% CIs and FDR-corrected p-values for each taxon x response.
#   9. TABLE 3 (manuscript): Cross-taxa meta-regressions testing the trophic
#      cascade hypothesis (e.g., does predator increase predict urchin decrease?).
#
# PRIMARY META-ANALYSIS STRUCTURE (Joint Model, No Outlier Removal):
#   yi ~ Taxa - 1, random = list(~1|MPA, ~1|Source)
#
#   In plain English: the mean effect size differs by taxon (fixed effect), and
#   effect sizes from the same MPA or the same monitoring program may be more
#   similar to each other than to random draws (random intercepts). The "- 1"
#   means each taxon gets its own coefficient directly (cell-means), rather than
#   one taxon being the baseline with others expressed as differences.
#
#   All taxa share common tau2(MPA) and tau2(Source), yielding more stable
#   estimates especially for small-k taxa (e.g., P. interruptus Bio k=2).
#   Statistical rationale: Borenstein & Higgins 2013, Cochrane Handbook Ch. 10,
#   Nakagawa & Santos 2012, Noble et al. 2017.
#
# SENSITIVITY: PER-TAXA MODEL STRUCTURE:
#   For each taxon-response separately:
#   yi ~ 1, random = ~1|MPA   (for k >= 5 and >= 3 MPAs)
#   yi ~ 1                     (simple RE for k < 5)
#   Shows robustness to modeling choice (per-taxa vs joint).
#
# WHY SEPARATE MODELS FOR BIOMASS AND DENSITY:
#   Biomass and density estimates for the same taxa/MPA/source are non-independent
#   (they come from the same surveys measuring the same organisms). Running them
#   together in one model would violate independence assumptions. Splitting into
#   separate biomass and density models avoids this problem.
#
# METHODS REFERENCE (from manuscript):
#   "We fit multilevel meta-analysis models with restricted maximum-likelihood
#    estimation (REML) using effect size estimates (ES) from pBACIPS and/or CI
#    analyses for each taxa using the rma.mv function in the metafor package.
#    Taxa was a fixed effect (moderator) in meta-analysis models with MPA as
#    a random effect."
#
# NOTE ON ASSUMPTIONS:
#   pBACIPS assumes parallel trends between MPA and reference sites in the
#   before-period. With typically short before-periods (3-5 years) and high
#   interannual variability, formal tests for parallel trends have low
#   statistical power. The before-period trends are implicitly assumed
#   parallel by the NLS model structure in 08_effect_sizes.R.
#
# INPUTS:
#   - SumStats.Final: Effect size estimates from 08_effect_sizes.R
#     (one row per MPA x taxon x response type, with columns Mean, SE, etc.)
#
# OUTPUTS (R objects kept in memory for downstream scripts):
#   - meta_biomass_full / meta_density_full: Primary joint rma.mv model objects
#   - meta_biomass / meta_density: Joint models after Cook's D outlier removal
#   - Table2: Summary dataframe (also exported as CSV)
#
# OUTPUTS (CSV files):
#   - tables/table_02_meta_analysis.csv                    Table 2: joint model estimates + heterogeneity
#   - tables/table_03_cross_taxa_meta_regression.csv       Table 3: trophic cascade meta-regressions
#   - tables/table_s_outlier_sensitivity.csv               Table S9: outlier sensitivity (4 methods)
#   - tables/table_s_kelp_leave1out.csv                    Table S8: leave-one-out M. pyrifera biomass
#   - tables/table_s_variance_components.csv               Table S2: variance components with CIs
#   - tables/table_s_source_sensitivity_models.csv         Table S3a: source RE model comparison (AIC/BIC)
#   - tables/table_s_source_sensitivity_coefficients.csv   Table S3b: source RE coefficient comparison
#   - outputs/filter_audit_meta_analysis.csv               Joint model filter audit (which rows were outliers)
#   - outputs/filter_audit_pertaxa_meta.csv                Per-taxa filter audit
#
# DEPENDENCIES:
#   Requires 00-08 scripts to be sourced first (loads SumStats.Final)
#
# AUTHORS: Emily Donham & Adrian Stier
# PROJECT: CA MPA Kelp Forest pBACIPS Analysis
# =============================================================================

####################################################################################################
## SECTION A: Prepare effect size data for meta-analysis ###########################################
####################################################################################################
#
# WHAT: Convert the individual MPA-level effect sizes from 08_effect_sizes.R into
#   the format required by metafor. The key transformation is computing the
#   sampling variance (vi = SE^2), which tells the meta-analysis how precisely
#   each effect size was estimated. Better-estimated effect sizes (lower SE)
#   get more weight in the overall average.
#
# WHY: metafor::rma.mv() requires a variance column (V) rather than SE.
#   We also split the data into separate biomass and density subsets because
#   these are non-independent (measured on the same organisms at the same sites).

# Compute sampling variance: vi = SE^2
# metafor weights each effect size by 1/vi, so more precise estimates count more
SumStats.Final$vi <- as.numeric(SumStats.Final$SE)^2

# Ensure numeric types (these may have been stored as character in the CSV from 08)
SumStats.Final$Mean <- as.numeric(SumStats.Final$Mean)
SumStats.Final$vi   <- as.numeric(SumStats.Final$vi)

# --- Data quality checks ---
# Flag any rows with non-positive SE or non-finite variance (would break the model)
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

# Flag extremely large effect sizes (|lnRR| > 5 means >150x difference MPA vs ref)
# These are not removed, just flagged -- they may be real (e.g., kelp at some sites)
extreme_es <- sum(abs(SumStats.Final$Mean) > 5, na.rm = TRUE)
if (extreme_es > 0) {
  warning("Found ", extreme_es, " effect sizes with |Mean| > 5 (may be outliers)")
}

# Convert grouping variables to factors (required by rma.mv formula interface)
SumStats.Final$Taxa <- factor(SumStats.Final$Taxa)
SumStats.Final$Source <- factor(SumStats.Final$Source)

# Split into biomass and density subsets
# These MUST be modeled separately because biomass and density are measured on
# the same organisms at the same sites -- they are not independent observations
biomass_data <- subset(SumStats.Final, Resp == "Bio")
density_data <- subset(SumStats.Final, Resp == "Den")

# Verify both subsets have data
stopifnot("No biomass data for meta-analysis" = nrow(biomass_data) > 0)
stopifnot("No density data for meta-analysis" = nrow(density_data) > 0)

cat("Biomass observations:", nrow(biomass_data), "\n")
cat("Density observations:", nrow(density_data), "\n")

####################################################################################################
## SECTION B: Joint Biomass Meta-Analysis (PRIMARY ANALYSIS) #######################################
####################################################################################################
#
# WHAT: Estimate the average MPA effect on BIOMASS for each of the 5 taxa, pooling
#   across all MPAs and monitoring programs. This is the primary analysis that
#   populates the biomass rows of manuscript Table 2.
#
# WHY A JOINT MODEL: By fitting all 5 taxa in one model, we "borrow strength"
#   across taxa -- the variance components (tau2_MPA and tau2_Source) are estimated
#   from all the data, not just from one taxon. This matters most for taxa with
#   very few effect sizes (e.g., lobster biomass has only k=2 MPAs).
#
# HOW TO READ THE OUTPUT: Each taxon's coefficient IS its estimated mean lnRR.
#   Positive = higher inside MPAs. The p-value tests whether lnRR differs from 0.
#   For example, lnRR = 0.61 means biomass is exp(0.61) = 1.84x higher inside MPAs.
#
# KNOWN LIMITATION: Biomass and density measurements come from the same organisms
#   at the same sites, so they are not independent. The (1|MPA) random effect
#   partially accounts for this, but the FDR correction across bio + den tests
#   may still be slightly anti-conservative. See 13_additional_analyses.R.

# --- Step 1: Fit the joint model on ALL biomass data (no outlier removal) ---
# This is the PRIMARY model. Cook's D outlier removal (Step 2-3) is for sensitivity only.
#
# Model formula explained:
#   yi = Mean:            effect size on log-ratio scale (lnRR)
#   V = vi:               sampling variance (SE^2) -- weights more precise estimates higher
#   mods = ~ Taxa - 1:    separate intercept for each taxon (cell-means parameterization)
#   random = ~1|MPA:      effect sizes from the same MPA may be correlated
#   random = ~1|Source:    effect sizes from the same monitoring program may be correlated
#   method = "REML":      restricted maximum likelihood (standard for inference)
#   test = "t":           use t-distribution for p-values (more conservative than z with small k)
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
    warning("Joint biomass meta-analysis model (primary) failed to converge: ", e$message)
    NULL
  }
)

# Initialize joint model objects to NULL (set below if model succeeds)
meta_biomass <- NULL
biomass_clean <- biomass_data
cooks_bio <- rep(0, nrow(biomass_data))
n_bio <- nrow(biomass_data)
cooks_threshold_bio <- 4 / n_bio
outliers_bio <- integer(0)
pseudo_I2_bio <- NA_real_

if (!is.null(meta_biomass_full)) {
  # Verify the model converged and produced valid estimates
  if (!is.null(meta_biomass_full$not.converged) && meta_biomass_full$not.converged) {
    warning("Biomass meta-analysis model did not converge")
  }
  stopifnot("Model coefficients contain NA" = !any(is.na(coef(meta_biomass_full))))
  stopifnot("Model coefficients not finite" = all(is.finite(coef(meta_biomass_full))))

  cat("\n--- Biomass model (full data) ---\n")
  print(summary(meta_biomass_full))

  # --- Step 2: Identify influential observations using Cook's distance ---
  # Cook's distance measures how much the model results would change if you removed
  # one observation. A high value means that single MPA-taxon combination is
  # disproportionately driving the overall result. The 4/n threshold is standard
  # in meta-analytic practice (Viechtbauer & Cheung 2010).
  cooks_bio <- cooks.distance(meta_biomass_full)
  p_bio <- length(coef(meta_biomass_full))

  # Flag observations exceeding the 4/n threshold as influential
  outliers_bio <- which(cooks_bio > cooks_threshold_bio)

  cat("\nBiomass Cook's distance threshold (4/n):", round(cooks_threshold_bio, 4), "\n")
  cat("Biomass outliers detected:", length(outliers_bio), "\n")
  if (length(outliers_bio) > 0) {
    cat("Outlier rows (Taxa, MPA, Mean):\n")
    print(biomass_data[outliers_bio, c("Taxa", "MPA", "Mean", "SE")])
  }

  # --- Step 3: Remove outliers and refit (SENSITIVITY, not primary) ---
  # The primary analysis uses meta_biomass_full (no removal). This cleaned model
  # (meta_biomass) is kept for the outlier sensitivity comparison in Table S9.
  if (length(outliers_bio) > 0) {
    pct_removed <- 100 * length(outliers_bio) / n_bio
    cat("Removing", length(outliers_bio), "outliers (", round(pct_removed, 1), "% of data)\n")
    if (pct_removed > 20) {
      warning("Removing >20% of biomass data as outliers - consider reviewing threshold")
    }
    biomass_clean <- biomass_data[-outliers_bio, ]
  } else {
    biomass_clean <- biomass_data
  }

  stopifnot("Insufficient biomass data after outlier removal" = nrow(biomass_clean) >= 10)

  # Refit the same joint model on the cleaned data
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
      warning("Biomass meta-analysis (clean) failed: ", e$message)
      NULL
    }
  )

  if (!is.null(meta_biomass)) {
    cat("\n--- Biomass model (outliers removed) ---\n")
    print(summary(meta_biomass))

    # --- Heterogeneity statistics ---
    # These describe how variable the MPA effects are across sites and programs.
    # tau2_MPA: variance in effect sizes attributable to differences among MPAs
    #   (e.g., some MPAs may have stronger effects due to enforcement, habitat, etc.)
    # tau2_Source: variance attributable to differences among monitoring programs
    #   (PISCO, KFM, LTER may estimate slightly different effect sizes due to methods)
    cat("\n--- Biomass Heterogeneity Statistics ---\n")
    cat("Between-MPA variance (tau²_MPA):", round(meta_biomass$sigma2[1], 4), "\n")
    cat("Between-Source variance (tau²_Source):", round(meta_biomass$sigma2[2], 4), "\n")

    # Pseudo-I2: what fraction of total variability is "real" heterogeneity vs.
    # just sampling noise? I2 > 75% = substantial heterogeneity, meaning MPAs
    # differ meaningfully in their effect sizes (expected for ecological reasons).
    # NOTE: Very high I2 (>90%) is normal here -- MPAs vary in size, age,
    # enforcement, habitat, and species composition. Table S2 gives the more
    # informative decomposition into MPA vs. Source components.
    total_hetero_bio <- sum(meta_biomass$sigma2)
    # Harmonic mean of sampling variances (Nakagawa & Santos 2012 recommendation
    # for multilevel models; more robust than arithmetic mean to outlier variances)
    typical_v_bio <- 1 / mean(1 / biomass_clean$vi)
    pseudo_I2_bio <- 100 * total_hetero_bio / (total_hetero_bio + typical_v_bio)
    cat("Pseudo-I² (total):", round(pseudo_I2_bio, 1), "%\n")
    cat("Interpretation: ",
        ifelse(pseudo_I2_bio < 25, "Low heterogeneity",
               ifelse(pseudo_I2_bio < 75, "Moderate heterogeneity", "High heterogeneity")), "\n")
  }
} else {
  cat("\n--- Joint biomass model (primary) skipped (model failed) ---\n")
  cat("WARNING: Primary analysis model failed. Per-taxa models will be used as fallback.\n")
}

####################################################################################################
## SECTION C: Joint Density Meta-Analysis (PRIMARY ANALYSIS) #######################################
####################################################################################################
#
# WHAT: Same analysis as Section B, but for DENSITY instead of biomass.
#   Estimates the average MPA effect on density for each of the 5 taxa.
#   Populates the density rows of manuscript Table 2.
#
# The model structure is identical to the biomass model:
#   yi ~ Taxa - 1, random = list(~1|MPA, ~1|Source)

# --- Step 1: Fit the joint model on ALL density data (no outlier removal) ---
meta_density_full <- tryCatch(
  rma.mv(
    yi     = Mean,
    V      = vi,
    mods   = ~ Taxa - 1,
    random = list(~ 1 | MPA, ~ 1 | Source),
    data   = density_data,
    method = "REML",
    test   = "t"
  ),
  error = function(e) {
    warning("Joint density meta-analysis model (primary) failed to converge: ", e$message)
    NULL
  }
)

# Initialize joint density model objects to NULL (set below if model succeeds)
meta_density <- NULL
density_clean <- density_data
cooks_den <- rep(0, nrow(density_data))
n_den <- nrow(density_data)
cooks_threshold_den <- 4 / n_den
outliers_den <- integer(0)
pseudo_I2_den <- NA_real_

if (!is.null(meta_density_full)) {
  # Verify convergence
  if (!is.null(meta_density_full$not.converged) && meta_density_full$not.converged) {
    warning("Density meta-analysis: optimizer reports non-convergence; results may be unreliable")
  }

  cat("\n--- Density model (full data) ---\n")
  print(summary(meta_density_full))

  # --- Step 2: Cook's distance outlier detection (same logic as biomass above) ---
  cooks_den <- cooks.distance(meta_density_full)
  outliers_den <- which(cooks_den > cooks_threshold_den)

  cat("\nDensity Cook's distance threshold (4/n):", round(cooks_threshold_den, 4), "\n")
  cat("Density outliers detected:", length(outliers_den), "\n")
  if (length(outliers_den) > 0) {
    cat("Outlier rows (Taxa, MPA, Mean):\n")
    print(density_data[outliers_den, c("Taxa", "MPA", "Mean", "SE")])
  }

  # --- Step 3: Remove outliers and refit (SENSITIVITY, not primary) ---
  if (length(outliers_den) > 0) {
    density_clean <- density_data[-outliers_den, ]
  } else {
    density_clean <- density_data
  }

  # Refit the same joint model on the cleaned data
  meta_density <- tryCatch(
    rma.mv(
      yi     = Mean,
      V      = vi,
      mods   = ~ Taxa - 1,
      random = list(~ 1 | MPA, ~ 1 | Source),
      data   = density_clean,
      method = "REML",
      test   = "t"
    ),
    error = function(e) {
      warning("Density meta-analysis (clean) failed: ", e$message)
      NULL
    }
  )

  if (!is.null(meta_density)) {
    cat("\n--- Density model (outliers removed) ---\n")
    print(summary(meta_density))

    # Heterogeneity statistics (same approach as biomass -- see Section B for explanation)
    cat("\n--- Density Heterogeneity Statistics ---\n")
    cat("Between-MPA variance (tau²_MPA):", round(meta_density$sigma2[1], 4), "\n")
    cat("Between-Source variance (tau²_Source):", round(meta_density$sigma2[2], 4), "\n")
    total_hetero_den <- sum(meta_density$sigma2)
    typical_v_den <- 1 / mean(1 / density_clean$vi)
    pseudo_I2_den <- 100 * total_hetero_den / (total_hetero_den + typical_v_den)
    cat("Pseudo-I² (total):", round(pseudo_I2_den, 1), "%\n")
    cat("Interpretation: ",
        ifelse(pseudo_I2_den < 25, "Low heterogeneity",
               ifelse(pseudo_I2_den < 75, "Moderate heterogeneity", "High heterogeneity")), "\n")
  }
} else {
  cat("\n--- Joint density model (primary) skipped (model failed) ---\n")
  cat("WARNING: Primary analysis model failed. Per-taxa models will be used as fallback.\n")
}

####################################################################################################
## SECTION D: Sensitivity -- Does the Source random effect matter? (Tables S3a, S3b) ###############
####################################################################################################
#
# WHAT: Compare the primary model (which includes both MPA and Source random effects)
#   to a simpler model with only the MPA random effect. If conclusions are the same
#   either way, we can be confident that differences among monitoring programs
#   (PISCO vs. KFM vs. LTER) are not driving the results.
#
# WHY THIS MATTERS: We only have 3-4 monitoring programs (Source levels), which is
#   below the recommended 5-6 minimum for reliable variance estimation. If the
#   Source random effect is poorly estimated, it could distort the taxa estimates.
#   This comparison shows whether that's a problem in practice.
#
# OUTPUT: Table S3a (AIC/BIC model comparison) and Table S3b (coefficient comparison)

cat("\n")
cat("====================================\n")
cat("SENSITIVITY ANALYSIS: Source Effect\n")
cat("====================================\n")

# Fit alternative models WITHOUT the Source random effect (MPA only)
meta_biomass_no_source <- if (!is.null(meta_biomass)) {
  tryCatch({
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
} else NULL

meta_density_no_source <- if (!is.null(meta_density)) {
  tryCatch({
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
} else NULL

# --- AIC/BIC model comparison ---
# Technical note: We must refit with ML (not REML) for valid AIC comparison.
# REML likelihoods cannot be compared across models with different random-effect
# structures. The REML models above are used for inference (coefficient estimates);
# these ML refits are used only for model selection (AIC/BIC).
cat("\n--- Model Comparison: With vs Without Source Random Effect ---\n")
cat("(Models refit with ML for valid AIC/BIC comparison)\n")

if (!is.null(meta_biomass) && !is.null(meta_biomass_no_source)) {
  # Refit both models with ML for AIC comparison
  meta_biomass_ml <- tryCatch(
    rma.mv(yi = Mean, V = vi, mods = ~ Taxa - 1,
           random = list(~ 1 | MPA, ~ 1 | Source),
           data = biomass_clean, method = "ML", test = "t"),
    error = function(e) { warning("ML refit (biomass with Source) failed: ", e$message); NULL })
  meta_biomass_no_source_ml <- tryCatch(
    rma.mv(yi = Mean, V = vi, mods = ~ Taxa - 1,
           random = list(~ 1 | MPA),
           data = biomass_clean, method = "ML", test = "t"),
    error = function(e) { warning("ML refit (biomass without Source) failed: ", e$message); NULL })

  if (!is.null(meta_biomass_ml) && !is.null(meta_biomass_no_source_ml)) {
    comparison_biomass <- data.frame(
      Model = c("With Source", "Without Source"),
      AIC = c(AIC(meta_biomass_ml), AIC(meta_biomass_no_source_ml)),
      BIC = c(BIC(meta_biomass_ml), BIC(meta_biomass_no_source_ml)),
      logLik = c(as.numeric(logLik(meta_biomass_ml)), as.numeric(logLik(meta_biomass_no_source_ml))),
      tau2_MPA = c(meta_biomass$sigma2[1], meta_biomass_no_source$sigma2[1]),
      tau2_Source = c(meta_biomass$sigma2[2], NA)
    )
    comparison_biomass$deltaAIC <- comparison_biomass$AIC - min(comparison_biomass$AIC)
    cat("\nBIOMASS Model Comparison:\n")
    print(comparison_biomass)
    cat("  -> Preferred model (lower AIC):", comparison_biomass$Model[which.min(comparison_biomass$AIC)], "\n")
  }
}

if (!is.null(meta_density) && !is.null(meta_density_no_source)) {
  # Refit both models with ML for AIC comparison
  meta_density_ml <- tryCatch(
    rma.mv(yi = Mean, V = vi, mods = ~ Taxa - 1,
           random = list(~ 1 | MPA, ~ 1 | Source),
           data = density_clean, method = "ML", test = "t"),
    error = function(e) { warning("ML refit (density with Source) failed: ", e$message); NULL })
  meta_density_no_source_ml <- tryCatch(
    rma.mv(yi = Mean, V = vi, mods = ~ Taxa - 1,
           random = list(~ 1 | MPA),
           data = density_clean, method = "ML", test = "t"),
    error = function(e) { warning("ML refit (density without Source) failed: ", e$message); NULL })

  if (!is.null(meta_density_ml) && !is.null(meta_density_no_source_ml)) {
    comparison_density <- data.frame(
      Model = c("With Source", "Without Source"),
      AIC = c(AIC(meta_density_ml), AIC(meta_density_no_source_ml)),
      BIC = c(BIC(meta_density_ml), BIC(meta_density_no_source_ml)),
      logLik = c(as.numeric(logLik(meta_density_ml)), as.numeric(logLik(meta_density_no_source_ml))),
      tau2_MPA = c(meta_density$sigma2[1], meta_density_no_source$sigma2[1]),
      tau2_Source = c(meta_density$sigma2[2], NA)
    )
    comparison_density$deltaAIC <- comparison_density$AIC - min(comparison_density$AIC)
    cat("\nDENSITY Model Comparison:\n")
    print(comparison_density)
    cat("  -> Preferred model (lower AIC):", comparison_density$Model[which.min(comparison_density$AIC)], "\n")
  }
}

# --- Coefficient comparison (from REML models -- these are the estimates we'd report) ---
# If the "with Source" and "without Source" columns are similar, the Source random
# effect is not meaningfully changing the taxa-level conclusions.
cat("\n--- Effect of Source on Coefficient Estimates ---\n")
if (!is.null(meta_biomass) && !is.null(meta_biomass_no_source)) {
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

if (!is.null(meta_density) && !is.null(meta_density_no_source)) {
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

# --- Export source sensitivity results to CSV ---

# Table S3a: Model comparison (AIC, BIC) for with/without Source random effect
source_model_comparison <- data.frame(
  Response = character(), Model = character(),
  AIC = numeric(), BIC = numeric(), logLik = numeric(),
  tau2_MPA = numeric(), tau2_Source = numeric(), deltaAIC = numeric(),
  stringsAsFactors = FALSE
)
if (exists("comparison_biomass")) {
  source_model_comparison <- rbind(source_model_comparison, comparison_biomass %>%
    dplyr::mutate(Response = "Biomass"))
}
if (exists("comparison_density")) {
  source_model_comparison <- rbind(source_model_comparison, comparison_density %>%
    dplyr::mutate(Response = "Density"))
}
if (nrow(source_model_comparison) > 0) {
  # Manuscript Table S3a
  write.csv(source_model_comparison,
            here::here("tables", "table_s_source_sensitivity_models.csv"), row.names = FALSE)
  cat("\nSource sensitivity model comparison exported to: tables/table_s_source_sensitivity_models.csv\n")
}

# Table S3b: Side-by-side coefficient estimates with/without Source random effect
source_coef_comparison <- data.frame(
  Response = character(), Taxa = character(),
  Est_with = numeric(), Est_without = numeric(),
  SE_with = numeric(), SE_without = numeric(),
  pval_with = numeric(), pval_without = numeric(),
  stringsAsFactors = FALSE
)
if (exists("coef_diff_bio")) {
  source_coef_comparison <- rbind(source_coef_comparison,
    coef_diff_bio %>% dplyr::mutate(Response = "Biomass"))
}
if (exists("coef_diff_den")) {
  source_coef_comparison <- rbind(source_coef_comparison,
    coef_diff_den %>% dplyr::mutate(Response = "Density"))
}
if (nrow(source_coef_comparison) > 0) {
  # Manuscript Table S3b
  write.csv(source_coef_comparison,
            here::here("tables", "table_s_source_sensitivity_coefficients.csv"), row.names = FALSE)
  cat("Source sensitivity coefficients exported to: tables/table_s_source_sensitivity_coefficients.csv\n")
}

####################################################################################################
## SECTION E: Variance Component Confidence Intervals (Table S2) ###################################
####################################################################################################
#
# WHAT: Report how much of the variation in MPA effect sizes comes from:
#   (a) real differences among MPAs (tau2_MPA), vs.
#   (b) differences among monitoring programs (tau2_Source), vs.
#   (c) sampling noise within each study.
#   Confidence intervals show how precisely we can estimate these components.
#
# WHY: Reviewers and readers want to know whether the high I2 (heterogeneity)
#   is driven by biological variation among MPAs (interesting) or by
#   methodological differences among PISCO/KFM/LTER (concerning).
#
# NOTE: With only 3-4 Source levels, the tau2_Source CI will be wide. That is
#   expected and should be reported transparently.

cat("\n")
cat("==========================================\n")
cat("VARIANCE COMPONENT CONFIDENCE INTERVALS\n")
cat("==========================================\n")
cat("(Using profile likelihood method from metafor::confint)\n")

# Profile likelihood CIs: the most reliable method for variance components
# (better than Wald CIs, which can go negative for variance parameters)
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

# --- Build a tidy summary table for Table S2 ---
# Technical note: metafor::confint() on rma.mv returns a nested list whose
# structure varies across metafor versions. The helper function below handles
# both the newer (numeric indexing) and older ("random" slot) formats, with
# a fallback to point estimates if neither works.

#' Extract variance component CIs from a metafor confint object
#' @param ci_obj   Output of confint(rma.mv model), or NULL
#' @param response_label  "Biomass" or "Density"
#' @param model_obj  The rma.mv model (used as fallback for point estimates)
#' @return List of data.frame rows, one per variance component
extract_varcomp_ci <- function(ci_obj, response_label, model_obj) {
  rows <- list()
  if (is.null(ci_obj)) return(rows)

  # metafor confint on rma.mv returns a list of per-component results
  # Try numeric indexing first (metafor >= 3.x)
  components <- c("MPA", "Source")
  for (i in seq_along(components)) {
    comp <- tryCatch({
      ci_comp <- ci_obj[[i]]
      if (!is.null(ci_comp) && "random" %in% names(ci_comp)) {
        # ci_comp$random is a matrix with row "sigma^2" or "tau^2"
        tau2_row <- ci_comp$random[1, ]  # First row is tau^2/sigma^2
        data.frame(
          Response = response_label,
          Component = components[i],
          Estimate = tau2_row["estimate"],
          CI_Lower = tau2_row["ci.lb"],
          CI_Upper = tau2_row["ci.ub"],
          stringsAsFactors = FALSE, row.names = NULL
        )
      } else {
        NULL
      }
    }, error = function(e) NULL)
    if (!is.null(comp)) rows[[length(rows) + 1]] <- comp
  }

  # Fallback: try the "random" slot (older metafor versions)
  if (length(rows) == 0 && "random" %in% names(ci_obj)) {
    for (i in seq_along(ci_obj$random)) {
      comp_name <- names(ci_obj$random)[i]
      if (is.null(comp_name)) comp_name <- paste0("Component_", i)
      rows[[length(rows) + 1]] <- data.frame(
        Response = response_label,
        Component = comp_name,
        Estimate = ci_obj$random[[i]]["estimate"],
        CI_Lower = ci_obj$random[[i]]["ci.lb"],
        CI_Upper = ci_obj$random[[i]]["ci.ub"],
        stringsAsFactors = FALSE, row.names = NULL
      )
    }
  }

  # If nothing worked, fall back to point estimates without CIs
  if (length(rows) == 0 && !is.null(model_obj)) {
    for (i in seq_along(model_obj$sigma2)) {
      rows[[length(rows) + 1]] <- data.frame(
        Response = response_label,
        Component = if (i <= length(components)) components[i] else paste0("Component_", i),
        Estimate = model_obj$sigma2[i],
        CI_Lower = NA_real_,
        CI_Upper = NA_real_,
        stringsAsFactors = FALSE, row.names = NULL
      )
    }
    if (length(rows) > 0) {
      warning(sprintf("Could not extract variance component CIs for %s; using point estimates only",
                       response_label))
    }
  }

  rows
}

tau2_rows <- c(
  extract_varcomp_ci(ci_biomass, "Biomass", meta_biomass),
  extract_varcomp_ci(ci_density, "Density", meta_density)
)

if (length(tau2_rows) > 0) {
  tau2_summary <- do.call(rbind, tau2_rows)
  rownames(tau2_summary) <- NULL

  # Add pseudo-I2 for each response (computed earlier for joint models)
  tau2_summary$Pseudo_I2 <- NA_real_
  tau2_summary$Pseudo_I2[tau2_summary$Response == "Biomass" & tau2_summary$Component == "MPA"] <- pseudo_I2_bio
  tau2_summary$Pseudo_I2[tau2_summary$Response == "Density" & tau2_summary$Component == "MPA"] <- pseudo_I2_den

  cat("\n--- Summary Table: Variance Components with CIs ---\n")
  print(tau2_summary)

  # Manuscript Table S2
  write.csv(tau2_summary, here::here("tables", "table_s_variance_components.csv"), row.names = FALSE)
  cat("\nVariance component summary exported to: tables/table_s_variance_components.csv\n")
} else if (!is.null(meta_biomass) || !is.null(meta_density)) {
  # Even if confint fails, export point estimates from joint models (if available)
  tau2_rows_fallback <- list()
  if (!is.null(meta_biomass)) {
    tau2_rows_fallback[[length(tau2_rows_fallback) + 1]] <- data.frame(
      Response = "Biomass", Component = "MPA",
      Estimate = meta_biomass$sigma2[1], CI_Lower = NA_real_, CI_Upper = NA_real_,
      Pseudo_I2 = pseudo_I2_bio, stringsAsFactors = FALSE)
    tau2_rows_fallback[[length(tau2_rows_fallback) + 1]] <- data.frame(
      Response = "Biomass", Component = "Source",
      Estimate = meta_biomass$sigma2[2], CI_Lower = NA_real_, CI_Upper = NA_real_,
      Pseudo_I2 = NA_real_, stringsAsFactors = FALSE)
  }
  if (!is.null(meta_density)) {
    tau2_rows_fallback[[length(tau2_rows_fallback) + 1]] <- data.frame(
      Response = "Density", Component = "MPA",
      Estimate = meta_density$sigma2[1], CI_Lower = NA_real_, CI_Upper = NA_real_,
      Pseudo_I2 = pseudo_I2_den, stringsAsFactors = FALSE)
    tau2_rows_fallback[[length(tau2_rows_fallback) + 1]] <- data.frame(
      Response = "Density", Component = "Source",
      Estimate = meta_density$sigma2[2], CI_Lower = NA_real_, CI_Upper = NA_real_,
      Pseudo_I2 = NA_real_, stringsAsFactors = FALSE)
  }
  tau2_summary <- do.call(rbind, tau2_rows_fallback)
  write.csv(tau2_summary, here::here("tables", "table_s_variance_components.csv"), row.names = FALSE)
  cat("\nVariance component summary (point estimates only) exported to: tables/table_s_variance_components.csv\n")
} else {
  cat("\nVariance component summary skipped (joint models not available)\n")
}

cat("\n")
cat("NOTE: The Source random effect has only 3-4 levels (PISCO, KFM, LTER, Landsat),\n")
cat("which is below the recommended minimum of 5-6 for reliable variance estimation.\n")
cat("Interpret tau²_Source estimates with caution; wide CIs are expected.\n")

####################################################################################################
## SECTION F: Per-Taxa Meta-Analysis (Sensitivity Check) ###########################################
####################################################################################################
#
# WHAT: Instead of fitting one big model with all 5 taxa (the joint model above),
#   fit a separate model for each taxon individually. This lets each taxon have
#   its own variance component, rather than sharing tau2 across all taxa.
#
# WHY: If the joint model and per-taxa models give similar answers, we can be
#   confident the results are not an artifact of the modeling choice. The joint
#   model is preferred (primary) because it borrows strength for small-k taxa,
#   but per-taxa models are the more traditional approach in ecology.
#
# MODEL SELECTION LOGIC:
#   k >= 5 and >= 3 MPAs: rma.mv() with MPA random effect (accounts for MPA clustering)
#   2 <= k < 5:           rma() simple random-effects model (not enough data for MPA RE)
#   k < 2:                report the single value descriptively (no model fitting possible)
#
# Cook's D outlier detection is also applied within each taxon as an additional
# sensitivity check (separate from the joint-model Cook's D in Sections B-C).

#' Fit per-taxa meta-analysis with within-taxon outlier detection
#'
#' Loops over each taxon in the data, fits the appropriate model based on sample
#' size, optionally removes Cook's D outliers, and returns tidy summary tables.
#'
#' @param data Data frame with columns: Taxa, MPA, Source, Mean, vi
#' @param response_label Character: "Biomass" or "Density"
#' @return List with: table (summary df), clean_data (filtered data), audit (all data with flags)
run_per_taxa_meta <- function(data, response_label) {
  taxa_list <- sort(unique(as.character(data$Taxa)))
  results_list <- list()
  results_full_list <- list()  # Full-model results (before Cook's D removal)
  clean_data_list <- list()
  full_data_list <- list()     # Full data per taxon (no outlier removal)
  audit_list <- list()

  for (taxon in taxa_list) {
    sub <- data[data$Taxa == taxon, ]
    k_input <- nrow(sub)

    cat(sprintf("\n  %s (%s): k_input = %d\n", taxon, response_label, k_input))

    # --- Handle k < 2: only one effect size, so no model can be fit ---
    # We just report the single value with a z-based CI. Interpret with caution.
    if (k_input < 2) {
      cat(sprintf("    -> k < 2: reporting descriptively\n"))
      mean_val <- as.numeric(sub$Mean[1])
      se_val <- as.numeric(sub$SE[1])
      z <- mean_val / se_val
      results_list[[taxon]] <- data.frame(
        Taxa = taxon, k = k_input,
        n_MPAs = length(unique(sub$MPA)),
        Estimate = round(mean_val, 4), SE = round(se_val, 4),
        tval = round(z, 4),
        pval = 2 * pnorm(-abs(z)),
        CI_lower = round(mean_val - 1.96 * se_val, 4),
        CI_upper = round(mean_val + 1.96 * se_val, 4),
        tau2 = NA_real_, I2 = NA_real_, QE = NA_real_, QEp = NA_real_,
        model_type = "descriptive (k<2)",
        Response = response_label, stringsAsFactors = FALSE
      )
      sub$Is_Outlier <- FALSE
      sub$Cooks_Distance <- NA_real_
      sub$Cooks_Threshold <- NA_real_
      sub$In_Final_Analysis <- TRUE
      sub$Exclusion_Reason <- "INCLUDED"
      audit_list[[taxon]] <- sub
      clean_data_list[[taxon]] <- sub
      full_data_list[[taxon]] <- sub
      results_full_list[[taxon]] <- results_list[[taxon]]  # Same as cleaned for k<2
      next
    }

    # --- Fit initial model (model type depends on sample size) ---
    # With enough data (k>=5, >=3 MPAs): multilevel model that accounts for MPA clustering
    # With less data (k 2-4): simple random-effects model (no MPA random effect)
    # If model fails: fall back to fixed-effect model as last resort
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

    # --- Extract results from the full model (before any outlier removal) ---
    # These go into the "no outlier removal" sensitivity row in the outlier table
    {
      coef_tbl_full <- coef(summary(model_full))
      stat_col_full <- if ("tval" %in% colnames(coef_tbl_full)) "tval" else "zval"
      if (inherits(model_full, "rma.mv")) {
        full_tau2 <- model_full$sigma2[1]
        full_type <- "multilevel RE (rma.mv)"
        full_QE <- if (!is.null(model_full$QE)) model_full$QE else NA_real_
        full_QEp <- if (!is.null(model_full$QEp)) model_full$QEp else NA_real_
        total_hetero_full <- sum(model_full$sigma2)
        typical_v_full <- 1 / mean(1 / sub$vi)
        full_I2 <- 100 * total_hetero_full / (total_hetero_full + typical_v_full)
      } else if (inherits(model_full, "rma.uni")) {
        full_tau2 <- model_full$tau2
        full_QE <- model_full$QE
        full_QEp <- model_full$QEp
        full_I2 <- model_full$I2
        full_type <- ifelse(model_full$method == "FE", "fixed effect (rma FE)", "simple RE (rma)")
      } else {
        full_tau2 <- NA_real_; full_QE <- NA_real_; full_QEp <- NA_real_
        full_I2 <- NA_real_; full_type <- "unknown"
      }
      results_full_list[[taxon]] <- data.frame(
        Taxa = taxon, k = model_full$k,
        n_MPAs = length(unique(sub$MPA)),
        Estimate = round(coef_tbl_full[1, "estimate"], 4),
        SE = round(coef_tbl_full[1, "se"], 4),
        tval = round(coef_tbl_full[1, stat_col_full], 4),
        pval = coef_tbl_full[1, "pval"],
        CI_lower = round(coef_tbl_full[1, "ci.lb"], 4),
        CI_upper = round(coef_tbl_full[1, "ci.ub"], 4),
        tau2 = round(full_tau2, 4),
        I2 = round(full_I2, 1),
        QE = round(full_QE, 4),
        QEp = full_QEp,
        model_type = full_type,
        Response = response_label, stringsAsFactors = FALSE
      )
      full_data_list[[taxon]] <- sub
    }

    # --- Per-taxa Cook's distance (within-taxon influence diagnostic) ---
    # Same approach as the joint model (Sections B-C), but applied within each taxon
    cooks_vals <- tryCatch(cooks.distance(model_full), error = function(e) rep(0, k_input))
    threshold <- 4 / k_input
    outlier_idx <- which(cooks_vals > threshold)

    cat(sprintf("    -> Cook's D threshold (4/%d): %.4f, outliers: %d\n",
                k_input, threshold, length(outlier_idx)))

    # Track audit info
    sub$Cooks_Distance <- cooks_vals
    sub$Cooks_Threshold <- threshold

    # --- Remove outliers and refit (only if we'd still have k >= 2 after removal) ---
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

    # --- Extract results from the final (possibly cleaned) model ---
    coef_tbl <- coef(summary(model))
    stat_col <- if ("tval" %in% colnames(coef_tbl)) "tval" else "zval"

    # Extract heterogeneity statistics (format differs between rma.mv and rma)
    if (inherits(model, "rma.mv")) {
      model_tau2 <- model$sigma2[1]
      model_type <- "multilevel RE (rma.mv)"
      # QE test of residual heterogeneity
      model_QE <- if (!is.null(model$QE)) model$QE else NA_real_
      model_QEp <- if (!is.null(model$QEp)) model$QEp else NA_real_
      # Pseudo-I2 for multilevel models
      total_hetero <- sum(model$sigma2)
      typical_v <- 1 / mean(1 / sub_clean$vi)
      model_I2 <- 100 * total_hetero / (total_hetero + typical_v)
    } else if (inherits(model, "rma.uni")) {
      model_tau2 <- model$tau2
      model_QE <- model$QE
      model_QEp <- model$QEp
      model_I2 <- model$I2
      model_type <- ifelse(model$method == "FE", "fixed effect (rma FE)",
                           "simple RE (rma)")
    } else {
      model_tau2 <- NA_real_; model_QE <- NA_real_; model_QEp <- NA_real_
      model_I2 <- NA_real_; model_type <- "unknown"
    }

    n_mpas <- length(unique(sub_clean$MPA))

    results_list[[taxon]] <- data.frame(
      Taxa = taxon, k = model$k,
      n_MPAs = n_mpas,
      Estimate = round(coef_tbl[1, "estimate"], 4),
      SE = round(coef_tbl[1, "se"], 4),
      tval = round(coef_tbl[1, stat_col], 4),
      pval = coef_tbl[1, "pval"],
      CI_lower = round(coef_tbl[1, "ci.lb"], 4),
      CI_upper = round(coef_tbl[1, "ci.ub"], 4),
      tau2 = round(model_tau2, 4),
      I2 = round(model_I2, 1),
      QE = round(model_QE, 4),
      QEp = model_QEp,
      model_type = model_type,
      Response = response_label, stringsAsFactors = FALSE
    )

    cat(sprintf("    -> Final k=%d, estimate=%.3f, SE=%.3f, p=%.4f\n",
                model$k, coef_tbl[1, "estimate"], coef_tbl[1, "se"], coef_tbl[1, "pval"]))
  }

  table_out <- do.call(rbind, results_list)
  rownames(table_out) <- NULL
  table_full_out <- do.call(rbind, results_full_list)
  rownames(table_full_out) <- NULL
  clean_out <- do.call(rbind, clean_data_list)
  full_out <- do.call(rbind, full_data_list)
  audit_out <- do.call(rbind, audit_list)

  list(table = table_out, table_full = table_full_out,
       clean_data = clean_out, full_data = full_out, audit = audit_out)
}

# --- Run per-taxa meta-analysis ---
cat("\n")
cat("============================================\n")
cat("PER-TAXA META-ANALYSIS (Sensitivity)\n")
cat("  Full dataset (no outlier removal)\n")
cat("  + Cook's D at 4/k\n")
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

# --- Outlier sensitivity table (Table S9): 4 approaches side by side ---
# This table lets readers verify that the main conclusions hold regardless
# of how outliers are handled. If all 4 approaches agree on the direction
# and significance of effects, the results are robust.
cat("\n")
cat("============================================\n")
cat("OUTLIER SENSITIVITY TABLE\n")
cat("============================================\n")

# Helper: extract a summary row from a joint rma.mv model for the sensitivity table
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

# Approach 1: Joint model, ALL data, no outliers removed (THIS IS THE PRIMARY ANALYSIS)
joint_full_bio <- extract_meta_table_helper(meta_biomass_full, biomass_data, "Biomass")
joint_full_den <- extract_meta_table_helper(meta_density_full, density_data, "Density")
joint_full_tbl <- rbind(joint_full_bio, joint_full_den)
joint_full_tbl$Method <- "Joint model, no removal (primary)"

# Approach 2: Per-taxa models, ALL data, no outliers removed
pertaxa_full_tbl <- rbind(pertaxa_biomass$table_full, pertaxa_density$table_full)
pertaxa_full_tbl_sens <- pertaxa_full_tbl[, c("Taxa", "k", "Estimate", "SE", "pval", "Response")]
pertaxa_full_tbl_sens$Method <- "Per-taxa models, no removal (sensitivity)"

# Approach 3: Per-taxa models WITH Cook's D outlier removal (within each taxon)
pertaxa_tbl <- rbind(pertaxa_biomass$table, pertaxa_density$table)
pertaxa_tbl_sens <- pertaxa_tbl[, c("Taxa", "k", "Estimate", "SE", "pval", "Response")]
pertaxa_tbl_sens$Method <- "Per-taxa Cook's D (4/k) (sensitivity)"

# Approach 4: Joint model WITH Cook's D outlier removal (legacy V5 approach)
joint_cooksd_tbl <- NULL
if (!is.null(meta_biomass) && !is.null(meta_density)) {
  joint_bio_tbl <- extract_meta_table_helper(meta_biomass, biomass_clean, "Biomass")
  joint_den_tbl <- extract_meta_table_helper(meta_density, density_clean, "Density")
  joint_cooksd_tbl <- rbind(joint_bio_tbl, joint_den_tbl)
  joint_cooksd_tbl$Method <- "Joint Cook's D (4/n) (legacy)"
}

sensitivity_table <- rbind(joint_full_tbl, pertaxa_full_tbl_sens,
                           pertaxa_tbl_sens, joint_cooksd_tbl)
sensitivity_table <- sensitivity_table[order(sensitivity_table$Response,
                                              sensitivity_table$Taxa,
                                              sensitivity_table$Method), ]

# Manuscript Table S9
write.csv(sensitivity_table, here::here("tables", "table_s_outlier_sensitivity.csv"), row.names = FALSE)
cat("Outlier sensitivity table saved to: tables/table_s_outlier_sensitivity.csv\n")
print(sensitivity_table)

# --- Leave-one-out analysis for kelp biomass (Table S8) ---
#
# WHAT: Drop each MPA one at a time and re-estimate the kelp biomass effect.
#   If the result flips from significant to non-significant (or vice versa)
#   when a single MPA is removed, the finding is fragile.
#
# WHY KELP SPECIFICALLY: The giant kelp biomass result is borderline significant
#   in the primary analysis. This test shows whether it depends on any single
#   MPA's extreme effect size (e.g., Scorpion SMR with very high lnRR).
cat("\n")
cat("============================================\n")
cat("LEAVE-ONE-OUT: M. pyrifera biomass\n")
cat("============================================\n")

kelp_clean <- pertaxa_biomass$full_data[pertaxa_biomass$full_data$Taxa == "M. pyrifera", ]
if (nrow(kelp_clean) >= 3) {
  # Technical note: metafor::leave1out() only works with rma() objects (not rma.mv),
  # so this uses a simple random-effects model without the MPA random effect.
  # Results may differ slightly from the multilevel primary model.
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

  # Manuscript Table S8
  write.csv(loo_df, here::here("tables", "table_s_kelp_leave1out.csv"), row.names = FALSE)
  cat("Leave-one-out results saved to: tables/table_s_kelp_leave1out.csv\n")
} else {
  cat("Insufficient kelp data for leave-one-out analysis\n")
}

####################################################################################################
## SECTION G: Build Table 2 -- The Main Results Table ##############################################
####################################################################################################
#
# WHAT: Assemble manuscript Table 2 from the PRIMARY joint models (meta_biomass_full
#   and meta_density_full -- no outlier removal). This is the central results table:
#   one row per taxon x response, showing the average MPA effect, 95% CI,
#   FDR-corrected p-value, and back-transformed response ratio (e.g., "1.84x").
#
# KEY COLUMNS IN THE OUTPUT:
#   Estimate:   mean lnRR (positive = higher inside MPAs)
#   CI_lower/upper: 95% confidence interval on lnRR
#   pval_fdr:   p-value corrected for multiple testing (Benjamini-Hochberg)
#   RR:         response ratio = exp(Estimate); RR=2 means 2x higher inside MPAs
#   Pct_Change: percent change = (RR - 1) * 100; e.g., +84% means 84% more inside MPAs

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
  # metafor column name is "zval" by default, but "tval" when test="t" was used
  stat_col <- if ("tval" %in% colnames(coef_table)) "tval" else "zval"

  # Strip the "Taxa" prefix from row names to get clean species names
  taxa_names <- gsub("^Taxa", "", rownames(coef_table))

  # k = number of effect sizes (MPA-level estimates) per taxon
  k_per_taxa <- vapply(taxa_names, function(taxon) {
    sum(data$Taxa == taxon, na.rm = TRUE)
  }, integer(1))

  # How many distinct MPAs contribute to each taxon's estimate
  n_mpas_per_taxa <- vapply(taxa_names, function(taxon) {
    length(unique(data$MPA[data$Taxa == taxon]))
  }, integer(1))

  # Heterogeneity statistics (shared across all taxa in the joint model)
  total_tau2 <- sum(model$sigma2)
  typical_v <- 1 / mean(1 / data$vi)
  pseudo_I2 <- round(100 * total_tau2 / (total_tau2 + typical_v), 1)
  model_QE <- if (!is.null(model$QE)) round(model$QE, 4) else NA_real_
  model_QEp <- if (!is.null(model$QEp)) model$QEp else NA_real_

  data.frame(
    Taxa       = taxa_names,
    k          = k_per_taxa,
    n_MPAs     = n_mpas_per_taxa,
    Estimate   = round(coef_table[, "estimate"], 4),
    SE         = round(coef_table[, "se"], 4),
    tval       = round(coef_table[, stat_col], 4),
    pval       = coef_table[, "pval"],  # Full precision for FDR correction
    CI_lower   = round(coef_table[, "ci.lb"], 4),
    CI_upper   = round(coef_table[, "ci.ub"], 4),
    tau2       = round(total_tau2, 4),
    I2         = pseudo_I2,
    QE         = model_QE,
    QEp        = model_QEp,
    model_type = "joint multilevel RE (rma.mv)",
    Response   = response,
    row.names  = NULL,
    stringsAsFactors = FALSE
  )
}

# --- PRIMARY TABLE 2: From the joint model with ALL data (no outlier removal) ---
# This is the table that goes in the manuscript. The joint model borrows strength
# across taxa for shared tau2 components (Borenstein & Higgins 2013; Noble et al. 2017).
Table2 <- rbind(
  extract_meta_table(meta_biomass_full, biomass_data, "Biomass"),
  extract_meta_table(meta_density_full, density_data, "Density")
)

# Alternative Table 2 versions (for sensitivity comparisons, not in main manuscript)
Table2_pertaxa <- rbind(pertaxa_biomass$table_full, pertaxa_density$table_full)     # Per-taxa, no removal
Table2_cooksd <- rbind(pertaxa_biomass$table, pertaxa_density$table)                # Per-taxa, Cook's D
Table2_joint_cooksd <- NULL                                                          # Joint, Cook's D (legacy)
if (!is.null(meta_biomass) && !is.null(meta_density)) {
  Table2_joint_cooksd <- rbind(
    extract_meta_table(meta_biomass, biomass_clean, "Biomass"),
    extract_meta_table(meta_density, density_clean, "Density")
  )
}

# Full data references for sample size reporting
biomass_full_pertaxa <- pertaxa_biomass$full_data
density_full_pertaxa <- pertaxa_density$full_data

# Report sample size summary
cat("\n--- Sample Size Summary (Full Dataset, No Outlier Removal) ---\n")
cat("Total biomass effect sizes:", nrow(biomass_full_pertaxa), "\n")
cat("  - Number of unique MPAs:", length(unique(biomass_full_pertaxa$MPA)), "\n")
cat("  - Number of unique data sources:", length(unique(biomass_full_pertaxa$Source)), "\n")
if ("N" %in% names(biomass_full_pertaxa)) {
  cat("  - Total underlying observations (N):", sum(as.numeric(biomass_full_pertaxa$N), na.rm = TRUE), "\n")
}
cat("Total density effect sizes:", nrow(density_full_pertaxa), "\n")
cat("  - Number of unique MPAs:", length(unique(density_full_pertaxa$MPA)), "\n")
cat("  - Number of unique data sources:", length(unique(density_full_pertaxa$Source)), "\n")
if ("N" %in% names(density_full_pertaxa)) {
  cat("  - Total underlying observations (N):", sum(as.numeric(density_full_pertaxa$N), na.rm = TRUE), "\n")
}

# Order taxa to match manuscript Table 2 layout
# (urchins first, then kelp, then predators -- matches trophic cascade narrative)
taxa_order <- c("S. purpuratus", "M. franciscanus", "M. pyrifera",
                "P. interruptus", "S. pulcher")
Table2$Taxa <- factor(Table2$Taxa, levels = taxa_order)
Table2 <- Table2[order(Table2$Response, Table2$Taxa), ]
rownames(Table2) <- NULL

# Flag taxa with very few effect sizes -- these estimates are less reliable
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

# Back-transform from lnRR to Response Ratio scale for the manuscript text
# lnRR = 0 -> RR = 1 -> no MPA effect
# lnRR > 0 -> RR > 1 -> higher inside MPAs (e.g., RR = 2 means 2x more)
# lnRR < 0 -> RR < 1 -> lower inside MPAs (e.g., RR = 0.5 means half as much)
Table2$RR        <- round(exp(Table2$Estimate), 3)
Table2$RR_lower  <- round(exp(Table2$CI_lower), 3)
Table2$RR_upper  <- round(exp(Table2$CI_upper), 3)
Table2$Pct_Change <- round((exp(Table2$Estimate) - 1) * 100, 1)

# Correct for multiple testing across all observed taxon-response coefficients
# (currently 9), using Benjamini-Hochberg FDR correction to control the false
# discovery rate.
# Use pval_fdr (not raw pval) for significance statements in the manuscript.
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

# --- Validate and export Table 2 ---
required_cols <- c("Taxa", "k", "n_MPAs", "Estimate", "SE", "tval", "pval",
                   "CI_lower", "CI_upper", "tau2", "I2", "QE", "QEp",
                   "model_type", "Response",
                   "RR", "RR_lower", "RR_upper", "Pct_Change", "pval_fdr")
# Include Flag column only when it has informative values (not all empty)
if (any(nchar(Table2$Flag) > 0, na.rm = TRUE)) {
  required_cols <- c(required_cols, "Flag")
}
missing_cols <- setdiff(required_cols, names(Table2))
if (length(missing_cols) > 0) {
  stop("Table2 is missing expected columns before export: ",
       paste(missing_cols, collapse = ", "))
}

# Manuscript Table 2
write.csv(Table2[, required_cols], here::here("tables", "table_02_meta_analysis.csv"), row.names = FALSE)
cat("\nTable 2 exported to:", here::here("tables", "table_02_meta_analysis.csv"), "\n")
cat("  Columns:", paste(required_cols, collapse = ", "), "\n")

####################################################################################################
## SECTION H: Filtering Audit -- Track exactly which observations enter the meta-analysis ##########
####################################################################################################
#
# WHAT: Create a row-by-row audit trail showing which effect sizes were included
#   in the joint meta-analysis and which were flagged as outliers (Cook's D).
#   Each row in the output CSV has the original data plus columns for:
#   Is_Outlier, Cooks_Distance, In_Final_Analysis, Exclusion_Reason.
#
# WHY: Transparency and reproducibility. Reviewers (and Emily) can verify exactly
#   which MPAs contributed to each taxa's estimate, and why any were excluded.
#   This is especially important for the lobster detail below, where k is very small.

cat("\n")
cat("==========================================\n")
cat("META-ANALYSIS FILTERING AUDIT\n")
cat("==========================================\n")

# --- Biomass audit ---
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

# --- Density audit ---
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

# --- Lobster detail ---
# P. interruptus has the fewest effect sizes of any taxon (often k=2 for biomass),
# so it's worth inspecting every single data point that goes into its estimate.
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
## SECTION I: Table 3 -- Cross-Taxa Meta-Regressions (Trophic Cascade Tests) ######################
####################################################################################################
#
# WHAT: Test whether the MPA effects on different trophic levels are correlated
#   in the direction predicted by a trophic cascade. For example:
#   - Do MPAs where purple urchins decreased more also show larger kelp increases?
#   - Do MPAs where lobsters increased more also show larger urchin decreases?
#
#   Each regression has one row per MPA (k is typically ~4), with the predictor
#   being one taxon's effect size and the response being another taxon's.
#
# WHY META-REGRESSION INSTEAD OF OLS:
#   Both the X and Y values are estimated with uncertainty (they are effect sizes
#   with SEs). Ordinary regression (lm) would ignore this uncertainty and bias
#   slopes toward zero. Meta-regression (metafor::rma with mods=) properly
#   weights by the response variable's sampling variance.
#
# IMPORTANT CAVEAT -- ATTENUATION BIAS:
#   Even with meta-regression, the predictor is measured with error, which
#   attenuates (shrinks) slopes toward zero. This means:
#   - Significant results are reliable (real signal overcame the bias)
#   - Non-significant results are ambiguous (could be no relationship, or
#     could be a real relationship masked by measurement error)
#   - With only k~4 MPAs per regression, statistical power is inherently limited
#
# ECOLOGICAL MODELS TESTED:
#   Model 1: S. purpuratus density  -> M. pyrifera biomass  (urchin grazing on kelp)
#   Model 2: M. franciscanus density -> M. pyrifera biomass (red urchin grazing on kelp)
#   Model 3: P. interruptus density -> S. purpuratus density (lobster predation on purple urchin)
#   Model 4: S. pulcher density -> S. purpuratus density    (sheephead predation on purple urchin)
#   Model 5: P. interruptus biomass -> S. purpuratus biomass (same as 3, but biomass)
#   Model 6: S. pulcher biomass -> S. purpuratus biomass    (same as 4, but biomass)

# --- Reshape data: one row per MPA with columns for each taxon's effect size ---
# We need each MPA's effect sizes for different taxa side-by-side in one row

# Check for duplicate MPA x Taxa_Resp combinations (will be averaged if found)
dup_check <- SumStats.Final %>%
  tidyr::unite("Taxa_Resp", Taxa, Resp, sep = "_", remove = FALSE) %>%
  dplyr::count(MPA, Taxa_Resp) %>%
  dplyr::filter(n > 1)
if (nrow(dup_check) > 0) {
  message("NOTE: ", nrow(dup_check), " MPA x Taxa_Resp duplicates will be averaged for Table 3 cross-taxa regressions")
}

# Pivot effect sizes to wide format: one column per taxon_response combination
es_wide_mean <- SumStats.Final %>%
  dplyr::select(Taxa, MPA, Mean, Resp) %>%          # Keep only needed columns
  dplyr::mutate(Mean = as.numeric(Mean)) %>%         # Ensure numeric
  tidyr::unite("Taxa_Resp", Taxa, Resp, sep = "_") %>%  # e.g., "S_purpuratus_Den"
  tidyr::pivot_wider(names_from = Taxa_Resp, values_from = Mean, values_fn = mean)
  # values_fn = mean: if an MPA has multiple sources, average their effect sizes

# Pivot sampling variances to wide format (needed to weight the meta-regression)
# When averaging n effect sizes: Var(mean) = sum(Var_i) / n^2 (assumes independence)
# NOTE: Most MPA-taxa-response combinations have only 1 source, so this rarely
# matters. For the few with 2+ sources, the overall uncertainty from k~4 MPAs
# dominates any variance underestimation from assumed independence.
es_wide_vi <- SumStats.Final %>%
  dplyr::select(Taxa, MPA, vi, Resp) %>%
  tidyr::unite("Taxa_Resp", Taxa, Resp, sep = "_vi_") %>%  # e.g., "M_pyrifera_vi_Bio"
  tidyr::pivot_wider(names_from = Taxa_Resp, values_from = vi,
                     values_fn = function(x) sum(x) / length(x)^2)

# Join effect sizes and their variances into one wide table
es_wide <- dplyr::left_join(es_wide_mean, es_wide_vi, by = "MPA")

# Clean column names: replace spaces/dots with underscores for formula compatibility
names(es_wide) <- gsub("[. ]+", "_", names(es_wide))

cat("\n============================\n")
cat("TABLE 3: Cross-taxa meta-regressions (metafor::rma)\n")
cat("============================\n")

#' Print and extract results from a single cross-taxa meta-regression
#'
#' For each trophic cascade test (e.g., "urchin density -> kelp biomass"),
#' this prints the key statistics to the console and returns a tidy data.frame
#' row for assembling Table 3.
#'
#' @param model An rma model object with one moderator
#' @param label Character description (e.g., "S. purpuratus density -> M. pyrifera biomass")
#' @param predictor Character: predictor taxon + response type
#' @param response Character: response taxon + response type
#' @return data.frame row with intercept, slope, CIs, QM test, tau2, I2, R2
extract_meta_reg <- function(model, label, predictor, response) {
  cat("\n---", label, "---\n")
  cat(sprintf("k = %d effect sizes\n", model$k))
  coef_tbl <- coef(summary(model))
  # Row 1 = intercept, Row 2 = moderator slope
  cat(sprintf("Intercept:  estimate = %.4f, SE = %.4f, p = %.4f\n",
              coef_tbl[1, "estimate"], coef_tbl[1, "se"], coef_tbl[1, "pval"]))
  cat(sprintf("Slope:      estimate = %.4f, SE = %.4f, p = %.4f\n",
              coef_tbl[2, "estimate"], coef_tbl[2, "se"], coef_tbl[2, "pval"]))
  cat(sprintf("  95%% CI for slope: [%.4f, %.4f]\n",
              coef_tbl[2, "ci.lb"], coef_tbl[2, "ci.ub"]))
  cat(sprintf("Test of moderator: QM(df=%d) = %.4f, p = %.4f\n",
              model$m, model$QM, model$QMp))
  cat(sprintf("Residual heterogeneity: QE(df=%d) = %.4f, p = %.4f\n",
              model$k - model$p, model$QE, model$QEp))
  cat(sprintf("tau² = %.4f, I² = %.1f%%\n", model$tau2, model$I2))
  cat(sprintf("R² (variance explained by moderator) = %.1f%%\n",
              ifelse(is.null(model$R2), NA, model$R2)))
  cat("\n")

  # Return structured row for CSV
  data.frame(
    Relationship = label,
    Predictor = predictor,
    Response_Var = response,
    k = model$k,
    Intercept_Est = round(coef_tbl[1, "estimate"], 4),
    Intercept_SE = round(coef_tbl[1, "se"], 4),
    Intercept_pval = coef_tbl[1, "pval"],
    Slope_Est = round(coef_tbl[2, "estimate"], 4),
    Slope_SE = round(coef_tbl[2, "se"], 4),
    Slope_pval = coef_tbl[2, "pval"],
    Slope_CI_lower = round(coef_tbl[2, "ci.lb"], 4),
    Slope_CI_upper = round(coef_tbl[2, "ci.ub"], 4),
    QM = round(model$QM, 4),
    QMp = model$QMp,
    QE = round(model$QE, 4),
    QEp = model$QEp,
    tau2 = round(model$tau2, 4),
    I2 = round(model$I2, 1),
    R2 = ifelse(is.null(model$R2), NA_real_, round(model$R2, 1)),
    stringsAsFactors = FALSE, row.names = NULL
  )
}

# Collector for all Table 3 meta-regression results
table3_rows <- list()

# Model 1: Does purple urchin density change predict kelp biomass change?
# Cascade prediction: negative slope (MPAs where urchins decreased more -> kelp increased more)
if (all(c("S_purpuratus_Den", "M_pyrifera_Bio") %in% names(es_wide))) {
  w_col <- "M_pyrifera_vi_Bio"
  req_cols <- c("S_purpuratus_Den", "M_pyrifera_Bio", w_col)
  if (w_col %in% names(es_wide)) {
    es_purp_macro <- es_wide[complete.cases(es_wide[, req_cols]), ]
    if (nrow(es_purp_macro) >= 3) {
      meta_purp_macro <- rma(yi = es_purp_macro$M_pyrifera_Bio,
                             vi = es_purp_macro[[w_col]],
                             mods = ~ S_purpuratus_Den,
                             data = es_purp_macro, method = "REML",
                             test = "knha")
      table3_rows[[length(table3_rows) + 1]] <- extract_meta_reg(
        meta_purp_macro,
        "S. purpuratus density -> M. pyrifera biomass",
        "S. purpuratus density", "M. pyrifera biomass")
    }
  }
}

# Model 2: Does red urchin density change predict kelp biomass change?
# Same cascade logic as Model 1, but for M. franciscanus instead of S. purpuratus
if (all(c("M_franciscanus_Den", "M_pyrifera_Bio") %in% names(es_wide))) {
  w_col <- "M_pyrifera_vi_Bio"
  req_cols <- c("M_franciscanus_Den", "M_pyrifera_Bio", w_col)
  if (w_col %in% names(es_wide)) {
    es_reds_macro <- es_wide[complete.cases(es_wide[, req_cols]), ]
    if (nrow(es_reds_macro) >= 3) {
      meta_reds_macro <- rma(yi = es_reds_macro$M_pyrifera_Bio,
                             vi = es_reds_macro[[w_col]],
                             mods = ~ M_franciscanus_Den,
                             data = es_reds_macro, method = "REML",
                             test = "knha")
      table3_rows[[length(table3_rows) + 1]] <- extract_meta_reg(
        meta_reds_macro,
        "M. franciscanus density -> M. pyrifera biomass",
        "M. franciscanus density", "M. pyrifera biomass")
    }
  }
}

# Model 3: Does lobster density increase predict purple urchin density decrease?
# Cascade prediction: negative slope (more lobsters -> fewer urchins via predation)
if (all(c("P_interruptus_Den", "S_purpuratus_Den") %in% names(es_wide))) {
  w_col <- "S_purpuratus_vi_Den"
  req_cols <- c("P_interruptus_Den", "S_purpuratus_Den", w_col)
  if (w_col %in% names(es_wide)) {
    es_lob_purp <- es_wide[complete.cases(es_wide[, req_cols]), ]
    if (nrow(es_lob_purp) >= 3) {
      meta_lob_purp <- rma(yi = es_lob_purp$S_purpuratus_Den,
                           vi = es_lob_purp[[w_col]],
                           mods = ~ P_interruptus_Den,
                           data = es_lob_purp, method = "REML",
                           test = "knha")
      table3_rows[[length(table3_rows) + 1]] <- extract_meta_reg(
        meta_lob_purp,
        "P. interruptus density -> S. purpuratus density",
        "P. interruptus density", "S. purpuratus density")
    }
  }
}

# Model 4: Does sheephead density increase predict purple urchin density decrease?
# Cascade prediction: negative slope (more sheephead -> fewer urchins via predation)
if (all(c("S_pulcher_Den", "S_purpuratus_Den") %in% names(es_wide))) {
  w_col <- "S_purpuratus_vi_Den"
  req_cols <- c("S_pulcher_Den", "S_purpuratus_Den", w_col)
  if (w_col %in% names(es_wide)) {
    es_sheep_purp <- es_wide[complete.cases(es_wide[, req_cols]), ]
    if (nrow(es_sheep_purp) >= 3) {
      meta_sheep_purp <- rma(yi = es_sheep_purp$S_purpuratus_Den,
                             vi = es_sheep_purp[[w_col]],
                             mods = ~ S_pulcher_Den,
                             data = es_sheep_purp, method = "REML",
                             test = "knha")
      table3_rows[[length(table3_rows) + 1]] <- extract_meta_reg(
        meta_sheep_purp,
        "S. pulcher density -> S. purpuratus density",
        "S. pulcher density", "S. purpuratus density")
    }
  }
}

# Model 5: Same as Model 3 but using biomass instead of density
# (lobster biomass increase -> purple urchin biomass decrease?)
if (all(c("P_interruptus_Bio", "S_purpuratus_Bio") %in% names(es_wide))) {
  w_col <- "S_purpuratus_vi_Bio"
  req_cols <- c("P_interruptus_Bio", "S_purpuratus_Bio", w_col)
  if (w_col %in% names(es_wide)) {
    es_lob_purp_bio <- es_wide[complete.cases(es_wide[, req_cols]), ]
    if (nrow(es_lob_purp_bio) >= 3) {
      meta_lob_purp_bio <- rma(yi = es_lob_purp_bio$S_purpuratus_Bio,
                               vi = es_lob_purp_bio[[w_col]],
                               mods = ~ P_interruptus_Bio,
                               data = es_lob_purp_bio, method = "REML",
                               test = "knha")
      table3_rows[[length(table3_rows) + 1]] <- extract_meta_reg(
        meta_lob_purp_bio,
        "P. interruptus biomass -> S. purpuratus biomass",
        "P. interruptus biomass", "S. purpuratus biomass")
    }
  }
}

# Model 6: Same as Model 4 but using biomass instead of density
# (sheephead biomass increase -> purple urchin biomass decrease?)
if (all(c("S_pulcher_Bio", "S_purpuratus_Bio") %in% names(es_wide))) {
  w_col <- "S_purpuratus_vi_Bio"
  req_cols <- c("S_pulcher_Bio", "S_purpuratus_Bio", w_col)
  if (w_col %in% names(es_wide)) {
    es_sheep_purp_bio <- es_wide[complete.cases(es_wide[, req_cols]), ]
    if (nrow(es_sheep_purp_bio) >= 3) {
      meta_sheep_purp_bio <- rma(yi = es_sheep_purp_bio$S_purpuratus_Bio,
                                 vi = es_sheep_purp_bio[[w_col]],
                                 mods = ~ S_pulcher_Bio,
                                 data = es_sheep_purp_bio, method = "REML",
                                 test = "knha")
      table3_rows[[length(table3_rows) + 1]] <- extract_meta_reg(
        meta_sheep_purp_bio,
        "S. pulcher biomass -> S. purpuratus biomass",
        "S. pulcher biomass", "S. purpuratus biomass")
    }
  }
}

# --- Export Table 3 to CSV ---
if (length(table3_rows) > 0) {
  Table3 <- do.call(rbind, table3_rows)
  rownames(Table3) <- NULL
  # Manuscript Table 3
  write.csv(Table3, here::here("tables", "table_03_cross_taxa_meta_regression.csv"), row.names = FALSE)
  cat("\nTable 3 (cross-taxa meta-regressions) exported to: tables/table_03_cross_taxa_meta_regression.csv\n")
  cat("  Models exported:", nrow(Table3), "\n")
  cat("  Columns:", paste(names(Table3), collapse = ", "), "\n")
} else {
  cat("\nWARNING: No Table 3 meta-regression models could be fit (insufficient cross-taxa data)\n")
}

####################################################################################################
## DONE ############################################################################################
####################################################################################################

cat("\nMeta-analysis complete.\n")
cat("\nR objects available for downstream scripts (10_temporal_analysis, 11_figures, etc.):\n")
cat("  meta_biomass_full / meta_density_full : Primary joint models (no outlier removal)\n")
cat("  meta_biomass / meta_density           : Joint models after Cook's D removal (sensitivity)\n")
cat("  Table2                                : Main results table (manuscript Table 2)\n")
cat("  Table2_pertaxa / Table2_cooksd        : Sensitivity versions of Table 2\n")
cat("\nExported CSV files (manuscript table mapping):\n")
cat("  1. tables/table_02_meta_analysis.csv                    -> Table 2 (main text)\n")
cat("  2. tables/table_03_cross_taxa_meta_regression.csv       -> Table 3 (main text)\n")
cat("  3. tables/table_s_variance_components.csv               -> Table S2\n")
cat("  4. tables/table_s_source_sensitivity_models.csv         -> Table S3a\n")
cat("  5. tables/table_s_source_sensitivity_coefficients.csv   -> Table S3b\n")
cat("  6. tables/table_s_kelp_leave1out.csv                    -> Table S8\n")
cat("  7. tables/table_s_outlier_sensitivity.csv               -> Table S9\n")
cat("  8. outputs/filter_audit_meta_analysis.csv               (audit trail, not in manuscript)\n")
cat("  9. outputs/filter_audit_pertaxa_meta.csv                (audit trail, not in manuscript)\n")
