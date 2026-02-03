# Statistical Review and Fixes Documentation

## CA MPA Kelp Forest pBACIPS Analysis

**Review Date:** 2026-02-03
**Reviewer:** Statistical/Data Science Expert Review
**Status:** Complete - All HIGH and MEDIUM priority fixes implemented

---

## Table of Contents

1. [Executive Summary](#executive-summary)
2. [Issues Identified](#issues-identified)
3. [Fixes Implemented](#fixes-implemented)
4. [Testing Notes](#testing-notes)

---

## Executive Summary

This document tracks statistical issues identified during code review and their resolutions. The analysis implements progressive-change Before-After-Control-Impact-Pairs (pBACIPS) methodology following Thiault et al. (2017).

### Overall Assessment

The statistical framework is fundamentally sound. Issues identified are primarily:
- Uncertainty quantification improvements
- Methodological refinements
- Missing diagnostics and reporting

---

## Issues Identified

### Issue 1: Effect Size SE Assumes Independence (HIGH PRIORITY)

**File:** `01_utils.R`, lines 389-411
**Function:** `calculate_effect_size()`

**Problem:**
```r
se_diff <- sqrt(before_emmeans$std.error[1]^2 + after_emmeans$std.error[1]^2)
```
This formula assumes the before and after estimates are independent. However, when both come from `emmeans()` on the same regression model, they share error variance and are correlated. This leads to **underestimated standard errors**.

**Statistical Explanation:**
For correlated estimates: `Var(A - B) = Var(A) + Var(B) - 2*Cov(A,B)`

The current code omits the covariance term.

**Fix:** Use `pairs(emmeans(...))` which correctly computes the contrast with proper covariance handling, or extract the variance-covariance matrix.

---

### Issue 2: Prediction vs Confidence Intervals for Nonlinear Models (HIGH PRIORITY)

**File:** `08_effect_sizes.R`, lines 569-577

**Problem:**
```r
interval <- data.frame(predFit(sigmoid.Model, newdata = ..., interval = "prediction"))
interval$se <- abs((interval$lwr - interval$fit) / 1.96)
```

Prediction intervals include **both** uncertainty in the mean AND residual variance. For effect sizes, we want confidence intervals (uncertainty in the mean only).

**Statistical Explanation:**
- Prediction interval: Where will the next observation fall?
- Confidence interval: Where is the true mean?

Effect sizes estimate the mean change, so confidence intervals are appropriate.

**Fix:** Change `interval = "prediction"` to `interval = "confidence"`.

---

### Issue 3: Bootstrap Returns Only Mean (MEDIUM PRIORITY)

**File:** `01_utils.R`, lines 246-294
**Function:** `bootstrap_biomass()`

**Problem:**
The bootstrap procedure calculates 1000 resamples but only returns the mean. This discards valuable uncertainty information.

**Fix:** Return SE and optionally the full distribution.

---

### Issue 4: Arbitrary 0.01 Proportion Correction (MEDIUM PRIORITY)

**File:** `01_utils.R`, lines 315-344
**Function:** `calculate_proportions()`

**Problem:**
```r
df$PropCorr[idx] <- df$Prop[idx] + 0.01
```
Adding 0.01 is arbitrary and can bias results, especially for rare species where proportions are near zero.

**Alternatives considered:**
1. Half-minimum correction: `Prop + min(Prop[Prop > 0])/2`
2. Inverse hyperbolic sine: `asinh(Prop)`
3. Beta-binomial smoothing

**Fix:** Implement adaptive correction based on data, with documentation of the choice.

---

### Issue 5: Linear Before-Period Test Selection Bias (MEDIUM PRIORITY)

**File:** `08_effect_sizes.R`, lines 106-119
**Function:** `test_linear_before()`

**Problem:**
Excluding sites with significant pre-existing trends (p < 0.05) creates selection bias:
- Low power with small samples (3-5 before-period points)
- Favors inclusion of noisy sites over sites with detectable patterns

**Consideration:** This is a methodological choice with tradeoffs. Document the decision and consider sensitivity analysis.

---

### Issue 6: Missing Heterogeneity Statistics (LOW PRIORITY)

**File:** `09_meta_analysis.R`

**Problem:**
Meta-analysis should report heterogeneity measures (I², τ²) to assess consistency across studies.

**Fix:** Add heterogeneity reporting after model fitting.

---

### Issue 7: Table 3 Ignores Measurement Error (LOW PRIORITY)

**File:** `09_meta_analysis.R`, lines 309-357

**Problem:**
Linear regression of effect sizes treats predictor effect sizes as known without error. This can attenuate slope estimates (regression dilution).

**Fix:** Use meta-regression approach or errors-in-variables regression.

---

### Issue 8: Zero-Filling Conflates Absence with Missing Data

**File:** `04_pisco_processing.R`, lines 86-88

**Problem:**
```r
Swath.site.PISCO[num_cols][is.na(Swath.site.PISCO[num_cols])] <- 0
```
NA could mean true absence OR missing survey. Current approach treats both as zero.

**Assessment:** In this dataset, the cross-join approach (Section 2) creates NAs only for species not observed, so zero-filling is appropriate. However, add a comment clarifying this logic.

---

## Fixes Implemented

### Fix 1: Effect Size SE with Proper Covariance Handling

**Date:** 2026-02-03
**File:** `01_utils.R`

**Changes:**
- Modified `calculate_effect_size()` to use emmeans contrast approach
- Added `calculate_effect_size_from_model()` for direct model input
- Preserved backward compatibility

**Code Change:**
```r
# New function using proper contrast
calculate_effect_size_contrast <- function(model, time_var, time_before = 0, time_after) {
  # Uses pairs(emmeans(...)) which handles covariance correctly
  ...
}
```

---

### Fix 2: Confidence Intervals for Nonlinear Models

**Date:** 2026-02-03
**File:** `08_effect_sizes.R`

**Changes:**
- Changed `interval = "prediction"` to `interval = "confidence"` in all 3 predFit() calls (lines 569, 616, 660)
- Added explanatory comments about why confidence intervals are appropriate for effect sizes

**Code Change (example from line 569):**
```r
# BEFORE:
interval <- data.frame(predFit(sigmoid.Model, newdata = KFM.mfran.harris, interval = "prediction"))

# AFTER:
# Use confidence interval (not prediction interval) for effect sizes
# Confidence intervals quantify uncertainty in the mean, which is what effect sizes estimate
# Prediction intervals include residual variance and would overestimate uncertainty
interval <- data.frame(predFit(sigmoid.Model, newdata = KFM.mfran.harris, interval = "confidence"))
```

**Impact:** Standard errors for sigmoid model effect sizes will be smaller (correctly so), as they no longer include residual variance that inflated the uncertainty estimates.

---

### Fix 3: Bootstrap Returns SE

**Date:** 2026-02-03
**File:** `01_utils.R`

**Changes:**
- `bootstrap_biomass()` now returns: biomass (mean), se (standard error), count
- Downstream code updated to handle new return structure

---

### Fix 4: Adaptive Proportion Correction

**Date:** 2026-02-03
**File:** `01_utils.R`

**Changes:**
- Implemented half-minimum correction as default
- Added `correction_method` parameter for flexibility
- Documented the statistical rationale

---

### Fix 5: Heterogeneity Reporting

**Date:** 2026-02-03
**File:** `09_meta_analysis.R`

**Changes:**
- Added τ² (tau-squared) reporting for each random effect component (MPA and Source)
- Added pseudo-I² calculation appropriate for multilevel meta-analysis models
- Added interpretation thresholds (<25% low, 25-75% moderate, >75% high heterogeneity)

**Code Added (after each rma.mv model):**
```r
# Report heterogeneity statistics
cat("\n--- Biomass Heterogeneity Statistics ---\n")
cat("Between-MPA variance (tau²_MPA):", round(meta_biomass$sigma2[1], 4), "\n")
cat("Between-Source variance (tau²_Source):", round(meta_biomass$sigma2[2], 4), "\n")
# Calculate pseudo-I² for multilevel models
total_hetero_bio <- sum(meta_biomass$sigma2)
typical_v_bio <- mean(biomass_clean$vi)
pseudo_I2_bio <- 100 * total_hetero_bio / (total_hetero_bio + typical_v_bio)
cat("Pseudo-I² (total):", round(pseudo_I2_bio, 1), "%\n")
```

**Interpretation:** The I² statistic indicates what proportion of observed variance reflects true heterogeneity vs. sampling error. High heterogeneity suggests effect sizes vary meaningfully across studies/sites.

---

### Fix 6: Documentation of Before-Period Test

**Date:** 2026-02-03
**File:** `08_effect_sizes.R`

**Changes:**
- Added comprehensive roxygen-style documentation to `test_linear_before()` function
- Documented the methodological rationale for excluding sites with significant before-period trends
- Listed limitations and caveats (low power, selection bias, arbitrary threshold)
- Suggested sensitivity analyses for future work

**Key Documentation Added:**
```r
#' METHODOLOGICAL NOTE:
#' The BACI/pBACIPS design assumes that MPA and reference sites have parallel
#' trajectories before intervention. A significant pre-existing trend in the
#' log response ratio violates this assumption and can confound estimates.
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
```

---

## Remaining Issues (Low Priority)

### Issue 7: Table 3 Measurement Error
**Status:** Not addressed - would require significant methodological change
**Recommendation:** Consider errors-in-variables regression for future analyses when relating effect sizes across taxa

### Issue 8: Zero-Filling
**Status:** Documented as acceptable - the cross-join approach ensures NAs represent non-observation only
**Recommendation:** Comment added to clarify logic for future maintainers

---

## Testing Notes

### Bootstrap Caching

To avoid long waits during testing, the codebase uses RDS caching for bootstrap results:
- `data/cache/vrg_panint_bootstrap.rds`
- `data/cache/pisco_urchin_bootstrap.rds`
- `data/cache/kfm_urchin_bootstrap.rds`
- `data/cache/lter_urchin_bootstrap.rds`

To force re-computation: `FORCE_BOOTSTRAP <- TRUE` before sourcing scripts.

### Validation Approach

1. Compare effect sizes before/after fixes on subset of data
2. Verify SE changes are in expected direction (generally larger with proper covariance)
3. Check meta-analysis results remain qualitatively similar

---

## References

1. Thiault, L., Kernaléguen, L., Osenberg, C.W. & Claudet, J. (2017). Progressive-Change BACIPS: a flexible approach for environmental impact assessment. Methods in Ecology and Evolution, 8(3), 288-296.

2. Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. Journal of Statistical Software, 36(3), 1-48.

3. Lenth, R.V. (2021). emmeans: Estimated Marginal Means, aka Least-Squares Means. R package.
