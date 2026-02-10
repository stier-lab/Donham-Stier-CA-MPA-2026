# Scope D: Meta-Analysis & Results Summary Audit

**Files audited:**
- `/Users/adrianstier/Donham-Stier-CA-MPA-2026/code/R/09_meta_analysis.R` (908 lines)
- `/Users/adrianstier/Donham-Stier-CA-MPA-2026/code/R/11_results_summary.R` (538 lines)

**Auditor:** Scope D agent
**Date:** 2026-02-09

---

## Issue Summary

| # | Severity | Location | Title |
|---|----------|----------|-------|
| 1 | **MODERATE** | 09, L535/L581 | FDR correction applied to rounded p-values |
| 2 | **MODERATE** | 09, L344-375 | Sensitivity analysis uses REML AIC/BIC for models with different random structures |
| 3 | **MODERATE** | 09, L755-759 | Table 3 meta-regressions use pre-outlier-removal data (SumStats.Final) |
| 4 | **MODERATE** | 09, L228-231 | Pseudo-I-squared formula uses mean(vi) instead of the correct weighted estimator |
| 5 | **MODERATE** | 09, L175/L257 | Cook's distance threshold uses hardcoded 4 instead of COOKS_DISTANCE_NUMERATOR constant |
| 6 | **MINOR** | 09, L458-489 | confint() extraction assumes $random structure that may not exist in metafor rma.mv |
| 7 | **MINOR** | 09, L49/L78 | Documentation says Table 2 is written to "project root" but code writes to data/ |
| 8 | **MINOR** | 09, L759/L765 | Table 3 pivot_wider uses values_fn = mean, silently averaging duplicate MPA-taxa entries |
| 9 | **MINOR** | 09, L811-898 | Table 3 meta-regressions use rma() (univariate) not rma.mv(), ignoring Source structure |
| 10 | **MINOR** | 11, L63/L67 | SIGNIFICANCE_ALPHA dependency not validated before use |
| 11 | **MINOR** | 11, L150 | Replicate-level significance inferred from CI overlap with zero, not from formal test |
| 12 | **MINOR** | 11, L117/L160 | Uses write_csv (readr) but 09 uses write.csv (base R) -- inconsistent row.names behavior |

---

## Detailed Findings

### Issue 1: FDR correction applied to rounded p-values [MODERATE]

**File:** `09_meta_analysis.R`, lines 535 and 581

**Problem:** In `extract_meta_table()` (line 535), p-values are rounded to 4 decimal places before being stored in `Table2$pval`:

```r
pval = round(coef_table[, "pval"], 4),  # line 535
```

Then at line 581, FDR correction is applied to these already-rounded values:

```r
Table2$pval_fdr <- p.adjust(Table2$pval, method = "fdr")  # line 581
```

Rounding p-values before FDR correction introduces quantization error. For example, a true p-value of 0.00005 would be rounded to 0.0001 before entering FDR correction, potentially affecting rank ordering and adjustment magnitude. With only 9 tests this is unlikely to change conclusions drastically, but it is methodologically incorrect. The FDR correction should be applied to the full-precision p-values, and rounding should happen only at the display/export stage.

**Impact:** Could alter FDR-corrected p-values, especially near significance boundaries. For this dataset with 9 tests and mostly well-separated p-values, the practical impact is likely small but not zero.

---

### Issue 2: Sensitivity analysis uses REML AIC/BIC for models with different random structures [MODERATE]

**File:** `09_meta_analysis.R`, lines 344-375

**Problem:** The sensitivity analysis compares models with vs. without the Source random effect using AIC and BIC. Both models are fit with `method = "REML"`. However, REML likelihoods (and thus REML-based AIC/BIC) are **not comparable** between models with different fixed-effects structures. While the fixed effects are the same here (`~ Taxa - 1`), the more subtle issue is that REML AIC/BIC for comparing random effects structures is acceptable ONLY when the fixed effects are identical -- which they are. So this is actually borderline acceptable.

However, the standard recommendation for comparing nested random effects structures is a likelihood ratio test using ML (not REML) estimation. The code does not perform an LRT, and the use of REML AIC/BIC for this comparison, while not strictly wrong when fixed effects match, is not the gold-standard approach. Using ML-based AIC or a formal LRT via `anova()` would be more rigorous.

**Impact:** The comparison results may be directionally correct but the AIC/BIC values are REML-based, which adds some uncertainty to the model selection conclusion. This is a defensible but suboptimal approach.

---

### Issue 3: Table 3 meta-regressions use pre-outlier-removal data [MODERATE]

**File:** `09_meta_analysis.R`, lines 755-759

**Problem:** The Table 3 cross-taxa meta-regressions (trophic cascade tests) pivot from `SumStats.Final` -- the original data that entered the meta-analysis, **not** the post-outlier-removal data (`biomass_clean`, `density_clean`). This means observations flagged as influential outliers by Cook's distance and removed from the Table 2 meta-analysis are still included in the Table 3 meta-regressions.

```r
es_wide_mean <- SumStats.Final %>%  # Uses ALL data, not outlier-cleaned data
  dplyr::select(Taxa, MPA, Mean, Resp) %>%
  ...
```

This is an inconsistency: Table 2 results reflect outlier-removed data, but Table 3 results include those outliers. Whether this is intentional (the meta-regressions are a separate analysis) or unintentional is ambiguous. If the outliers were removed for valid statistical reasons, they should arguably be excluded from all downstream analyses.

**Impact:** Trophic cascade relationship estimates may be influenced by the same outliers that were deemed problematic for the main meta-analysis.

---

### Issue 4: Pseudo-I-squared formula uses simple mean of vi [MODERATE]

**File:** `09_meta_analysis.R`, lines 228-231, 292-295

**Problem:** The pseudo-I-squared calculation uses `mean(biomass_clean$vi)` as the "typical within-study variance":

```r
typical_v_bio <- mean(biomass_clean$vi)
pseudo_I2_bio <- 100 * total_hetero_bio / (total_hetero_bio + typical_v_bio)
```

For multilevel meta-analysis models (`rma.mv`), the standard I-squared formula from Higgins & Thompson does not directly apply. The code correctly labels this as "Pseudo-I-squared," but the formula `tau2 / (tau2 + mean(vi))` is a rough approximation. The more appropriate approach for `rma.mv` models would be to use the formula from Nakagawa & Santos (2012) or the one implemented in `metafor::fitstats()`, which uses a weighted average of vi (weighting by the inverse of vi, i.e., the precision-weighted typical variance).

Using the unweighted mean of vi gives outsized influence to imprecise studies with large sampling variances, inflating the denominator and potentially underestimating I-squared.

**Impact:** The reported pseudo-I-squared values may be biased, though the direction and magnitude depend on the variance distribution. The interpretation thresholds (low/moderate/high) are arbitrary guidelines anyway, so this is unlikely to change qualitative conclusions but is technically imprecise.

---

### Issue 5: Hardcoded Cook's distance numerator [MODERATE]

**File:** `09_meta_analysis.R`, lines 175, 257

**Problem:** The script uses hardcoded `4` in the Cook's distance threshold calculation:

```r
cooks_threshold_bio <- 4 / n_bio   # line 175
cooks_threshold_den <- 4 / n_den   # line 257
```

The project already defines `COOKS_DISTANCE_NUMERATOR <- 4` in `00c_analysis_constants.R` (line 129) specifically for this purpose. The hardcoded value works correctly but defeats the purpose of the centralized constants file. If someone changes the constant, the meta-analysis script would not reflect the change.

**Impact:** No functional bug currently, but violates the project's convention of centralizing magic numbers. Creates a maintenance risk.

---

### Issue 6: confint() extraction assumes $random list structure [MINOR]

**File:** `09_meta_analysis.R`, lines 458-489

**Problem:** The code checks for `"random" %in% names(ci_biomass)` before extracting variance component CIs. However, for `rma.mv` objects, `confint()` returns a list structure that depends on the metafor version. In recent versions of metafor (>= 3.0), `confint()` on an `rma.mv` object returns a list indexed numerically (e.g., `ci_biomass[[1]]` for the first random effect), not with a `$random` component.

The code has a comment acknowledging this: "structure depends on metafor version" (line 459). The `if ("random" %in% names(...))` guard means that if the structure does not have a `$random` field, the code silently skips extraction and produces an empty `tau2_summary` dataframe. This is a graceful failure, but it means the variance component CI table may silently not be generated.

**Impact:** The variance component CI summary table (`table_s_variance_components.csv`) may be empty in certain metafor versions. The `tryCatch` wrapper on `confint()` itself (lines 424-436) handles errors, but the structural mismatch would cause silent data loss.

---

### Issue 7: Documentation says Table 2 written to "project root" [MINOR]

**File:** `09_meta_analysis.R`, lines 49, 78

**Problem:** The header comments state:
- Line 49: `table_02_meta_analysis.csv: Written to project root`
- Line 78: `Exported to project root as "table_02_meta_analysis.csv"`

But the actual code (line 596) writes to `here::here("data", "table_02_meta_analysis.csv")`, which is the `data/` subdirectory. This matches the CLAUDE.md convention but contradicts the header documentation.

**Impact:** Documentation inconsistency only. No functional issue.

---

### Issue 8: pivot_wider silently averages duplicate MPA-taxa entries [MINOR]

**File:** `09_meta_analysis.R`, lines 759, 765

**Problem:** When creating the wide-format data for Table 3 meta-regressions:

```r
tidyr::pivot_wider(names_from = Taxa_Resp, values_from = Mean, values_fn = mean)
```

The `values_fn = mean` argument silently averages effect sizes when multiple entries exist for the same MPA-taxa-response combination (e.g., when both PISCO and KFM have data for the same MPA). This is a reasonable approach but:

1. It gives equal weight to effect sizes regardless of their precision (SE)
2. It silently collapses data without informing the user how many entries were averaged
3. The corresponding SE averaging (`values_fn = mean` for SE) is incorrect -- averaging SEs is not the proper way to combine them. The correct approach would be inverse-variance weighting or propagating uncertainty.

**Impact:** The averaged SEs used as weights in the Table 3 meta-regressions (`vi = SE^2`) are incorrect when multiple sources are averaged. This affects the precision of trophic cascade relationship estimates but does not bias the point estimates.

---

### Issue 9: Table 3 meta-regressions use rma() instead of rma.mv() [MINOR]

**File:** `09_meta_analysis.R`, lines 811-898

**Problem:** The Table 3 trophic cascade meta-regressions use `metafor::rma()` (standard random-effects meta-regression) rather than `rma.mv()`. This means they do not include the Source random effect that the main models use. The code comments (lines 745-748) justify using `rma()` for its ability to weight by sampling variance, but this loses the multilevel structure.

Since the data has been collapsed to MPA-level (one row per MPA), there is only one observation per MPA, making the MPA random effect unnecessary. However, the Source random effect is still relevant if multiple sources contributed to the averaged MPA-level estimate.

**Impact:** Minor inconsistency in modeling approach between Table 2 (multilevel) and Table 3 (standard). Given the small k values and the MPA-level aggregation, this is defensible but should be noted.

---

### Issue 10: SIGNIFICANCE_ALPHA used without existence check [MINOR]

**File:** `11_results_summary.R`, lines 63, 67, 191

**Problem:** The `extract_model_summary()` function and markdown generation use `SIGNIFICANCE_ALPHA` without checking whether it exists. If `00c_analysis_constants.R` was not sourced (e.g., running `11_results_summary.R` in isolation or after an error in the pipeline), this would cause an error.

The script does check for `meta_biomass`, `meta_density`, `Table2`, and `SumStats.Final` with `exists()` guards, but `SIGNIFICANCE_ALPHA` has no such guard.

**Impact:** Would cause a runtime error if constants are not loaded. In normal pipeline execution this is never triggered.

---

### Issue 11: Replicate-level significance based on CI overlap, not formal test [MINOR]

**File:** `11_results_summary.R`, line 150

**Problem:** Replicate-level significance is determined by:

```r
Significant = (CI_Lower > 0 | CI_Upper < 0)
```

This tests whether the confidence interval excludes zero. However, the CI column from `SumStats.Final` contains a mix of t-distribution-based CIs (from `qt(0.975, df)`) and z-distribution-based CIs (from `1.96 * SE`), depending on the model type used in `08_effect_sizes.R`. For pBACIPS models using `emmeans::pairs()`, the CI is t-based. For NLS models (sigmoid, asymptotic), the CI uses z = 1.96.

This inconsistency means the "significance" determination uses different Type I error rates across different effect sizes. For small-sample cases, the z-based CI is anticonservative (narrower than it should be).

**Impact:** Some replicate-level significance calls in the appendix table may be slightly inconsistent, but this affects only the summary display, not the formal meta-analysis results.

---

### Issue 12: Inconsistent CSV writing functions [MINOR]

**File:** `11_results_summary.R`, lines 117, 160 vs. `09_meta_analysis.R`, line 596

**Problem:** `11_results_summary.R` uses `readr::write_csv()` (lines 117, 160) while `09_meta_analysis.R` uses base `write.csv()` (line 596). The key difference is that `write.csv()` includes row names by default (the code explicitly sets `row.names = FALSE`), while `write_csv()` never includes row names. If anyone removes the `row.names = FALSE` from the `write.csv()` calls, the CSV format would break.

**Impact:** No functional issue currently, but inconsistent coding style across the two scripts.

---

## Items Verified as Correct

1. **Meta-analysis model specification:** `rma.mv()` is correctly specified with `yi = Mean`, `V = vi`, `mods = ~ Taxa - 1` (cell-means parameterization), `random = list(~ 1 | MPA, ~ 1 | Source)` (crossed random effects), `method = "REML"`, and `test = "t"`. This is appropriate for the study design.

2. **Biomass/density separation:** Running separate models for biomass and density is correctly justified by the non-independence of these measures from the same surveys.

3. **Cook's distance workflow:** Outlier detection is correctly applied in the right order: (1) fit initial model, (2) compute Cook's distance, (3) remove outliers, (4) refit. The 4/n threshold is standard.

4. **FDR correction scope:** The Benjamini-Hochberg correction is applied across all 9 taxa-response tests in Table 2 jointly (not separately by response type). This is appropriate because the tests share common random effects and the family-wise interpretation spans the full table.

5. **Low-k flagging:** Taxa with k < 5 are correctly flagged as "preliminary (k<5)" in Table 2 (line 572).

6. **Table 2 generation:** The `extract_meta_table()` function correctly pulls estimates, SEs, test statistics, p-values, and CIs from the model coefficient table and matches them to the clean (outlier-removed) data.

7. **Results summary extraction:** `11_results_summary.R` dynamically extracts all statistics from model objects -- no hardcoded result values were found. All numbers in the markdown summary are generated programmatically.

8. **Heterogeneity tau-squared:** The `sigma2` components from `rma.mv()` are correctly reported as tau-squared values for each random effect (MPA and Source).

9. **Density model correctly omits M. pyrifera:** The RESULTS_SUMMARY.md shows 4 taxa in the density table (no kelp density, which makes biological sense -- kelp biomass is the relevant metric). This is not a bug; it reflects that M. pyrifera density is not measured/analyzed.

---

## Recommendations (Priority Order)

1. **Fix FDR rounding (Issue 1):** Store full-precision p-values for FDR correction, round only for display.
2. **Use ML for sensitivity analysis comparison (Issue 2):** Refit models with `method = "ML"` for the AIC/BIC comparison, or add a formal LRT.
3. **Clarify outlier treatment in Table 3 (Issue 3):** Either use cleaned data or explicitly document that Table 3 uses full data by design.
4. **Use constant for Cook's threshold (Issue 5):** Replace `4` with `COOKS_DISTANCE_NUMERATOR`.
5. **Fix SE averaging in Table 3 pivot (Issue 8):** Use inverse-variance weighted average for SEs.
