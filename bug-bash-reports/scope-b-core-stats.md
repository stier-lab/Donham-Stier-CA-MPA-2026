# Scope B: Core Statistics Audit Report

**Files audited:**
- `/Users/adrianstier/Donham-Stier-CA-MPA-2026/code/R/02_pBACIPS_function.R` (1,804 lines)
- `/Users/adrianstier/Donham-Stier-CA-MPA-2026/code/R/08_effect_sizes.R` (1,710 lines)

**Audit date:** 2026-02-09
**Auditor:** Scope B agent (Core Statistics)

---

## Executive Summary

The core statistics module is generally well-engineered with extensive fallback chains for NLS convergence, proper covariance-aware SE estimation via `emmeans::pairs()`, and thorough input validation. However, several issues were identified ranging from a critical variable-naming mismatch that likely causes `predFit` failures on sigmoid models, to moderate concerns about AICc comparison fairness when fallback models replace the intended nonlinear forms, and a semantic ambiguity in `time.model` coding that merges the first "After" year with all "Before" years.

**Summary counts:** 3 CRITICAL, 7 MODERATE, 5 MINOR

---

## CRITICAL Issues

### C1. `predFit` newdata variable name mismatch for sigmoid models

**Location:** `08_effect_sizes.R` lines 1077-1084, 1161-1168, 1240-1247
**Also relevant:** `02_pBACIPS_function.R` lines 1488-1492 (`mySIGfun_standalone`)

**Description:**
The sigmoid standalone function fits the model using a local variable `time_offset` (defined as `time.model + 0.01`):
```r
time_offset <- time.model + 0.01
foSIG <- delta ~ (M * (time_offset / L)^K) / (1 + (time_offset / L)^K) + B
```

The NLS model formula therefore expects a variable called `time_offset` in any `newdata`. However, `08_effect_sizes.R` creates prediction data with columns named `time.model` and `time.true`:
```r
pred_data <- data.frame(
  time.model = c(0, EFFECT_SIZE_TIME_YEARS),
  time.true = c(...)
)
interval <- data.frame(predFit(sigmoid.Model, newdata = pred_data, interval = "confidence"))
```

`predFit` from the `investr` package requires that the `newdata` contain a column matching the predictor variable name used in the model formula. Since the model was fit with `time_offset`, it needs a `time_offset` column, not `time.model`. This would cause `predFit` to either error (caught by `tryCatch` returning `NULL`) or, if it falls through to a different code path, produce incorrect predictions.

**Impact:** If `predFit` fails silently (returns NULL), the sigmoid effect size row is simply skipped, and the analysis may be missing effect sizes for Harris Point M. franciscanus, Scorpion M. pyrifera, and Gull Island P. interruptus. If these are critical data points, results could be incomplete. Moreover, even if `predFit` succeeds (e.g., if the model was returned from a fallback that uses different variable names), the predictions at `time.model = 0` would not include the 0.01 offset, meaning the "before" prediction is evaluated at `time_offset = 0` rather than `0.01`, introducing a small numerical discrepancy.

**Additionally:** When the sigmoid model falls through to a fallback (piecewise linear, quadratic, GAM), those fallbacks use `time` as the predictor, not `time_offset` and not `time.model`. This compounds the mismatch.

### C2. `warnOnly = TRUE` allows non-converged NLS models to propagate

**Location:** `02_pBACIPS_function.R` -- 12 occurrences (lines 324, 371, 508, 584, 926, 973, 1103, 1154, 1401, 1446, 1547, 1599)

**Description:**
The `nls.control(warnOnly = TRUE)` setting tells R to return the model object even when the algorithm has not converged, issuing only a warning rather than an error. These non-converged models:
- Have unreliable parameter estimates
- Have unreliable residual variance (which feeds into AICc)
- Can "win" AICc competition against legitimately converged models

There is no post-fit convergence check (e.g., checking `model$convInfo$isConv` or the convergence tolerance). The code only checks that coefficients are finite and within bounds, but a non-converged model can have finite coefficients that are simply at whatever values the optimizer stopped at.

**Impact:** A non-converged model could be selected as the "best" model by AICc, producing incorrect effect size estimates. Since these models are then used for prediction (via `predFit`) and SE estimation, the resulting effect sizes and their uncertainties could be arbitrarily wrong.

### C3. SumStats type coercion via `c()` for sigmoid rows silently converts numeric to character

**Location:** `08_effect_sizes.R` lines 1106-1108, 1188-1190, 1267-1269

**Description:**
The sigmoid effect size rows are added using:
```r
SumStats[nrow(SumStats) + 1, ] <- c("M. franciscanus", "Harris Point SMR", mean_es, pSE, pSD, pCI,
                                      "Sigmoid", "KFM", "Bio", "Y", "Y", "pBACIPS", "N", n_obs_original,
                                      NA, NA)
```

The `c()` function coerces all elements to a single type. Since string elements are present, ALL values (including `mean_es`, `pSE`, `pSD`, `pCI`, `n_obs_original`) are coerced to character strings. While the later `as.numeric()` conversion at lines 1336-1341 recovers the numeric values, this intermediate character representation means:
- If `mean_es` is a floating point number like `0.123456789`, the string representation may lose precision
- The row assignment `SumStats[nrow(SumStats) + 1, ]` can silently convert the entire dataframe to character if the existing columns had been typed differently
- This is fragile and error-prone compared to the `data.frame()` approach used by the helper functions

**Impact:** Potential precision loss in effect sizes for the three sigmoid-fitted MPA-taxa combinations. The helper functions (`add_step_effect_size`, `add_linear_effect_size`) correctly use `data.frame()` construction and `rbind()`, which preserves types. The sigmoid code path uses a different, less robust pattern.

---

## MODERATE Issues

### M1. `time.model` coding: first "After" year gets value 0, same as all "Before" years

**Location:** `08_effect_sizes.R` line 309

**Description:**
```r
dat$time.model[idx] <- c(rep(0, n_before), seq(0, n_after - 1))
```
The "Before" period gets all zeros, and the "After" period starts at 0 and counts up. This means the first year of the "After" period has `time.model = 0`, which is identical to the "Before" period values.

**Implications:**
1. In the step model (`delta ~ period`), this is irrelevant since `period` is a factor ("Before"/"After").
2. In the linear model (`delta ~ time.model`), the first After observation has the same predictor value (0) as all Before observations. This means the model cannot distinguish the first After year from the Before period based on `time.model` alone.
3. In the sigmoid/asymptotic models, the first After data point effectively contributes to the "baseline" estimate, not the "change" estimate.
4. The `ProgressiveChangeBACIPS` function uses `time.model == 0` to identify Before period (line 865: `period <- ifelse(time.model == 0, "Before", "After")`), so the first After year is actually classified as "Before" in the step model's period factor.

**Impact:** The first After observation is misclassified as "Before" in the step model. This dilutes the estimated MPA effect by mixing one post-MPA observation into the control group. For pBACIPS with the linear/asymptotic/sigmoid models, the first After year is treated as baseline, which underestimates early MPA effects. This is a systematic bias across all pBACIPS analyses.

### M2. AICc comparison is unfair when fallback models replace nonlinear models

**Location:** `02_pBACIPS_function.R` lines 1022-1047, 1203-1228 (fallback chains)

**Description:**
When the asymptotic NLS model fails to converge, the code substitutes fallback models in this order: stable parameterization, self-starting, piecewise linear, quadratic, GAM. Similarly for sigmoid. These fallback models are then entered into AICc competition as if they were the original model class.

The problem is:
1. A piecewise linear model (3 parameters) competing as the "asymptotic" slot has different model complexity than a true Michaelis-Menten (3 parameters, but nonlinear)
2. A GAM model's effective degrees of freedom are not directly comparable via AICc to parametric models
3. A quadratic model (3 parameters) is fundamentally different from an asymptotic model

The AICc comparison assumes all models are fit to the same data with comparable likelihoods. Mixing model classes invalidates this assumption.

**Impact:** Model selection could be biased toward or against the "asymptotic" or "sigmoid" slots depending on which fallback model was substituted. This doesn't break the pipeline but could lead to incorrect model selection, which affects which effect size estimation path is used.

### M3. NLS SE estimation for sigmoid assumes z-distribution (1.96) instead of t-distribution

**Location:** `08_effect_sizes.R` lines 1092-1098, 1176-1180, 1255-1260

**Description:**
For sigmoid models, the SE is back-calculated from `predFit` confidence intervals:
```r
se_before <- abs((interval$lwr[1] - interval$fit[1]) / 1.96)
se_after <- abs((interval$lwr[2] - interval$fit[2]) / 1.96)
pSE <- sqrt(se_before^2 + se_after^2)
```
This divides by 1.96 (z-critical), but `predFit` with `interval = "confidence"` uses a t-distribution with `n - p` degrees of freedom. For small samples (e.g., n=10, p=4, df=6), `qt(0.975, 6) = 2.447`, which means the true SE is `CI_width / (2 * 2.447)`, not `CI_width / (2 * 1.96)`. Dividing by 1.96 when the actual critical value is 2.447 would **overestimate** the SE by a factor of ~1.25.

Additionally, the independence assumption (`pSE = sqrt(se_before^2 + se_after^2)`) ignores the covariance between predictions at t=0 and t=11 from the same model, which inflates the SE further.

**Impact:** The SE for sigmoid effect sizes is systematically overestimated, making these estimates appear less precise than they actually are. In the meta-analysis (which weights by `1/SE^2`), these observations receive less weight than they should. The code comments acknowledge both biases and note they "approximately cancel," but for small samples the t-distribution bias alone can be substantial (25%+ inflation). This is conservative but not methodologically correct.

### M4. `run_pbacips_loop` p-value extraction assumes specific coefficient counts from fallback models

**Location:** `08_effect_sizes.R` lines 422-431

**Description:**
```r
psi <- if (!is.null(result$sigmoid)) summary(result$sigmoid)$coefficients[, 4] else rep(NA, 4)
pa <- if (!is.null(result$asymptotic)) summary(result$asymptotic)$coefficients[, 4] else rep(NA, 3)
c(ps, pl, psi, pa)
```
The code expects the sigmoid model to have exactly 4 coefficients and the asymptotic to have exactly 3. But when fallback models are used:
- Piecewise linear: 3 coefficients (intercept, after, time_after) -- matches asymptotic but not sigmoid
- Quadratic: 3 coefficients -- matches asymptotic but not sigmoid
- GAM: `summary()` returns a different structure entirely, and `$coefficients[, 4]` may not exist or may have different dimensions
- Self-starting models (SSasymp, SSlogis): 3 parameters each, not 4

When `psi` extracts a 3-element vector instead of 4, or `pa` gets a different count, the total length of `pvals` will not be 9. The `rbind` into `pvals.poly` (initialized with 9 columns) would then trigger the "numbers of columns of arguments do not match" warning.

**Impact:** This is the **root cause of the previously flagged rbind column mismatch warning**. When fallback models are substituted, the p-value vector has the wrong length, causing R to recycle values or produce a warning. However, `pvals.poly` is not used downstream for effect size calculations (it's returned from `run_pbacips_loop` but never referenced again), so this is a diagnostic-only issue that does not affect final results.

### M5. `Norm.poly` initialized with 9 columns but receives 1 value per rbind

**Location:** `08_effect_sizes.R` lines 391, 441

**Description:**
```r
Norm.poly <- data.frame(matrix(ncol = 9, nrow = 0))
# ...
Norm.poly <- rbind(Norm.poly, norm_p)
```
`Norm.poly` is initialized as a 9-column dataframe, but each iteration appends a single scalar (`norm_p`). R's `rbind` will recycle the scalar across all 9 columns, producing a row of 9 identical values. Like `pvals.poly`, this is not used downstream but indicates sloppy initialization.

**Impact:** No effect on results since `Norm.poly` is only returned from the function and never consumed. But it indicates the dataframe initialization was copied from the `pvals.poly` pattern without adjusting the column count.

### M6. Inconsistent data source used for KFM urchin analyses: `All.RR.sub` vs `All.RR.sub.trans`

**Location:** `08_effect_sizes.R` lines 973, 990, 1000, 1030 vs 1117, 1199, 1288

**Description:**
KFM urchin density and biomass analyses (S. purpuratus and M. franciscanus) use `All.RR.sub`:
```r
KFM.purps.den <- subset(All.RR.sub, source == "KFM" & y == "Strongylocentrotus purpuratus" & resp == "Den")
KFM.reds.den <- subset(All.RR.sub, source == "KFM" & y == "Mesocentrotus franciscanus" & resp == "Den")
```

But KFM kelp, lobster, and sheephead analyses use `All.RR.sub.trans`:
```r
KFM.macro <- subset(All.RR.sub.trans, source == "KFM" & y == "Macrocystis pyrifera" & resp == "Bio")
KFM.lob <- subset(All.RR.sub.trans, source == "KFM" & y == "Panulirus interruptus")
```

All LTER and PISCO analyses use `All.RR.sub.trans`. The difference between these two datasets (which is created in `07_combine_data.R`) likely involves some transformation or filtering. Using `All.RR.sub` for urchins may be intentional (e.g., urchin data doesn't need the transformation) or may be an oversight.

**Impact:** If `All.RR.sub.trans` applies a transformation that the urchin data should also receive, the KFM urchin effect sizes would be calculated on a different scale than other taxa. Without seeing `07_combine_data.R`, this cannot be definitively classified as a bug vs. intentional, but the inconsistency warrants investigation.

### M7. `preprocess_for_nls` winsorizes data but winsorized values are not used in all code paths

**Location:** `02_pBACIPS_function.R` lines 172-209, 906-907, 1077-1078

**Description:**
The `preprocess_for_nls` function winsorizes extreme delta values (caps at 3 SD) and returns processed data. The `nls.lm` code path uses `processed$delta` and `processed$time`, but the `nls()` port algorithm code path uses the original `delta` and `time.model` variables from the enclosing scope:
```r
# Uses original data (not preprocessed):
nls(foAsy, start = parStart, algorithm = "port", ...)
# Uses preprocessed data:
nls.lm(par = parStart, fn = residFun, observed = processed$delta, time.model = processed$time, ...)
```

This means the port algorithm attempts are fit to the original (non-winsorized) data, while the LM fallback uses winsorized data. If the port algorithm succeeds on the original data, the model will be influenced by extreme values that the preprocessing was designed to mitigate.

**Impact:** Inconsistency in outlier handling between fitting strategies. The effect depends on how many extreme values exist in the data. In practice, the port algorithm may succeed more often precisely because it doesn't winsorize, but the model quality may be lower.

---

## MINOR Issues

### m1. `df` column in `AICc.test` is calculated as `length(coef(m)) + 1`, which is not standard for all model classes

**Location:** `02_pBACIPS_function.R` lines 1256-1258

**Description:**
For `aov` models, `length(coef(m))` returns the number of group means (not the number of estimated parameters). For GAM fallbacks, `coef()` returns parametric coefficients only, not the effective degrees of freedom of smooth terms. This `df` column is only used for display (not for AICc calculation, which is handled by `AICc(m)`), but it could be misleading in diagnostic output.

### m2. Global state mutation in `log_model_fit` and `MODEL_FIT_LOG`

**Location:** `02_pBACIPS_function.R` lines 73-99

**Description:**
The `MODEL_FIT_LOG` is a global list modified via `<<-`. If the pipeline is sourced multiple times in the same session, the log accumulates entries from all runs without clearing. The `clear_model_fit_log()` function exists but is not called by `run_all.R`.

### m3. `run_dharma_diagnostics` is defined in both `02_pBACIPS_function.R` (line 647) and `08_effect_sizes.R` (line 102) with different signatures

**Location:** `02_pBACIPS_function.R` line 647, `08_effect_sizes.R` line 102

**Description:**
Two functions with the same name but different signatures:
- `02_pBACIPS_function.R`: `run_dharma_diagnostics(model, plot = FALSE)` -- returns a list
- `08_effect_sizes.R`: `run_dharma_diagnostics(model, taxa, mpa, source, model_type)` -- returns a dataframe

Since `08_effect_sizes.R` sources `02_pBACIPS_function.R` at line 55, the version in `02_pBACIPS_function.R` is loaded first, then overwritten by the version in `08_effect_sizes.R`. This is not a bug per se (the one used in `08_effect_sizes.R` is the correct one for that context), but it means the version defined in `02_pBACIPS_function.R` is never actually used.

### m4. `test_linear_before` function uses `return(NULL)` inside tryCatch which returns from tryCatch, not from the outer function

**Location:** `08_effect_sizes.R` lines 357-368

**Description:**
```r
result <- tryCatch({
  lmBefore <- lm(lnDiff ~ time.true, data = temp)
  coef_summary <- summary(lmBefore)$coefficients
  if (nrow(coef_summary) < 2) {
    warning("test_linear_before: no slope for ", j)
    return(NULL)  # This return is from the tryCatch expression, not from the for loop body
  }
  coef_summary[, 4]
}, error = function(e) { ... })
```
The `return(NULL)` inside `tryCatch({...})` exits the expression block and causes `result` to be `NULL`. The subsequent `if (!is.null(result))` check handles this correctly. So this is functionally correct but semantically confusing -- the comment says "return" but it only exits the tryCatch block.

### m5. Missing explicit `investr` library load

**Location:** `08_effect_sizes.R` uses `predFit` at lines 1084, 1168, 1247

**Description:**
The `predFit` function from the `investr` package is used without an explicit `library(investr)` or `investr::predFit()` call in `08_effect_sizes.R`. It relies on `investr` being loaded by `00_libraries.R` (confirmed at line 95 of that file). While this works, it creates a hidden dependency that could break if the library loading order changes.

---

## Bootstrap Implementation Review

**Location:** `01_utils.R` (function `bootstrap_biomass`), `04_pisco_processing.R` (usage)

The bootstrap implementation in `01_utils.R` (`bootstrap_biomass` function) is straightforward and correctly implemented:

1. **Resampling:** Draws `n` individuals (where `n` = observed count) from the size frequency population with replacement, repeated `n_resamples` (default 1000) times. This is a standard nonparametric bootstrap.

2. **Size population construction:** Creates a vector `a` by expanding the size frequency table (repeating each size by its count). This correctly reconstructs the empirical distribution.

3. **Biomass conversion:** Applies the species-specific biomass function to each resampled individual, then sums per resample. This is correct for estimating total biomass.

4. **Output:** Returns both the mean and SE (standard deviation of bootstrap sums). The SE is correctly calculated as `sd(bootstrap_sums)`, which represents uncertainty due to unknown individual sizes.

5. **Edge cases:** Handles `n = 0` (returns 0 biomass), `n = NA` (returns NA), and empty size frequency data (returns NA). These are all appropriate.

6. **Reproducibility:** Random seeds are set before bootstrap loops in the processing scripts (e.g., `set.seed(12345)` in `04_pisco_processing.R`).

**One note:** The bootstrap SE represents uncertainty from the size-to-biomass conversion only. It does NOT represent sampling uncertainty in the count itself (i.e., whether the count is a good estimate of true abundance). This is an appropriate design choice given that counts are treated as fixed observations.

---

## Cook's Distance Implementation Review

Cook's distance outlier detection is implemented in `09_meta_analysis.R` (not in the two audited files), but the constant `COOKS_DISTANCE_NUMERATOR = 4` is defined in `00c_analysis_constants.R`. The implementation in the meta-analysis uses the standard `4/n` threshold. This appears correct. Note that the meta-analysis uses `cooks.distance()` from the `metafor` package's `rma.mv` objects, which is the appropriate function for multilevel meta-regression models.

---

## Root Cause Analysis: rbind Column Mismatch (9 vs 10 columns)

**Root cause identified:** See issue M4 above.

The `pvals.poly` dataframe is initialized with 9 columns (line 390). The p-value extraction code (lines 422-427) constructs a vector by concatenating:
- `ps`: 1 value (step model F-test p-value)
- `pl`: 2 values (linear model intercept + slope p-values)
- `psi`: 4 values (sigmoid M, B, L, K p-values) -- BUT fallback models may return 3 values
- `pa`: 3 values (asymptotic M, B, L p-values) -- BUT fallback models may return different counts

Total expected: 1 + 2 + 4 + 3 = 10 values, but the dataframe has 9 columns.

Wait -- re-examining: `colnames(pvals.poly)` assigns 9 names: "Step", "Linear", "Asymptotic M", "Asymptotic B", "Asymptotic L", "Sigmoid M", "Sigmoid B", "Sigmoid L", "Sigmoid K". That's 9 column names. But the p-value vector `c(ps, pl, psi, pa)` has 1 + 2 + 4 + 3 = 10 elements.

**The `pl` vector has 2 elements** (intercept and slope p-values), but only 1 column name ("Linear") is assigned. This means the second element of `pl` (the slope p-value) is placed into the "Asymptotic M" column, shifting all subsequent values by one position.

**This is the root cause:** The column naming is off by one. The "Linear" column actually contains the intercept p-value, and the first "Asymptotic" column contains the linear slope p-value. All sigmoid p-values are shifted into the asymptotic columns, etc.

However, since `pvals.poly` is never used downstream for any calculation (it's only returned from `run_pbacips_loop` and the return value is assigned to a variable like `pbacips.reds.den` which is later `rm()`'d), this misalignment has no impact on results.

The actual rbind warning occurs because the initialization has 9 columns but the vector being rbind'd has 10 elements. R handles this by recycling/truncating with a warning.

---

## Summary Table

| ID | Severity | Description | File | Lines |
|----|----------|-------------|------|-------|
| C1 | CRITICAL | predFit newdata uses wrong variable name for sigmoid models | 08_effect_sizes.R | 1077-1084 |
| C2 | CRITICAL | warnOnly=TRUE allows non-converged NLS models | 02_pBACIPS_function.R | 12 locations |
| C3 | CRITICAL | c() coerces numeric to character in sigmoid SumStats rows | 08_effect_sizes.R | 1106-1108 |
| M1 | MODERATE | First After year gets time.model=0, same as Before period | 08_effect_sizes.R | 309 |
| M2 | MODERATE | AICc comparison unfair with fallback model substitution | 02_pBACIPS_function.R | 1022-1047 |
| M3 | MODERATE | NLS SE uses z-distribution (1.96) instead of t-distribution | 08_effect_sizes.R | 1092-1098 |
| M4 | MODERATE | p-value extraction assumes specific coefficient counts from fallbacks | 08_effect_sizes.R | 422-431 |
| M5 | MODERATE | Norm.poly initialized with wrong column count | 08_effect_sizes.R | 391 |
| M6 | MODERATE | Inconsistent All.RR.sub vs All.RR.sub.trans for KFM urchins | 08_effect_sizes.R | 973-1030 |
| M7 | MODERATE | Winsorization only applied in nls.lm path, not port algorithm path | 02_pBACIPS_function.R | 906-928 |
| m1 | MINOR | df column calculation incorrect for aov/GAM models | 02_pBACIPS_function.R | 1256-1258 |
| m2 | MINOR | Global MODEL_FIT_LOG accumulates across sessions | 02_pBACIPS_function.R | 73-99 |
| m3 | MINOR | run_dharma_diagnostics defined twice with different signatures | Both files | 647 / 102 |
| m4 | MINOR | Confusing return(NULL) inside tryCatch expression block | 08_effect_sizes.R | 357-368 |
| m5 | MINOR | Missing explicit investr library reference | 08_effect_sizes.R | 1084 |
