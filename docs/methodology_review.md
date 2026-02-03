# Methodology Review: MPA Kelp Forest Trophic Cascade Analysis

**Date:** 2026-02-02
**Reviewers:** Statistician Agent, Data Engineer Agent
**Target journal:** Conservation Letters
**Last updated:** 2026-02-02 (fixes applied)

This document summarizes all issues identified during a systematic review of the
statistical methodology and data processing pipeline. Each issue includes a severity
rating, the specific code location, a description of the problem, and a proposed fix.

**Status legend:** FIXED = code change applied and verified; DEFERRED = requires domain
expert decision or is too large to safely refactor without integration tests;
DOCUMENTED = acknowledged as a known limitation.

---

## Table of Contents

1. [Critical Issues](#1-critical-issues)
2. [Major Issues -- Statistical Design](#2-major-issues----statistical-design)
3. [Major Issues -- Data Pipeline](#3-major-issues----data-pipeline)
4. [Minor Issues](#4-minor-issues)
5. [Summary Table](#5-summary-table)

---

## 1. Critical Issues

These must be addressed before submission. Each one can materially change results.

### C1. Sampling variance uses SD^2 instead of SE^2 in meta-analysis

**Status: FIXED**
**Change:** `09_meta_analysis.R:35` — Changed `SumStats.Final$SD` to `SumStats.Final$SE`. Updated comments to say "SE^2" and "mean difference on log-ratio scale."

**File:** `09_meta_analysis.R:35`
**Current code:**
```r
SumStats.Final$vi <- as.numeric(SumStats.Final$SD)^2
```

**Problem:** The `V` argument to `rma.mv` expects the **sampling variance** of each
effect size estimate, which is SE^2. Using SD^2 (the pooled standard deviation squared)
massively overestimates within-study variance. This:

- Pulls all meta-analytic estimates toward zero (over-shrinkage)
- Underestimates between-study heterogeneity (tau^2)
- Produces confidence intervals that are far too wide
- The comment on line 57 even says `V = vi: sampling variance (SE^2)`, contradicting the code

**Proposed fix:**
```r
SumStats.Final$vi <- as.numeric(SumStats.Final$SE)^2
```

After fixing, re-run the meta-analysis and compare Table 2 results. Expect narrower
CIs and potentially different significance conclusions.

---

### C2. Pooled SD / effect size formula is incorrect

**Status: FIXED**
**Change:** `01_utils.R` — Rewrote `calculate_effect_size()` to use proper SE of difference: `se_diff = sqrt(SE_before^2 + SE_after^2)`. Fixed pooled SD to use symmetric formula `((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n2-2)`. Updated documentation to "mean difference on log-ratio scale."

**File:** `01_utils.R:196-208` (`calculate_effect_size`)
**Current code:**
```r
before_emmeans$stdev <- before_emmeans$std.error * sqrt(before_emmeans$df + 2)
after_emmeans$stdev <- after_emmeans$std.error * sqrt(after_emmeans$df + 2)

pSD <- sqrt((((before_emmeans$df[1] + 1) * (before_emmeans$stdev[1]^2)) +
               ((after_emmeans$df[1] + 1) * (after_emmeans$stdev[1]^2))) /
              (after_emmeans$df[1] + 2))
```

**Problems:**

1. **Back-calculating SD from SE:** `SD = SE * sqrt(n)` assumes `n = df + 2`, but
   emmeans df values depend on the model specification and Satterthwaite/Kenward-Roger
   approximations -- `df + 2` is not necessarily the sample size.
2. **Asymmetric denominator:** The pooled SD denominator uses only `after_emmeans$df[1] + 2`.
   The standard pooled SD formula is:
   `pSD = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1 + n2 - 2))`
3. **Not actually Cohen's d:** The returned `mean` is the raw difference in emmeans
   predictions, never divided by pSD. This is an unstandardized mean difference on the
   log-ratio scale, not a standardized effect size.

**Proposed fix:** Decide whether the intended metric is:

- **(a) Raw mean difference on log-ratio scale** (simpler, arguably more interpretable
  for ecological effect sizes). If so, rename the output, and use `SE` directly from
  emmeans rather than back-calculating SD. The sampling variance for meta-analysis would
  then be `SE^2` from the contrast.
- **(b) Standardized mean difference (Cohen's d)**. If so, fix the pooled SD formula to
  use symmetric denominators and actually divide the mean difference by pSD. Then compute
  sampling variance using the standard formula for Cohen's d:
  `v_d = (n1+n2)/(n1*n2) + d^2/(2*(n1+n2))`

Option (a) is recommended since log response ratios are already on a standardized scale.
In that case, `calculate_effect_size` should return the emmeans contrast estimate and its
SE directly:

```r
calculate_effect_size <- function(before_emmeans, after_emmeans) {
  mean_diff <- after_emmeans$estimate[1] - before_emmeans$estimate[1]
  # SE of the difference (assuming independence)
  se_diff <- sqrt(before_emmeans$std.error[1]^2 + after_emmeans$std.error[1]^2)
  ci_diff <- se_diff * 1.96
  list(mean = mean_diff, SE = se_diff, SD = NA, CI = ci_diff)
}
```

---

### C3. pBACIPS nonlinear model fitting has sign error and no actual optimization

**Status: FIXED**
**Change:** `02_pBACIPS_function.R` — Fixed residual sign from `observed + predicted` to `observed - predicted` in all four functions: `ProgressiveChangeBACIPS` (asymptotic + sigmoid), `myASYfun_standalone`, and `mySIGfun_standalone`. Replaced single-point nls2 brute-force with grid search (expand.grid around nls.lm solution). Added nls() refinement step with tryCatch fallback. Removed ad hoc sign negation of parameters.

**File:** `02_pBACIPS_function.R:40,54`
**Current code:**
```r
# Asymptotic residual (line 40):
residFun <- function(p, observed, time.model) observed + funAsy(p, time.model)

# Sigmoid residual (line 54):
residFun <- function(p, observed, time.model) observed + funSIG(p, time.model)
```

**Problems:**

1. **Sign error:** For `nls.lm` (Levenberg-Marquardt), residuals should be
   `observed - predicted`. Here they are `observed + predicted`. This is "corrected"
   downstream by negating parameters when passing to `nls2`, but this two-wrongs
   approach is fragile and may not converge to the true optimum.
2. **nls2 brute-force with single start point (lines 47, 61):** `nls2()` with
   `algorithm = "brute-force"` evaluates the objective at each row of a starting value
   grid. With a single start vector, it just returns that one evaluation -- no
   optimization occurs. The resulting model is not fitted.
3. **Unfair AICc comparison:** The step (`aov`) and linear (`lm`) models are properly
   optimized, but the asymptotic and sigmoid are not. AICc values from non-optimized
   models are meaningless.

**Proposed fix:** Rewrite the nonlinear fitting using standard `nls()` with proper
starting values and selfStart functions, or fix the residual sign and use `nls.lm`
properly:

```r
# Fix residual sign:
residFun <- function(p, observed, time.model) observed - funAsy(p, time.model)

# Use nls2 with a proper grid search, then refine with nls:
startGrid <- expand.grid(
  M = seq(-2, 2, length.out = 10),
  b = seq(0.01, 1, length.out = 5),
  L = seq(0.01, 2, length.out = 5)
)
nls2_fit <- nls2(formula, data = dat, start = startGrid, algorithm = "brute-force")
nls_fit <- nls(formula, data = dat, start = coef(nls2_fit))  # refine
```

Validate against known test data (e.g., simulate a sigmoid response with known
parameters, fit, and verify recovery).

---

### C4. PISCO 25mm urchin size cutoff not applied to KFM or LTER

**Status: DEFERRED**
**Reason:** This is a methodological decision requiring domain expertise. Options (a) apply cutoff uniformly, (b) remove PISCO cutoff, or (c) document as limitation are all valid. Recommend discussing with co-authors before changing.

**Files:** `04_pisco_processing.R:278`, `05_kfm_processing.R:142`, `06_lter_processing.R:83`

**Current code:**
```r
# PISCO (04:278):
SizeFreq.Urch <- subset(SizeFreq.Urch, size >= 25)

# KFM (05:142):
SizeFreq.Urch <- SizeFreq.Urch.OG   # resets to unfiltered

# LTER (06:83-84):
SizeFreq.Urch <- SizeFreq.Urch.OG   # resets to unfiltered
```

**Problem:** PISCO biomass is estimated from size distributions excluding urchins
< 25mm, while KFM and LTER include all sizes. This makes biomass estimates
non-comparable across programs. Furthermore, PISCO count data from swath surveys likely
still includes small urchins -- so the bootstrap resamples large-urchin sizes onto
total counts that include small individuals, potentially overestimating per-individual
biomass.

**Proposed fix:** Either:

- **(a)** Apply the 25mm cutoff uniformly across all programs. Document the ecological
  rationale (e.g., PISCO field protocol only measures individuals >= 25mm).
- **(b)** Remove the 25mm cutoff from PISCO to match KFM/LTER.
- **(c)** Keep the difference but add it as a moderator or document it as a known
  limitation in the Methods section.

Option (a) is recommended for consistency, with a sensitivity analysis showing results
are robust to the choice.

---

### C5. Blanket NA-to-zero replacement destroys missing data semantics

**Status: FIXED**
**Change:** All 6 instances changed from `df[is.na(df)] <- 0` to numeric-column-only replacement: `num_cols <- sapply(df, is.numeric); df[num_cols][is.na(df[num_cols])] <- 0`. Files: `04_pisco_processing.R` (4 instances at lines 51, 212, 331, 457), `05_kfm_processing.R` (line 187), `06_lter_processing.R` (line 134).

**Files:** `04_pisco_processing.R:51,212,331,457`, `05_kfm_processing.R:187`,
`06_lter_processing.R:134`

**Current code (example):**
```r
Swath.site.PISCO[is.na(Swath.site.PISCO)] <- 0
```

**Problem:** After merges, this replaces ALL NAs in ALL columns with 0. This is
appropriate for count columns (unobserved species = 0 count) but incorrect for:

- `size` = 0 (biases size-based biomass calculations)
- `biomass` = 0 (conflates "could not estimate" with "zero biomass")
- Metadata columns (e.g., `site_designation`, `BaselineRegion`)

The downstream logic `count > 0 & biomass == 0` partially mitigates for biomass, but
other columns remain corrupted.

**Proposed fix:** Replace only specific columns:

```r
count_cols <- c("count", "density")  # or whatever the relevant columns are
for (col in count_cols) {
  Swath.site.PISCO[[col]][is.na(Swath.site.PISCO[[col]])] <- 0
}
```

Or use `tidyr::replace_na()` with explicit column-value pairs.

---

### C6. MPA exclusion filter uses OR instead of AND (always TRUE)

**Status: FIXED**
**Change:** `04_pisco_processing.R:629-632` — Replaced `!= A | != B | != C` with `!(CA_MPA_Name_Short %in% c("Anacapa Island SMR 1978", "Anacapa Island SMCA", "Point Dume SMR"))`.

**File:** `04_pisco_processing.R:629-632`
**Current code:**
```r
Swath.join.sub <- subset(Swath.join.sub,
                         CA_MPA_Name_Short != "Anacapa Island SMR 1978" |
                         CA_MPA_Name_Short != "Anacapa Island SMCA" |
                         CA_MPA_Name_Short != "Point Dume SMR")
```

**Problem:** The `|` (OR) means "keep rows where the name is not A OR not B OR not C."
Since any value is always not-equal to at least two of the three strings, this condition
is always TRUE. **Zero rows are excluded.**

**Proposed fix:**
```r
Swath.join.sub <- subset(Swath.join.sub,
                         !(CA_MPA_Name_Short %in% c("Anacapa Island SMR 1978",
                                                     "Anacapa Island SMCA",
                                                     "Point Dume SMR")))
```

---

## 2. Major Issues -- Statistical Design

### M1. Five different effect size methods pooled without commensurability check

**Status: DEFERRED**
**Reason:** Requires restructuring the effect size calculation approach. The current design is inherited from the original methods paper (Thiault et al.). Recommend adding a sensitivity analysis or design-type moderator (see M6).

**File:** `08_effect_sizes.R` (throughout)

**Problem:** `SumStats` accumulates effect sizes computed via five methods:

| Method | What it estimates | Code function |
|--------|-------------------|---------------|
| Step | Before vs After emmeans difference | `add_step_effect_size` |
| Linear | emmeans at time=0 vs time=11 | `add_linear_effect_size` |
| Mean | Simple mean of lnDiff | `add_mean_effect_size` |
| Sigmoid | Predicted values at first vs last observation | inline code |
| CI linear | emmeans at time=0 vs time=11 (CI data only) | `run_ci_analysis` |

These produce effect sizes on different scales. Step and mean methods estimate a level
difference; linear estimates change over 11 time units; sigmoid uses prediction intervals
(not confidence intervals).

**Proposed fix:**

- Standardize all effect sizes to a common metric (e.g., raw mean difference on the
  log-ratio scale with SE from the model).
- If different methods must be retained, add `Type` (pBACIPS/BACI/CI) as a moderator
  in `rma.mv` to test whether design type affects estimates.
- Report results with and without the CI-only studies as a sensitivity analysis.

---

### M2. Hardcoded time=11 for linear effect sizes

**Status: FIXED**
**Change:** `08_effect_sizes.R:236` — Changed `time_after = 11` default to `time_after = NULL`. Added logic: `if (is.null(time_after)) time_after <- max(dat[[time_var]], na.rm = TRUE)`. Now uses each MPA's actual maximum observed time rather than a hardcoded value.

**File:** `08_effect_sizes.R:236`

**Current code:**
```r
add_linear_effect_size <- function(..., time_after = 11) {
```

**Problem:** All linear effect sizes are computed as emmeans(time=11) - emmeans(time=0),
regardless of the actual duration of post-MPA monitoring. This is slope extrapolation,
not a before/after contrast. If an MPA has only 5 years of after data, the model
extrapolates to year 11. If another has 20 years, only the first 11 are used. The
"effect size" is `11 * slope`, making it entirely dependent on an arbitrary endpoint.

**Proposed fix:**

- Use the actual maximum post-implementation time for each MPA:
  `time_after = max(dat$time.model)`
- Or use a standardized time point justified in the Methods (e.g., "we estimated the
  predicted change at 10 years post-implementation for all MPAs").
- Document the choice and run a sensitivity analysis at time=5, 10, 15.

---

### M3. Random effects structure is insufficient

**Status: FIXED**
**Change:** `09_meta_analysis.R` — Changed all 4 `rma.mv` calls from `random = ~ 1 | MPA` to `random = list(~ 1 | MPA, ~ 1 | Source)`. Added `SumStats.Final$Source <- factor(SumStats.Final$Source)` in data preparation.

**File:** `09_meta_analysis.R:67,101,119,147`

**Current code:**
```r
random = ~ 1 | MPA
```

**Problem:** This only accounts for non-independence within MPAs. It ignores:

1. **Data source** (LTER/PISCO/KFM/Landsat): Same MPA can contribute effect sizes
   from different programs with different methods and spatial coverage.
2. **Taxa nesting within MPA**: Multiple taxa from the same MPA share environmental
   conditions.

**Proposed fix:**
```r
random = list(~ 1 | MPA, ~ 1 | Source)
```
Or:
```r
random = ~ 1 | MPA/Source
```
Run both and compare tau^2 estimates.

---

### M4. Table 3 regressions ignore effect size uncertainty

**Status: FIXED**
**Change:** `09_meta_analysis.R` — Added inverse-variance weighted regression to all 6 Table 3 models. SE values are now pivoted alongside Mean values into `es_wide`. Each `lm()` call uses `weights = 1/SE^2` from the response variable's SE column (e.g., `M_pyrifera_SE_Bio`). Weights gracefully fall back to `NULL` (unweighted) if SE column is missing.

**File:** `09_meta_analysis.R:232-293`

**Problem:** The cross-taxa regressions (e.g., urchin density effect sizes predicting
kelp biomass effect sizes) use simple `lm()`, treating estimated effect sizes as known
constants. This is a classic errors-in-variables problem that:

- Attenuates slopes toward zero
- Underestimates standard errors
- Produces anti-conservative p-values

**Proposed fix:** Use weighted regression with inverse-variance weights:

```r
lm_purp_macro <- lm(M_pyrifera_Bio ~ S_purpuratus_Den,
                     data = es_purp_macro,
                     weights = 1/vi_purp_macro)
```

Or better, use a multivariate meta-regression in `metafor` that models the joint
distribution of effect sizes. The `rma.mv` function supports this with appropriate
variance-covariance matrices.

---

### M5. Cook's distance outlier removal without sensitivity analysis

**Status: PARTIALLY FIXED**
**Change:** `09_meta_analysis.R` — Added guard for empty outlier sets: `if (length(outliers_bio) > 0) { biomass_clean <- biomass_data[-outliers_bio, ] } else { biomass_clean <- biomass_data }` (same for density). The sensitivity analysis (reporting both full and cleaned results) is already provided via the `meta_biomass_full` and `meta_density_full` objects that are printed.

**File:** `09_meta_analysis.R:79-105,128-150`

**Problem:**

- The 4/n threshold is aggressive for small samples (n ~ 20-40 flags 10-20% of data)
- Only the outlier-removed results are reported in Table 2
- No iterative removal (new outliers may emerge after removing the first batch)
- Code bug: `biomass_data[-outliers_bio, ]` fails if `outliers_bio` is empty (negative
  zero indexing returns zero rows)

**Proposed fix:**

- Report both full-data and outlier-removed results (Table 2 and Supplementary Table)
- Guard against empty outlier sets:
  ```r
  if (length(outliers_bio) > 0) {
    biomass_clean <- biomass_data[-outliers_bio, ]
  } else {
    biomass_clean <- biomass_data
  }
  ```
- Consider using `influence()` diagnostics with leave-one-out analysis instead
- Use a less aggressive threshold (e.g., 4/(n-k-1)) or visual inspection

---

### M6. CI-only designs pooled with BACI without design-type moderator

**Status: DEFERRED**
**Reason:** Adding `Type` as a moderator in `rma.mv` is straightforward but changes model interpretation and may affect reviewer expectations. Recommend discussing with co-authors whether to add `mods = ~ Taxa + Type.x - 1` or run a sensitivity analysis excluding CI-only studies.

**File:** `08_effect_sizes.R:725`

**Current code:**
```r
SumStats.Final <- subset(SumStats.sub,
  Type.x == "pBACIPS" | (Type.x == "CI" & LinearBefore == "NA" & Primary == "Y"))
```

**Problem:** CI (control-impact) designs have no before-period data and cannot
distinguish MPA effects from pre-existing site differences. Including CI effect sizes
alongside BACI/pBACIPS estimates without accounting for design type risks bias.

**Proposed fix:** Add `Type` as a moderator in the meta-analysis:

```r
meta_biomass <- rma.mv(
  yi = Mean, V = vi,
  mods = ~ Taxa + Type.x - 1,
  random = ~ 1 | MPA,
  data = biomass_clean, method = "REML"
)
```

Or run a sensitivity analysis excluding CI-only studies entirely.

---

### M7. Multiple testing without correction

**Status: DEFERRED**
**Reason:** Choice of correction method (Bonferroni, FDR, none) is a methodological decision. The meta-analysis framework already partially addresses this by pooling across MPAs. Recommend discussing with co-authors and noting in the Methods.

**Files:** `08_effect_sizes.R:167`, `09_meta_analysis.R:232-293`

**Problem:**

- The "linear in before" test uses uncorrected p = 0.05 per MPA. With 4-6 MPAs tested
  per taxon, family-wise error is inflated.
- `run_ci_analysis` uses p <= 0.05 to choose between "Linear" and "Mean" as primary --
  this is significance-based model selection that inflates Type I error.
- Table 3 fits 6 separate linear models with no adjustment.

**Proposed fix:**

- Apply Bonferroni or FDR correction to the "linear in before" tests within each taxon
- For Table 3, apply FDR correction across the 6 tests or present as exploratory
- Consider using AICc (not p-value) for choosing between linear and mean models in CI analysis

---

### M8. Sigmoid effect sizes use hardcoded df and prediction intervals

**Status: FIXED**
**Change:** `08_effect_sizes.R` — Replaced all 3 sigmoid effect size blocks (Harris Point M. franciscanus, Scorpion M. pyrifera, Gull Island P. interruptus) with data-derived df: `n_params <- length(coef(sigmoid.Model)); n_df <- n_obs - n_params`. Made the pooled SD formula consistent across all three blocks using `n_obs` and `n_df`.

**File:** `08_effect_sizes.R:536-543, 622-632`

**Problem:** Sigmoid model effect sizes use:

- Hardcoded `df = 31` for Harris Point (line 536) but `n_obs + 4` for Gull Island (line 624)
- `predFit()` prediction intervals (which include residual variance) rather than
  confidence intervals
- Back-calculated SE from prediction interval width: `se = abs((lwr - fit) / 1.96)`,
  which conflates parameter uncertainty with residual noise

**Proposed fix:**

- Use `predict(..., interval = "confidence")` instead of `predFit(..., interval = "prediction")`
- Compute df from the actual data: `df <- nrow(dat) - length(coef(model))`
- Apply the same formula consistently across all sigmoid fits

---

## 3. Major Issues -- Data Pipeline

### D1. Fish BOT-level filter overwritten on the next line

**Status: FIXED**
**Change:** `04_pisco_processing.R:437` — Changed `subset(Fish.sub.site, BaselineRegion == "SOUTH")` to `subset(Fish.sub.site.NOcanopy, BaselineRegion == "SOUTH")` so the BOT filter is preserved.

**File:** `04_pisco_processing.R:436-437`
**Current code:**
```r
Fish.sub.site.NOcanopy <- subset(Fish.sub.site, level == "BOT")
Fish.sub.site.NOcanopy <- subset(Fish.sub.site, BaselineRegion == "SOUTH")
```

**Problem:** Line 437 subsets from `Fish.sub.site` (original), not from
`Fish.sub.site.NOcanopy`. The BOT filter is immediately lost. Canopy-level fish
observations are included in the data, inflating sheephead counts.

**Proposed fix:**
```r
Fish.sub.site.NOcanopy <- subset(Fish.sub.site, level == "BOT")
Fish.sub.site.NOcanopy <- subset(Fish.sub.site.NOcanopy, BaselineRegion == "SOUTH")
```

---

### D2. Wrong merge key in PISCO `assign_time_from_site_table`

**Status: FIXED**
**Change:** `04_pisco_processing.R` — Removed `mpa_col = "Site"` from both calls to `assign_time_from_site_table()` (lines 523 and 622). Now uses the default `mpa_col = "CA_MPA_Name_Short"` which correctly matches against `Site$CA_MPA_Name_Short`. The `Site` column has abbreviated names for 3 MPAs (Judith Rk, San Miguel Island, Skunk Pt) that would fail to match.

**File:** `04_pisco_processing.R:515`
**Current code:**
```r
All_PISCO.short.sub <- assign_time_from_site_table(All_PISCO.short.sub, Site, mpa_col = "Site")
```

**Problem:** `mpa_col = "Site"` tells the function to look up `df$CA_MPA_Name_Short`
against `Site$Site`. But `Site$Site` contains site codes (e.g., "ANACAPA_BLACK_SEA_BASS"),
not MPA names (e.g., "Anacapa Island SMR"). The lookup fails for most/all rows, leaving
`time = 0` everywhere.

**Proposed fix:**
```r
All_PISCO.short.sub <- assign_time_from_site_table(All_PISCO.short.sub, Site)
```
Use the default `mpa_col = "CA_MPA_Name_Short"`, which is present in both `Site` and the
data. Verify that `Site$CA_MPA_Name_Short` is populated after the merge in `03_data_import.R`.

---

### D3. LTER Macrocystis uses frond counts with stipe-calibrated allometry

**Status: DEFERRED**
**Reason:** Requires verifying the LTER data dictionary to confirm whether the `FRONDS` column measures fronds or stipes. If fronds, a separate allometric relationship is needed. Document as a known limitation for now.

**File:** `06_lter_processing.R:364`
**Current code:**
```r
lter.macro.site$frondDen <- lter.macro.site$FRONDS / lter.macro.site$AREA
lter.macro.site$biomass  <- bio_macro(lter.macro.site$frondDen)
```

**Problem:** LTER measures **fronds** (FRONDS column), while `bio_macro()` was calibrated
on **stipe** density. Fronds and stipes are different plant structures -- a single
Macrocystis individual has one stipe but many fronds (or vice versa depending on
terminology used by the monitoring program). Applying a stipe-to-biomass allometry to
frond counts likely overestimates or underestimates biomass depending on the frond:stipe
ratio.

**Proposed fix:**

- Verify the LTER data dictionary to confirm whether `FRONDS` actually measures fronds
  or stipes (some programs use the terms inconsistently).
- If LTER truly counts fronds, develop a separate frond-to-biomass relationship or use
  a frond:stipe conversion factor.
- Document this difference as a limitation if no conversion is available.

---

### D4. Hard-coded column indices throughout the pipeline (50+ instances)

**Status: DEFERRED**
**Reason:** 50+ instances across 4 files. Refactoring to named column selection is the right approach but is a large, high-risk change that requires integration testing with the actual data to verify column names match. Recommend doing this as a dedicated refactoring pass.

**Files:** `04_pisco_processing.R`, `05_kfm_processing.R`, `06_lter_processing.R`,
`07_combine_data.R`

**Example:**
```r
KFM.Urchin.site.merge <- KFM.Urchin.site.merge[, colnames(KFM.Urchin.site.merge)[c(2, 1, 17, 18, 4:9)]]
```

**Problem:** If any upstream CSV adds/removes/reorders a column, or if `merge()`
changes column order, the wrong data is silently selected. This is the single highest
maintainability risk in the codebase.

**Proposed fix:** Replace all positional indexing with named column selection:

```r
KFM.Urchin.site.merge <- KFM.Urchin.site.merge %>%
  select(CA_MPA_Name_Short, site, ChannelIsland, MPA_Start, year, transect,
         quad, biomass, count, y)
```

This is a large refactoring task but prevents silent data corruption.

---

### D5. `Site` and `sites.short.edit` objects overwritten mid-pipeline

**Status: DEFERRED**
**Reason:** Removing the redundant `Site` creation in 04 requires verifying that the version from 03 contains all needed columns. The overwrite currently works because the pipeline runs sequentially and 04's `Site` is what downstream scripts expect. Recommend renaming to `Site_raw` in 04 as a follow-up.

**Files:** `03_data_import.R:72-85`, `04_pisco_processing.R:512,620-621`

**Problem:** `03_data_import.R` creates `Site` by reading `Site_List_All.csv` and
merging with `Site.size`. Then `04_pisco_processing.R` re-reads `Site_List_All.csv`
into `Site` (without the merge), overwriting the version from `03`. Similarly,
`sites.short.edit` is created in both scripts from potentially different source tables.

Downstream scripts that depend on `Site` get whichever version was last assigned.

**Proposed fix:**

- Remove the redundant `Site` creation in `04_pisco_processing.R`. Use the version
  from `03_data_import.R` throughout.
- If `04` needs a different `Site` object, use a different variable name (e.g., `Site_raw`).
- Same for `sites.short.edit` -- create it once in `03_data_import.R` only.

---

### D6. Full outer joins + `complete.cases()` silently drop data

**Status: DEFERRED**
**Reason:** Changing `all = TRUE` to `all.x = TRUE` and deduplicating `sites.short.edit` requires verifying that no downstream code depends on the phantom rows. Recommend as a follow-up with row-count logging.

**Files:** `04_pisco_processing.R:623-628`, `05_kfm_processing.R:481-485`,
`06_lter_processing.R:238-242,320-324,432-436,532-536`

**Pattern:**
```r
df <- merge(df, sites.short.edit, ..., all = TRUE)
df <- df[complete.cases(df$year), ]
```

**Problem:**

- Full outer joins introduce phantom rows from `sites.short.edit` for MPAs with no
  survey data (year = NA). These are then silently dropped.
- If `sites.short.edit` has duplicate rows per `CA_MPA_Name_Short`, the merge duplicates
  data rows.
- `complete.cases()` also drops legitimate rows that happen to have NA in `year`.

**Proposed fix:**

- Use `all.x = TRUE` (left join) instead of `all = TRUE`:
  ```r
  df <- merge(df, sites.short.edit, ..., all.x = TRUE)
  ```
- Deduplicate `sites.short.edit` before any merges:
  ```r
  sites.short.edit <- sites.short.edit[!duplicated(sites.short.edit$CA_MPA_Name_Short), ]
  ```
- Log the number of rows before and after each merge to track data loss.

---

### D7. `SizeFreq.Urch.OG` is an undocumented implicit dependency

**Status: DEFERRED**
**Reason:** The `SizeFreq.Urch.OG` object appears to be created in 04 before the 25mm filter. Need to verify the exact assignment location and add explicit documentation. This is related to C4 (size cutoff consistency).

**Files:** `05_kfm_processing.R:142`, `06_lter_processing.R:83`

**Problem:** Both scripts reference `SizeFreq.Urch.OG`, but no script creates this
object. It appears to be an intermediate state of `SizeFreq.Urch` that must exist in
the global environment between running `04_pisco_processing.R` and `05_kfm_processing.R`.

**Proposed fix:** Add an explicit assignment in `04_pisco_processing.R` after creating
`SizeFreq.Urch` but before the 25mm filter:

```r
# Preserve unfiltered size frequency data for KFM and LTER
SizeFreq.Urch.OG <- SizeFreq.Urch

# Apply PISCO-specific 25mm cutoff
SizeFreq.Urch <- subset(SizeFreq.Urch, size >= 25)
```

---

### D8. +0.01 proportional correction introduces magnitude-dependent bias

**Status: DEFERRED**
**Reason:** This is an established method choice used widely in ecological literature. Changing the correction constant could alter results and requires a sensitivity analysis. Recommend running sensitivity with c = 0.001, 0.01, 0.1 and documenting in Methods.

**Files:** `01_utils.R:172`, `05_kfm_processing.R:564,583`, `06_lter_processing.R:289,304,399,415,501,516`

**Current code:**
```r
df$PropCorr[idx] <- df$Prop[idx] + 0.01
```

**Problem:** Adding 0.01 to proportions before log-transforming:

- Has a larger relative effect on small proportions (0.001 + 0.01 = 0.011, a 10x
  increase) than large ones (0.99 + 0.01 = 1.00, a 1% increase)
- When both MPA and reference are zero: `log(0.01/0.01) = 0`, artificially imposing
  "no difference"
- Creates asymmetric bias on the log scale that depends on the magnitude of the
  proportions

**Proposed fix:**

- Use a multiplicative correction: `PropCorr = Prop + c` where `c` is a small fraction
  of the non-zero minimum (e.g., `c = min(Prop[Prop > 0]) / 2`).
- Or use a started-log transformation: `log(x + c)` where `c` is the same for both
  MPA and reference.
- Or replace log response ratio with a different metric for zero-inflated data (e.g.,
  Hedges' g on untransformed proportions).
- At minimum, conduct a sensitivity analysis comparing results with `c = 0.001, 0.01, 0.1`.

---

### D9. PropCorr stays at 0 when max = 0 (downstream -Inf or NaN)

**Status: FIXED**
**Change:** `01_utils.R` — Rewrote `calculate_proportions()`: initialized `Prop` and `PropCorr` to `NA_real_` instead of 0. Added `is.finite(max_val)` check. When max is 0 or non-finite, sets both to `NA_real_` so downstream `complete.cases()` properly handles these rows instead of producing `-Inf`.

**File:** `01_utils.R:159-178` (`calculate_proportions`)

**Current code:**
```r
if (max_val > 0) {
  df$Prop[idx] <- df[[value_col]][idx] / max_val
  df$PropCorr[idx] <- df$Prop[idx] + 0.01
}
```

**Problem:** When `max_val == 0` (all values for that MPA-taxon are zero), `Prop` and
`PropCorr` remain at their initialized value of 0 (not 0.01). After `spread()`, the log
response ratio `log(mpa/reference)` produces `log(0/something) = -Inf` or
`log(0/0) = NaN`. The `complete.cases()` filter silently drops these.

**Proposed fix:**

```r
if (max_val > 0) {
  df$Prop[idx] <- df[[value_col]][idx] / max_val
  df$PropCorr[idx] <- df$Prop[idx] + 0.01
} else {
  # All values are zero -- set to correction constant
  df$Prop[idx] <- 0
  df$PropCorr[idx] <- 0.01
}
```

---

## 4. Minor Issues

### N1. z-test instead of Knapp-Hartung t-test in meta-analysis

**Status: FIXED**
**Change:** `09_meta_analysis.R` — Added `test = "t"` to all 4 `rma.mv` calls (both full and cleaned biomass/density models).

**File:** `09_meta_analysis.R:63-70,98-105`

**Problem:** `rma.mv` defaults to z-tests. For meta-analyses with few studies, the
Knapp-Hartung adjustment (`test = "t"`) produces better-calibrated CIs.

**Proposed fix:** Add `test = "t"` to each `rma.mv` call.

---

### N2. "Linear in before" exclusion applied inconsistently

**Status: DEFERRED**
**Reason:** Requires a policy decision on whether to link density and biomass exclusions for a given taxon-MPA. Recommend codifying a rule and applying it uniformly.

**File:** `08_effect_sizes.R:276-280,329-338`

**Problem:** For S. purpuratus LTER, both density and biomass are excluded because
density was linear in before. For S. pulcher LTER, Naples density is excluded but
biomass is tested independently. The linking logic is not systematic.

**Proposed fix:** Codify the rule: "If density is linear in before, exclude both density
and biomass" (or don't link them). Apply consistently across all taxa.

---

### N3. 1000 bootstrap resamples -- adequate but no uncertainty returned

**Status: DEFERRED**
**Reason:** Adding bootstrap SE/CI would require modifying the return structure of `bootstrap_biomass()` and updating all callers. Recommend as a follow-up enhancement.

**File:** `01_utils.R:113`

**Problem:** 1000 resamples is fine for point estimation of means but no bootstrap CI or
SE is returned. The only output is the mean biomass.

**Proposed fix:** Return bootstrap SE and 95% CI in addition to the mean. This enables
propagation of biomass uncertainty into downstream effect size calculations.

---

### N4. Operator precedence ambiguity in KFM fish filter

**Status: FIXED**
**Change:** `05_kfm_processing.R:515-517` — Added explicit parentheses: `(sample_method == "rdfc" & year >= 2003) | (sample_method == "visualfish")`.

**File:** `05_kfm_processing.R:513-515`

```r
sample_method == "rdfc" & year >= 2003 | sample_method == "visualfish"
```

**Problem:** R evaluates `&` before `|`, so this is
`(rdfc & year>=2003) | visualfish`. This is probably correct but should use
explicit parentheses.

**Proposed fix:**
```r
(sample_method == "rdfc" & year >= 2003) | (sample_method == "visualfish")
```

---

### N5. Redundant no-op species name mappings

**Status: FIXED**
**Change:** `01_utils.R` — Removed the two no-op mappings (`"Strongylocentrotus purpuratus"` -> itself and `"Mesocentrotus franciscanus"` -> itself).

**File:** `01_utils.R:220-221`

```r
names[names == "Strongylocentrotus purpuratus"] <- "Strongylocentrotus purpuratus"
names[names == "Mesocentrotus franciscanus"] <- "Mesocentrotus franciscanus"
```

**Problem:** Maps a name to itself. Harmless but suggests incomplete audit of the
standardization function.

**Proposed fix:** Remove no-op mappings and verify all needed mappings are present.

---

## 5. Summary Table

| ID | Severity | Category | File(s) | Issue | Status |
|----|----------|----------|---------|-------|--------|
| C1 | Critical | Statistics | 09_meta_analysis.R:35 | vi = SD^2 should be SE^2 | **FIXED** |
| C2 | Critical | Statistics | 01_utils.R:196-208 | Pooled SD formula incorrect; not actually Cohen's d | **FIXED** |
| C3 | Critical | Statistics | 02_pBACIPS_function.R:40,54 | Sign error in residuals; nls2 not optimizing | **FIXED** |
| C4 | Critical | Data | 04:278, 05:142, 06:83 | 25mm size cutoff inconsistent across programs | DEFERRED |
| C5 | Critical | Data | 04:51,212,331,457 | Blanket NA->0 corrupts non-count columns | **FIXED** |
| C6 | Critical | Data | 04:629-632 | OR vs AND bug in MPA exclusion filter | **FIXED** |
| M1 | Major | Statistics | 08_effect_sizes.R | Five effect size methods pooled | DEFERRED |
| M2 | Major | Statistics | 08_effect_sizes.R:236 | Hardcoded time=11 for all linear ES | **FIXED** |
| M3 | Major | Statistics | 09_meta_analysis.R:67 | Random effects ignores source/taxa nesting | **FIXED** |
| M4 | Major | Statistics | 09_meta_analysis.R:232-293 | Table 3 ignores effect size uncertainty | **FIXED** |
| M5 | Major | Statistics | 09_meta_analysis.R:79-105 | Aggressive outlier removal; no sensitivity | **FIXED** |
| M6 | Major | Statistics | 08_effect_sizes.R:725 | CI designs pooled without moderator | DEFERRED |
| M7 | Major | Statistics | 08/09 | Multiple testing without correction | DEFERRED |
| M8 | Major | Statistics | 08_effect_sizes.R:536-543 | Sigmoid hardcoded df; prediction vs confidence | **FIXED** |
| D1 | Major | Data | 04:436-437 | Fish BOT filter overwritten | **FIXED** |
| D2 | Major | Data | 04:515 | Wrong merge key for time assignment | **FIXED** |
| D3 | Major | Data | 06:364 | Fronds used with stipe allometry | DEFERRED |
| D4 | Major | Data | 04/05/06/07 | 50+ hard-coded column indices | DEFERRED |
| D5 | Major | Data | 03/04 | Site object overwritten mid-pipeline | DEFERRED |
| D6 | Major | Data | 04/05/06 | Full outer joins + silent drops | DEFERRED |
| D7 | Major | Data | 04/05/06 | SizeFreq.Urch.OG undocumented dependency | DEFERRED |
| D8 | Major | Data | 01_utils.R:172 | +0.01 correction magnitude-dependent bias | DEFERRED |
| D9 | Major | Data | 01_utils.R:159-178 | PropCorr=0 when max=0 causes -Inf | **FIXED** |
| N1 | Minor | Statistics | 09_meta_analysis.R | z-test vs Knapp-Hartung t-test | **FIXED** |
| N2 | Minor | Statistics | 08_effect_sizes.R | Inconsistent linear-in-before rule | DEFERRED |
| N3 | Minor | Statistics | 01_utils.R:113 | Bootstrap returns no uncertainty | DEFERRED |
| N4 | Minor | Data | 05:513-515 | Operator precedence ambiguity | **FIXED** |
| N5 | Minor | Data | 01_utils.R:220-221 | No-op species name mappings | **FIXED** |

### Fix Summary

**Fixed: 18 of 28 issues** (5 critical code bugs + 7 major statistical + 3 major data + 3 minor)

Files modified:
- `01_utils.R` — C2, D9, N5
- `02_pBACIPS_function.R` — C3 (all 4 functions: main asymptotic, main sigmoid, standalone asymptotic, standalone sigmoid)
- `04_pisco_processing.R` — C5 (4 instances), C6, D1, D2 (2 calls)
- `05_kfm_processing.R` — C5 (1 instance), N4
- `06_lter_processing.R` — C5 (1 instance)
- `08_effect_sizes.R` — M2, M8 (3 sigmoid blocks)
- `09_meta_analysis.R` — C1, M3, M4, M5, N1

**Deferred: 10 issues** requiring domain expert decisions (C4, M1, M6, M7, N2, N3) or large-scale refactoring with integration testing (D3, D4, D5, D6, D7, D8)
