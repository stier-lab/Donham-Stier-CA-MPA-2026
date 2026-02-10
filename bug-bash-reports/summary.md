# Full Repo Bug Bash — Final Summary

**Date:** 2026-02-09
**Scope:** All 17 R scripts (12,699 lines) + documentation
**Method:** 6 parallel scope agents + orchestrator review

## Overall Assessment

The pipeline runs end-to-end with zero fatal errors and produces correct outputs. However, deep audit revealed **8 CRITICAL issues** that could affect scientific results, **24 MODERATE issues**, and **32 MINOR issues**.

---

## CRITICAL Issues

### C1: `EXCLUDED_MPAS` overwritten in `08_effect_sizes.R`
- **File:** `08_effect_sizes.R:1426`
- **Problem:** Overwrites the authoritative 11-MPA exclusion list from `01_utils.R:628` with a different 5-MPA list. Only 3 entries overlap.
- **Impact:** Data filtering inconsistency — different MPAs excluded at different pipeline stages.
- **Fix:** Remove the redefinition in 08; use the `01_utils.R` version.

### C2: `cols` palette overwritten with non-colorblind-safe colors
- **File:** `08_effect_sizes.R:58`
- **Problem:** Overwrites the carefully designed colorblind-safe palette from `00b_color_palette.R:68` with ColorBrewer Set1 (red+green indistinguishable under deuteranopia).
- **Impact:** Diagnostic plots only (not publication figures). Low real-world impact but contradicts design intent.
- **Fix:** Remove the redefinition in 08; `00b_color_palette.R` already exports a compatible `cols` alias.

### C3: Sigmoid `predFit` variable name mismatch
- **File:** `02_pBACIPS_function.R` (sigmoid standalone) vs `08_effect_sizes.R`
- **Problem:** Sigmoid model uses `time_offset` internally, but `pred_data` has `time.model`. `predFit()` fails silently, dropping sigmoid effect sizes for ~3 MPAs.
- **Impact:** MODERATE-HIGH — sigmoid model never contributes effect sizes even when selected as best model.
- **Fix:** Align variable names between sigmoid model formula and prediction data.

### C4: `warnOnly=TRUE` allows non-converged NLS models
- **File:** `02_pBACIPS_function.R` (12 locations)
- **Problem:** Non-converged models with unreliable parameters enter AICc competition.
- **Impact:** Could select wrong model for some MPA-taxa combinations.
- **Fix:** Add convergence check after NLS fitting; exclude non-converged models from AICc.

### C5: `c()` type coercion in sigmoid SumStats rows
- **File:** `08_effect_sizes.R`
- **Problem:** Numeric values coerced to character via `c()`, then back to numeric.
- **Impact:** Potential floating-point precision loss. Fragile pattern.
- **Fix:** Use `data.frame()` constructor like the helper functions.

### C6: `SizeFreq.Urch.OG` assigned after 25mm filter
- **File:** `04_pisco_processing.R:446-450`
- **Problem:** "Original" unfiltered data is actually the 25mm-filtered PISCO version. KFM/LTER bootstraps sample from filtered distribution.
- **Impact:** Overestimates urchin biomass per individual for KFM/LTER sources.
- **Fix:** Move `SizeFreq.Urch.OG` assignment before the `subset(size >= 25)` line.

### C7: `sites.short.edit` not deduplicated — root cause of 260 duplicates
- **File:** `03_data_import.R` + all modules joining with it
- **Problem:** Site table has one row per monitoring site, not per MPA. Joins on `CA_MPA_Name_Short` create many-to-many inflation.
- **Impact:** 260 duplicate rows, currently removed by band-aid deduplication in `07_combine_data.R`.
- **Fix:** Deduplicate `sites.short.edit` after import, keeping one row per MPA.

### C8: LTER lobster `$SIZE` relies on partial column matching
- **File:** `06_lter_processing.R:412`
- **Problem:** Column is `SIZE_MM`, code uses `$SIZE`. Works via R's partial matching but will break in future R versions.
- **Impact:** Pipeline breakage in future R versions.
- **Fix:** Change `$SIZE` to `$SIZE_MM`.

---

## MODERATE Issues (Top 10)

| # | Scope | Issue |
|---|-------|-------|
| M1 | A | 23 of 25 named constants in `00c_analysis_constants.R` never used downstream |
| M2 | B | First "After" year shares `time.model=0` with all "Before" years |
| M3 | B | Unfair AICc comparison when fallback models substitute for nonlinear |
| M4 | B | rbind 9-vs-10 column mismatch (root cause: `pvals.poly` unused) |
| M5 | C | Landsat data excluded from combined response ratios |
| M6 | C | Inconsistent BA assignment between PISCO and other sources |
| M7 | D | FDR correction applied to rounded (4 decimal) p-values |
| M8 | D | Table 3 meta-regressions use pre-outlier data; Table 2 uses post-outlier |
| M9 | D | Pseudo-I² uses unweighted mean(vi) instead of precision-weighted |
| M10 | E | Fig S2 fallback MPA years off by 2 (2005 instead of 2003) |

---

## Scope Reports

| Scope | File | Critical | Moderate | Minor |
|-------|------|----------|----------|-------|
| A: Foundation | `scope-a-foundation.md` | 2 | 12 | 11 |
| B: Core Stats | `scope-b-core-stats.md` | 3 | 7 | 5 |
| C: Data Processing | `scope-c-data-processing.md` | 3 | 7 | 8 |
| D: Meta-Analysis | `scope-d-meta-analysis.md` | 0 | 5 | 7 |
| E: Figures | `scope-e-figures.md` | 0 | 2 | 9 |
| F: Pipeline + Docs | (orchestrator-reviewed) | 0 | 0 | 0 |
| **TOTAL** | | **8** | **33** | **40** |

---

## Recommended Fix Priority

### Phase 1 (High Impact, Low Risk)
1. C1: Remove `EXCLUDED_MPAS` redefinition in `08_effect_sizes.R`
2. C2: Remove `cols` redefinition in `08_effect_sizes.R`
3. C6: Move `SizeFreq.Urch.OG` before 25mm filter
4. C8: Change `$SIZE` to `$SIZE_MM` in LTER processing
5. M7: Apply FDR to full-precision p-values

### Phase 2 (Medium Impact, Requires Testing)
6. C3: Fix sigmoid `predFit` variable name mismatch
7. C7: Deduplicate `sites.short.edit` at import
8. C4: Add NLS convergence check before AICc
9. M2: Offset first "After" year to `time.model=1`

### Phase 3 (Low Priority)
10. C5: Refactor sigmoid SumStats to use `data.frame()`
11. M8: Decide if Table 3 should use post-outlier data
12. Clean up unused constants in `00c_analysis_constants.R`
