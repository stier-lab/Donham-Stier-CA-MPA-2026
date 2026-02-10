# Scope A: Foundation Module Audit Report

**Auditor:** Claude Opus 4.6 (Scope Agent A)
**Date:** 2026-02-09
**Files audited:**
- `/Users/adrianstier/Donham-Stier-CA-MPA-2026/code/R/00_libraries.R` (189 lines)
- `/Users/adrianstier/Donham-Stier-CA-MPA-2026/code/R/00b_color_palette.R` (517 lines)
- `/Users/adrianstier/Donham-Stier-CA-MPA-2026/code/R/00c_analysis_constants.R` (147 lines)
- `/Users/adrianstier/Donham-Stier-CA-MPA-2026/code/R/01_utils.R` (874 lines)

---

## CRITICAL Issues (Breaks pipeline or produces wrong results)

### C1. `EXCLUDED_MPAS` silently overwritten in `08_effect_sizes.R`

**Location:** `01_utils.R` line 628 defines `EXCLUDED_MPAS` with 11 entries. `08_effect_sizes.R` line 1426 redefines it with only 5 entries.

**Problem:** The foundation module defines `EXCLUDED_MPAS` as the authoritative exclusion list containing 11 MPAs (including "Crystal Cove SMCA", "Laguna Beach SMR", "South La Jolla SMR", "Vandenberg SMR", "N/A", "Anacapa Island SMR 1978", "Anacapa Island SMCA", "Carrington Point SMR"). However, `08_effect_sizes.R` overwrites this with a completely different, shorter list of 5 MPAs (adding "San Miguel Island SC" and "Judith Rk SMR" which are NOT in the original list).

This means:
- The data filtering in `04_pisco_processing.R` (line 720) uses the original 11-MPA list.
- The FilterAudit in `08_effect_sizes.R` (line 1429) uses the replacement 5-MPA list.
- The two exclusion lists are **inconsistent** -- different MPAs are excluded at different pipeline stages.
- "Crystal Cove SMCA", "Laguna Beach SMR", "South La Jolla SMR", "Vandenberg SMR", and "Carrington Point SMR" are excluded from raw data processing but **not** from the FilterAudit transparency check, so the audit will not flag them as excluded.
- Conversely, "San Miguel Island SC" and "Judith Rk SMR" appear only in the `08_effect_sizes.R` list but never in the foundation list, so they are excluded from the audit but never from the raw data processing stage.

**Impact:** The FilterAudit (a statistical transparency tool) does not accurately reflect the actual exclusion criteria applied upstream. This could mislead reviewers about which MPAs were removed and why. Additionally, if `EXCLUDED_MPAS` is referenced after `08_effect_sizes.R` runs (e.g., in `09_meta_analysis.R` or `10_figures.R`), the overwritten value would be used.

**Severity:** CRITICAL -- potential for inconsistent data filtering and misleading transparency reporting.

### C2. `cols` color palette overwritten in `08_effect_sizes.R` with different colors

**Location:** `00b_color_palette.R` line 68 defines `cols <- col_taxa_short` (the colorblind-safe palette). `08_effect_sizes.R` line 58 redefines it with completely different hex codes:
- Foundation: `Sheephead = "#2A7B8E"`, `Kelp = "#4A7C59"`, `Lobs = "#D4933B"`, `Purps = "#7B68A6"`, `Reds = "#B85A4C"`
- 08_effect_sizes.R: `Sheephead = "#377eb8"`, `Kelp = "#4daf4a"`, `Lobs = "#ff7f00"`, `Purps = "#984ea3"`, `Reds = "#e41a1c"`

**Problem:** The replacement colors are from the ColorBrewer `Set1` palette, which is NOT the colorblind-safe palette designed in `00b_color_palette.R`. The `#e41a1c` red and `#4daf4a` green will be indistinguishable under deuteranopia (red-green color blindness). The entire purpose of the custom palette in `00b_color_palette.R` was to avoid this exact problem.

**Impact:** Any diagnostic plots generated during effect size calculation will use the non-colorblind-safe palette. While these may be internal-only plots, the overwrite also means any code running after `08_effect_sizes.R` that references `cols` will get the wrong colors.

**Severity:** CRITICAL -- defeats the colorblind-safety design goal and creates palette inconsistency.

---

## MODERATE Issues (Could affect results or maintainability)

### M1. Survey area constants defined but hardcoded magic numbers used downstream

**Location:** `00c_analysis_constants.R` lines 25-34 define `PISCO_SWATH_AREA_M2 <- 60`, `KFM_QUAD_AREA_M2 <- 2`, `LTER_LOBSTER_PLOT_AREA_M2 <- 300`, `LTER_FISH_SURVEY_AREA_M2 <- 80`.

**Problem:** These constants are NEVER used anywhere in the codebase. Instead, the magic numbers are hardcoded directly:
- `04_pisco_processing.R` lines 185, 258, 396, 559, 673, 695, 821: `/ 60`
- `05_kfm_processing.R` lines 311-312: `/ 2`
- `06_lter_processing.R` lines 479-480: `/ 300`
- `06_lter_processing.R` lines 680-681: `/ 80`

**Impact:** The constants exist to centralize these values, but the downstream code ignores them. If a survey area value ever needs to change, editing the constant would have no effect -- the hardcoded values would need to be found and updated individually, with risk of missing some.

**Severity:** MODERATE -- no current bug, but a significant maintainability hazard. The constant definitions give false confidence that changing them would propagate downstream.

### M2. Sheephead allometric constants defined but hardcoded magic numbers used

**Location:** `00c_analysis_constants.R` lines 140-141 define `SHEEPHEAD_BIOMASS_ALLOMETRIC_A <- 0.0144` and `SHEEPHEAD_BIOMASS_ALLOMETRIC_B <- 3.04`.

**Problem:** These constants are NEVER referenced downstream. Instead, `06_lter_processing.R` line 621 hardcodes: `lter.fish.sub.site$biomass <- 0.0144 * (lter.fish.sub.site$SIZE ^ 3.04)`. Same issue as M1 -- the constant exists but is not used.

**Severity:** MODERATE -- same maintainability concern as M1.

### M3. Bootstrap seed constants and `N_BOOTSTRAP_RESAMPLES` defined but not used

**Location:** `00c_analysis_constants.R` lines 86-94 define `N_BOOTSTRAP_RESAMPLES <- 1000` and five `SEED_*` constants.

**Problem:**
- `N_BOOTSTRAP_RESAMPLES` is never referenced. The default `n_resamples = 1000` is hardcoded in `bootstrap_biomass()` (line 269 of `01_utils.R`) and also hardcoded in `06_lter_processing.R` line 214.
- The `SEED_*` constants are never referenced by name. Instead, the same numeric values are hardcoded directly: `04_pisco_processing.R` uses `set.seed(12345)` (line 311) and `set.seed(12346)` (line 471); `05_kfm_processing.R` uses `set.seed(12347)` (line 223) and `set.seed(12348)` (line 395).
- The bootstrap scripts do NOT use `run_cached_bootstrap()` which would properly use the seed constants. They roll their own caching logic with inline `set.seed()` calls.

**Impact:** The `run_cached_bootstrap()` utility function was designed to centralize caching and seeding logic, but it is never called. The seeds happen to match the constants, but this is fragile -- editing the constants would not change the actual seeds used.

**Severity:** MODERATE -- same pattern as M1/M2, but multiplied across more constants.

### M4. Many statistical constants in `00c_analysis_constants.R` are never referenced

**Location:** `00c_analysis_constants.R` Section 5 (lines 98-129).

**Problem:** The following constants are defined but never used anywhere in the codebase:
- `Z_CRITICAL_95_CI <- 1.96` -- `1.96` is hardcoded in 8+ locations across `08_effect_sizes.R` and `10_figures.R`
- `MIN_OBS_SHAPIRO_TEST <- 3` -- never referenced
- `MAX_OBS_SHAPIRO_TEST <- 5000` -- never referenced
- `DHARMA_N_SIMULATIONS <- 250` -- never referenced
- `MAX_HETEROSCEDASTICITY_CORRELATION_LM <- 0.5` -- never referenced
- `MAX_HETEROSCEDASTICITY_CORRELATION_NLS <- 0.6` -- never referenced
- `MAX_OUTLIERS_NLS <- 1` -- never referenced
- `NLS_NORMALITY_THRESHOLD <- 0.01` -- never referenced
- `COOKS_DISTANCE_NUMERATOR <- 4` -- never referenced

**Impact:** These constants were created to centralize diagnostic thresholds, but `08_effect_sizes.R` (which is where these thresholds would be used) does not reference them. The entire Section 5 of `00c_analysis_constants.R` appears to be aspirational infrastructure that was never wired up.

**Severity:** MODERATE -- no current bug, but the constants file overstates what is actually centralized.

### M5. Temporal constants `SURVEY_MONTH_START/END` not used

**Location:** `00c_analysis_constants.R` lines 55-56.

**Problem:** `SURVEY_MONTH_START <- 5` and `SURVEY_MONTH_END <- 10` are defined but never referenced. Instead, `04_pisco_processing.R` line 113 hardcodes: `month >= 5 & month <= 10`.

**Severity:** MODERATE -- same pattern as above.

### M6. `MINIMUM_YEARS_OF_DATA` constant not used

**Location:** `00c_analysis_constants.R` line 68.

**Problem:** `MINIMUM_YEARS_OF_DATA <- 5` is defined but never referenced. `04_pisco_processing.R` line 704 hardcodes `dplyr::filter(n_years >= 5)`.

**Severity:** MODERATE -- same pattern.

### M7. `PISCO_URCHIN_MIN_SIZE_MM` constant not used

**Location:** `00c_analysis_constants.R` line 45.

**Problem:** `PISCO_URCHIN_MIN_SIZE_MM <- 25` is defined but never referenced. `04_pisco_processing.R` line 446 hardcodes: `subset(SizeFreq.Urch, size >= 25)`.

**Severity:** MODERATE -- same pattern.

### M8. `UCSB_LOBSTER_SIZE_START_YEAR` and `VRG_LOBSTER_SIZE_START_YEAR` not used

**Location:** `00c_analysis_constants.R` lines 60-61.

**Problem:** These constants (`2010` and `2011`) are never referenced. `04_pisco_processing.R` hardcodes `year > 2009` (line 207, equivalent to >= 2010) and `year > 2010` (line 277, equivalent to >= 2011).

Note: The constants use `>=` semantics (start year) while the code uses `>` (year minus one). This is a subtle mismatch -- if someone wired up the constants, they'd need to use `year >= UCSB_LOBSTER_SIZE_START_YEAR` not `year > UCSB_LOBSTER_SIZE_START_YEAR` to get equivalent behavior, but the off-by-one is not documented.

**Severity:** MODERATE -- unused constant with semantics mismatch.

### M9. `KFM_RDFC_SURVEY_START_YEAR` not used

**Location:** `00c_analysis_constants.R` line 64.

**Problem:** `KFM_RDFC_SURVEY_START_YEAR <- 2003` is defined but never referenced. `05_kfm_processing.R` line 665 hardcodes: `year >= 2003`.

**Severity:** MODERATE -- same pattern.

### M10. `assign_time_from_site_table` generates false positive warnings

**Location:** `01_utils.R` lines 103-109.

**Problem:** The warning logic checks `df$time == 0` to find "unmatched" MPAs, but legitimate MPAs with no post-implementation data would also have `time = 0` for all observations. This means the warning would incorrectly flag properly matched MPAs that simply have no "After" period data.

```r
unmatched <- unique(df$CA_MPA_Name_Short[df$time == 0])
```

This will report MPAs where ALL observations occurred before or during the implementation year, even if they were correctly matched in the site table.

**Impact:** Nuisance warnings that could mask real problems (cry-wolf effect), or could cause users to dismiss legitimate warnings about truly unmatched MPAs.

**Severity:** MODERATE -- misleading diagnostics.

### M11. `calculate_proportions` adaptive correction adds correction to ALL values including non-zero

**Location:** `01_utils.R` line 389.

**Problem:** The adaptive zero-correction adds the correction value to ALL proportions, not just the zeros:
```r
df$PropCorr[idx] <- df$Prop[idx] + correction
```

This means a non-zero proportion of, say, 0.5 becomes 0.5 + correction, while the maximum proportion of 1.0 becomes 1.0 + correction. This shifts the entire distribution upward, which is appropriate for the zero-correction purpose (ensuring no zero values before taking logs), but the comment "Use half of the minimum non-zero proportion" and the Aitchison (1986) citation suggest a multiplicative replacement was intended for the zeros only.

In the "fixed" method, the same issue exists: all values get +0.01, not just the zeros.

**Note:** This is likely intentional (adding a constant preserves ratios between MPA and reference when both get the same additive correction), but the documentation implies it is a zero-replacement strategy from compositional data analysis, which would typically only replace the zeros. The distinction matters because the additive approach slightly shrinks the range of the proportions.

**Severity:** MODERATE -- the behavior is consistent but the documentation is misleading about what it does and why.

### M12. `bootstrap_biomass` builds size population vector via repeated `c()` concatenation

**Location:** `01_utils.R` lines 282-286.

**Problem:**
```r
a <- numeric(0)
for (j in seq_along(t2)) {
  reps <- rep(size_freq_table$size[t2[j]], size_freq_table$count[t2[j]])
  a <- c(a, reps)
}
```

Growing a vector inside a loop with `c()` is a classic R performance anti-pattern. Each iteration copies the entire vector. For large size frequency tables, this could be very slow. Since the bootstrap is called many thousands of times across the pipeline, this could significantly affect total runtime.

A vectorized approach like `rep(size_freq_table$size[t2], size_freq_table$count[t2])` would be much faster.

**Severity:** MODERATE -- performance issue that compounds across many bootstrap iterations.

---

## MINOR Issues (Cosmetic, documentation, or low-risk)

### m1. `wesanderson` package loaded only for `pal` variable that is never used

**Location:** `00_libraries.R` line 113 loads `library(wesanderson)`. Line 183 creates `pal <- wes_palette("Zissou1", 25, type = "continuous")`.

**Problem:** The `pal` variable is never referenced anywhere in the codebase. The color palette system in `00b_color_palette.R` completely replaces it. Loading `wesanderson` adds an unnecessary package dependency.

**Severity:** MINOR -- unnecessary dependency, no functional impact.

### m2. `validated_merge` and `run_cached_bootstrap` utility functions are never called

**Location:** `01_utils.R` lines 831 and 786.

**Problem:** These two utility functions were created but are never used anywhere in the pipeline. The bootstrap scripts roll their own caching logic, and merge operations use `merge()` directly without validation.

**Severity:** MINOR -- dead code. No functional impact, but adds maintenance burden.

### m3. `shape_source_open`, `fill_ba`, and `alpha_ba` palette variables never used outside `00b_color_palette.R`

**Location:** `00b_color_palette.R` lines 197, 171, 176.

**Problem:** `shape_source_open` (open-point variants for Before period) is defined but never used outside the palette preview function. `fill_ba` is only used in the palette preview. `alpha_ba` is never used at all.

**Severity:** MINOR -- defined infrastructure that is not yet utilized in figures.

### m4. `col_type` palette only used in palette preview

**Location:** `00b_color_palette.R` line 218.

**Problem:** The analysis type colors (`pBACIPS`, `BACI`, `CI`) are defined but never used in any actual figure. Only used in the palette preview function.

**Severity:** MINOR -- no impact.

### m5. `label_source` only used inside `scale_shape_source`

**Location:** `00b_color_palette.R` line 205.

**Problem:** `label_source` is only used as the `labels` argument inside the `scale_shape_source()` convenience function. It is not independently referenced. This is fine, but means the labels "KFM (NPS)", "SBC LTER" etc. are only applied when using the convenience function, not when manually using `scale_shape_manual()`.

**Severity:** MINOR -- working as designed, but note for consistency.

### m6. `theme_mpa()` uses `theme_bw` base while CLAUDE.md specifies `theme_minimal` or `theme_classic`

**Location:** `00b_color_palette.R` line 246.

**Problem:** The global CLAUDE.md instructions say "Start from `theme_minimal(base_size = 11)` or `theme_classic(base_size = 11)`, never the default grey theme." The `theme_mpa()` function starts from `theme_bw()` with `base_size = 9`. While `theme_bw()` is not the "default grey theme" and the custom theme overrides most elements anyway, this is technically a deviation from the project's own guidelines.

**Severity:** MINOR -- `theme_mpa()` effectively creates a custom theme that achieves the same goals (clean, minimal, no gridlines), so this is a style inconsistency rather than a functional issue.

### m7. Comment on line 156 of `00_libraries.R` says `sf` is "not currently used" but it IS loaded

**Location:** `00_libraries.R` lines 156 and 123.

**Problem:** Line 156 says "# sf: Simple Features for spatial data (not currently used)" in the Miscellaneous Utilities section, but `library(sf)` is explicitly loaded on line 123 in the Mapping and Spatial section. The commented-out `library(sf)` on line 163 (with "Not currently used" note) further confuses matters. This is a documentation error -- `sf` IS used and IS loaded.

**Severity:** MINOR -- confusing documentation, no functional impact.

### m8. `Rmisc` package loads `plyr` functions, partially undermining the dplyr masking strategy

**Location:** `00_libraries.R` line 160.

**Problem:** The script carefully loads `plyr` first (line 31), then `tidyverse` (line 41), then re-loads `dplyr` (line 46) to ensure dplyr's `summarise()` takes precedence. However, `library(Rmisc)` on line 160 depends on and loads `plyr` internally. If `Rmisc` triggers any re-attachment of `plyr` functions, it could alter the masking order. In practice, since `plyr` is already loaded, the re-load is usually a no-op, but the intent is fragile. The `summarySE()` function from `Rmisc` itself uses `plyr::ddply()` internally, which can interact poorly with grouped tibbles.

**Severity:** MINOR -- the current loading order likely works, but the interaction is fragile and poorly documented.

### m9. `calculate_effect_size` function still has unreachable code after `stop()`

**Location:** `01_utils.R` lines 486-518.

**Problem:** The function calls `stop()` on line 487, so lines 497-518 (the actual computation) can never execute. This dead code after the `stop()` is confusing -- the function should either be removed entirely or the dead code should be removed, leaving only the deprecation notice.

**Severity:** MINOR -- dead code, no functional impact since `stop()` prevents execution.

### m10. `patchwork` and `scales` packages missing from `00_libraries.R`

**Location:** `00_libraries.R`.

**Problem:** `10_figures.R` requires both `patchwork` (line 485) and `scales` (line 154) packages, but neither is loaded in `00_libraries.R`. They are loaded ad-hoc within `10_figures.R`. While this works (since `10_figures.R` handles its own loading), it breaks the centralized dependency management pattern established by `00_libraries.R`.

**Severity:** MINOR -- works fine, but inconsistent with the centralization pattern.

### m11. `col_source_map` uses key "NPS-KFM" while `shape_source` and `label_source` use "KFM"

**Location:** `00b_color_palette.R` lines 145-149 vs lines 189-194.

**Problem:** Map marker colors use the key `"NPS-KFM"` but shape and label lookups use `"KFM"`. If a figure tries to map both color and shape by the same source variable, the keys will not match. The data would need to use one naming convention or the other.

**Severity:** MINOR -- different lookup tables can have different keys, but this creates a potential mismatch if both are used on the same data.

---

## Summary Table

| ID | Severity | Category | Description |
|----|----------|----------|-------------|
| C1 | CRITICAL | Cross-file | `EXCLUDED_MPAS` overwritten in 08_effect_sizes.R with different list |
| C2 | CRITICAL | Cross-file | `cols` palette overwritten with non-colorblind-safe colors |
| M1 | MODERATE | Constants | `PISCO_SWATH_AREA_M2` etc. defined but hardcoded `/60` used |
| M2 | MODERATE | Constants | Sheephead allometric constants defined but hardcoded |
| M3 | MODERATE | Constants | Bootstrap seeds/count defined but hardcoded |
| M4 | MODERATE | Constants | 9 statistical threshold constants never referenced |
| M5 | MODERATE | Constants | Survey month constants not used |
| M6 | MODERATE | Constants | `MINIMUM_YEARS_OF_DATA` not used |
| M7 | MODERATE | Constants | `PISCO_URCHIN_MIN_SIZE_MM` not used |
| M8 | MODERATE | Constants | Lobster start year constants not used + semantics mismatch |
| M9 | MODERATE | Constants | `KFM_RDFC_SURVEY_START_YEAR` not used |
| M10 | MODERATE | Logic | `assign_time_from_site_table` false positive warnings |
| M11 | MODERATE | Logic | `calculate_proportions` adds correction to all values, not just zeros |
| M12 | MODERATE | Performance | `bootstrap_biomass` uses `c()` concatenation in loop |
| m1 | MINOR | Dependencies | `wesanderson` loaded but `pal` never used |
| m2 | MINOR | Dead code | `validated_merge` and `run_cached_bootstrap` never called |
| m3 | MINOR | Dead code | `shape_source_open`, `fill_ba`, `alpha_ba` unused |
| m4 | MINOR | Dead code | `col_type` only used in preview |
| m5 | MINOR | Scope | `label_source` only used via convenience function |
| m6 | MINOR | Style | `theme_mpa()` uses `theme_bw`, not `theme_minimal`/`theme_classic` |
| m7 | MINOR | Docs | `sf` documented as "not currently used" but IS loaded |
| m8 | MINOR | Dependencies | `Rmisc` re-loads `plyr`, fragile masking |
| m9 | MINOR | Dead code | Unreachable code after `stop()` in deprecated function |
| m10 | MINOR | Dependencies | `patchwork` and `scales` missing from central library loader |
| m11 | MINOR | Naming | `col_source_map` uses "NPS-KFM" key vs "KFM" elsewhere |

---

## Overall Assessment

The foundation module is well-documented and thoughtfully designed, with proper separation of concerns across the four files. The color palette system (`00b_color_palette.R`) is particularly well-crafted with detailed colorblind-safety rationale.

However, there is a **systemic problem**: `00c_analysis_constants.R` defines 25+ named constants, but the vast majority are never actually referenced by downstream scripts. The downstream code continues to use hardcoded magic numbers. This creates a false sense of centralization -- editing the constants file would have no effect on the actual analysis. Of the constants defined:
- **Used downstream:** `EFFECT_SIZE_TIME_YEARS`, `SIGNIFICANCE_ALPHA` (only in `11_results_summary.R`)
- **Not used downstream:** All 4 survey area constants, all 5 bootstrap seeds, `N_BOOTSTRAP_RESAMPLES`, all 9 statistical thresholds, both sheephead allometric constants, both survey month constants, both lobster start year constants, `KFM_RDFC_SURVEY_START_YEAR`, `MINIMUM_YEARS_OF_DATA`, `PISCO_URCHIN_MIN_SIZE_MM`

The two CRITICAL issues (C1 and C2) involve downstream scripts overwriting foundation-level objects (`EXCLUDED_MPAS` and `cols`) with different values, creating inconsistency in the pipeline's data filtering and color scheme.
