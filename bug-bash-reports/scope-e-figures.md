# Scope E: Figures Module Audit Report

**File:** `/Users/adrianstier/Donham-Stier-CA-MPA-2026/code/R/10_figures.R` (2,538 lines)
**Auditor:** Claude Opus 4.6
**Date:** 2026-02-09

---

## Executive Summary

The figures module is well-structured with robust error handling, comprehensive input validation, and consistent styling. I identified **2 moderate issues**, **8 minor issues**, and **0 critical (pipeline-breaking) issues**. The code successfully defends against common failure modes (missing packages, missing data objects, cairo_pdf failures). The primary concerns are: (1) stale PISCO short codes in Figure 4 raw data filters that will never match after upstream standardization, and (2) incorrect fallback MPA implementation years in Figure S2.

---

## 1. save_fig() Function (Lines 277-380)

### Findings

**MINOR-01: Patchwork detection heuristic may miss some patchwork objects**
- Line 282: `is_patchwork <- inherits(plot, "patchwork") || (inherits(plot, "gg") && !is.null(plot$patches))`
- The `$patches` field is an internal implementation detail of patchwork. Some patchwork compositions (e.g., using `wrap_plots()`) may not expose `$patches` on the outer object. In practice this only affects Strategy 1 selection; the fallback strategies (2-4) will still work.
- **Impact:** Low. The fallback strategies handle this gracefully.

**No issues found with:**
- cairo_pdf failure handling: Four strategies attempted in order (patchwork-specific pdf(), cairo_pdf via ggsave, standard pdf via ggsave, raw pdf() device). Each has error catching.
- Output file verification: Both PDF and PNG existence and non-zero size are verified (lines 355-356).
- Old file removal before save to prevent stale artifacts (line 285).
- PNG is treated as critical (line 366-368), PDF failure is a warning only (lines 369-376). This is a reasonable design decision.

---

## 2. should_render() Guards (Lines 264-273)

### All 10 figure sections are properly guarded:

| Figure | Guard | Line |
|--------|-------|------|
| Fig 1 | `should_render("fig01")` | 487 |
| Fig 2 | `should_render("fig02")` | 926 |
| Fig 3 | `should_render("fig03")` | 1244 |
| Fig 4 | `should_render("fig04")` | 1348 |
| Fig S1 | `should_render("fig_s01")` | 1137 |
| Fig S2 | `should_render("fig_s02")` | 1731 |
| Fig S3 | `should_render("fig_s03")` | 1942 |
| Fig S4 | `should_render("fig_s04")` | 2133 |
| Fig S5 | `should_render("fig_s05")` | 2222 |
| Fig S6 | `should_render("fig_s06")` | 2338 |

### Shared variables defined outside guards:

The following are correctly defined at top-level scope (outside any guard), ensuring they are available regardless of which figures are rendered:

- `taxa_levels` (line 79)
- `source_levels` (line 81)
- `excluded_mpas` (lines 104-106)
- `trophic_assignment` (lines 109-120)
- `taxa_col` (set at lines 254-260)
- `save_fig()`, `standardize_status()`, `shorten_mpa_name()` (helper functions)
- Dimension constants `FIG2_DIMS` through `FIG_S2_DIMS` (lines 94-98)

**MINOR-02: Dimension constants for S3, S4, S5 defined inside guards**
- `FIG_S3_DIMS` (line 1945, inside `fig_s03` guard)
- `FIG_S4_DIMS` (line 2136, inside `fig_s04` guard)
- `FIG_S5_DIMS` (line 2225, inside `fig_s05` guard)
- These are only used within their respective guards, so this is functionally correct. However, it breaks the convention established at lines 87-98 where `FIG2_DIMS` through `FIG_S2_DIMS` are defined in a central "dimension constants" section. If another section ever needed to reference these dimensions (e.g., a summary table of figure sizes), it would fail.
- **Impact:** Negligible. Cosmetic inconsistency only.

---

## 3. Data Dependencies and Column Name Mismatches

### Figure 2 (Lines 926-1130)

**No issues found.** The code correctly:
- Filters `All.Resp.sub` using `source == "KFM"`, `CA_MPA_Name_Short == "Scorpion SMR"`, and `taxon_name` (the correct column for this dataframe).
- Filters `All.RR.sub.trans` using `y` (the correct taxa column for this dataframe).
- References `lnDiff` (correct column name in `All.RR.sub.trans`).
- Calls `standardize_status()` before grouping (line 957), fixing a previously documented bug.

### Figure 3 (Lines 1244-1339)

**No issues found.** The code correctly:
- Uses `Table2` with columns `Taxa`, `Response`, `Estimate`, `CI_lower`, `CI_upper` (validated at lines 225-229).
- Uses `SumStats.Final` with columns `Taxa`, `MPA`, `Mean`, `SE`, `Resp` (validated at lines 235-239).
- Converts `Resp` ("Den"/"Bio") to `Response` ("Density"/"Biomass") at line 1263 for consistency with `Table2`.

### Figure 4 (Lines 1348-1724)

**MODERATE-01: Stale PISCO short codes in raw data filters that will never match**
- Lines 1481, 1489, 1527, 1535: The raw data filters for panels B and D include PISCO short codes `"STRPURAD"` and `"MACPYRAD"` in the `y %in% c(...)` filter:
  ```r
  y %in% c("Strongylocentrotus purpuratus", "S. purpuratus", "STRPURAD")
  y %in% c("Macrocystis pyrifera", "M. pyrifera", "MACPYRAD")
  ```
- After upstream standardization in `07_combine_data.R` (line 142: `standardize_species_names()`), these PISCO codes are converted to full scientific names. The codes `"STRPURAD"` and `"MACPYRAD"` should never appear in `All.RR.sub.trans$y` at this point.
- Similarly, abbreviated names like `"S. purpuratus"` and `"M. pyrifera"` should never appear in `All.RR.sub.trans$y` since standardization converts to full names only (e.g., "Strongylocentrotus purpuratus").
- **Impact:** These dead filter values do not cause errors (they simply match zero rows), but they give a false impression that the code handles un-standardized data. They are misleading to future maintainers. Panels A and C (lines 1454-1461, 1502-1509) do not include PISCO short codes, creating an inconsistency within the same figure.
- **Data correctness risk:** None. The filters are inclusive (OR-based), and the full scientific names will always match. No data is lost.

**MINOR-03: Inconsistent taxa filtering between Figure 4 panels**
- Panels A and B (biomass): Include both urchin species (`S. purpuratus` + `M. franciscanus`).
- Panels C and D (density): Include only `S. purpuratus` for urchin density.
- This appears intentional (red urchin density may not be available from all sources), but is not documented in any comment.

### Figure S1 (Lines 1137-1237)

**No issues found.** Correctly uses `SumStats.Final` columns `Taxa`, `MPA`, `Mean`, `SE`, `CI`, `Source`, `Resp`, `AnalysisType`.

### Figure S2 (Lines 1731-1931)

**MODERATE-02: Incorrect fallback MPA implementation years**
- Line 1740: `MPA_Start = c(2012, 2005, 2005)` for `c("Naples SMCA", "Scorpion SMR", "Anacapa Island SMR 2003")`
- The fallback year for Scorpion SMR and Anacapa Island SMR 2003 is `2005`, but both Channel Islands MPAs were implemented in `2003` (as confirmed by `08_effect_sizes.R` line 968 and the Figure 2 fallback at line 934).
- This fallback is used only if the `Site` table does not contain these MPAs (which is unlikely in normal operation, but the fallback exists for robustness).
- **Impact:** If the Site table lookup fails for these MPAs, the MPA implementation line in the time series panels would be plotted 2 years too late, misrepresenting the intervention timing. The post-MPA LOESS smoothers would also start 2 years too late (line 1841: `filter(d, year >= mpa_start_yr)`).

### Figure S3 (Lines 1942-2125)

**No issues found.** Correctly uses:
- `All.RR.sub.trans` with `BA == "After"`, `time >= 0`
- `trophic_assignment[.data[[taxa_col]]]` where `taxa_col = "y"` -- maps full scientific names to trophic levels
- All referenced columns (`lnDiff`, `BA`, `time`, `CA_MPA_Name_Short`) exist in `All.RR.sub.trans`

### Figure S4 (Lines 2128-2214)

**No issues found.** Similar data access pattern to S3, with correct column references.

### Figure S5 (Lines 2217-2328)

**MINOR-04: Model colors include "Sigmoid" but Model column may rarely contain it**
- Line 2238-2244: `col_model` defines colors for "Step", "Linear", "Asymptotic", "Sigmoid", "Mean".
- In `08_effect_sizes.R`, "Sigmoid" values are only assigned for specific KFM macrocystis cases (lines 1107, 1189, 1268). If sigmoid fits fail (which the code handles at line 1150), no sigmoid entries exist.
- **Impact:** None. The missing model type simply won't appear in the stacked bar. The color definition is harmless.

**MINOR-05: Variance component panel depends on `meta_biomass` and `meta_density` objects**
- Lines 2268-2298: Panel S5B only renders if both `meta_biomass` and `meta_density` exist.
- These are model objects from `09_meta_analysis.R`. They should always exist if the pipeline runs correctly. The `exists()` check provides correct fallback behavior.
- The fallback (lines 2322-2325) renders only panel A without a guide area, which is correct.

### Figure S6 (Lines 2331-2526)

**No issues found.** Correctly handles the loop over taxa with proper data filtering and fallback response type selection.

---

## 4. Axis Limits and Scales

### Data Clipping Assessment

**Figure S1 (Forest plot):**
- Line 1209: `scale_x_continuous(limits = c(-x_limit, x_limit), breaks = x_breaks)` where `x_limit = max(abs(x_range_all)) * 1.1`. The 10% padding prevents clipping.

**Figure 3 (Mean effects):**
- Line 1322: `coord_cartesian(ylim = c(fig3_y_min - 0.9, fig3_y_max + 0.5))`. Uses `coord_cartesian` (not `scale_y_continuous limits`), so no data points are excluded from statistical calculations. The 0.9 and 0.5 unit padding are hardcoded rather than proportional, but are generous enough for typical lnRR ranges.
- **MINOR-06: Asymmetric padding in Fig 3 y-limits.** The bottom padding (0.9) includes room for n= labels, and the top padding (0.5) is just for breathing room. The label position at `fig3_y_min - fig3_pad` (where `fig3_pad = 0.4`) is within the bottom padding of 0.9, so labels will not be clipped.

**Figure 4 (Trophic scatter):**
- Line 1650: `coord_fixed(ratio = 1, xlim = lim_x, ylim = lim_y, expand = FALSE)`. The symmetric limits (`lim_val`) are calculated as `max(max_abs_x, max_abs_y)` with 15% padding (line 1607-1609). No data clipping.

**Figure S2 (Time series):**
- Line 1858: `scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 1))`. This is a hard limit. Any `lnDiff` values outside [-2, 2] will be silently dropped by `scale_y_continuous`.
- **MINOR-07: Hard y-axis limits in Fig S2 may clip extreme values.** If any lnRR values exceed +/-2, they will be removed before rendering (ggplot2 drops out-of-range points with `scale_*` limits, unlike `coord_cartesian`). The code should use `coord_cartesian(ylim = c(-2, 2))` instead to avoid data loss.

**Figure S4 (Heatmap):**
- Line 2191: `limits = c(-3, 3), oob = scales::squish`. Values outside [-3, 3] are squished to the boundary rather than dropped. This is correct.

**Figure S6 (Appendix):**
- Line 2473: `coord_cartesian(ylim = c(-y_lim, y_lim))` -- correctly uses `coord_cartesian`, so no data clipping.

---

## 5. Factor Levels Consistency

**taxa_levels** (line 79): `c("S. purpuratus", "M. franciscanus", "M. pyrifera", "P. interruptus", "S. pulcher")`
- These match the abbreviated names used in `SumStats.Final$Taxa` and `Table2$Taxa`.
- These match the `col_taxa` names from `00b_color_palette.R`.
- Consistent throughout Figures 3, S1, S5.

**source_levels** (line 81): `c("KFM", "LTER", "PISCO", "Landsat")`
- Match the `shape_source` names from `00b_color_palette.R`.
- Used correctly in Figures S1, S2.

**MINOR-08: Figure S2 factor levels include both full and abbreviated taxa names**
- Lines 1761-1767: Sets `y` factor levels to include both "Strongylocentrotus purpuratus" and "S. purpuratus" (and similarly for all five species).
- Since `All.RR.sub.trans$y` should only contain full names after standardization, the abbreviated names in the factor levels will never be observed. This is harmless but misleading.
- The code does later recode (lines 1774-1785) to abbreviated names via `forcats::fct_recode`, so this is more of a defensive measure.

---

## 6. Expression Rendering

All `expression()` calls in axis labels were reviewed:

| Location | Expression | Renders correctly? |
|----------|------------|-------------------|
| Line 844 | `expression(Kelp~biomass~(g~m^{-2}))` | Yes |
| Line 1048 | `expression(Density~(ind~m^{-2}))` | Yes |
| Line 2101 | `expression("Effect trajectory slope (lnRR" ~ yr^{-1} * ")")` | Yes |
| Line 2289 | `expression(tau^2)` | Yes |

**No issues found.** All expressions use standard `plotmath` syntax with `~` for spacing and `^{}` for superscripts.

---

## 7. Patchwork Assembly

### Figure 1 (Map + panels)
- Lines 890-902: `(panels_row) / main_map` assembled correctly with `plot_layout(heights = c(2.4, 0.6))`.
- `guides = "collect"` used on panels_row (line 891) to merge time series legends.
- `& theme(...)` applied to panels_row (lines 892-897) to enforce L-shaped axes.

### Figure 2 (Data processing pipeline)
- Lines 1112-1126: 4 panels assembled with `plot_layout(ncol = 4, guides = "collect")`.
- `& theme(...)` correctly propagates axis styling.

### Figure 4 (Trophic cascade)
- Lines 1712-1721: `(panel_A | panel_B) / (panel_C | panel_D)` with `plot_annotation(tag_levels = "a")`.
- Correct use of lowercase "a" tag level with parenthetical formatting.

### Figure S2 (All taxa time series)
- Lines 1905-1914: Three panels stacked vertically with `guide_area()` for explicit legend placement.
- `heights = c(1, 1, 1, 0.16)` correctly allocates space for three equal panels and a thin guide area.

### Figure S3 (Temporal dynamics)
- Lines 2110-2118: `panel_S3A / (panel_S3B | panel_S3C)` with `plot_annotation(tag_levels = "a")`.
- Correct 3-panel layout.

### Figure S5 (Statistical transparency)
- Lines 2300-2320: `panel_S5A / panel_S5B + patchwork::guide_area()`.
- `heights = c(1, 0.8, 0.12)` for three elements (A, B, guide_area).
- The `&` operator at line 2305 propagates theme to all elements.

**No patchwork assembly issues found.**

---

## 8. Dimension Constants

| Figure | Constant | Width (cm) | Height (cm) | Notes |
|--------|----------|-----------|------------|-------|
| Fig 1 | `FIG1_DIMS` | 17 | 12 | Double-column max (correct) |
| Fig 2 | `FIG2_DIMS` | 17 | 7.5 | Double-column (correct for 4-panel row) |
| Fig 3 | `FIG3_DIMS` | 14 | 10 | Single-column+ (reasonable) |
| Fig 4 | `FIG4_DIMS` | 17 | 17 | Square 2x2 panel (correct for `coord_fixed`) |
| Fig S1 | `FIG_S1_DIMS` | 17.8 | 22 | CL max width (correct) |
| Fig S2 | `FIG_S2_DIMS` | 17.8 | 22 | CL max width (correct) |
| Fig S3 | `FIG_S3_DIMS` | 17 | 17 | Square (reasonable for 3-panel layout) |
| Fig S4 | `FIG_S4_DIMS` | 17 | 14 | Standard (correct) |
| Fig S5 | `FIG_S5_DIMS` | 17 | 11 | Standard (correct) |
| Fig S6 | Dynamic | min(17, n_cols*4) | max(10, n_rows*4) | Scales with data (correct) |

All widths are within the Conservation Letters double-column maximum of 178mm (~17.8cm). Single-column figures (Fig 3 at 14cm) are within the 85mm single-column limit equivalent. No issues found.

---

## 9. Dead Code

**Line 1933-1935: Comment referencing removed Figure 5**
```r
# (Figure 5 removed -- trophic cascade now shown in Figure 4 with dual layers)
# The following block is commented out but preserved for reference.
```
There is no actual commented-out code block below this comment. The comment is orphaned -- it references code that was apparently already deleted but the comment was left behind.
- **MINOR-09: Orphaned comment about removed Figure 5** (lines 1933-1935). No actual dead code, just a stale comment.

**No `if(FALSE)` blocks found.** No large commented-out code sections.

---

## Summary of All Issues

### CRITICAL (0)
None found. The pipeline will not break due to any issues in this file.

### MODERATE (2)

| ID | Description | Lines | Impact |
|----|-------------|-------|--------|
| MOD-01 | Stale PISCO short codes (`STRPURAD`, `MACPYRAD`) and abbreviated names in Figure 4 raw data filters will never match post-standardization | 1481, 1489, 1527, 1535 | No data loss (inclusive OR filter), but misleading to maintainers and inconsistent with panels A/C which omit these codes |
| MOD-02 | Incorrect fallback MPA implementation years for Scorpion SMR and Anacapa Island SMR 2003 in Figure S2 (2005 instead of 2003) | 1740 | If Site table lookup fails, MPA implementation lines and LOESS smoothers will start 2 years too late |

### MINOR (9)

| ID | Description | Lines | Impact |
|----|-------------|-------|--------|
| MIN-01 | Patchwork detection heuristic relies on internal `$patches` field | 282 | Low; fallback strategies compensate |
| MIN-02 | Dimension constants for S3/S4/S5 defined inside guards, breaking convention | 1945, 2136, 2225 | Cosmetic inconsistency only |
| MIN-03 | Inconsistent urchin species filtering between biomass and density panels in Fig 4 (undocumented) | 1423-1441 | Likely intentional but not commented |
| MIN-04 | `col_model` includes "Sigmoid" which may rarely appear in data | 2241 | Harmless unused color |
| MIN-05 | Variance component panel silently missing if meta objects don't exist | 2268 | Correct fallback behavior, just undocumented to user |
| MIN-06 | Asymmetric padding in Fig 3 y-limits (0.9 bottom vs 0.5 top) is hardcoded | 1322 | Works for typical lnRR ranges but fragile |
| MIN-07 | Fig S2 uses `scale_y_continuous(limits=...)` instead of `coord_cartesian()`, silently dropping extreme lnRR values | 1858 | Potential data loss for lnRR outside [-2, 2] |
| MIN-08 | Fig S2 factor levels include both full and abbreviated taxa names (redundant) | 1761-1767 | Harmless but misleading |
| MIN-09 | Orphaned comment about removed Figure 5 with no accompanying dead code | 1933-1935 | Cosmetic |

---

## Positive Observations

1. **Robust save_fig()**: Four fallback strategies for PDF generation with file existence verification is excellent defensive coding.
2. **Comprehensive input validation**: All required data objects, columns, and palette objects are checked before any figure construction begins.
3. **Selective rendering**: The `should_render()` / `RENDER_FIGURES` system allows targeted figure regeneration, which is valuable during iterative development.
4. **Consistent theming**: `theme_mpa()` and the scale helper functions from `00b_color_palette.R` ensure visual consistency across all figures.
5. **L-shaped axes enforcement**: The `& theme(panel.border = element_blank(), axis.line = element_blank(), axis.line.x.bottom = ..., axis.line.y.left = ...)` pattern is consistently applied via patchwork's `&` operator.
6. **Status standardization**: `standardize_status()` is called before grouping in Figure 2, preventing the documented "separate groups for control/reference/outside" bug.
7. **Dynamic axis limits**: Most figures calculate axis limits from data rather than using hardcoded values (exception: Fig S2 y-limits).
