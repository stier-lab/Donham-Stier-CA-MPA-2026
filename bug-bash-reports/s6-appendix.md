# Bug Fix Report: S6 Appendix Figure X-Axis Label Crowding

**Date:** 2026-02-09
**File:** `/Users/adrianstier/Donham-Stier-CA-MPA-2026/code/R/10_figures.R`
**Issue:** X-axis year labels crowding and overlapping in narrow facet panels of S6 appendix figures
**Status:** FIXED

---

## Problem Description

The Figure S6 appendix (site-level log response ratio time series) consists of 5 multi-panel figures, one per taxa:
- `fig_s06_appendix_pinterruptus.pdf`
- `fig_s06_appendix_spulcher.pdf`
- `fig_s06_appendix_spurpuratus.pdf`
- `fig_s06_appendix_mfranciscanus.pdf`
- `fig_s06_appendix_mpyrifera.pdf`

Each figure displays time series for multiple MPAs in a 4-column facet grid. The time series spans approximately 43 years (1980-2023), which creates a wide data range but physically narrow facet panels (width = total_width / 4).

**Root Cause:**
At line 2413 of `10_figures.R`, the code computed year axis breaks as:
```r
by_val <- if (year_span > 30) 10 else 5
```

For the ~43-year span, this produces `by_val = 10`, yielding 6 breaks: **1980, 1990, 2000, 2010, 2020, 2030**

In narrow facet panels (each ~170mm / 4 ≈ 42.5mm wide), 6 four-digit year labels render too densely, causing overlap and illegibility.

---

## Solution Applied

**Changed line 2413 from:**
```r
by_val <- if (year_span > 30) 10 else 5
```

**To:**
```r
by_val <- if (year_span > 30) 20 else 10
```

### Impact

This change reduces year breaks for long time series (>30-year span) from 10-year to 20-year intervals, producing 3 breaks instead of 6:
- **1980, 2000, 2020** (vs. 1980, 1990, 2000, 2010, 2020, 2030)

For shorter time series (≤30-year span), breaks decrease from 5-year to 10-year intervals.

### Rationale

1. **Spacing improvement:** Three breaks fit comfortably in narrow facet columns without overlap
2. **Readability:** 20-year intervals are still frequent enough to enable temporal interpretation
3. **Consistency:** Same method used in other figure types (e.g., Figure S2 panels) with appropriate interval scaling
4. **Small-multiple principle:** Readers can still compare temporal patterns across 20-year windows across MPAs in the same facet grid

---

## Verification

### Code Location
- **File:** `/Users/adrianstier/Donham-Stier-CA-MPA-2026/code/R/10_figures.R`
- **Lines:** 2408-2416 (S6 figure generation block)
- **Function scope:** Within `if (should_render("fig_s06"))` block (lines 2306-2494)

### Figure Scope
This fix applies to all 5 S6 appendix figures generated in the loop at line 2346:
```r
for (taxa_i in taxa_plot_order) {  # loops through 5 taxa
  # ... code that uses by_val ...
}
```

### Testing Notes
No other figures are affected by this change. Line 1814 (Figure S2) intentionally uses the earlier parameter settings and was not modified, as S2 panels are wider (not faceted into 4 columns) and do not exhibit label crowding.

---

## Files Modified

| File | Lines | Change |
|------|-------|--------|
| `code/R/10_figures.R` | 2413 | `10 else 5` → `20 else 10` |

---

## Impact on Manuscript

- **Affected figures:** Figure S6 (appendix site-level time series, 5 PDFs)
- **User-facing change:** Cleaner, more readable x-axis labels in narrow facets
- **No impact on:** Data, statistics, or analysis results
- **No impact on:** Other main-text or supplemental figures

---

## Commit Context

This fix is part of the s6-appendix visual refinement sprint, addressing readability issues identified during figure review for Conservation Letters submission.
