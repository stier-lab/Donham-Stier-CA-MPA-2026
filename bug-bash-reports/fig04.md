# Figure 4 Bug Report: Inconsistent Axis Labels in Trophic Cascade Scatterplots

## File Modified

`/Users/adrianstier/Donham-Stier-CA-MPA-2026/code/R/10_figures.R`

## Issue

Figure 4 is a 2x2 panel figure showing trophic cascade relationships:
- Top row (panels a, b): biomass effects
- Bottom row (panels c, d): density effects

The top-row panels (a and b) used generic taxon labels (e.g., "Predator biomass effect (lnRR)", "Urchin biomass effect (lnRR)"). However, the bottom-row panels (c and d) mixed generic and species-specific labels with italicized binomial names. Specifically:

- **Panel C y-axis** used `expression(italic("S. purpuratus") * " density effect (lnRR)")` instead of `"Urchin density effect (lnRR)"`
- **Panel D x-axis** used `expression(italic("S. purpuratus") * " density effect (lnRR)")` instead of `"Urchin density effect (lnRR)"`

This created a visual inconsistency: the same conceptual axis (urchin effect) was labeled differently depending on which panel it appeared in and whether it represented biomass or density. The inconsistency was flagged by review agents.

## Root Cause

The density panels were constructed at a different time than the biomass panels, and the developer used species-specific italic labels (`expression(italic("S. purpuratus") * ...)`) for the urchin axes in the density row, while the biomass row used plain generic labels. The predator x-axis in Panel C was already generic ("Predator density effect (lnRR)"), making the inconsistency partial even within the density row.

## Fix Applied

Changed the two species-specific axis labels in the density panels to match the generic naming convention used in the biomass panels.

### Panel C (line 1682)

**Before:**
```r
expression(italic("S. purpuratus") * " density effect (lnRR)")
```

**After:**
```r
"Urchin density effect (lnRR)"
```

### Panel D (line 1691)

**Before:**
```r
expression(italic("S. purpuratus") * " density effect (lnRR)")
```

**After:**
```r
"Urchin density effect (lnRR)"
```

## Result

All four panels now use consistent generic taxon labels:

| Panel | x-axis | y-axis |
|-------|--------|--------|
| (a) Biomass: Predator vs Urchin | Predator biomass effect (lnRR) | Urchin biomass effect (lnRR) |
| (b) Biomass: Urchin vs Kelp | Urchin biomass effect (lnRR) | Kelp biomass effect (lnRR) |
| (c) Density: Predator vs Urchin | Predator density effect (lnRR) | Urchin density effect (lnRR) |
| (d) Density: Urchin vs Kelp | Urchin density effect (lnRR) | Kelp biomass effect (lnRR) |

The naming convention is now uniform: `{Generic taxon} {response metric} effect (lnRR)` across all panels and axes.

## Verification

The fix requires re-rendering Figure 4 by running the `fig04` section of `10_figures.R`. No other figures or code are affected.
