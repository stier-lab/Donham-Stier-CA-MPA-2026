# Bug-Bash Report: Figures S1, S3, S4

**File modified:** `/Users/adrianstier/Donham-Stier-CA-MPA-2026/code/R/10_figures.R`
**Date:** 2026-02-09

---

## Figure S1 -- Forest Plot: Missing x-axis tick labels on top-row panels

### Problem

The forest plot uses `facet_wrap(~ Taxa, ncol = 2, scales = "free_y")` which, by default in ggplot2, only renders x-axis tick labels on the bottom row of faceted panels. The top-row panels (S. purpuratus, M. franciscanus) lacked x-axis tick labels entirely, making it impossible for readers to gauge effect-size magnitudes in those panels without visually tracing down to the bottom row.

### Root cause

`facet_wrap()` in ggplot2 defaults to drawing axis ticks and labels only on the outermost panels (bottom for x-axis, left for y-axis). Since `scales = "free_y"` was set (not `"free"` or `"free_x"`), the x-axis is shared across all panels -- but the labels were still suppressed on interior rows by default layout logic.

### Fix (line 1199)

Added `axes = "all_x"` to the `facet_wrap()` call:

```r
facet_wrap(~ Taxa, ncol = 2, scales = "free_y", axes = "all_x") +
```

The `axes` parameter (available since ggplot2 3.5.0; this project uses 4.0.1) controls whether axis ticks and labels are drawn on non-edge panels. `"all_x"` draws x-axis ticks and labels on every panel row.

### Scope

Single-line change. No behavioral side effects -- the x-axis scale and breaks remain the same (`scale_x_continuous(limits = c(-x_limit, x_limit), breaks = x_breaks)`), they are simply now visible in all panels.

---

## Figure S3 -- Temporal Dynamics: Predators/Kelp color indistinguishable

### Problem

The trophic-level color mapping uses:
- Predators: `#5C8A70` (muted sage green)
- Kelp: `#4A7C59` (forest green, sourced from `col_taxa["M. pyrifera"]`)

These two greens differ by only ~18 units in hue and ~10 units in lightness (CIELAB), making them very difficult to distinguish, particularly for readers with deuteranopia (red-green color blindness). Both lines also used the same solid linetype, providing no secondary visual cue.

### Root cause

The `trophic_colors` named vector assigned two similar green hues to different trophic levels. Changing colors would be a larger design decision (affects palette consistency across figures), but adding linetype differentiation is a minimal, non-breaking way to resolve the ambiguity.

### Fix (lines 1974-1996, 2035-2041)

1. Defined a `trophic_linetypes` named vector just inside the `if` block (line 1976):

```r
trophic_linetypes <- c("Predators" = "solid", "Urchins" = "solid", "Kelp" = "dashed")
```

2. Added `linetype = Trophic_Level` to `aes()` in Panel A (line 1980-1981):

```r
aes(x = time, y = mean_lnRR, color = Trophic_Level,
    fill = Trophic_Level, linetype = Trophic_Level)
```

3. Added `scale_linetype_manual()` to Panel A (line 1996):

```r
scale_linetype_manual(values = trophic_linetypes, name = "Trophic Level", drop = FALSE) +
```

4. Added the same `linetype` aesthetic and scale to Panel B (lines 2035-2041):

```r
aes(x = time, y = cumulative_mean, color = Trophic_Level,
    linetype = Trophic_Level)
...
scale_linetype_manual(values = trophic_linetypes, name = "Trophic Level") +
```

### Design choice

Kelp gets a dashed line because it is the indirect/bottom trophic level (the cascade endpoint), while Predators and Urchins (the direct actors) use solid lines. This semantic encoding reinforces the trophic-cascade narrative: direct effects are solid, indirect effects are dashed.

### Legend note

Because `color`, `fill`, and `linetype` all share the same `name = "Trophic Level"`, ggplot2 will merge them into a single legend entry per level. The legend keys will show both the color and the linetype, improving accessibility.

### Scope

Panel C (slope boxplot) is not affected -- it uses points and pointranges, not lines, so linetype is not relevant there.

---

## Figure S4 -- Space-Time Heatmap: Ambiguous white cells

### Problem

The diverging color scale uses `mid = "white"` at `midpoint = 0`, meaning cells with a near-zero effect size appear white. Cells with no data (missing MPA-year combinations) also appear as empty/white background, making it impossible to distinguish "no data" from "effect near zero."

### Root cause

The plot only draws `geom_tile()` for MPA-year combinations that exist in `heatmap_data`. Missing combinations produce no tile at all, leaving the panel background showing through -- which is also white (from `theme_mpa`).

### Fix (lines 2166-2203)

1. Built a complete MPA x year background grid (lines 2166-2174):

```r
all_mpas <- levels(heatmap_data$CA_MPA_Name_Short)
all_times <- seq(min(heatmap_data$time), max(heatmap_data$time))
bg_grid <- tidyr::expand_grid(
  CA_MPA_Name_Short = factor(all_mpas, levels = all_mpas),
  time = all_times
)
```

2. Added a background `geom_tile` layer with `fill = "grey90"` painted first (line 2178):

```r
geom_tile(data = bg_grid, fill = "grey90", color = "white", linewidth = 0.3) +
```

3. The existing data `geom_tile` (now with explicit `aes(fill = mean_lnRR)`) paints over the grey background wherever data exists (line 2180).

4. Set `na.value = "grey90"` in `scale_fill_gradient2()` for consistency (line 2189).

5. Added a caption annotation `"Grey cells = no data available"` (line 2195) and styled it at the bottom-left (line 2202):

```r
plot.caption = element_text(size = 7, color = "grey40", hjust = 0)
```

### Visual result

- Cells with data: colored on the green-white-purple diverging scale (near-zero effects appear white)
- Cells without data: light grey (`grey90`)
- The grey is visually distinct from the white midpoint, resolving the ambiguity

### Scope

The background grid only covers the range `min(time)` to `max(time)` for the MPAs present in the data, matching the natural extent of the heatmap. No extra whitespace or phantom rows are introduced.

---

## Summary of changes

| Figure | Line(s) | Change | Risk |
|--------|---------|--------|------|
| S1 | 1199 | Added `axes = "all_x"` to `facet_wrap()` | None -- purely additive display change |
| S3 (A) | 1974-1996 | Added `trophic_linetypes` vector, `linetype` aes, `scale_linetype_manual()` | None -- additive visual encoding, legend merges automatically |
| S3 (B) | 2035-2041 | Added `linetype` aes and `scale_linetype_manual()` | None -- same pattern as Panel A |
| S4 | 2166-2203 | Added grey background grid, `na.value`, caption | None -- background layer is painted first, data layer overpaints |

All fixes are additive and non-breaking. No existing data processing, statistical calculations, or other figures are affected.
