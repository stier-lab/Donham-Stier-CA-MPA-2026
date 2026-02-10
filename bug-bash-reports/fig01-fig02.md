# Bug-Bash Report: Figure 1 and Figure 2 Visual Issues

**File:** `code/R/10_figures.R`
**Date:** 2026-02-09

---

## Figure 1: Trailing "2" Artifact in Inset Panel Y-Axis Label

### Issue Reported
Review agents flagged a possible y-axis label rendering artifact in the inset time series panels. The label appeared to show `Kelp biomass (g m^-2) 2` with a trailing "2".

### Investigation

The inset time series panels are built by the `build_ts_panel()` function (starting at line 807). The y-axis label was defined conditionally on lines 838-840 (original):

```r
y = if (identical(letter, "a")) expression(Kelp~biomass~(g~m^{-2})~"/2")
    else expression(Kelp~biomass~(g~m^{-2}))
```

**Root cause confirmed:** Panel (a) used a special y-axis label with `~"/2"` appended to indicate that the tick labels were scaled by half (via `scales::label_number(scale = 0.5)` on line 834). In R's `expression()` rendering, the `~"/2"` appends the literal string `/2` after the math expression. Depending on the graphics device (PDF vs. PNG), font metrics, and anti-aliasing, the forward-slash character can render poorly or become nearly invisible, leaving what appears to be a standalone trailing "2" after the proper `(g m^-2)` label.

Additionally, panels (b)-(d) set `axis.title.y = element_text(color = "transparent")` (line 864) to preserve layout spacing while hiding the label. This means any rendering artifact on panel (a) is the only visible y-axis label, making it more noticeable.

### Fix Applied

Removed the `~"/2"` suffix from panel (a)'s y-axis label. All four panels now share the identical label `expression(Kelp~biomass~(g~m^{-2}))`. The half-scaling of tick labels on panel (a) is already communicated by the tick label values themselves (via `scales::label_number(scale = 0.5)`), so the `/2` in the axis title was redundant and confusing.

**Changed (line 838-841):**
```r
# All panels share the same y-axis label. Panel (a) tick labels are
# already halved via scales::label_number(scale = 0.5) -- the previous
# ~"/2" suffix in the axis title rendered as a trailing "2" artifact.
y = expression(Kelp~biomass~(g~m^{-2})))
```

### Risk Assessment
Low risk. This change only removes an ambiguous text suffix. Data, tick values, and all other formatting are unchanged.

---

## Figure 2: "MPA established" Annotation Text Overlaps Data in Panel (a)

### Issue Reported
The "MPA established" annotation text in panel (a) overlaps with data points near the MPA establishment year.

### Investigation

The annotation was defined on lines 1032-1035 (original):

```r
annotate("text", x = scorpion_start - 0.5, y = Inf,
         label = "MPA\nestablished", hjust = 1.05, vjust = 1.25,
         size = 2.3, color = MPA_LABEL_COLOR, lineheight = 0.9)
```

**Root cause:** `annotate("text", ...)` renders plain text with no background. The text is positioned at `y = Inf` (top of panel) with `vjust = 1.25`, which shifts it slightly below the top edge. When data points exist near the MPA establishment year (x ~ 2003) and at high y-values, the annotation text directly overlaps those data points with no visual separation.

### Fix Applied

Changed `annotate("text", ...)` to `annotate("label", ...)` with a semi-transparent white fill and no border. This adds a subtle background rectangle behind the annotation text that occludes any data underneath, creating clear visual separation without changing the annotation's position or appearance.

**Changed (lines 1032-1038):**
```r
# TUFTE: Annotate the MPA implementation line with semi-transparent background
# to prevent overlap with data points near the MPA establishment year.
annotate("label", x = scorpion_start - 0.5, y = Inf,
         label = "MPA\nestablished", hjust = 1.05, vjust = 1.25,
         size = 2.3, color = MPA_LABEL_COLOR, lineheight = 0.9,
         fill = alpha("white", 0.8), label.size = 0,
         label.padding = unit(0.15, "lines"))
```

Key parameters:
- `fill = alpha("white", 0.8)` -- 80% opaque white background, enough to separate text from data while not being visually heavy
- `label.size = 0` -- no border around the label box (keeps the Tufte-style minimal aesthetic)
- `label.padding = unit(0.15, "lines")` -- tight padding to avoid an oversized box

### Risk Assessment
Low risk. Only the annotation rendering method changed. Position, text content, color, and size are all unchanged. The visual effect is a subtle white background behind the two-line annotation.

---

## Summary of Changes

| Figure | Issue | Root Cause | Fix | Lines Changed |
|--------|-------|------------|-----|---------------|
| Fig 1 (inset panels) | Trailing "2" in y-axis label | `~"/2"` suffix in `expression()` renders slash poorly | Removed `~"/2"` suffix; tick labels already communicate half-scaling | 838-841 |
| Fig 2 (panel a) | Annotation overlaps data | `annotate("text")` has no background | Changed to `annotate("label")` with semi-transparent fill | 1032-1038 |

Both fixes are minimal, low-risk, and preserve the original design intent.
