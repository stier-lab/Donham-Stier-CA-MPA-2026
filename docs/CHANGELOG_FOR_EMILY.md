# Pipeline Changes for Manuscript Revision

**Prepared by:** Adrian Stier
**Last Updated:** 2026-02-13
**Purpose:** Summary of all changes between the original analysis (pBACIPS_PISCO_V10.R) and the current modular pipeline

---

## Action Items for Emily

| Priority | Action | Time |
|----------|--------|------|
| **1** | Review new Table 2 (below) — 5 of 9 effects are now significant | 15 min |
| **2** | Update Abstract: kelp is now **+257%** (was ~266% in V5 draft, then ~80% in joint model) | 5 min |
| **3** | Update Results section with new Table 2 values and significance statements | 30 min |
| **4** | Review new Figure 0 (trophic cascade concept diagram) — include in main text? | 5 min |
| **5** | Decide on Figure 5 placement: Results or Discussion? | 5 min |
| **6** | Review Methods text for per-taxa meta-analysis description (see below) | 10 min |

---

## The Big Change: Per-Taxa Meta-Analysis

The most important methodological change is switching from a **joint multi-taxa model** to **per-taxa meta-analysis**.

### What changed

- **Old (joint model):** One `rma.mv(yi ~ Taxa - 1, ...)` model fit to all species simultaneously. Cook's D outliers detected at a global 4/n threshold (n = all observations across all taxa).
- **New (per-taxa):** Separate `rma.mv()` models fit for each taxon-response combination. Cook's D outliers detected within each taxon at 4/k threshold (k = number of effect sizes for that taxon).

### Why this matters

The joint model flagged **62% of observations** as outliers because effect sizes from one taxon (e.g., kelp with large positive effects) looked like "outliers" relative to other taxa (e.g., urchins with negative effects). This is statistically indefensible — the whole point of fitting taxa as moderators is that they differ.

The per-taxa approach asks the correct question: "Is this kelp effect size unusual *relative to other kelp effect sizes*?" This retains far more data and produces more powerful, defensible results.

### Results: 5 of 9 effects now significant (was 2 of 9)

The three effects that gained significance were always present in the data — they were masked by inappropriate outlier removal in the joint model.

---

## Current Table 2: Meta-Analysis Results

These are the values in `data/table_02_meta_analysis.csv`. All p-values are FDR-corrected.

### Biomass

| Taxa | Common Name | k | Estimate | SE | RR | % Change | p (FDR) | Sig? |
|------|-------------|---|----------|----|----|----------|---------|------|
| S. purpuratus | Purple urchin | 10 | -0.876 | 0.130 | 0.42 | **-58%** | <0.001 | YES |
| M. franciscanus | Red urchin | 12 | +0.243 | 0.284 | 1.27 | +27% | 0.411 | no |
| M. pyrifera | Giant kelp | 29 | **+1.272** | 0.259 | **3.57** | **+257%** | **<0.001** | **YES** |
| P. interruptus | Spiny lobster | 6 | +0.940 | 0.378 | 2.56 | +156% | 0.083 | no |
| S. pulcher | Sheephead | 17 | **+0.829** | 0.154 | **2.29** | **+129%** | **<0.001** | **YES** |

### Density

| Taxa | Common Name | k | Estimate | SE | RR | % Change | p (FDR) | Sig? |
|------|-------------|---|----------|----|----|----------|---------|------|
| S. purpuratus | Purple urchin | 12 | **-1.450** | 0.281 | **0.23** | **-77%** | **<0.001** | **YES** |
| M. franciscanus | Red urchin | 12 | -0.428 | 0.378 | 0.65 | -35% | 0.317 | no |
| P. interruptus | Spiny lobster | 15 | **+1.307** | 0.337 | **3.70** | **+270%** | **0.003** | **YES** |
| S. pulcher | Sheephead | 18 | +0.235 | 0.123 | 1.27 | +27% | 0.095 | no |

### Key narrative points

1. **Kelp biomass is 3.6x higher inside MPAs (+257%)** — highly significant, robust to leave-one-out (29/29 permutations significant)
2. **Lobster density is 3.7x higher (+270%)** — now significant with per-taxa analysis
3. **Purple urchin biomass is 58% lower, density 77% lower** — both highly significant
4. **Sheephead biomass is 2.3x higher (+129%)** — significant
5. **Red urchin effects are non-significant** in both biomass and density
6. **Trophic cascade evidence (Fig 4d):** Urchin density predicts kelp biomass (beta = -0.51, p = 0.021, R-squared = 84%)

---

## New Supplemental Tables

| File | Description |
|------|-------------|
| `data/table_s_outlier_sensitivity.csv` | Compares 3 outlier approaches: no removal, per-taxa Cook's D (primary), joint Cook's D |
| `data/table_s_kelp_leave1out.csv` | Leave-one-out analysis for kelp biomass — significant in all 29/29 permutations |
| `data/table_s_ar1_sensitivity.csv` | AR1 temporal autocorrelation sensitivity check |
| `data/table_s_cascade_consistency.csv` | Per-MPA cascade consistency scores |
| `data/table_s_moderator_meta_regression.csv` | Protection level (SMR vs SMCA) and location moderators |
| `data/table_s_temporal_meta_regression.csv` | Species-level temporal slopes |

---

## Figure Inventory

### Main Text (6 figures)

| Figure | File | Description |
|--------|------|-------------|
| **Fig 0** | `fig_00_trophic_concept` | NEW: Conceptual diagram of trophic cascade hypothesis (MPA → Predators → Urchins → Kelp) |
| **Fig 1** | `fig_01_mpa_map` | MPA map with bathymetry + 4 kelp time series insets. Site labels now have white halos for readability; inset headers show MPA type (SMCA/SMR) |
| **Fig 2** | `fig_02_data_processing` | Pipeline illustration (raw → standardized → lnRR). Now labeled "Example: M. pyrifera at Scorpion SMR" |
| **Fig 3** | `fig_03_mean_effects` | Meta-analytic effect sizes by taxa. Now shows common names + significance stars (FDR-corrected) |
| **Fig 4** | `fig_04_trophic_scatter` | 4-panel trophic cascade scatter. Solid lines = significant, dashed = non-significant. Panel (d) is the key finding (R-squared = 84%) |
| **Fig 5** | `fig_05_recovery_curves` | Recovery trajectories over time (5 species). Spaghetti lines subdued; kelp panel y-axis tightened |

### Supplemental (9+ figures)

| Figure | File | Description |
|--------|------|-------------|
| **S1** | `fig_s01_forest_plot` | Forest plot: effect sizes by MPA for each taxa |
| **S2** | `fig_s02_all_taxa_timeseries` | Refactored to 5x3 facet grid (species rows x site columns) — no more overplotting |
| **S3** | `fig_s03_recovery_curves` | GAM recovery curves. Now has legend identifying "GAM smooth (95% CI)" vs "Individual MPA" |
| **S4** | `fig_s04_cascade_phase` | Species-pair phase portraits |
| **S5** | `fig_s05_triptych_heatmap` | Space-time heatmaps |
| **S6** | `fig_s06_slope_comparison` | Per-MPA slope comparison |
| **S7** | `fig_s07_statistical_transparency` | Model selection + variance components. Whitespace reduced; value labels added to tau-squared bars |
| **S8** | `fig_s08_appendix_*` | Site-level lnRR time series (5 files, one per taxa) |
| **S9** | `fig_s09_moderator_comparisons` | SMR vs SMCA and mainland vs Channel Islands comparisons |

---

## Other Methodological Changes

These were implemented before the per-taxa switch and remain in the current pipeline:

### Random effects structure
- **Old:** `random = ~1|MPA`
- **New:** `random = list(~1|MPA, ~1|Source)` (crossed random effects for MPA and data source)
- 52% of MPAs are sampled by multiple sources, justifying the crossed structure

### Zero-correction method
- **Old:** Fixed `+0.01` added to all zero proportions
- **New:** Adaptive constant (half the minimum non-zero proportion)

### Effect size SE for NLS models
- Now uses delta method with full variance-covariance matrix instead of assuming independence

### Fallback model exclusion
- Non-NLS fallback models are tagged and excluded from AICc competition

### Back-transformed response ratios
- Figures 3, S1, S9 now display on the response ratio (RR) scale
- Table 2 CSV includes RR, RR_lower, RR_upper, Pct_Change columns

---

## Suggested Methods Text (Per-Taxa Meta-Analysis)

> "We estimated the overall effect of MPA protection on each taxon-response combination using separate multilevel meta-analytic models (metafor::rma.mv, Viechtbauer 2010). For taxa with k >= 5 effect sizes from >= 3 MPAs, models included MPA as a random effect to account for spatial non-independence. For taxa with fewer effect sizes, simple random-effects models were used. Within each taxon, influential observations were identified using Cook's distance at a 4/k threshold and removed prior to final model fitting. P-values were corrected for multiple comparisons using the false discovery rate (FDR) method across all 9 taxon-response tests."

---

## Code Structure

The original monolithic script has been refactored into numbered modules:

```
00_libraries.R           Package dependencies
00b_color_palette.R      Colors, shapes, theme_mpa()
00c_analysis_constants.R Named constants, exclusion lists
01_utils.R               Utility functions, save_fig()
02_pBACIPS_function.R    Core pBACIPS methodology
03_data_import.R         Import size frequency data
04_pisco_processing.R    PISCO data processing
05_kfm_processing.R      KFM/NPS data processing
06_ltr_processing.R      LTER data processing
06b_landsat_processing.R Satellite kelp canopy data
07_combine_data.R        Combine all sources
08_effect_sizes.R        Calculate effect sizes (pBACIPS)
09_meta_analysis.R       Per-taxa meta-analysis + sensitivity tables
10_temporal_analysis.R   Temporal dynamics (Figs S3-S6)
11_figures.R             Main text + supplemental figures
12_results_summary.R     Results CSVs and markdown summary
13_additional_analyses.R Moderator analyses (Fig S9)
run_all.R                Pipeline orchestration (~2.2 min)
```

Run the full pipeline with: `source(here::here("code", "R", "run_all.R"))`

---

## Audit Trail

Detailed filtering records are in `outputs/`:
- `filter_audit_effect_sizes.csv` — Every observation's inclusion/exclusion reason at the effect size stage
- `filter_audit_meta_analysis.csv` — Joint model outlier detection (sensitivity comparison)
- `filter_audit_pertaxa_meta.csv` — Per-taxa outlier detection (primary analysis)
- `filter_summary_by_taxa.csv` — Summary counts by taxa
- `data_flow_summary.csv` — k-values through each pipeline stage
