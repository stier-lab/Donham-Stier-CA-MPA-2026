# CA MPA Kelp Forest Analysis - Project Guidelines

## Project Overview

This project analyzes the effects of California Marine Protected Areas (MPAs) on kelp forest ecosystems using the progressive-change Before-After-Control-Impact-Pairs (pBACIPS) methodology. The analysis combines data from multiple monitoring programs (PISCO, KFM/MBON, LTER) and produces publication-quality figures and statistical summaries.

**Target Journal:** Conservation Letters
**Manuscript Title:** Restoration of trophic cascades within kelp forests following the establishment of a network of marine protected areas
**Authors:** Emily Donham (lead) & Adrian Stier
**Status:** Analysis pipeline complete; manuscript revision in progress (Emily leading rewrite)

## Conservation Letters Figure & Design Specifications

Sourced from [Wiley Electronic Artwork Guidelines](https://authorservices.wiley.com/asset/photos/electronic_artwork_guidelines.pdf) and [Conservation Letters Author Guidelines](https://conbio.onlinelibrary.wiley.com/hub/journal/1755263x/homepage/forauthors.html).

### Dimensions
- **Single-column (small):** 80 mm width (quarter-page figures)
- **Double-column (large):** 180 mm width (half-page or full-page figures)
- **Minimum pixel width:** 1800 px at any size
- Figures must be created between 80-180 mm width; figures outside this range will be modified during composition (may decrease quality)

### Resolution
- **Line art** (graphs, scatter plots, flowcharts): **600-1000 DPI** preferred (must be legible at 80 mm / 1800 px unmagnified)
- **Images** (photographs, maps): **300 DPI** minimum
- If a figure contains both line art and images, follow the line art (higher) resolution guidelines

### File Formats
| Stage | Line Art | Images |
|-------|----------|--------|
| Peer review | EPS, PDF (preferred); any standard type acceptable | TIFF, PNG, EPS (preferred) |
| Post-acceptance | EPS, PDF (preferred) | TIFF, PNG, EPS (preferred) |

- When in doubt, submit a PDF
- 1 figure per file preferred; name files `Figure_1.tiff`, `Figure_2.pdf`, etc.

### Typography & Line Weights
- **Minimum font size:** 8 pt at final printed size — all words and symbols must be large enough for easy reading
- **Font:** Sans-serif preferred (Helvetica, Arial); use consistent font throughout all figures
- **Line thickness:** No less than 0.5 pt (0.2 mm); evenly balanced stroke weights
- Larger fonts improve readability, especially for line art at 80 mm width

### File Size
- **Individual figures:** < 10 MB each
- **Complete submission (zipped):** < 500 MB total
- Remove excess white space around figures to reduce file size

### Legends & Captions
- Use Arabic numerals in order of appearance in text
- Place figure legends in a **separate section in the manuscript, after references** (preferred) — acceptable anywhere that clearly indicates which figure it explains
- Legends should succinctly describe content and explain all abbreviations/symbols
- All figures must be cited in main text in the order they appear

### Conservation Letters-Specific Requirements
- **Error bars required:** All figures showing statistical data must display error bars (95% CIs preferred); authors must state what error bars represent
- **95% CIs in text and figures:** Any article reporting p-values must also report 95% CIs
- **Max figures+tables:** ~8 combined in main text (additional in Supporting Information)
- **Supporting Information:** Submit in separate files
- **Cover/social media image:** Provide an image for social media promotion; if selected for cover, 300 DPI at reproduction size in CMYK required

### Submission Checklist (from Wiley)
- [ ] All figures included as separate files or in a single PDF/Word document
- [ ] Each figure has an accompanying legend explaining content, abbreviations, and symbols
- [ ] All figures cited in main text in numbered order
- [ ] All text/symbols large enough for easy reading (min 8 pt at final size)
- [ ] Figures saved in preferred file types (EPS/PDF for line art, TIFF/PNG for images)
- [ ] Each figure file < 10 MB
- [ ] Figures created between 80-180 mm width at 300-600 DPI
- [ ] Figure files named with figure number only (e.g., `Figure_1.tiff`)

## Directory Structure

```
Donham-Stier-CA-MPA-2026/
├── code/
│   └── R/                    # R scripts (numbered 00-13)
├── data/
│   ├── harmonized/           # Analysis-ready CSVs from data-processing repo (~1 MB, tracked in git)
│   ├── cache/                # Bootstrap and intermediate results (.rds)
│   ├── MPA/                  # MPA boundary shapefiles (for Fig 1 map only)
│   └── (PISCO/, MBON/, LTER/, LANDSAT/ no longer needed — raw data lives in data repo)
├── docs/                     # Documentation and manuscript
├── tables/                   # Generated manuscript tables (CSV)
├── plots/                    # Generated figures (PDF + PNG at 600 DPI)
├── outputs/                  # Filter audits, data flow summaries, replicate effects
├── logs/                     # Pipeline execution logs
├── renv.lock                 # R package version lockfile (for reproducibility)
└── MPA_Kelp.Rproj            # RStudio project file
```

### Sibling Data-Processing Repo

Raw data processing (scripts 03-07) lives in a separate repo:

```
~/Donham-Stier-CA-MPA-Data-2026/    (data processing)
  data/ (raw monitoring files, 1.3 GB)
    → scripts 03-07
    → output/harmonized/ (4 CSVs, ~1 MB)

~/Donham-Stier-CA-MPA-2026/          (this repo — analysis)
  data/harmonized/ (4 CSVs, tracked in git)
    → scripts 08-13
    → figures, tables, meta-analysis
```

## Manuscript Figure Mapping

### Main Text Figures

| MS Figure | Description | Code Output File | Script |
|-----------|-------------|------------------|--------|
| Figure 1 | Map of MPAs with Channel Islands | `fig_01_mpa_map.pdf` | 11_figures.R |
| Figure 2 | Trophic cascade case studies: 3×3 grid (predators/urchins/kelp rows × 3 Channel Islands sites), before/after MPA with linear trends | `fig_02_cascade_case_studies.pdf` | 11_figures.R |
| Figure 3 | Meta-analytic mean effect sizes by taxa (RR-scaled axis) | `fig_03_mean_effects.pdf` | 11_figures.R |
| Figure 4 | Recovery trajectories: 3×2 trophic grid (predators/herbivores/producer rows, 2 species per row), biomass lmer predictions with 95% CI | `fig_04_recovery_curves.pdf` | 11_figures.R |
| Figure 5 (candidate, new) | Heatwave × MPA cascade: per-taxon lnRR (emmeans ± 95% CI) across before/during/after the 2014-16 MHW, RR-scaled axis, Southern California. Companion to Kumagai et al. 2024. | `fig_heatwave_cascade.pdf` | 14_heatwave_analysis.R |
| Figure S (candidate, new) | Cascade regressions across MPA-years (kelp~urchin, urchin~lobster), density lnRR; cf. Kumagai 2024 Fig 7 | `fig_heatwave_cascade_regression.pdf` | 14_heatwave_analysis.R |
| Figure S (candidate, new) | Analytical multiverse: method bridge (one analytical flip per step, Kumagai-style → ours) + cross-substrate data effect; after-before lnRR change per taxon, RR-scaled | `fig_s_methods_multiverse.pdf` | 15_methods_comparison.R |
| Figure S (candidate, new) | MPA effect vs environmental gradients: per-SD meta-regression slopes (forest, by taxon) for MHW intensity, cold-spell intensity, latitude, log MPA size (Kumagai PCA covariates, here modeled) | `fig_s_env_moderators.pdf` | 16_environmental_moderators.R |
| Figure S (candidate, new) | Reproduced Eisaguirre 2020 Fig 3 (urchin~sunflower star, pre-disease) and Fig 4 (urchin~sheephead by protection, post-disease) | `fig_s_eisaguirre_fig3_urchin_seastar.pdf`, `fig_s_eisaguirre_fig4_urchin_sheephead.pdf` | 17_eisaguirre_reproduction.R |
| Figure S (candidate, new) | Predicting MPA effectiveness: random-forest variable importance + observed-vs-LOO-predicted kelp lnRR (LOO-CV R² negative) | `fig_s_effectiveness_predictors.pdf` | 18_mpa_effectiveness_predictors.R |

Environmental-moderator supplement (16_environmental_moderators.R) tests whether the MPA effect varies across the environmental gradients Kumagai et al. (2024) MAPPED in their PCA (Fig 3) but never put in a model — closing the one covariate gap surfaced by script 15. Per-MPA covariates obtainable without fabrication: MHW cumulative intensity during 2014-16 (thermal stress; their 1-km raster at MPA centroids), cold-spell cumulative intensity (thermal regime/upwelling; their 1-km grid, nearest cell), latitude, log MPA size. Per taxon, random-effects meta-regression (`metafor::rma`, **Knapp-Hartung** small-sample t-test — essential at k=9-15; density for animals, biomass for kelp), BH-FDR across 20 tests. NOT available per-MPA (documented, not fabricated): wave exposure (`hsmax`; their per-site NetCDF not in mirror — only a single regional SBC wave series exists), nitrate, depth, human gravity. Tables: `table_s_mpa_env_covariates.csv`, `table_s_env_moderators.csv`; section + legend in `docs/methods_comparison_supplement.md`.
KEY: **clean null — no moderator survives FDR.** Strongest = red urchin effect declines with MHW intensity (−1.84 lnRR/SD, p=0.003, FDR=0.052; suggestive only). **Giant kelp effect is NOT modulated by any gradient** (all p≥0.17) → resilience is not an artifact of thermal/biogeographic reserve placement. Practical upshot: omitting Kumagai's environmental covariates does not appear to bias the headline conclusions. (Without Knapp-Hartung the same k=9-15 regressions give implausible p down to 1e-18 — overconfident IV weighting.)

Methods-comparison supplement (15_methods_comparison.R) is a formal **analytical multiverse** quantifying how the method of analysis (vs the data) drives variation in MPA × heatwave effect sizes, contrasting our paired/proportion/lnRR + meta-analysis against Kumagai's pooled/raw-density GLMM. It walks a 6-waypoint **method bridge** on a common dataset, flipping ONE choice per step (scale → distribution/link → reference/pairing → AR1 autocorrelation → IV weighting), and separately runs the pooled waypoints on BOTH datasets to isolate the **data effect**. Common currency = after-before change in ln(MPA/Reference) per taxon. Tables: `table_s_methods_crosswalk.csv` (qualitative axis-by-axis), `table_s_methods_multiverse.csv` (44 cells), `table_s_methods_decomposition.csv` (per-axis mean/max |Δln| + significance flips). Draft: `docs/methods_comparison_supplement.md`.
KEY: **giant kelp resilience is method-invariant** (significant at every waypoint, both datasets); **lobster** needs pairing to appear (p=0.86 pooled-raw → p≤0.008 paired); **urchin inference is set by the AR1 choice** — point estimate barely moves (mean |Δln|=0.14) but significance flips (purple after WP3 p=1e-6 → WP4/AR1 p=0.12 → WP5/naive-meta p=1e-22), reproducing script 14's caveat; and the **single largest lever is the dataset itself** (mean |Δln|=0.63, max 1.83, all 5 taxa flip significance) — data > any one method axis.

Heatwave models (14_heatwave_analysis.R) use glmmTMB with site + Source random intercepts and an **AR1(year|MPA)** residual structure as the PRIMARY spec — annual lnRR is strongly temporally autocorrelated (lag-1 ACF 0.2-0.6; `table_heatwave_diagnostics.csv`), and naive lme4 (no AR1) inflates p by many orders of magnitude (shown as `p_naive` in `table_heatwave_contrasts.csv`). Tables: `table_heatwave_period_effects.csv` (RR by period), `table_heatwave_contrasts.csv` (AR1 + naive), `table_heatwave_diagnostics.csv` (DHARMa + ACF), `table_heatwave_cascade_regression.csv`, `table_heatwave_per_mpa_exposure.csv`, `table_heatwave_pbacips_integrated.csv`, `table_heatwave_sensitivity.csv`, `table_heatwave_crosscheck_kumagai.csv`.
KEY (rigorous, AR1): only **giant kelp** is robust — far more abundant inside MPAs during/after the MHW (p<1e-4) with a heatwave-specific effect beyond recovery (p=0.03). Lobster higher inside after (p=0.005) but recovery-driven. Urchin period-suppression and the lobster->urchin link are NOT robust under AR1 (purple after p=0.19; urchin~lobster p=0.91) — the cascade-mediation claim is only weakly supported once autocorrelation is handled. Draft: `docs/heatwave_section_draft.md`; legends: `docs/heatwave_figure_legends.md`.

**Figure 1 Details:**
- Base map: Southern California coastline with Channel Islands
- Site markers: Shape indicates data source (square=NPS-KFM, circle=LTER, triangle=PISCO)

**Figure 2 Details:**
- 3 rows (Predators, Urchins, Producers) × 3 columns (Scorpion SMR, Gull Island SMR, SB Island SMR)
- Before-period data shown with lighter alpha; after-period with full weight + linear trend
- Vertical dotted line at time=0 (MPA implementation)
- Y-axis shared per trophic row, RR-scaled

### Supplemental Figures

**Note:** The SI document renumbered figures sequentially after dropping S3, S7, and S10. The disk filenames retain their original numbering. The table below maps SI figure numbers to their disk files.

| SI Figure | Description | Disk File | Script |
|-----------|-------------|-----------|--------|
| Figure S1 | Data processing pipeline (raw → proportion → lnRR) | `fig_s01_data_processing.pdf` | 11_figures.R |
| Figure S2 | Forest plot: effect sizes by MPA for each taxa (RR-scaled axis) | `fig_s02_forest_plot.pdf` | 11_figures.R |
| Figure S3 | GAM recovery curves with MPA spaghetti + lmer overlay (t≤15) | `fig_s04_recovery_curves.pdf` | 10_temporal_analysis.R |
| Figure S4 | Species-pair phase portraits (network-wide mean lnRR) | `fig_s05_cascade_phase.pdf` | 10_temporal_analysis.R |
| Figure S5 | Per-MPA slope comparison and cascade consistency | `fig_s06_slope_comparison.pdf` | 10_temporal_analysis.R |
| Figure S6 | Combined moderator comparisons: SMR vs SMCA, Islands vs mainland | `fig_s09_moderator_comparisons.pdf` | 13_additional_analyses.R |
| Figure S7a-e | Site-level appendix: lnRR time series per taxa (5 files) | `fig_s08_appendix_*.pdf` | 11_figures.R |
| Figure S8 | DHARMa residual diagnostics for NLS effect-size models | `fig_s11_dharma_diagnostics.pdf` | 11_figures.R |
| Figure S9 | Meta-analysis funnel plots + Egger's test (publication bias) | `fig_s12_funnel_plots.pdf` | 11_figures.R |
| Figure S10 | lmer temporal model residual diagnostics (4-panel) | `fig_s13_lmer_residuals.pdf` | 11_figures.R |
| Figure S11 | NLS model selection frequencies & DW autocorrelation | `fig_s14_model_selection.pdf` | 11_figures.R |
| Figure S12 | Cook's distance & outlier sensitivity summary | `fig_s15_sensitivity_summary.pdf` | 11_figures.R |

Dropped from SI (not included in manuscript; code still renderable via `should_render()`):
- `fig_s03_all_taxa_timeseries` (was old S3)
- `fig_s07_statistical_transparency` (was old S7)
- `fig_s10_recovery_bio_den` (was old S10)

### Tables

| MS Table | Description | Code Output File | Script |
|----------|-------------|------------------|--------|
| Table 1 | Average density/biomass by taxa and source | `tables/average_responses.csv` | data-processing repo (07) |
| Table 2 | Meta-analysis summary statistics (with heterogeneity) | `tables/table_02_meta_analysis.csv` | 09_meta_analysis.R |
| Table 3 | Cross-taxa meta-regressions (trophic cascade tests) | `tables/table_03_cross_taxa_meta_regression.csv` | 09_meta_analysis.R |
| Table S1 | Per-MPA data availability + taxa matrix (merged) | `tables/table_s_mpa_data_availability.csv` + `tables/table_s_mpa_taxa_matrix.csv` | 12_results_summary.R |
| Table S1b | Sample sizes by taxa × response × source | `tables/table_s_sample_sizes.csv` | 12_results_summary.R |
| Table S2 | Variance components with CIs | `tables/table_s_variance_components.csv` | 09_meta_analysis.R |
| Table S3a | Source RE sensitivity: model comparison | `tables/table_s_source_sensitivity_models.csv` | 09_meta_analysis.R |
| Table S3b | Source RE sensitivity: coefficient comparison | `tables/table_s_source_sensitivity_coefficients.csv` | 09_meta_analysis.R |
| Table S4 | Temporal meta-regression (stacked: pooled + bio + den) | `tables/table_s_temporal_meta_regression.csv` + `tables/table_s_temporal_meta_regression_bio.csv` + `tables/table_s_temporal_meta_regression_den.csv` | 10_temporal_analysis.R |
| Table S5 | Cascade consistency scores by MPA | `tables/table_s_cascade_consistency.csv` | 10_temporal_analysis.R |
| Table S6 | Moderator meta-regression + diagnostics (merged) | `tables/table_s_moderator_meta_regression.csv` + `tables/table_s_moderator_diagnostics.csv` | 13_additional_analyses.R |
| Table S7 | Temporal diagnostics: GAM linearity + AR1 sensitivity (merged) | `tables/table_s_gam_linearity_diagnostics.csv` + `tables/table_s_ar1_sensitivity.csv` | 10_temporal_analysis.R |
| Table S8 | Kelp leave-one-out sensitivity | `tables/table_s_kelp_leave1out.csv` | 09_meta_analysis.R |
| Table S9 | Outlier sensitivity (4 methods) | `tables/table_s_outlier_sensitivity.csv` | 09_meta_analysis.R |

Tables dropped from SI (data still in CSVs but not shown as numbered SI tables):
- Old Table S1c (`table_s_analysis_summary.csv`) — absorbed into S1 overview text
- Old moderator subgroup table (`table_s_moderator_subgroup_effects.csv`) — dropped

Additional outputs (not numbered as manuscript tables):

| File | Description | Script |
|------|-------------|--------|
| `tables/table_s_lmer_model_fit.csv` | lmer model fit statistics (AIC, BIC, R-squared, variance) | 10_temporal_analysis.R |
| `tables/table_s_species_slopes.csv` | Species-level derived slopes from lmer (pooled, bio, den) | 10_temporal_analysis.R |
| `tables/table_s_slope_summary_by_species.csv` | Slope summary by species (direction, significance) | 10_temporal_analysis.R |
| `tables/table_s_per_mpa_slopes.csv` | Per-MPA slope estimates by species | 10_temporal_analysis.R |
| `tables/table_s_phase_portrait_correlations.csv` | Species-pair phase portrait correlations | 10_temporal_analysis.R |
| `tables/table_s_phase_portrait_yearly_means.csv` | Yearly mean lnRR by species for phase portraits | 10_temporal_analysis.R |
| `tables/table_data_provenance_rr.csv` | Response ratio data provenance (MPAs, years per source) | data-processing repo (07) |
| `tables/table_data_provenance_raw.csv` | Raw data provenance (observations per source/taxon) | data-processing repo (07) |
| `tables/table_sample_sizes_per_taxa.csv` | Sample size breakdown per taxa (input counts) | 08_effect_sizes.R |
| `tables/table_model_selection.csv` | NLS model type selection per MPA/taxa | 08_effect_sizes.R |
| `data/sumstats_final.csv` | Full SumStats.Final export (all 146 effect sizes) | 08_effect_sizes.R |
| `outputs/table_cascade_analysis.csv` | Cascade meta-regression results | 11_figures.R |
| `outputs/model_results_summary.csv` | Meta-analysis results in machine-readable format | 12_results_summary.R |
| `outputs/replicate_effects.csv` | All MPA-taxa-response replicate effect sizes | 12_results_summary.R |
| `outputs/data_filtering_steps.csv` | Row counts at each pipeline filtering step | data-processing repo (07) |
| `outputs/data_flow_summary.csv` | k-values through each pipeline stage | 12_results_summary.R |
| `outputs/filter_summary_by_taxa.csv` | Filtering summary counts by taxa | 08_effect_sizes.R |
| `tables/table_s_lmer_prediction_params.csv` | lmer prediction parameters for recovery curves | 10_temporal_analysis.R |
| `outputs/filter_audit_effect_sizes.csv` | Per-observation inclusion/exclusion audit | 08_effect_sizes.R |
| `outputs/filter_audit_meta_analysis.csv` | Joint model outlier detection audit | 09_meta_analysis.R |
| `outputs/filter_audit_pertaxa_meta.csv` | Per-taxa outlier detection audit | 09_meta_analysis.R |
| `outputs/lmer_residual_diagnostics.csv` | lmer residuals/fitted values for diagnostic plots | 10_temporal_analysis.R |
| `outputs/lmer_ranef_diagnostics.csv` | lmer random effect BLUPs for diagnostic plots | 10_temporal_analysis.R |
| `outputs/model_fallback_audit.csv` | NLS fallback model audit (which MPAs used fallback) | 08_effect_sizes.R |
| `data/model_diagnostics.csv` | DHARMa + NLS residual diagnostics per model | 08_effect_sizes.R |

## File Naming Conventions

### General Principles
- **Use lowercase snake_case** for all generated files
- **Use descriptive names** that indicate content
- **Include version/figure numbers** with zero-padding (e.g., `01`, `02`)
- **Avoid spaces and special characters** in filenames

### Figures (`plots/`)

**Main text figures:** `fig_{NN}_{descriptive_name}.{pdf|png}`
**Supplemental figures:** `fig_s{NN}_{descriptive_name}.{pdf|png}`

Examples:
- `fig_01_mpa_map.pdf` (Main text Figure 1)
- `fig_s01_data_processing.pdf` (Supplemental Figure S1)

Guidelines:
- Always export both PDF (vector) and PNG (600 DPI raster)
- Use zero-padded two-digit figure numbers
- Main text: `fig_01`, `fig_02`, etc.
- Supplemental: `fig_s01`, `fig_s02`, etc.

### Table Outputs (`tables/`)

Pattern: `{descriptive_name}.csv` or `table_{NN}_{description}.csv`

Examples:
- `average_responses.csv`
- `table_02_meta_analysis.csv`

Guidelines:
- Use snake_case for all output files
- Tables referenced in manuscript use `table_{NN}_{description}.csv`
- Summary statistics use descriptive names without numbers

### Cache Files (`data/cache/`)

Pattern: `{source}_{taxon}_{operation}.rds`

Examples:
- `pisco_urchin_bootstrap.rds`
- `kfm_urchin_bootstrap.rds`
- `lter_urchin_bootstrap.rds`
- `vrg_panint_bootstrap.rds`

### R Scripts (`code/R/`)

Pattern: `{NN}_{descriptive_name}.R` or `{NN}{letter}_{descriptive_name}.R`

Pipeline order:
```
00_download_data.R       - Download data from Dryad (placeholder — needs DOI after upload)
00_libraries.R           - Package dependencies
00b_color_palette.R      - Color scheme and ggplot theme
00c_analysis_constants.R - Named constants (survey areas, size thresholds, exclusions)
01_utils.R               - Utility functions (effect size, figure rendering)
02_pBACIPS_function.R    - Core statistical function
03_load_harmonized_data.R - Load harmonized CSVs (replaces old scripts 03-07)
08_effect_sizes.R        - Calculate effect sizes
09_meta_analysis.R       - Multilevel meta-analysis (joint model primary, per-taxa sensitivity)
10_temporal_analysis.R   - Temporal dynamics (SI Figs S3-S5, temporal tables)
11_figures.R             - Publication figures (Figs 1-4, S1-S2, S7a-e, S8-S12 diagnostics)
13_additional_analyses.R - Moderator analyses (SI Fig S6, moderator tables)
14_heatwave_analysis.R   - MPA cascade x 2014-16 marine heatwave (Southern CA; Kumagai 2024 comparison)
15_methods_comparison.R  - Analytical multiverse: how method of analysis drives effect-size variation vs Kumagai (SI supplement)
16_environmental_moderators.R - Env covariates (Kumagai PCA: MHW/cold-spell/latitude/MPA size) as meta-regression moderators (SI)
17_eisaguirre_reproduction.R - Reproduce Eisaguirre 2020 (paired NCI cascade) from raw PISCO (SI; standalone)
18_mpa_effectiveness_predictors.R - Predict among-MPA variation in effectiveness (env+design+trophic; LOO-CV) (SI; standalone)
12_results_summary.R     - Generate results CSVs and markdown summary
run_all.R                - Pipeline orchestration
run_figures_only.R       - Fast figure regeneration (~17s all, ~4s single)
```

Note: Data processing scripts (03_data_import, 04_pisco_processing, 05_kfm_processing,
06_lter_processing, 06b_landsat_processing, 07_combine_data) now live in the sibling
data-processing repo (`Donham-Stier-CA-MPA-Data-2026`).

## Code Style Guidelines

### File Paths
- **Always use `here::here()`** for file paths
- Never use absolute paths or `setwd()`
- Example: `here::here("data", "cache", "pisco_urchin_bootstrap.rds")`

### Output Files
- Table outputs go to `tables/`
- Figures go to `plots/`
- Cache/intermediate files go to `data/cache/`
- Raw/input data stays in `data/`
- Never write files to project root

### Statistical Methods
- Document methodological choices with inline comments
- Report uncertainty (SE, CI) for all estimates
- Use `emmeans::pairs()` for contrasts (proper covariance handling)
- Use confidence intervals (not prediction intervals) for effect sizes
- Report heterogeneity statistics (I², τ²) for meta-analyses
- **Primary meta-analysis**: Joint multilevel `rma.mv` with `~ Taxa - 1` moderator and `random = list(~1|MPA, ~1|Source)`. Per-taxa separate models are sensitivity only.
- **P-value correction**: Benjamini-Hochberg FDR applied to all meta-analytic p-values
- NLS effect size SEs use delta method (`nls_difference_se()` in `08_effect_sizes.R`)
- Non-NLS fallback models are excluded from AICc competition
- Temporal meta-regression includes AR1 sensitivity analysis (Section C2 in `10_temporal_analysis.R`)

### R Code Style
- Use tidyverse style (pipes, dplyr verbs)
- Comment non-obvious operations
- Use roxygen-style documentation for functions
- Keep functions in `01_utils.R` or `02_pBACIPS_function.R`

## Color Palette & Figure Conventions

All colors, shapes, and themes are centralized in `code/R/00b_color_palette.R`. Never hardcode hex values in figure scripts — use the named palette vectors and convenience scales defined there.

### Canonical Trophic Order

All figures use a consistent species ordering reflecting the trophic cascade (top predators first, primary producer last):

| Position | Abbreviated | Full Name |
|----------|-------------|-----------|
| 1 (top)  | P. interruptus | Panulirus interruptus |
| 2        | S. pulcher | Semicossyphus pulcher |
| 3        | S. purpuratus | Strongylocentrotus purpuratus |
| 4        | M. franciscanus | Mesocentrotus franciscanus |
| 5 (bottom) | M. pyrifera | Macrocystis pyrifera |

This order is enforced by:
- `taxa_levels` in `11_figures.R` (abbreviated names)
- `species_order_full` in `10_temporal_analysis.R` and `13_additional_analyses.R` (full names)

When using `coord_flip()` with horizontal bar charts, apply `rev(levels)` so the visual top-to-bottom matches trophic order.

### Named Palette Vectors

| Vector | Purpose | Used in |
|--------|---------|---------|
| `col_taxa` | 5 species-specific colors | Figs 2, 3, 4, S2–S5, S7a-e |
| `col_response` | Density (`Den`) vs Biomass (`Bio`) | Fig 3, S2 |
| `col_site` | Inside/Outside MPA status | appendix panels (SI Fig S7) |
| `col_model` | NLS model types (blue/grey family) | diagnostics only (SI Fig S7 dropped) |
| `col_variance` | Random effect variance components | diagnostics only (SI Fig S7 dropped) |
| `col_effect` | Positive/negative direction (Okabe-Ito) | Fig S1 conceptual diagram |
| `col_cascade` | Cascade consistency TRUE/FALSE | SI Fig S5 panel b |
| `col_map` | Land/ocean/coastline | Fig 1 map |
| `col_ba` | Before/After periods | Time series panels |
| `col_type` | Analysis type (pBACIPS/BACI/CI) | Diagnostics |

Key design principle: `col_model`, `col_variance`, `col_effect`, and `col_cascade` are all visually distinct from `col_taxa` to prevent confusion between biological and methodological dimensions.

### Convenience Scale Functions

Drop-in ggplot2 scale layers defined in `00b_color_palette.R`:

```r
scale_color_taxa()       # col_taxa values
scale_fill_taxa()
scale_color_response()   # col_response values
scale_fill_response()
scale_color_site()       # col_site_all values (handles Inside/Outside, MPA/Reference, etc.)
scale_fill_site()
scale_shape_source()     # shape_source + label_source
```

### Response Ratio (RR) Axis Scale

Data is stored on the lnRR (log response ratio) scale internally. For display, most figures show back-transformed RR tick labels using helpers from `01_utils.R`:

```r
scale_y_rr(y_lo, y_hi, name = "Response ratio (MPA / Reference)")
scale_x_rr(x_lo, x_hi, name = "Response ratio (MPA / Reference)")
```

These select tick positions from a pool of human-readable RR values (0.01, 0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8, 20, 100, ...) that fall within the visible data range. RR = 1 (lnRR = 0) means no MPA effect.

**Exceptions** (remain on raw lnRR scale): Fig S1 panel c (conceptual diagram), SI Fig S4 phase portraits (inverted axes).

### Theme

All figures use `theme_mpa()` (defined in `00b_color_palette.R`):
- Based on `theme_bw()` with `base_size = 9`, Helvetica font
- L-shaped axes (bottom + left only), no gridlines
- Clean strip text (italic, no background)
- Tight margins for journal layout

## Data Setup

This analysis repo consumes harmonized CSVs produced by the data-processing repo. The harmonized CSVs (~1 MB) are tracked in git under `data/harmonized/`, so **no additional data setup is needed to run the analysis pipeline**.

### Default: Just run the pipeline

The harmonized CSVs are already in the repo. Run the analysis directly:

```r
source(here::here("code", "R", "run_all.R"))
```

### Option 1: Symlink MPA shapefiles from Google Drive (for Fig 1 map)

Only `data/MPA/` is needed from raw data (for the map figure). Other raw data directories are no longer needed.

```bash
GDRIVE="/Users/$(whoami)/Library/CloudStorage/GoogleDrive-astier@ucsb.edu/My Drive/Stier Lab/People/Emily Donham/Projects/Kelp MPA/data"
PROJECT="path/to/Donham-Stier-CA-MPA-2026/data"

ln -s "$GDRIVE/MPA" "$PROJECT/MPA"
```

### Option 2: Figures only (no raw data needed)

Cached intermediate results are tracked in git (`data/cache/*.rds`, ~8 MB total). To regenerate all figures without running the analysis:

```r
source(here::here("code", "R", "run_figures_only.R"))
```

### Harmonized data files (`data/harmonized/`)

These CSVs are the handoff between the data-processing and analysis repos:

| File | Rows | Description |
|------|------|-------------|
| `harmonized_response_ratios.csv` | ~3,400 | Log response ratios (All.RR.sub.trans) |
| `harmonized_raw_responses.csv` | ~6,200 | Raw density/biomass measurements (All.Resp.sub) |
| `harmonized_landsat_rr.csv` | ~720 | Satellite kelp response ratios (Landsat.RR) |
| `harmonized_site_metadata.csv` | 34 | MPA site metadata (Site) |

See `DATA_DICTIONARY.md` in the data-processing repo for full column documentation.

### Cache files (`data/cache/`)

| File | Contents |
|------|----------|
| `*_bootstrap.rds` | Bootstrap biomass estimates (from data repo) |
| `figures_snapshot.rds` | Snapshot for run_figures_only.R |
| `presplit_reference.rds` | Pre-split reference checksums |
| `*.rds` (map data) | Bathymetry, hillshade, MPA boundaries |

## Running the Analysis

```r
# Full pipeline (harmonized CSVs are tracked in git — no raw data needed)
source(here::here("code", "R", "run_all.R"))

# Figures only (uses cached snapshot, even faster)
source(here::here("code", "R", "run_figures_only.R"))
```

## Key Datasets

| Object | Description | Loaded/Created in |
|--------|-------------|-------------------|
| `Site` | MPA metadata and implementation dates | 03_load_harmonized_data.R |
| `All.RR.sub.trans` | Combined log response ratios | 03_load_harmonized_data.R |
| `All.Resp.sub` | Combined raw response data | 03_load_harmonized_data.R |
| `Landsat.RR` | Satellite kelp response ratios | 03_load_harmonized_data.R |
| `SumStats.Final` | Effect sizes for all MPA-taxa combinations | 08_effect_sizes.R |
| `Table2` | Meta-analysis summary by taxa (from joint model) | 09_meta_analysis.R |
| `meta_biomass_full` | Joint rma.mv model for biomass (primary, no outlier removal) | 09_meta_analysis.R |
| `meta_density_full` | Joint rma.mv model for density (primary, no outlier removal) | 09_meta_analysis.R |

## Adding New Analyses

When adding new figures or outputs:

1. **Main text figures**: `fig_{NN}_{name}.{pdf|png}` - increment from existing
2. **Supplemental figures**: `fig_s{NN}_{name}.{pdf|png}` - use s prefix
3. **Table outputs**: Use snake_case, write to `tables/`
4. **Update this file** with the manuscript mapping

## Documentation

- `docs/MPA_Kelp_MS_V5.pdf` - Archived original manuscript draft (reference for revision)
- `docs/ANALYSIS_REVISIONS.md` - Comprehensive comparison of original archived code vs current pipeline (12 sections)
- `docs/analysis_revisions.html` - Styled HTML rendering of ANALYSIS_REVISIONS.md (with sidebar TOC)
- `docs/BEFORE_AFTER_COMPARISON.md` - Line-by-line manuscript revision guide (old vs new values)
- `docs/email_draft_for_emily.md` - Summary email for Emily with reading order and action items
- `docs/RESULTS_SUMMARY.md` - Auto-generated results summary from 12_results_summary.R
- `docs/results_report.html` - Interactive HTML report (generate by knitting `code/R/results_report.Rmd`)
- `docs/supporting_information.Rmd` - Supporting Information R Markdown source
- `docs/supporting_information.html` - Rendered Supporting Information HTML
- `docs/si_style.css` - Supporting Information custom CSS
- `docs/si_script.html` - Supporting Information JavaScript (lightbox, progress bar)
- `docs/revisions_style.css` - Analysis revisions HTML styling
- `docs/revisions_script.html` - Analysis revisions HTML JavaScript (scroll spy, progress bar)
- This file (`CLAUDE.md`) - Project conventions for AI assistants

## Related Work / External Comparison (in progress)

We are running a **compare-and-contrast** against **Kumagai et al. 2024, *Global Change Biology*** — "Marine Protected Areas That Preserve Trophic Cascades Promote Resilience of Kelp Forests to Marine Heatwaves." It is a parallel study on an overlapping system and dataset (same PISCO/MLPA subtidal surveys + the same T. Bell SBC LTER Landsat product) that reaches a convergent conclusion: Southern California no-take MPAs preserve the lobster/sheephead → urchin → kelp cascade, conferring resilience to the 2014–2016 marine heatwave (MHW), with no effect in Central California (sea otters protected statewide; sunflower stars lost to SSWD).

Goals of the comparison: (1) confirm our PISCO numbers match theirs at shared sites; (2) contrast their approach (GLMM + permutation resistance/recovery framed around the MHW; Hamilton & Caselle 2015 sheephead allometry) against ours (pBACIPS before/after MPA establishment + multilevel meta-analysis); (3) scope **our own MHW analysis** — we currently include no climate/heatwave covariate.

**Formal method-vs-data decomposition (done):** `15_methods_comparison.R` quantifies how much of the difference between the two studies is method vs data via an analytical multiverse — a 6-waypoint method bridge (one analytical flip per step) on a common dataset, plus the pooled waypoints run on both datasets. See the "Methods-comparison supplement" note under the figure-mapping section above and `docs/methods_comparison_supplement.md`. Headline: kelp resilience is method-invariant; urchin inference hinges on the AR1 autocorrelation choice; the dataset itself is the single largest lever.

**Environmental-covariate gap (done):** `16_environmental_moderators.R` tests the environmental predictors Kumagai mapped in their PCA but never modeled (MHW & cold-spell intensity, latitude, MPA size — wave/nitrate/depth unavailable per-MPA without fabrication) as meta-regression moderators of our per-MPA effect sizes (Knapp-Hartung). Clean null: nothing survives FDR; giant-kelp resilience is not modulated by any gradient; omitting these covariates does not appear to bias the headline conclusions. See the "Environmental-moderator supplement" note above.

**Predicting MPA effectiveness (done):** `18_mpa_effectiveness_predictors.R` assembles the fullest per-MPA predictor set (environmental: MHW/cold-spell/latitude; design: log area, SMR/SMCA, island, establishment year; trophic: per-MPA lobster/sheephead/urchin lnRR) and asks how much among-MPA variation in effectiveness (giant-kelp biomass lnRR; purple-urchin suppression) we can PREDICT, via univariate IV screen (Knapp-Hartung+FDR), multi-predictor meta-regression, random-forest OOB R², and **leave-one-out CV R²** (the honest arbiter). RESULT: essentially unpredictable. 0/15 predictors survive FDR; kelp in-sample pseudo-R²=0.65 is overfitting (7 preds, k=12) — RF OOB R²=−0.53 and **LOO-CV R²=−0.40** (negative = worse than the mean); urchin worse (LOO −4.6). Among-MPA effectiveness variation is real but idiosyncratic / not captured by available covariates at k~10-17. Practical upshot: report the robust AVERAGE effect + its heterogeneity (τ²/I²); do NOT claim to predict which MPAs are most effective. Tables `table_s_effectiveness_predictors.csv`, `table_s_effectiveness_models.csv`; fig `fig_s_effectiveness_predictors`; doc `docs/mpa_effectiveness_predictors.md`.

**Eisaguirre et al. 2020 reproduction (done):** `17_eisaguirre_reproduction.R` rebuilds Eisaguirre et al. (2020, Ecology 101(5):e02993, paired NCI MPA cascade, predator size structure) FROM the raw PISCO fish+swath data we hold (they archived no data/code; same PISCO program). Standalone (reads `~/Donham-Stier-CA-MPA-Data-2026/data/PISCO/`; NOT in run_all). Western NCI no-take SMRs + references, 2010-2017, pre(2010-12)/post(2014-16). RESULT: design-level results reproduce — post-disease top model = sheephead(dens×len)+protection (exact), urchins higher outside MPAs (MPA −1.35 log, p<0.001; 582 vs 327/60m²), sheephead larger inside (TL 32 vs 29), size quartiles 19/26/39 vs their 24.6/31.5/38.2, bivariate sea-star↓urchin (−0.22), sheephead↓urchin (−0.23), urchin↓kelp (−0.40) all match. BUT the PARTIAL coefficients do NOT reproduce: sea-star partial flips positive (between-reef gradient absorbed by reef RE), and the size-structured sheephead suppression dissolves under sheephead-protection collinearity (the same collinearity they flagged for lobster). Same lesson as scripts 14-16: structure robust, partial mechanism attribution specification-sensitive. Substitutions documented (14 reefs not their unnamed 7; in-situ temp not MUR SST). Tables `table_s_eisaguirre_reproduction.csv` (reproduced vs published), `table_s_eisaguirre_models.csv` (AICc); figs `fig_s_eisaguirre_fig3/fig4`; doc `docs/eisaguirre_reproduction.md`.

- Their code + data: cloned to `~/kumagai2024-comparison/` (full GitHub clone + Zenodo CC-BY-4.0 snapshot). See `~/kumagai2024-comparison/PROVENANCE.md` for sources/integrity. Key file: `repo/Processed_data/MLPA_data_summarized_wo_siteblocks.csv` (their processed subtidal data) + `repo/Processed_data/SST/MHW_cummulative_intensity_1km.tif` (their MHW-intensity raster).
- MHW definitions / SBC LTER temperature context: `~/sbc-oceanography/` (`R/31_marine_heatwaves.R`, OISST-based MHW figure, mooring/satellite/reanalysis temperature data).
- Kumagai PDF: `references/pdfs/Kumagai2024.pdf`.

## Known Limitations & Methodological Notes

These limitations are documented in code comments throughout the pipeline and should be acknowledged in the manuscript discussion section.

### Statistical
- **Bio/Den non-independence**: Biomass and density are measured on the same individuals. The `(1|MPA)` random effect partially accounts for this, but residual correlation may inflate effective sample sizes. Documented in `09_meta_analysis.R`.
- **Cross-taxa attenuation bias**: Cross-taxa meta-regressions (Table 3) use estimated effect sizes as moderators, introducing errors-in-variables bias that attenuates slopes toward zero. With k=4, power is limited. Documented in `09_meta_analysis.R`.
- **pBACIPS parallel trends**: The before-period is typically short (3-5 years), limiting power to formally test the parallel trends assumption.
- **Pseudo-I2 interpretation**: Very high I2 (>90%) reflects ecological diversity of MPAs, not model deficiency. Variance components (Table S2) provide more informative decomposition.
- **Moderator power**: SMR vs SMCA and Islands vs Mainland comparisons are exploratory due to small, unbalanced subgroup sizes.

### Data Processing
- **Proportion-based lnRR**: Response ratios are computed on proportions (standardized by site total), not raw densities. This is non-standard but necessary to harmonize across monitoring programs with different transect areas.
- **Cross-program biomass bootstrap**: KFM and LTER urchin biomass is estimated by bootstrapping from PISCO size-frequency distributions, assuming similar size structure across programs.
- **25mm urchin threshold**: Applied only to PISCO data; KFM/LTER don't collect size data so all observed urchins are included.
- **bio_macro stipe/frond note**: The kelp biomass conversion may receive frond density (LTER) vs stipe density (PISCO). Documented in data repo `01_utils.R` pending domain verification.

### Ecological
- **Extreme kelp effect sizes**: Some MPA-taxa combinations (e.g., Scorpion SMR M. pyrifera lnRR ~ 9) reflect near-zero reference-site kelp, not measurement error. Handled by inverse-variance weighting.
- **No climate covariates**: ENSO, marine heatwaves, and sea star wasting disease are potential confounds not formally included as covariates. Should be acknowledged in Discussion.
- **SSWD confound**: Sea Star Wasting Disease (2013-2014) affected urchin populations differentially inside vs outside MPAs during the study period.
- **Fishing displacement**: Potential concentration of fishing effort near MPA boundaries is not modeled.
- **M. franciscanus mixed response**: Red urchin biomass is significantly higher inside MPAs (+82%, FDR = 0.004) but density is non-significant (-33%, FDR = 0.365). This is consistent with prior findings: Malakhoff & Miller (2021, *Proc. R. Soc. B* 288: 20203061) found red urchin biomass nearly quadrupled (+397%) inside Channel Islands reserves, driven by release from fishing pressure; Teck et al. (2017, *Biol. Conserv.* 210: 321–330) showed red urchins were larger inside MPAs with greater adult biomass and reproductive biomass density. The most parsimonious explanation is a direct fishing effect (reduced harvest → larger individuals → higher biomass per capita) overlaid on indirect cascade effects. This complicates the simple trophic cascade narrative but is ecologically coherent and well-precedented.

### Reproducibility
- **renv.lock**: Package versions are captured in `renv.lock` at the repo root.
- **Dryad DOI**: `00_download_data.R` has a placeholder DOI (`10.5061/dryad.XXXXXXX`); update after Dryad upload.
- **DW test seeded**: Durbin-Watson permutation test is seeded for reproducibility.
- **Effect size extraction time**: Main analysis extracts at a fixed time point; sensitivity to alternative extraction times (t=8, 10, 14) is documented as a recommended future analysis.

## Contact

Authors: Emily Donham & Adrian Stier
Project: Conservation Letters manuscript on MPA effects on kelp forest trophic cascades

---

## File Ownership (parallel work)
- `code/` — R analysis scripts, each can be edited independently
- `data/` — READ ONLY
- `docs/` — documentation, independent
- `tables/` — generated manuscript tables, safe to regenerate
- `outputs/` — generated results, safe to regenerate
- `plots/` — generated figures, safe to regenerate
- `logs/` — processing logs, reference only
