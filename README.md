# Restoration of Trophic Cascades in California MPA Kelp Forests

[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue.svg)](https://www.r-project.org/)
[![Status](https://img.shields.io/badge/Status-In%20Preparation-yellow.svg)]()

> **Manuscript Title:** Restoration of trophic cascades within kelp forests following the establishment of a network of marine protected areas

---

> **For co-authors — the short version.** You don't need to run anything to write the
> paper. Figures are in `plots/`, tables in `tables/`, and plain-language results
> summaries in `docs/` (start with `docs/RESULTS_SUMMARY.md`). The main paper analyses
> are scripts `08`–`13` plus the in-pipeline resilience modules `14/15/16/19/21/22/23`;
> Figure 5 comes from script `22`. Script `24` checks that the resilience tables,
> figures, and headline values still match the manuscript and supplement. A plain tour
> is in `docs/RESILIENCE_SYNTHESIS.md`, and the audit trail is in
> `docs/RESILIENCE_PIPELINE_INTEGRATION.md`. The `CLAUDE.md` file in the root is a technical config
> for an AI coding assistant — you can ignore it. Questions on any of it — just ask Adrian.

---

## Overview

This repository contains the complete analysis pipeline for evaluating the effects of California's Marine Protected Area (MPA) network on kelp forest ecosystems. Using the **progressive-change Before-After-Control-Impact-Pairs (pBACIPS)** methodology, we quantify how MPA protection has influenced trophic cascades—from predator recovery through urchin suppression to kelp restoration.

This is the public analysis/code repository. Current Word manuscripts, cover
letters, Supporting Information construction, citation-audit working files, and
publisher PDF source copies live in the private sibling manuscript repository,
`Donham-Stier-CA-MPA-MS-2026`.

### Key Findings

The analysis examines:
- **Predator recovery** (lobster, sheephead) inside vs. outside MPAs
- **Urchin population changes** in response to predator recovery
- **Kelp forest biomass** trajectories following MPA establishment
- **Trophic cascade strength** across the California MPA network

---

## Authors

**Emily Donham**<sup>1</sup> & **Adrian Stier**<sup>1</sup>

<sup>1</sup> Department of Ecology, Evolution, and Marine Biology, University of California, Santa Barbara

---

## Repository Structure

This is the **analysis repo**. It consumes harmonized CSVs produced by the sibling data-processing repo ([Donham-Stier-CA-MPA-Data-2026](https://github.com/stier-lab/Donham-Stier-CA-MPA-Data-2026)). The harmonized CSVs (~1 MB) are tracked in git, so no raw monitoring data is needed to run the main analysis.

```
Donham-Stier-CA-MPA-2026/
│
├── code/
│   └── R/                            # Analysis scripts
│       ├── 00_libraries.R            # Package dependencies
│       ├── 00b_color_palette.R       # Color scheme & ggplot theme
│       ├── 00c_analysis_constants.R  # Named constants & site exclusions
│       ├── 01_utils.R                # Utility functions
│       ├── 02_pBACIPS_function.R     # Core pBACIPS methodology
│       ├── 03_load_harmonized_data.R # Load harmonized CSVs
│       ├── 08_effect_sizes.R         # Calculate effect sizes
│       ├── 09_meta_analysis.R        # Multilevel meta-analysis
│       ├── 10_temporal_analysis.R    # Temporal dynamics (SI Figs S3-S5)
│       ├── 11_figures.R              # Publication figures (Figs 1-4, S1-S2, S7a-e)
│       ├── 12_results_summary.R      # Generate results summaries
│       ├── 13_additional_analyses.R  # Moderator analyses (SI Fig S6)
│       ├── 14-23_*.R                 # Resilience / heatwave robustness suite
│       ├── resilience_modules.R      # Shared 14-23 module membership
│       ├── 24_resilience_pipeline_check.R # Paper/SI resilience concordance gate
│       ├── run_all.R                 # Pipeline orchestration
│       ├── run_resilience.R          # Full resilience suite
│       └── run_figures_only.R        # Fast figure regeneration (~17s)
│
├── data/
│   ├── harmonized/                   # Analysis-ready CSVs (tracked in git)
│   └── cache/                        # Bootstrap results & intermediates
│
├── tables/                           # Generated manuscript tables (CSV)
├── plots/                            # Generated figures (PDF & PNG at 600 DPI)
├── outputs/                          # Filter audits, data flow summaries
├── docs/                             # Documentation & manuscript drafts
└── README.md                         # This file
```

---

## Data Sources

This analysis integrates long-term monitoring data from three major programs:

| Source | Full Name | Coverage | Data Type |
|--------|-----------|----------|-----------|
| **PISCO** | Partnership for Interdisciplinary Studies of Coastal Oceans | Statewide | Fish, invertebrate, algae surveys |
| **KFM/MBON** | Kelp Forest Monitoring / Marine Biodiversity Observation Network | Channel Islands | Long-term ecological monitoring |
| **LTER** | Santa Barbara Coastal Long-Term Ecological Research | Santa Barbara Channel | Comprehensive ecosystem data |
| **Landsat** | NASA/USGS Satellite Imagery | Statewide | Kelp canopy cover (remote sensing) |

### Data

**No raw data is needed to run this analysis.** The harmonized CSVs in `data/harmonized/` are tracked in git (~1 MB total). Raw monitoring data (~1.3 GB) lives in the sibling [data-processing repo](https://github.com/stier-lab/Donham-Stier-CA-MPA-Data-2026).

---

## Setup Instructions

### Prerequisites

- **R** (>= 4.0)
- **RStudio** (recommended)

### 1. Clone the Repository

```bash
git clone https://github.com/stier-lab/sbc-2026-donham-kelp-mpa-cascade.git Donham-Stier-CA-MPA-2026
cd Donham-Stier-CA-MPA-2026
```

### 2. Install R Dependencies

```r
# Recommended: restore exact package versions from lockfile
renv::restore()

# Alternative: install packages manually
source("code/R/00_libraries.R")
```

Package versions are captured in `renv.lock` for reproducibility. Key packages include:
`tidyverse`, `here`, `metafor`, `emmeans`, `lme4`, `cowplot`, `patchwork`, `sf`

### 3. Run the Analysis

```r
# Full pipeline (~2.3 min in the latest complete run)
source(here::here("code", "R", "run_all.R"))

# Or just regenerate figures (~17s)
source(here::here("code", "R", "run_figures_only.R"))
```

No raw monitoring data setup is required for the main pipeline. The paper-grade resilience concordance check expects the local Kumagai mirror at `~/kumagai2024-comparison` so scripts `15` and `16` can regenerate the method-vs-data comparison and cold-spell moderator rows used in the manuscript/SI audit trail. The full resilience suite can also use the sibling raw-data repo for the standalone Eisaguirre and SSWD modules.

---

## Running the Analysis

### Full Pipeline

```r
# Run complete analysis (~2.3 min in the latest complete run)
source(here::here("code", "R", "run_all.R"))
```

The pipeline loads harmonized CSVs, calculates effect sizes, runs meta-analysis, generates figures, runs the in-pipeline resilience subset, checks resilience concordance against the manuscript/SI claims, and produces results summaries. Progress is logged to the console and a timestamped file in `logs/`.

### Individual Scripts

```r
# Fast figure regeneration (~17s all figures, ~4s single figure)
source(here::here("code", "R", "run_figures_only.R"))

# Render only specific figures
RENDER_FIGURES <- c("fig03", "fig04")
source(here::here("code", "R", "run_figures_only.R"))
```

---

## Analysis Pipeline

```
  ┌─────────────────────────────────┐
  │  data/harmonized/ (4 CSVs)      │
  │  Produced by data-processing    │
  │  repo; tracked in git           │
  └───────────────┬─────────────────┘
                  ▼
  ┌─────────────────────────────────┐
  │ 03_load_harmonized_data.R       │
  │ Load & validate CSVs            │
  └───────────────┬─────────────────┘
                  ▼
  ┌─────────────────────────────────┐
  │ 08_effect_sizes.R               │
  │ pBACIPS effect sizes            │
  │ 4 candidate models per MPA      │
  └───────────────┬─────────────────┘
                  ▼
  ┌─────────────────────────────────┐
  │ 09_meta_analysis.R              │
  │ Multilevel meta-analysis        │
  │ (metafor::rma.mv)               │
  └───────────────┬─────────────────┘
                  ▼
  ┌─────────────────────────────────┐
  │ 10_temporal_analysis.R          │
  │ Temporal dynamics, SI Figs S3-S5│
  └───────────────┬─────────────────┘
                  ▼
  ┌─────────────────────────────────┐
  │ 11_figures.R                    │
  │ Publication figures             │
  │ Figs 1-4, S1-S2, S7a-e         │
  └───────────────┬─────────────────┘
                  ▼
  ┌─────────────────────────────────┐
  │ 13_additional_analyses.R        │
  │ Moderator analyses, SI Fig S6   │
  └───────────────┬─────────────────┘
                  ▼
  ┌─────────────────────────────────┐
  │ 14/15/16/19/21/22/23            │
  │ Resilience + heatwave robustness│
  └───────────────┬─────────────────┘
                  ▼
  ┌─────────────────────────────────┐
  │ 24_resilience_pipeline_check.R  │
  │ Figure 5 + SI concordance gate  │
  └───────────────┬─────────────────┘
                  ▼
  ┌─────────────────────────────────┐
  │ 12_results_summary.R            │
  │ Auto-generate results CSVs      │
  └─────────────────────────────────┘
```

---

## Outputs

### Main Text Figures

| Figure | Description | Output File |
|--------|-------------|-------------|
| **Figure 1** | Map of MPAs with Channel Islands | `plots/fig_01_mpa_map.pdf` |
| **Figure 2** | Trophic cascade case studies: 3×3 grid (predators/urchins/kelp × 3 sites), before/after MPA | `plots/fig_02_cascade_case_studies.pdf` |
| **Figure 3** | Meta-analytic mean effect sizes by taxa, RR-scaled axis | `plots/fig_03_mean_effects.pdf` |
| **Figure 4** | Recovery trajectories: 3×2 trophic grid (biomass), lmer predictions with 95% CI | `plots/fig_04_recovery_curves.pdf` |
| **Figure 5** | Giant-kelp resistance and recovery through the 2014-16 marine heatwave | `plots/fig_kelp_resilience.pdf` |

### Supplemental Figures

| Figure | Description | Output File |
|--------|-------------|-------------|
| **Figure S1** | Data processing pipeline (raw → standardized → lnRR) | `plots/fig_s01_data_processing.pdf` |
| **Figure S2** | Forest plot: effect sizes by MPA, RR-scaled axis | `plots/fig_s02_forest_plot.pdf` |
| **Figure S3** | Species-level GAM recovery curves with lmer overlay (t≤15) | `plots/fig_s04_recovery_curves.pdf` |
| **Figure S4** | Species-pair phase portraits of trophic cascade | `plots/fig_s05_cascade_phase.pdf` |
| **Figure S5** | Per-MPA slope comparison and cascade consistency | `plots/fig_s06_slope_comparison.pdf` |
| **Figure S6** | Combined moderator comparisons: SMR vs SMCA, Islands vs mainland | `plots/fig_s09_moderator_comparisons.pdf` |
| **Figure S7a-e** | Site-level appendix (5 taxa-specific panels) | `plots/fig_s08_appendix_*.pdf` |
| **Figure S8** | DHARMa residual diagnostics for NLS models | `plots/fig_s11_dharma_diagnostics.pdf` |
| **Figure S9** | Funnel plots + Egger's test for publication bias | `plots/fig_s12_funnel_plots.pdf` |
| **Figure S10** | lmer residual diagnostics (4-panel) | `plots/fig_s13_lmer_residuals.pdf` |
| **Figure S11** | NLS model selection frequency + DW statistics | `plots/fig_s14_model_selection.pdf` |
| **Figure S12** | Cook's D values + 4-method sensitivity comparison | `plots/fig_s15_sensitivity_summary.pdf` |
| **Figure S13** | Outlier removal rationale: 4-panel visualization | `plots/fig_s16_outlier_rationale.pdf` |
| **Figure S14** | Temporal trajectories by outlier status (lnRR) | `plots/fig_s17_temporal_outlier_trajectories.pdf` |
| **Figure S15** | Raw data trajectories by outlier status | `plots/fig_s18_raw_trajectories_outlier_status.pdf` |

Note: SI renumbered figures but disk filenames retain original numbering.

### Resilience / Robustness Figures

The in-pipeline resilience suite also produces SI and internal robustness figures: `fig_heatwave_cascade`, `fig_heatwave_cascade_regression`, `fig_heatwave_replication`, `fig_resistance_recovery`, `fig_kelp_resilience_paired`, `fig_ecological_memory`, `fig_resilience_stability`, `fig_s_methods_multiverse`, and `fig_s_env_moderators`. The full standalone suite can additionally produce `fig_compound_disturbance`, `fig_s_effectiveness_predictors`, and the Eisaguirre reproduction panels when external raw-data/comparison inputs are available. See `docs/RESILIENCE_SYNTHESIS.md` for manuscript status.

### Tables

| Table | Description | Output File |
|-------|-------------|-------------|
| **Table 1** | Average density/biomass by taxa | `tables/average_responses.csv` (data repo) |
| **Table 2** | Meta-analysis summary statistics | `tables/table_02_meta_analysis.csv` |
| **Table 3** | Cross-taxa trophic cascade meta-regressions | `tables/table_03_cross_taxa_meta_regression.csv` |
| **Table S1** | MPA data availability + taxa coverage (merged) | `tables/table_s_mpa_data_availability.csv` |
| **Table S1b** | Sample sizes by taxa × response × source | `tables/table_s_sample_sizes.csv` |
| **Table S2** | Variance components with CIs | `tables/table_s_variance_components.csv` |
| **Table S3a/b** | Source random effect sensitivity | `tables/table_s_source_sensitivity_*.csv` |
| **Table S4** | Temporal meta-regression (pooled + bio + den) | `tables/table_s_temporal_meta_regression.csv` |
| **Table S5** | Cascade consistency scores by MPA | `tables/table_s_cascade_consistency.csv` |
| **Table S6** | Moderator meta-regression + diagnostics | `tables/table_s_moderator_meta_regression.csv` |
| **Table S7** | Temporal diagnostics: GAM linearity + AR1 | `tables/table_s_gam_linearity_diagnostics.csv` |
| **Table S8** | Kelp leave-one-out sensitivity | `tables/table_s_kelp_leave1out.csv` |
| **Table S9** | Outlier sensitivity (4 methods) | `tables/table_s_outlier_sensitivity.csv` |

Resilience outputs checked by `code/R/24_resilience_pipeline_check.R`:

| Output | Description | Script |
|--------|-------------|--------|
| `tables/table_resistance_recovery.csv` | Main Figure 5 resistance/recovery values for kelp and cascade taxa | 22 |
| `tables/table_resistance_recovery_pairs.csv` + `tables/table_resistance_recovery_pair_counts.csv` | Pair-level support for the 9/10 and 8/10 reserve-count statements | 22 |
| `tables/table_resistance_recovery_sensitivity.csv` | Baseline-window, paired-test, CI, and leave-one-reserve-out checks for the kelp headline | 22 |
| `tables/table_heatwave_replication.csv` | Two-heatwave repeatability support | 19 |
| `tables/table_s_methods_multiverse.csv` + `tables/table_s_methods_decomposition.csv` | Kumagai-style-to-ours method bridge and data-effect support | 15 |
| `tables/table_s_env_moderators.csv` | Tests of environmental covariates mapped by Kumagai et al. | 16 |
| `tables/table_resilience_stability.csv` + `tables/table_ecological_memory.csv` | Stability and reserve-level memory support | 21, 23 |

---

## Statistical Methods

### pBACIPS Methodology

We implement the **progressive-change Before-After-Control-Impact-Pairs** approach following [Thiault et al. (2017)](https://doi.org/10.1111/2041-210X.12655). This method:

1. Pairs each MPA site with a nearby reference site
2. Calculates the log response ratio (ln RR) between paired sites
3. Fits 4 candidate models and selects best via AICc:
   - **Step**: Abrupt shift after MPA implementation
   - **Linear**: Gradual constant-rate recovery
   - **Asymptotic**: Fast initial change that saturates
   - **Sigmoid**: Delayed onset followed by rapid change

### Meta-Analysis

Effect sizes are synthesized using multilevel meta-analysis (`metafor::rma.mv`) with:
- Taxa as fixed effect moderator
- Random effects for MPA and data source
- Heterogeneity assessment (I², τ²)
- Separate models for biomass and density responses

### Quality Control

- Before-period linearity tests to verify parallel trends assumption
- Bootstrap uncertainty quantification for biomass estimates
- Cross-validation of effect sizes across data sources

---

## Key Dependencies

| Package | Purpose |
|---------|---------|
| `tidyverse` | Data wrangling and visualization |
| `here` | Reproducible file paths |
| `metafor` | Meta-analysis models |
| `emmeans` | Estimated marginal means |
| `lme4` | Mixed-effects models |
| `investr` | Nonlinear model inference |
| `cowplot` | Figure composition |
| `patchwork` | Multi-panel figure assembly |
| `sf` / `rnaturalearth` | Geographic mapping |

---

## Documentation

- **`docs/ANALYSIS_REVISIONS.md`** - Comprehensive comparison of V5 archived code vs current pipeline
- **`docs/RESULTS_SUMMARY.md`** - Auto-generated results summary with all meta-analysis values
- **`docs/RESILIENCE_SYNTHESIS.md`** - Plain-language synthesis of scripts 14-23
- **`docs/RESILIENCE_PIPELINE_INTEGRATION.md`** - Paper/SI resilience claim-to-output audit trail
- **`docs/README.md`** - Public documentation map and notes on removed manuscript/PDF artifacts

Generated browser HTML, manuscript drafts, and source-PDF libraries are
intentionally not tracked in this public analysis repo. Regenerate paper-facing
manuscript products from the private manuscript repo.

### Resilience analysis suite

Scripts 14–23 form an integrated **resilience module**. The main pipeline runs
`14/15/16/19/21/22/23`, then `24_resilience_pipeline_check.R` verifies that the
generated tables and figures still match the manuscript and Supporting Information.
Run `code/R/run_resilience.R` for the broader standalone suite; it adds the external
input modules `17/18/20` when the raw PISCO or comparison inputs are available.

The suite tests resilience across facets: core heatwave response (14), repeatability
across two heatwaves (19), resistance/recovery on the state variable (22), ecological
memory/reserve-consistency (23), temporal stability (21), method-vs-data robustness
(15), environmental moderators (16), generality across stressor types / sea-star
wasting (20), predictability of effectiveness (18), and cross-study reproduction of
Eisaguirre 2020 (17). **Theme:** robust, repeatable giant-kelp resilience to marine
heatwaves — kelp *grows* inside reserves while it declines outside (resistance 2.18x
vs 0.94x, recovery to 2.38x vs 0.83x), the *same* reserves respond in both heatwaves,
and the result is method-invariant, not modulated by environmental gradients, and not
driven by sea-star wasting. The less-resolved questions are consistent across the suite:
the urchin-mediated mechanism, predator diversity, and *predicting which reserves are
most effective* are sensitive to modeling choices or not yet explained (reserve
effectiveness is repeatable across events but not predictable from the covariates we
have), and protection does not clearly reduce year-to-year variability (the resilience
shows up in the mean state, not the variance). These are honest limitations to note,
not problems with the core result.

### Related work (compare & contrast)

We are comparing this analysis against **Kumagai et al. 2024** (*Global Change Biology*, [doi:10.1111/gcb.17620](https://doi.org/10.1111/gcb.17620)) — *"Marine Protected Areas That Preserve Trophic Cascades Promote Resilience of Kelp Forests to Marine Heatwaves."* It is a parallel study on an overlapping system/dataset (same PISCO/MLPA subtidal + T. Bell SBC LTER Landsat) with a convergent conclusion. We are using it to (1) cross-check our PISCO numbers, (2) contrast their GLMM/permutation, heatwave-resilience approach with our pBACIPS + meta-analysis, and (3) scope our own marine-heatwave (MHW) analysis. Their code + data are mirrored locally (GitHub clone + Zenodo CC-BY-4.0 snapshot); see `CLAUDE.md` → "Related Work / External Comparison" for details.

A formal **method-vs-data decomposition** (`code/R/15_methods_comparison.R`, SI supplement `docs/methods_comparison_supplement.md`) runs an analytical multiverse — a 6-waypoint method bridge that flips one analytical choice at a time from the Kumagai-style endpoint to ours on a common dataset, plus the pooled waypoints run on both datasets — to measure how much of the between-study variation in MPA × heatwave effect sizes is driven by method versus data. Headline: giant-kelp resilience is method-invariant; the urchin inference hinges on whether temporal autocorrelation is modelled; and the dataset itself is the single largest lever.

A companion **environmental-moderator analysis** (`code/R/16_environmental_moderators.R`) tests the **full set** of predictors Kumagai et al. mapped in their PCA but never modeled — marine-heatwave and cold-spell intensity, latitude, MPA size, **wave exposure**, **nitrate (nutrients)**, and **human gravity** — as meta-regression moderators of our per-MPA effect sizes. Wave exposure (HSMAX), nitrate, temperature, npp, and depth are derived for **all 34 MPAs including the Channel Islands** from Bell's (2023) island-inclusive kelp-canopy-environment product (EDI `knb-lter-sbc.162`, the per-pixel NetCDF behind Kumagai's `hsmax` and Wanner et al. 2024's island analysis), and human gravity (a population-pressure index) from Kumagai's per-kelp-patch grid (`code/R/extract_kelp_env_covariates.R` → `data/per_mpa_kelp_env.csv`). The covariates recover the expected gradients (wave: San Miguel 4.9 m → Catalina lee 0.65 m; nitrate: cold upwelling 3.0 µM → warm San Diego 0.39 µM; gravity: remote islands 0 → urban mainland ~22,000). The result is a clean null across all 35 tests: no moderator survives FDR correction, and giant-kelp resilience is not modulated by any environmental gradient, indicating that omitting these covariates does not bias the headline conclusions.

We also **reproduce Eisaguirre et al. (2020, *Ecology*)** — the closest paired-design Channel-Islands analogue of our cascade claim — from the raw PISCO monitoring data we hold, since they archived no data or code (`code/R/17_eisaguirre_reproduction.R`, doc `docs/eisaguirre_reproduction.md`). The design-level results reproduce cleanly (post-disease top model identical; urchins higher outside MPAs; sheephead larger inside; bivariate sea-star→urchin, sheephead→urchin, and urchin→kelp signs all match), but the *partial* size-structured suppression does not survive once protection is in the model (sheephead–protection collinearity) — the same lesson as the Kumagai comparison: structural conclusions are robust, fine-grained partial-mechanism attribution is specification-sensitive.

We also test whether a **different stressor type** — the 2013–14 sea-star wasting disease (SSWD) — drove the cascade response (`code/R/20_compound_disturbance.R`, doc `docs/compound_disturbance.md`). The key check is the distribution of the sunflower star (*Pycnopodia*, the one sea star that is a major urchin predator): it is a cold-water species and was **patchy in Southern California — present at only ~1/3 of reserves (cold/northern/island) and essentially absent at the rest**. Because that pre-SSWD abundance varies among reserves while the heatwave hit them all, we can ask whether the cascade response scales with sunflower-star loss. It does not: kelp recovery is the same (+1.5 lnRR) whether or not a reserve ever had sunflower stars (interaction p=0.78). So in our system the cascade resilience is attributable to the **heatwave plus fishing protection, not to sea-star loss** — the keystone-removal mechanism Eisaguirre documented at the cold Channel Islands does not generalize where the sunflower star was absent. This removes a confound and reinforces the heatwave interpretation.

Because our series runs to 2023, we can also test whether the resilience **repeats across two marine heatwaves** (`code/R/19_heatwave_replication.R`, doc `docs/heatwave_replication.md`) — the 2014–16 event and a second 2018–20 event ("Blob 2.0") that the single-event analysis had folded into a flat "after" baseline. Giant-kelp resilience repeats: it is strongly elevated inside reserves in *both* heatwaves (even more so in the second) and tracks heatwave intensity beyond ongoing recovery — the first demonstration that the response is not a one-off, which a single-event design cannot establish. Predator and urchin components remain recovery-driven or marginal, as in the single-event analysis.

Finally, we test whether the among-MPA variation in **effectiveness** is *predictable* (`code/R/18_mpa_effectiveness_predictors.R`, doc `docs/mpa_effectiveness_predictors.md`): assembling environmental, reserve-design, and trophic predictors per MPA and validating with leave-one-out CV. It is not — no predictor survives FDR, and out-of-sample R² is negative (kelp LOO-CV R² = −0.40), so the variation is real but idiosyncratic and not captured by available covariates at this sample size. The defensible statement is the robust *average* MPA effect plus its meta-analytic heterogeneity, not a prediction of which reserves are most effective.

---

## Known Limitations

Key methodological caveats documented throughout the pipeline (see `docs/ANALYSIS_REVISIONS.md` for full details):

- **Bio/Den non-independence**: Biomass and density are measured on the same individuals; the `(1|MPA)` random effect partially accounts for this but residual correlation may inflate effective sample sizes
- **Cross-taxa attenuation bias**: Cross-taxa meta-regressions (Table 3) use estimated effect sizes as moderators, introducing errors-in-variables bias
- **Proportion-based lnRR**: Response ratios are computed on proportions (standardized by site total), not raw densities -- necessary to harmonize across monitoring programs with different transect areas
- **Cross-program biomass bootstrap**: KFM and LTER urchin biomass is estimated by bootstrapping from PISCO size-frequency distributions
- **Climate/stressor attribution**: The core pBACIPS analysis is paired but observational; heatwave, environmental-moderator, and sea-star-wasting checks are now included in scripts 14–23, but ENSO and other overlapping regional drivers cannot be fully isolated
- **Extreme kelp effect sizes**: Some MPA-taxa combinations reflect near-zero reference-site kelp, handled by inverse-variance weighting
- **Reproducibility**: Package versions captured in `renv.lock`; Dryad DOI is currently a placeholder (`10.5061/dryad.XXXXXXX`)

---
## Citation

If you use this code or data, please cite:

> Donham, E. & Stier, A. (2026). Restoration of trophic cascades within kelp forests following the establishment of a network of marine protected areas. Manuscript in preparation.

Working title for the current Journal of Applied Ecology draft: *Partial
trophic-cascade recovery and kelp resilience in a southern California marine
protected-area network*.

---

## License

This project is shared for academic collaboration. Please contact the authors before using for other purposes.

---

## Acknowledgments

We thank the Partnership for Interdisciplinary Studies of Coastal Oceans (PISCO), the National Park Service Kelp Forest Monitoring Program, and the Santa Barbara Coastal LTER for data access. Funding was provided by [funding sources].

---

<p align="center">
  <i>University of California, Santa Barbara</i><br>
  <i>Department of Ecology, Evolution, and Marine Biology</i>
</p>
