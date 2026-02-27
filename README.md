# Restoration of Trophic Cascades in California MPA Kelp Forests

[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue.svg)](https://www.r-project.org/)
[![Target Journal](https://img.shields.io/badge/Target-Conservation%20Letters-green.svg)](https://conbio.onlinelibrary.wiley.com/journal/1755263x)
[![Status](https://img.shields.io/badge/Status-In%20Preparation-yellow.svg)]()

> **Manuscript Title:** Restoration of trophic cascades within kelp forests following the establishment of a network of marine protected areas

---

## Overview

This repository contains the complete analysis pipeline for evaluating the effects of California's Marine Protected Area (MPA) network on kelp forest ecosystems. Using the **progressive-change Before-After-Control-Impact-Pairs (pBACIPS)** methodology, we quantify how MPA protection has influenced trophic cascades—from predator recovery through urchin suppression to kelp restoration.

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

This is the **analysis repo**. It consumes harmonized CSVs produced by the sibling data-processing repo ([Donham-Stier-CA-MPA-Data-2026](https://github.com/stier-lab/Donham-Stier-CA-MPA-Data-2026)). The harmonized CSVs (~1 MB) are tracked in git, so no raw data is needed to run the analysis.

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
│       ├── run_all.R                 # Pipeline orchestration
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
git clone https://github.com/stier-lab/Donham-Stier-CA-MPA-2026.git
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
# Full pipeline (~0.6 min)
source(here::here("code", "R", "run_all.R"))

# Or just regenerate figures (~17s)
source(here::here("code", "R", "run_figures_only.R"))
```

No symlinks or raw data setup required.

---

## Running the Analysis

### Full Pipeline

```r
# Run complete analysis (~0.6 min)
source(here::here("code", "R", "run_all.R"))
```

The pipeline loads harmonized CSVs, calculates effect sizes, runs meta-analysis, generates figures, and produces results summaries. Progress is logged to the console.

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

Note: SI renumbered figures but disk filenames retain original numbering.

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
- **`docs/MPA_Kelp_MS_V5.pdf`** - Current manuscript draft
- **`docs/RESULTS_SUMMARY.md`** - Auto-generated results summary with all meta-analysis values
- **`docs/results_report.html`** - Interactive HTML report of all pipeline outputs
- **`docs/supporting_information.Rmd`** - Supporting Information source (R Markdown)
- **`docs/supporting_information.html`** - Rendered Supporting Information
- **`docs/si_style.css`** - Supporting Information custom styles
- **`docs/si_script.html`** - Supporting Information JavaScript (lightbox, progress bar)

---

## Known Limitations

Key methodological caveats documented throughout the pipeline (see `docs/ANALYSIS_REVISIONS.md` for full details):

- **Bio/Den non-independence**: Biomass and density are measured on the same individuals; the `(1|MPA)` random effect partially accounts for this but residual correlation may inflate effective sample sizes
- **Cross-taxa attenuation bias**: Cross-taxa meta-regressions (Table 3) use estimated effect sizes as moderators, introducing errors-in-variables bias
- **Proportion-based lnRR**: Response ratios are computed on proportions (standardized by site total), not raw densities -- necessary to harmonize across monitoring programs with different transect areas
- **Cross-program biomass bootstrap**: KFM and LTER urchin biomass is estimated by bootstrapping from PISCO size-frequency distributions
- **No climate covariates**: ENSO, marine heatwaves, and sea star wasting disease are potential confounds not formally included as covariates
- **Extreme kelp effect sizes**: Some MPA-taxa combinations reflect near-zero reference-site kelp, handled by inverse-variance weighting
- **Reproducibility**: Package versions captured in `renv.lock`; Dryad DOI is currently a placeholder (`10.5061/dryad.XXXXXXX`)

---
## Citation

If you use this code or data, please cite:

> Donham, E. & Stier, A. (2026). Restoration of trophic cascades within kelp forests following the establishment of a network of marine protected areas. *Conservation Letters* (in preparation).

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
