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

```
Donham-Stier-CA-MPA-2026/
│
├── code/
│   └── R/                          # Analysis scripts (numbered 00-11)
│       ├── 00_libraries.R          # Package dependencies
│       ├── 00b_color_palette.R     # Color scheme & ggplot theme
│       ├── 00c_analysis_constants.R# Named constants & site exclusions
│       ├── 01_utils.R              # Utility functions
│       ├── 02_pBACIPS_function.R   # Core pBACIPS methodology
│       ├── 03_data_import.R        # Import & clean data
│       ├── 04_pisco_processing.R   # PISCO data processing
│       ├── 05_kfm_processing.R     # KFM/NPS data processing
│       ├── 06_lter_processing.R    # LTER data processing
│       ├── 06b_landsat_processing.R# Satellite kelp canopy data
│       ├── 07_combine_data.R       # Merge all data sources
│       ├── 08_effect_sizes.R       # Calculate effect sizes
│       ├── 09_meta_analysis.R      # Multilevel meta-analysis
│       ├── 10_figures.R            # Publication figures
│       ├── 11_results_summary.R    # Generate results summaries
│       └── run_all.R               # Pipeline orchestration
│
├── data/
│   ├── cache/                      # Bootstrap results & intermediates
│   ├── MBON/ → [Google Drive]      # KFM/NPS monitoring data
│   ├── PISCO/ → [Google Drive]     # PISCO monitoring data
│   ├── LTER/ → [Google Drive]      # LTER monitoring data
│   ├── LANDSAT/ → [Google Drive]   # Satellite kelp canopy data
│   └── ALL_sizefreq_2024.csv →     # Size frequency data
│
├── plots/                          # Generated figures (PDF & PNG)
├── docs/                           # Documentation & manuscript drafts
├── CLAUDE.md                       # AI assistant guidelines
└── README.md                       # This file
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

### Data Access

Due to file sizes (~1.2 GB total), raw data is stored in Google Drive and symlinked to the project directory. See [Setup Instructions](#setup-instructions) below.

---

## Setup Instructions

### Prerequisites

- **R** (>= 4.0)
- **RStudio** (recommended)
- **Google Drive** desktop app with access to shared data folder

### 1. Clone the Repository

```bash
git clone https://github.com/your-username/Donham-Stier-CA-MPA-2026.git
cd Donham-Stier-CA-MPA-2026
```

### 2. Install R Dependencies

```r
# Open R/RStudio and run:
source("code/R/00_libraries.R")
```

This will check for and install required packages including:
`tidyverse`, `here`, `metafor`, `emmeans`, `lme4`, `cowplot`, `maps`, `mapdata`

### 3. Set Up Data Symlinks

Large data files are stored in Google Drive. Create symlinks to connect them:

```bash
# Set paths (adjust for your system)
GDRIVE="/Users/$(whoami)/Library/CloudStorage/GoogleDrive-astier@ucsb.edu/My Drive/Stier Lab/People/Emily Donham/Projects/Kelp MPA/data"
PROJECT="/path/to/Donham-Stier-CA-MPA-2026/data"

# Create symlinks
ln -s "$GDRIVE/MBON" "$PROJECT/MBON"
ln -s "$GDRIVE/PISCO" "$PROJECT/PISCO"
ln -s "$GDRIVE/LTER" "$PROJECT/LTER"
ln -s "$GDRIVE/LANDSAT" "$PROJECT/LANDSAT"
ln -s "$GDRIVE/ALL_sizefreq_2024.csv" "$PROJECT/ALL_sizefreq_2024.csv"
```

### 4. Verify Setup

```r
# Check that data files are accessible
list.files(here::here("data", "PISCO"))
list.files(here::here("data", "MBON"))
```

---

## Running the Analysis

### Full Pipeline

```r
# Run complete analysis from data import through figures
source(here::here("code", "R", "run_all.R"))
```

The pipeline executes scripts in sequence with error handling and checkpoint validation. Progress is logged to the console.

### Individual Scripts

Run specific analysis steps:

```r
# Just regenerate figures
source(here::here("code", "R", "10_figures.R"))

# Just run meta-analysis
source(here::here("code", "R", "09_meta_analysis.R"))
```

### Bootstrap Options

Bootstrap calculations are cached to speed up subsequent runs:

```r
# Use cached bootstrap results (default)
source(here::here("code", "R", "run_all.R"))

# Force recalculation of all bootstraps
FORCE_BOOTSTRAP <- TRUE
source(here::here("code", "R", "run_all.R"))
```

---

## Analysis Pipeline

```
┌─────────────────────────────────────────────────────────────────┐
│                        DATA IMPORT                               │
│  03_data_import.R: Site metadata, size frequency data           │
└─────────────────────────────────┬───────────────────────────────┘
                                  │
         ┌────────────────────────┼────────────────────────┬────────┐
         ▼                        ▼                        ▼        ▼
┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐ │
│ 04_pisco        │    │ 05_kfm          │    │ 06_lter         │ │
│ processing.R    │    │ processing.R    │    │ processing.R    │ │
│                 │    │                 │    │                 │ │
│ PISCO swath &   │    │ KFM monitoring  │    │ LTER time       │ │
│ fish surveys    │    │ data            │    │ series          │ │
└────────┬────────┘    └────────┬────────┘    └────────┬────────┘ │
         │                      │                      │           │
         │                      │                      │           ▼
         │                      │                      │  ┌─────────────────┐
         │                      │                      │  │ 06b_landsat     │
         │                      │                      │  │ processing.R    │
         │                      │                      │  │                 │
         │                      │                      │  │ Satellite kelp  │
         │                      │                      │  │ canopy data     │
         │                      │                      │  └────────┬────────┘
         │                      │                      │           │
         └──────────────────────┼──────────────────────┴───────────┘
                                ▼
                  ┌─────────────────────────┐
                  │ 07_combine_data.R       │
                  │                         │
                  │ Merge all sources       │
                  │ Calculate log RR        │
                  └────────────┬────────────┘
                               ▼
                  ┌─────────────────────────┐
                  │ 08_effect_sizes.R       │
                  │                         │
                  │ pBACIPS effect sizes    │
                  │ 4 candidate models      │
                  └────────────┬────────────┘
                               ▼
                  ┌─────────────────────────┐
                  │ 09_meta_analysis.R      │
                  │                         │
                  │ Multilevel meta-analysis│
                  │ (metafor::rma.mv)       │
                  └────────────┬────────────┘
                               ▼
                  ┌─────────────────────────┐
                  │ 10_figures.R            │
                  │                         │
                  │ Publication figures     │
                  │ Main text + supplement  │
                  └────────────┬────────────┘
                               ▼
                  ┌─────────────────────────┐
                  │ 11_results_summary.R    │
                  │                         │
                  │ Auto-generate results   │
                  │ CSVs and markdown       │
                  └─────────────────────────┘
```

---

## Outputs

### Main Text Figures

| Figure | Description | Output File |
|--------|-------------|-------------|
| **Figure 1** | Map of MPAs with Channel Islands + kelp time series insets | `plots/fig_01_mpa_map.pdf` |
| **Figure 2** | Data processing pipeline visualization | `plots/fig_02_data_processing.pdf` |
| **Figure 3** | Mean effect sizes by taxa (meta-analysis) | `plots/fig_03_mean_effects.pdf` |
| **Figure 4** | Urchin vs kelp effect size relationship | `plots/fig_04_urchin_kelp_scatter.pdf` |

### Supplemental Figures

| Figure | Description | Output File |
|--------|-------------|-------------|
| **Figure S1** | Forest plot: effect sizes by MPA | `plots/fig_s01_forest_plot.pdf` |
| **Figure S2** | All taxa time series at example MPAs | `plots/fig_s02_all_taxa_timeseries.pdf` |

### Tables

| Table | Description | Output File |
|-------|-------------|-------------|
| **Table 1** | Average density/biomass by taxa | `data/average_responses.csv` |
| **Table 2** | Meta-analysis summary statistics | `data/table_02_meta_analysis.csv` |

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
| `maps` / `mapdata` | Geographic mapping |

---

## Documentation

- **`CLAUDE.md`** - Project conventions for AI assistants
- **`docs/STATISTICAL_REVIEW.md`** - Statistical methodology and fixes
- **`docs/methodology_review.md`** - Additional methodology notes
- **`docs/MPA_Kelp_MS_V5.pdf`** - Current manuscript draft

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
