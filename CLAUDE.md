# CA MPA Kelp Forest Analysis - Project Guidelines

## Project Overview

This project analyzes the effects of California Marine Protected Areas (MPAs) on kelp forest ecosystems using the progressive-change Before-After-Control-Impact-Pairs (pBACIPS) methodology. The analysis combines data from multiple monitoring programs (PISCO, KFM/MBON, LTER) and produces publication-quality figures and statistical summaries.

**Target Journal:** Conservation Letters
**Manuscript Title:** Restoration of trophic cascades within kelp forests following the establishment of a network of marine protected areas

## Directory Structure

```
Donham-Stier-CA-MPA-2026/
├── code/
│   └── R/                    # R scripts (numbered 00-10)
├── data/
│   ├── cache/                # Bootstrap and intermediate results (.rds)
│   ├── LANDSAT/              # Satellite kelp canopy data
│   ├── LTER/                 # LTER monitoring data
│   ├── MBON/                 # KFM/MBON monitoring data
│   └── PISCO/                # PISCO monitoring data
├── docs/                     # Documentation and manuscript
├── plots/                    # Generated figures
├── logs/                     # Pipeline execution logs
└── outputs/                  # Additional analysis outputs
```

## Manuscript Figure Mapping

### Main Text Figures

| MS Figure | Description | Code Output File | Script |
|-----------|-------------|------------------|--------|
| Figure 1 | Map of MPAs with Channel Islands + inset kelp time series | `fig_01_mpa_map_composite.pdf` | 10_figures.R |
| Figure 2 | Data processing pipeline (raw → proportion → lnRR) | `fig_02_data_processing.pdf` | 10_figures.R |
| Figure 3 | Mean effect sizes by taxa (meta-analysis) | `fig_03_mean_effects.pdf` | 10_figures.R |
| Figure 4 | Urchin vs kelp lnRR scatterplot | `fig_04_urchin_kelp_scatter.pdf` | 10_figures.R |

**Figure 1 Details:**
- Base map: Southern California coastline with Channel Islands
- Site markers: Shape indicates data source (square=NPS-KFM, circle=LTER, triangle=PISCO)
- Inset panels (a-f): Kelp biomass time series at Campus Point, Point Vicente, Harris Point, South Point, Gull Island, Santa Barbara Island
- Also outputs `fig_01_map_only.pdf` (map without insets)

### Supplemental Figures

| MS Figure | Description | Code Output File | Script |
|-----------|-------------|------------------|--------|
| Figure S1 | Forest plot: effect sizes by MPA for each taxa | `fig_s01_forest_plot.pdf` | 10_figures.R |
| Figure S2 | All taxa time series at example MPAs | `fig_s02_all_taxa_timeseries.pdf` | 10_figures.R |

### Tables

| MS Table | Description | Code Output File | Script |
|----------|-------------|------------------|--------|
| Table 1 | Average density/biomass by taxa and source | `average_responses.csv` | 07_combine_data.R |
| Table 2 | Meta-analysis summary statistics | `table_02_meta_analysis.csv` | 09_meta_analysis.R |
| Table S1 | Data availability matrix | (in manuscript) | N/A |

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
- `fig_01_kelp_timeseries.pdf` (Main text Figure 1)
- `fig_s01_forest_plot.pdf` (Supplemental Figure S1)

Guidelines:
- Always export both PDF (vector) and PNG (300 DPI raster)
- Use zero-padded two-digit figure numbers
- Main text: `fig_01`, `fig_02`, etc.
- Supplemental: `fig_s01`, `fig_s02`, etc.

### Data Outputs (`data/`)

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
00_libraries.R        - Package dependencies
00b_color_palette.R   - Color scheme and ggplot theme
01_utils.R            - Utility functions
02_pBACIPS_function.R - Core statistical function
03_data_import.R      - Import size frequency data
04_pisco_processing.R - PISCO data processing
05_kfm_processing.R   - KFM/NPS data processing
06_lter_processing.R  - LTER data processing
06b_landsat_processing.R - Landsat satellite data
07_combine_data.R     - Combine all sources
08_effect_sizes.R     - Calculate effect sizes
09_meta_analysis.R    - Multilevel meta-analysis
10_figures.R          - Publication figures (including Fig 1 map)
run_all.R             - Pipeline orchestration
```

## Code Style Guidelines

### File Paths
- **Always use `here::here()`** for file paths
- Never use absolute paths or `setwd()`
- Example: `here::here("data", "cache", "pisco_urchin_bootstrap.rds")`

### Output Files
- Data outputs go to `data/`
- Figures go to `plots/`
- Cache/intermediate files go to `data/cache/`
- Never write files to project root

### Statistical Methods
- Document methodological choices with inline comments
- Report uncertainty (SE, CI) for all estimates
- Use `emmeans::pairs()` for contrasts (proper covariance handling)
- Use confidence intervals (not prediction intervals) for effect sizes
- Report heterogeneity statistics (I², τ²) for meta-analyses

### R Code Style
- Use tidyverse style (pipes, dplyr verbs)
- Comment non-obvious operations
- Use roxygen-style documentation for functions
- Keep functions in `01_utils.R` or `02_pBACIPS_function.R`

## Running the Analysis

```r
# Full pipeline
source(here::here("code", "R", "run_all.R"))

# Skip bootstrap (use cached results)
# Ensure cache files exist, then run normally

# Force bootstrap recalculation
FORCE_BOOTSTRAP <- TRUE
source(here::here("code", "R", "run_all.R"))
```

## Key Datasets

| Object | Description | Created in |
|--------|-------------|------------|
| `Site` | MPA metadata and implementation dates | 03_data_import.R |
| `All.RR.sub.trans` | Combined log response ratios | 07_combine_data.R |
| `All.Resp.sub` | Combined raw response data | 07_combine_data.R |
| `SumStats.Final` | Effect sizes for all MPA-taxa combinations | 08_effect_sizes.R |
| `Table2` | Meta-analysis summary by taxa | 09_meta_analysis.R |

## Adding New Analyses

When adding new figures or outputs:

1. **Main text figures**: `fig_{NN}_{name}.{pdf|png}` - increment from existing
2. **Supplemental figures**: `fig_s{NN}_{name}.{pdf|png}` - use s prefix
3. **Data outputs**: Use snake_case, write to `data/`
4. **Update this file** with the manuscript mapping

## Documentation

- `docs/MPA_Kelp_MS_V5.pdf` - Current manuscript draft
- `docs/STATISTICAL_REVIEW.md` - Statistical methodology and fixes
- This file (`CLAUDE.md`) - Project conventions for AI assistants

## Contact

Authors: Emily Donham & Adrian Stier
Project: Conservation Letters manuscript on MPA effects on kelp forest trophic cascades
