# CA MPA Kelp Forest Analysis - Project Guidelines

## Project Overview

This project analyzes the effects of California Marine Protected Areas (MPAs) on kelp forest ecosystems using the progressive-change Before-After-Control-Impact-Pairs (pBACIPS) methodology. The analysis combines data from multiple monitoring programs (PISCO, KFM/MBON, LTER) and produces publication-quality figures and statistical summaries.

**Target Journal:** Conservation Letters
**Manuscript Title:** Restoration of trophic cascades within kelp forests following the establishment of a network of marine protected areas

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
│   ├── cache/                # Bootstrap and intermediate results (.rds)
│   ├── LANDSAT/              # Satellite kelp canopy data
│   ├── LTER/                 # LTER monitoring data
│   ├── MBON/                 # KFM/MBON monitoring data
│   └── PISCO/                # PISCO monitoring data
├── docs/                     # Documentation and manuscript
├── plots/                    # Generated figures (PDF + PNG at 600 DPI)
├── outputs/                  # Filter audits, data flow summaries, replicate effects
├── logs/                     # Pipeline execution logs
└── .agent-review/            # Multi-agent review synthesis
```

## Manuscript Figure Mapping

### Main Text Figures

| MS Figure | Description | Code Output File | Script |
|-----------|-------------|------------------|--------|
| Figure 1 | Map of MPAs with Channel Islands + inset kelp time series | `fig_01_mpa_map.pdf` | 11_figures.R |
| Figure 2 | Data processing pipeline (raw → proportion → lnRR) | `fig_02_data_processing.pdf` | 11_figures.R |
| Figure 3 | Meta-analytic mean effect sizes by taxa (RR-scaled axis) | `fig_03_mean_effects.pdf` | 11_figures.R |
| Figure 4 | Trophic cascade scatter: 4-panel (a) predator biomass vs urchin biomass, (b) urchin biomass vs kelp biomass, (c) predator density vs urchin density, (d) urchin density vs kelp biomass | `fig_04_trophic_scatter.pdf` | 11_figures.R |
| Figure 5 | Recovery trajectories: linear trends of lnRR over time for all 5 species, with lmer slopes and t=11 marker | `fig_05_recovery_curves.pdf` | 11_figures.R |

**Figure 1 Details:**
- Base map: Southern California coastline with Channel Islands
- Site markers: Shape indicates data source (square=NPS-KFM, circle=LTER, triangle=PISCO)
- Inset panels (a-d): Kelp biomass time series at Campus Point, Harris Point, South Point, Santa Barbara Island

### Supplemental Figures

| MS Figure | Description | Code Output File | Script |
|-----------|-------------|------------------|--------|
| Figure S1 | Forest plot: effect sizes by MPA for each taxa (RR-scaled axis) | `fig_s01_forest_plot.pdf` | 11_figures.R |
| Figure S2 | All taxa time series at example MPAs | `fig_s02_all_taxa_timeseries.pdf` | 11_figures.R |
| Figure S3 | Species-level GAM recovery curves with MPA spaghetti | `fig_s03_recovery_curves.pdf` | 10_temporal_analysis.R |
| Figure S4 | Species-pair phase portraits of trophic cascade | `fig_s04_cascade_phase.pdf` | 10_temporal_analysis.R |
| Figure S5 | Species-level space-time heatmaps (5 panels) | `fig_s05_triptych_heatmap.pdf` | 10_temporal_analysis.R |
| Figure S6 | Per-MPA slope comparison and cascade consistency | `fig_s06_slope_comparison.pdf` | 10_temporal_analysis.R |
| Figure S7 | Model selection distribution and variance components | `fig_s07_statistical_transparency.pdf` | 11_figures.R |
| Figure S8 | Site-level appendix: lnRR time series per taxa (5 files) | `fig_s08_appendix_*.pdf` | 11_figures.R |
| Figure S9 | Combined moderator comparisons: (a) SMR vs SMCA protection level, (b) Channel Islands vs mainland | `fig_s09_moderator_comparisons.pdf` | 13_additional_analyses.R |

### Tables

| MS Table | Description | Code Output File | Script |
|----------|-------------|------------------|--------|
| Table 1 | Average density/biomass by taxa and source | `average_responses.csv` | 07_combine_data.R |
| Table 2 | Meta-analysis summary statistics | `table_02_meta_analysis.csv` | 09_meta_analysis.R |
| Table S1 | Data availability matrix | (in manuscript) | N/A |
| Table S2 | Temporal meta-regression coefficients | `table_s_temporal_meta_regression.csv` | 10_temporal_analysis.R |
| Table S3 | Cascade consistency scores by MPA | `table_s_cascade_consistency.csv` | 10_temporal_analysis.R |
| Table S4 | Moderator meta-regression (type, size, location) | `table_s_moderator_meta_regression.csv` | 13_additional_analyses.R |
| Table S5 | AR1 sensitivity: Durbin-Watson diagnostics | `table_s_ar1_sensitivity.csv` | 10_temporal_analysis.R |

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
- `fig_s01_forest_plot.pdf` (Supplemental Figure S1)

Guidelines:
- Always export both PDF (vector) and PNG (600 DPI raster)
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
00_libraries.R           - Package dependencies
00b_color_palette.R      - Color scheme and ggplot theme
00c_analysis_constants.R - Named constants (survey areas, size thresholds, exclusions)
01_utils.R               - Utility functions
02_pBACIPS_function.R    - Core statistical function
03_data_import.R         - Import size frequency data
04_pisco_processing.R    - PISCO data processing
05_kfm_processing.R      - KFM/NPS data processing
06_lter_processing.R     - LTER data processing
06b_landsat_processing.R - Landsat satellite data
07_combine_data.R        - Combine all sources
08_effect_sizes.R        - Calculate effect sizes
09_meta_analysis.R       - Multilevel meta-analysis
10_temporal_analysis.R   - Temporal dynamics appendix (Figs S3-S6)
11_figures.R             - Publication figures (Figs 1-5, S1-S2, S7-S8)
13_additional_analyses.R - Moderator analyses (Fig S9, moderator table)
12_results_summary.R     - Generate results CSVs and markdown summary
run_all.R                - Pipeline orchestration
run_figures_only.R       - Fast figure regeneration (~17s all, ~4s single)
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
- NLS effect size SEs use delta method (`nls_difference_se()` in `08_effect_sizes.R`)
- Non-NLS fallback models are excluded from AICc competition
- Temporal meta-regression includes AR1 sensitivity analysis (Section C2 in `10_temporal_analysis.R`)

### R Code Style
- Use tidyverse style (pipes, dplyr verbs)
- Comment non-obvious operations
- Use roxygen-style documentation for functions
- Keep functions in `01_utils.R` or `02_pBACIPS_function.R`

## Data Setup

Large data files are stored in Google Drive and symlinked to the project. This allows:
- Syncing across computers via Google Drive
- Keeping the git repo small
- Sharing data with collaborators

**Google Drive location:**
```
/Users/[username]/Library/CloudStorage/GoogleDrive-astier@ucsb.edu/My Drive/Stier Lab/People/Emily Donham/Projects/Kelp MPA/data/
```

**To set up on a new machine:**
```bash
GDRIVE="/Users/$(whoami)/Library/CloudStorage/GoogleDrive-astier@ucsb.edu/My Drive/Stier Lab/People/Emily Donham/Projects/Kelp MPA/data"
PROJECT="path/to/Donham-Stier-CA-MPA-2026/data"

ln -s "$GDRIVE/MBON" "$PROJECT/MBON"
ln -s "$GDRIVE/PISCO" "$PROJECT/PISCO"
ln -s "$GDRIVE/LTER" "$PROJECT/LTER"
ln -s "$GDRIVE/LANDSAT" "$PROJECT/LANDSAT"
ln -s "$GDRIVE/MPA" "$PROJECT/MPA"
ln -s "$GDRIVE/ALL_sizefreq_2024.csv" "$PROJECT/ALL_sizefreq_2024.csv"
```

**Data files in Google Drive:**
- `MBON/` - KFM/NPS monitoring data (~1.1 GB)
- `PISCO/` - PISCO monitoring data (~113 MB)
- `LTER/` - LTER monitoring data (~27 MB)
- `LANDSAT/` - Satellite kelp canopy data
- `MPA/` - California MPA boundary shapefiles (~1.5 MB)
- `ALL_sizefreq_2024.csv` - Size frequency data (~37 MB)

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
- `docs/methodology_review.md` - Historical issue tracking and fixes (28 issues, 19 fixed)
- `docs/CHANGELOG_FOR_EMILY.md` - Collaborator changelog with action items
- `docs/ADDITIONAL_ANALYSES_PLAN.md` - Plan and rationale for supplemental analyses
- `.agent-review/SYNTHESIS.md` - Multi-agent review synthesis with prioritized action items
- This file (`CLAUDE.md`) - Project conventions for AI assistants

## Contact

Authors: Emily Donham & Adrian Stier
Project: Conservation Letters manuscript on MPA effects on kelp forest trophic cascades

---

## File Ownership (parallel work)
- `code/` — R analysis scripts, each can be edited independently
- `data/` — READ ONLY
- `docs/` — documentation, independent
- `outputs/` — generated results, safe to regenerate
- `plots/` — generated figures, safe to regenerate
- `logs/` — processing logs, reference only
