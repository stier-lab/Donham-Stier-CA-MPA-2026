# Full Repo Bug Bash Plan

## Scope
Systematic end-to-end audit of all 17 R scripts (12,699 lines) plus documentation.

## Module Map

| Scope | Scripts | Lines | Focus |
|-------|---------|-------|-------|
| **A: Foundation** | 00_libraries, 00b_color_palette, 00c_analysis_constants, 01_utils | 1,927 | Package deps, palette consistency, constants, utility functions |
| **B: Core Stats** | 02_pBACIPS_function, 08_effect_sizes | 3,514 | pBACIPS methodology, model fitting, effect size calculation |
| **C: Data Processing** | 03_data_import, 04_pisco, 05_kfm, 06_lter, 06b_landsat, 07_combine | 3,112 | Data ingestion, processing, merging |
| **D: Meta-Analysis + Summary** | 09_meta_analysis, 11_results_summary | 1,446 | Meta-regression, FDR, output tables, results markdown |
| **E: Figures** | 10_figures | 2,538 | All 14 publication figures |
| **F: Pipeline + Docs** | run_all, run_figures_only, CLAUDE.md, README.md, plots/README.md, CHANGELOG | 361+docs | Orchestration, documentation accuracy |

## Wave 1: Parallel Scope Agents (6 agents)

Each agent: read all files in scope, check for bugs/issues, fix what's safe, report.

## Wave 2: Review Agent
- Run full pipeline via `run_all.R`
- Cross-module integration validation
- Final summary
