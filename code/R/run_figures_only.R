#' ---
#' title: "Fast Figure Regeneration"
#' description: "Regenerate figures without re-running the full pipeline"
#' ---
#'
#' Loads a cached snapshot of pipeline data objects and sources only the
#' minimal scripts needed by 11_figures.R. Requires a prior full pipeline
#' run to create the snapshot.
#'
#' Usage:
#'   # Render all figures (~17 seconds vs ~2.3 minutes for full pipeline)
#'   source(here::here("code", "R", "run_figures_only.R"))
#'
#'   # Render only specific figures
#'   RENDER_FIGURES <- c("fig02", "fig03")
#'   source(here::here("code", "R", "run_figures_only.R"))
#'
#'   # Available figure names (should_render keys):
#'   #   fig01, fig02, fig03, fig04 (from 11_figures.R)
#'   #   fig_s01, fig_s02, fig_s08 (from 11_figures.R)
#'   #   fig_s04, fig_s05, fig_s06 (from 10_temporal_analysis.R, = SI Figs S3-S5)
#'   #   fig_s09 (from 13_additional_analyses.R, = SI Fig S6)
#'   #   fig_s11 (from 11_figures.R, DHARMa model diagnostics, = SI Fig S8)
#'   #   fig_s12 (from 11_figures.R, meta-analysis funnel plots + Egger's test, = SI Fig S9)
#'   #   fig_s13 (from 11_figures.R, lmer residual diagnostics, = SI Fig S10)
#'   #   fig_s14 (from 11_figures.R, NLS model selection & DW autocorrelation, = SI Fig S11)
#'   #   fig_s15 (from 11_figures.R, sensitivity & robustness summary, = SI Fig S12)
#'   #   fig_s16 (from 11_figures.R, outlier removal rationale, 4-panel, = SI Fig S13)
#'   #   fig_s17 (from 11_figures.R, temporal outlier trajectories, = SI Fig S14)
#'   #   fig_s18 (from 11_figures.R, raw data trajectories by outlier status, = SI Fig S15)
#'   #   fig_s03, fig_s07, fig_s10 (from 11_figures.R, dropped from SI but still renderable)

# Preserve RENDER_FIGURES across rm() using options (survives workspace clear)
if (exists("RENDER_FIGURES", envir = .GlobalEnv)) {
  options(.render_figures_saved = get("RENDER_FIGURES", envir = .GlobalEnv))
}

rm(list = ls())
t0 <- Sys.time()

# Restore RENDER_FIGURES so 11_figures.R can use it for selective rendering
if (!is.null(getOption(".render_figures_saved"))) {
  RENDER_FIGURES <- getOption(".render_figures_saved")
  options(.render_figures_saved = NULL)
}

cat("========================================================================\n")
cat("  Fast Figure Regeneration (run_figures_only.R)\n")
cat("========================================================================\n\n")

# --- 1. Load minimal dependency scripts ---
cat("Loading libraries...\n")
source(here::here("code", "R", "00_libraries.R"))

cat("Loading color palette and theme...\n")
source(here::here("code", "R", "00b_color_palette.R"))

cat("Loading analysis constants...\n")
source(here::here("code", "R", "00c_analysis_constants.R"))

cat("Loading utility functions...\n")
source(here::here("code", "R", "01_utils.R"))

# --- 2. Load figures snapshot ---
snapshot_path <- here::here("data", "cache", "figures_snapshot.rds")
if (!file.exists(snapshot_path)) {
  stop("Figures snapshot not found at: ", snapshot_path,
       "\nRun the full pipeline first: source(here::here('code', 'R', 'run_all.R'))",
       call. = FALSE)
}

snapshot <- readRDS(snapshot_path)
cat("Loaded snapshot from: ", format(snapshot$snapshot_time, "%Y-%m-%d %H:%M:%S"), "\n")

# Check staleness (use file modification time for reliability. Embedded snapshot_time
# may not reflect the actual file write time if the snapshot was copied or restored)
snapshot_age <- difftime(Sys.time(), file.mtime(snapshot_path), units = "hours")
if (as.numeric(snapshot_age) > 24) {
  warning("Snapshot is ", round(as.numeric(snapshot_age), 1),
          " hours old. Consider re-running the full pipeline.",
          call. = FALSE)
}

# Unpack required data objects into global environment
required_objs <- c("All.RR.sub.trans", "All.Resp.sub", "SumStats.Final", "Table2", "Site")
for (obj_name in required_objs) {
  if (is.null(snapshot[[obj_name]])) {
    stop("Snapshot missing required object: ", obj_name,
         "\nRe-run the full pipeline to regenerate the snapshot.", call. = FALSE)
  }
  assign(obj_name, snapshot[[obj_name]], envir = .GlobalEnv)
}

# Unpack optional objects
if (!is.null(snapshot$meta_biomass_full)) assign("meta_biomass_full", snapshot$meta_biomass_full, envir = .GlobalEnv)
if (!is.null(snapshot$meta_density_full)) assign("meta_density_full", snapshot$meta_density_full, envir = .GlobalEnv)
if (!is.null(snapshot$meta_biomass)) assign("meta_biomass", snapshot$meta_biomass, envir = .GlobalEnv)
if (!is.null(snapshot$meta_density)) assign("meta_density", snapshot$meta_density, envir = .GlobalEnv)
if (!is.null(snapshot$Landsat.RR)) assign("Landsat.RR", snapshot$Landsat.RR, envir = .GlobalEnv)
if (!is.null(snapshot$pertaxa_biomass)) assign("pertaxa_biomass", snapshot$pertaxa_biomass, envir = .GlobalEnv)
if (!is.null(snapshot$pertaxa_density)) assign("pertaxa_density", snapshot$pertaxa_density, envir = .GlobalEnv)

cat("Data objects loaded:", paste(required_objs, collapse = ", "), "\n\n")

# --- 3. Temporal analysis (fig_s04-s06 on disk = SI Figs S3-S5; + prediction params for Fig 4) ---
# Must run before 11_figures.R so that table_s_lmer_prediction_params_*.csv exists
source(here::here("code", "R", "10_temporal_analysis.R"))

# --- 3b. Generate main and supplemental figures ---
source(here::here("code", "R", "11_figures.R"))

# --- 4. Additional analyses (S9 + moderator table) ---
source(here::here("code", "R", "13_additional_analyses.R"))

elapsed <- round(difftime(Sys.time(), t0, units = "secs"), 1)
cat("\n========================================================================\n")
cat("  Figures-only complete in ", elapsed, " seconds\n")
cat("========================================================================\n")
