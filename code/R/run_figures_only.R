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
#'   RENDER_FIGURES <- c("fig03", "fig04")
#'   source(here::here("code", "R", "run_figures_only.R"))
#'
#'   # Available figure names:
#'   #   fig01, fig02, fig03, fig04, fig05 (from 11_figures.R)
#'   #   fig_s01, fig_s02 (from 11_figures.R)
#'   #   fig_s03, fig_s04, fig_s05, fig_s06 (from 10_temporal_analysis.R)
#'   #   fig_s07, fig_s08 (from 11_figures.R)
#'   #   fig_s09 (from 13_additional_analyses.R)

# Preserve RENDER_FIGURES across rm() using options (survives workspace clear)
if (exists("RENDER_FIGURES", envir = .GlobalEnv)) {
  options(.render_figures_saved = get("RENDER_FIGURES", envir = .GlobalEnv))
}

rm(list = ls())
t0 <- Sys.time()

# Restore RENDER_FIGURES so 10_figures.R can use it for selective rendering
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

# --- 2. Load figures snapshot ---
snapshot_path <- here::here("data", "cache", "figures_snapshot.rds")
if (!file.exists(snapshot_path)) {
  stop("Figures snapshot not found at: ", snapshot_path,
       "\nRun the full pipeline first: source(here::here('code', 'R', 'run_all.R'))",
       call. = FALSE)
}

snapshot <- readRDS(snapshot_path)
cat("Loaded snapshot from: ", format(snapshot$snapshot_time, "%Y-%m-%d %H:%M:%S"), "\n")

# Check staleness
snapshot_age <- difftime(Sys.time(), snapshot$snapshot_time, units = "hours")
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

# Unpack optional meta-analysis model objects (used by Fig S5)
if (!is.null(snapshot$meta_biomass)) assign("meta_biomass", snapshot$meta_biomass, envir = .GlobalEnv)
if (!is.null(snapshot$meta_density)) assign("meta_density", snapshot$meta_density, envir = .GlobalEnv)

cat("Data objects loaded:", paste(required_objs, collapse = ", "), "\n\n")

# --- 3. Generate figures ---
source(here::here("code", "R", "11_figures.R"))

# --- 4. Additional analyses (S9 + moderator table) ---
source(here::here("code", "R", "13_additional_analyses.R"))

elapsed <- round(difftime(Sys.time(), t0, units = "secs"), 1)
cat("\n========================================================================\n")
cat("  Figures-only complete in ", elapsed, " seconds\n")
cat("========================================================================\n")
