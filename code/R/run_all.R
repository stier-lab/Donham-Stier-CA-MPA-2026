#' ---
#' title: "Master Orchestration Script"
#' description: "Sources all analysis modules in order for the CA MPA kelp forest pBACIPS project"
#' author: "Emily Donham & Adrian Stier"
#' project: "Donham-Stier-CA-MPA-2026"
#' ---
#'
#' HOW TO RUN THE FULL ANALYSIS:
#'   Open R in the project directory, then:
#'     source(here::here("code", "R", "run_all.R"))
#'   This takes ~2 minutes and produces all tables, figures, and summary docs.
#'
#' HOW TO REGENERATE JUST THE FIGURES (faster, ~17 seconds):
#'     source(here::here("code", "R", "run_figures_only.R"))
#'
#' WHAT EACH SCRIPT DOES (in order):
#'   00_libraries.R         - Load R packages (tidyverse, metafor, lme4, etc.)
#'   00b_color_palette.R    - Colors, shapes, and ggplot theme for all figures
#'   00c_analysis_constants.R - Named constants (survey areas, size thresholds, etc.)
#'   01_utils.R             - Shared helper functions + excluded MPA list
#'   02_pBACIPS_function.R  - The core pBACIPS method (step/linear/asymptotic/sigmoid)
#'   03_load_harmonized_data.R - Load the 4 harmonized CSVs into R objects
#'   08_effect_sizes.R      - Calculate effect sizes for each MPA x taxa x response
#'                            -> produces SumStats.Final (142 effect sizes in the latest run)
#'   09_meta_analysis.R     - Multilevel meta-analysis: Table 2, Table 3, variance components
#'   10_temporal_analysis.R - Recovery trajectories, phase portraits, cascade consistency
#'   11_figures.R           - All main text (Figs 1-4) and most SI figures
#'   13_additional_analyses.R - Moderator comparisons (SMR vs SMCA, Islands vs Mainland)
#'   RESILIENCE MODULE (scripts 14-23; membership defined in resilience_modules.R) -
#'     the in-pipeline subset runs here, the full suite via run_resilience.R:
#'       14_heatwave_analysis.R        - MPA cascade x 2014-16 marine heatwave (Kumagai 2024 comparison)
#'       15_methods_comparison.R       - Analytical multiverse: method vs data drives effect sizes (SI)
#'       16_environmental_moderators.R - Env covariates (Kumagai PCA) as meta-regression moderators (SI)
#'       19_heatwave_replication.R     - Does resilience repeat across TWO heatwaves? (SI)
#'       21_resilience_stability.R     - Temporal stability (CV) inside vs outside (SI)
#'       22_resistance_recovery.R      - Resistance/recovery on the state variable (main-text kelp fig)
#'       23_ecological_memory.R        - Do the same reserves respond in both heatwaves? (SI)
#'       24_resilience_pipeline_check.R - Concordance gate for paper/SI resilience claims
#'     (17/18/20 need raw PISCO / the Kumagai mirror; run via run_resilience.R only.)
#'   12_results_summary.R   - Summary CSVs, RESULTS_SUMMARY.md, data flow audit
#'
#' KEY OUTPUTS:
#'   plots/          - All figures as PDF + PNG
#'   tables/         - All manuscript tables as CSV
#'   data/sumstats_final.csv - The 142 effect sizes used in meta-analysis
#'   docs/RESULTS_SUMMARY.md - Auto-generated plain-English results summary

####################################################################################################
## Setup ###########################################################################################
####################################################################################################

# Clear workspace
rm(list = ls())

# Record pipeline start time
pipeline_start <- Sys.time()

# Ensure here is available for path resolution (needed early for log file path)
if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required. Install it with install.packages('here').")
}

####################################################################################################
## File-based logging setup ########################################################################
####################################################################################################

# Create logs directory if it doesn't exist
log_dir <- here::here("logs")
if (!dir.exists(log_dir)) {
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
}

# Create timestamped log file
log_filename <- paste0("pipeline_", format(pipeline_start, "%Y%m%d_%H%M%S"), ".log")
log_filepath <- file.path(log_dir, log_filename)

# Helper function to write to both console and log file
log_message <- function(..., sep = "", append = TRUE) {
  msg <- paste(..., sep = sep)
  cat(msg)
  cat(msg, file = log_filepath, append = append)
}

# Initialize log file with header
log_message("", append = FALSE)  # Create/overwrite file
log_message("========================================================================\n")
log_message("  CA MPA Kelp Forest pBACIPS Analysis Pipeline - Log File\n")
log_message("  Started: ", format(pipeline_start, "%Y-%m-%d %H:%M:%S"), "\n")
log_message("  Log file: ", log_filepath, "\n")
log_message("========================================================================\n")
log_message("\n")

# Log system information
log_message("System Information:\n")
log_message("  R version: ", R.version.string, "\n")
log_message("  Platform:  ", R.version$platform, "\n")
log_message("  Working dir: ", getwd(), "\n")
log_message("  Project root: ", here::here(), "\n")
log_message("\n")

cat("\n")
cat("========================================================================\n")
cat("  CA MPA Kelp Forest pBACIPS Analysis Pipeline\n")
cat("  Started: ", format(pipeline_start, "%Y-%m-%d %H:%M:%S"), "\n")
cat("  Log file: ", log_filepath, "\n")
cat("========================================================================\n")
cat("\n")

####################################################################################################
## Helper: source a module with error handling #####################################################
####################################################################################################

source_module <- function(filename, label) {
  start_msg <- paste0("---- [", label, "] Sourcing ", filename, " ...\n")
  cat(start_msg)
  log_message(start_msg)

  t0 <- Sys.time()
  tryCatch(
    {
      source(here::here("code", "R", filename), local = FALSE)
      elapsed <- round(difftime(Sys.time(), t0, units = "secs"), 1)
      end_msg <- paste0("---- [", label, "] ", filename, " completed successfully (",
                        elapsed, " s)\n\n")
      cat(end_msg)
      log_message(end_msg)
    },
    error = function(e) {
      err_msg <- paste0("\n****** ERROR in module [", label, "] ", filename, " ******\n",
                        "Message: ", conditionMessage(e), "\n",
                        "Pipeline halted.\n")
      cat(err_msg)
      log_message(err_msg)
      stop("Failed at module: ", filename, " -- ", conditionMessage(e), call. = FALSE)
    }
  )
}

####################################################################################################
## Helper: verify that expected objects exist in the global environment #############################
####################################################################################################

check_objects <- function(obj_names, module_label) {
  missing <- obj_names[!vapply(obj_names, exists, logical(1), envir = globalenv())]
  if (length(missing) > 0) {
    msg <- paste0(
      "Checkpoint failed after [", module_label, "]. ",
      "Missing expected objects: ", paste(missing, collapse = ", ")
    )
    log_message(msg, "\n")
    stop(msg, call. = FALSE)
  }
  check_msg <- paste0("  Checkpoint OK: ", paste(obj_names, collapse = ", "), " found.\n\n")
  cat(check_msg)
  log_message(check_msg)
}

####################################################################################################
## Run pipeline ####################################################################################
####################################################################################################

# --- 00: Package dependencies --------------------------------------------------------------------
source_module("00_libraries.R", "00")

# --- 00b: Color palette and theme ----------------------------------------------------------------
source_module("00b_color_palette.R", "00b")

# --- 00c: Analysis constants (magic numbers) -----------------------------------------------------
source_module("00c_analysis_constants.R", "00c")

# --- 01: Utility functions ------------------------------------------------------------------------
source_module("01_utils.R", "01")

# --- 02: pBACIPS function definition --------------------------------------------------------------
source_module("02_pBACIPS_function.R", "02")

# --- 03: Load harmonized data (from data-processing repo or local data/harmonized/) ---------------
source_module("03_load_harmonized_data.R", "03")
check_objects(c("All.RR.sub.trans", "All.Resp.sub", "Site", "Landsat.RR"), "03")

# --- 08: Effect sizes -----------------------------------------------------------------------------
source_module("08_effect_sizes.R", "08")
check_objects(c("SumStats", "SumStats.Final"), "08")

# --- 09: Meta-analysis ---------------------------------------------------------------------------
source_module("09_meta_analysis.R", "09")
check_objects(c("Table2", "meta_biomass_full", "meta_density_full"), "09")

# --- Save figures snapshot (enables run_figures_only.R) ----------------------------------------
snapshot_path <- here::here("data", "cache", "figures_snapshot.rds")
figures_snapshot <- list(
  All.RR.sub.trans  = All.RR.sub.trans,
  All.Resp.sub      = All.Resp.sub,
  SumStats.Final    = SumStats.Final,
  Table2            = Table2,
  Site              = Site,
  Landsat.RR        = if (exists("Landsat.RR")) Landsat.RR else NULL,
  meta_biomass_full = if (exists("meta_biomass_full")) meta_biomass_full else NULL,
  meta_density_full = if (exists("meta_density_full")) meta_density_full else NULL,
  meta_biomass      = if (exists("meta_biomass")) meta_biomass else NULL,
  meta_density      = if (exists("meta_density")) meta_density else NULL,
  pertaxa_biomass   = if (exists("pertaxa_biomass")) pertaxa_biomass else NULL,
  pertaxa_density   = if (exists("pertaxa_density")) pertaxa_density else NULL,
  snapshot_time     = Sys.time()
)
saveRDS(figures_snapshot, snapshot_path)
snap_msg <- paste0("  Figures snapshot saved to: ", snapshot_path, "\n\n")
cat(snap_msg)
log_message(snap_msg)

# --- 10: Temporal dynamics appendix --------------------------------------------------------
source_module("10_temporal_analysis.R", "10")

# --- 11: Figures ----------------------------------------------------------------------------------
source_module("11_figures.R", "11")

# --- 13: Additional analyses (moderator comparisons, supplemental fig S9) ---
source_module("13_additional_analyses.R", "13")

# --- RESILIENCE MODULE (scripts 14-23) -----------------------------------------------------------
# Heatwave resilience + robustness suite. Module membership is defined once in
# resilience_modules.R (shared with run_resilience.R). The in-pipeline subset runs here;
# scripts 15/16 use the local Kumagai mirror for the full comparator/cold-spell checks.
# The full suite adds 17/18/20, which need raw PISCO or other external comparison inputs.
source(here::here("code", "R", "resilience_modules.R"))
for (resil_module in RESILIENCE_MODULES_IN_PIPELINE) {
  source_module(resil_module, sub("_.*", "", resil_module))
}
source_module("24_resilience_pipeline_check.R", "24")

# --- 12: Results Summary --------------------------------------------------------------------------
source_module("12_results_summary.R", "12")

####################################################################################################
## Pipeline summary ################################################################################
####################################################################################################

pipeline_end <- Sys.time()
elapsed_total <- difftime(pipeline_end, pipeline_start, units = "mins")

# Pipeline completion message
completion_msg <- paste0(
  "\n",
  "========================================================================\n",
  "  Pipeline complete\n",
  "  Finished: ", format(pipeline_end, "%Y-%m-%d %H:%M:%S"), "\n",
  "  Total elapsed time: ", round(as.numeric(elapsed_total), 2), " minutes\n",
  "========================================================================\n",
  "\n"
)
cat(completion_msg)
log_message(completion_msg)

# Print dimensions of key objects
summary_header <- "Summary of key objects:\n-----------------------------------------------------------------------\n"
cat(summary_header)
log_message(summary_header)

summary_objects <- list(
  "Site"             = "Site",
  "All.RR.sub.trans" = "All.RR.sub.trans",
  "All.Resp.sub"     = "All.Resp.sub",
  "Landsat.RR"       = "Landsat.RR",
  "SumStats"         = "SumStats",
  "SumStats.Final"   = "SumStats.Final",
  "Table2"           = "Table2"
)

for (obj_label in names(summary_objects)) {
  obj_name <- summary_objects[[obj_label]]
  if (exists(obj_name, envir = globalenv())) {
    obj <- get(obj_name, envir = globalenv())
    if (is.data.frame(obj)) {
      obj_msg <- sprintf("  %-20s %d rows x %d cols\n", obj_label, nrow(obj), ncol(obj))
    } else {
      obj_msg <- sprintf("  %-20s class: %s\n", obj_label, paste(class(obj), collapse = ", "))
    }
    cat(obj_msg)
    log_message(obj_msg)
  }
}

summary_footer <- paste0(
  "-----------------------------------------------------------------------\n",
  "Log file saved to: ", log_filepath, "\n",
  "Done.\n"
)
cat(summary_footer)
log_message(summary_footer)
