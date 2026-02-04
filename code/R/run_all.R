#' ---
#' title: "Master Orchestration Script"
#' description: "Sources all analysis modules in order for the CA MPA kelp forest pBACIPS project"
#' author: "Emily Donham & Adrian Stier"
#' project: "Donham-Stier-CA-MPA-2026"
#' ---
#'
#' This script runs the full analysis pipeline from data import through publication
#' figures. Each module is sourced in sequence with error handling and checkpoint
#' validation of key intermediate objects.
#'
#' Usage:
#'   source("code/R/run_all.R")
#'   # or from the project root:
#'   source(here::here("code", "R", "run_all.R"))
#'
#' Pipeline:
#'   00  - Load packages and set global options
#'   00b - Define unified color palette and ggplot theme
#'   01  - Utility functions (biomass conversions, species standardization)
#'   02  - Define pBACIPS function
#'   03  - Import size frequency data and site metadata
#'   04  - PISCO data wrangling
#'   05  - KFM / NPS data processing
#'   06  - LTER data processing
#'   06b - Landsat remote sensing data
#'   07  - Combine all data sources
#'   08  - Calculate effect sizes (pBACIPS)
#'   09  - Multilevel meta-analysis
#'   10  - Publication figures

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

# --- 03: Data import ------------------------------------------------------------------------------
source_module("03_data_import.R", "03")
check_objects(c("Site"), "03")

# --- 04: PISCO processing -------------------------------------------------------------------------
source_module("04_pisco_processing.R", "04")
check_objects(c("Swath.join.sub"), "04")

# --- 05: KFM processing --------------------------------------------------------------------------
source_module("05_kfm_processing.R", "05")
check_objects(c("KFM.join.ave"), "05")

# --- 06: LTER processing -------------------------------------------------------------------------
source_module("06_lter_processing.R", "06")
check_objects(c("LTER.join.ave"), "06")

# --- 06b: Landsat remote sensing ------------------------------------------------------------------
source_module("06b_landsat_processing.R", "06b")

# --- 07: Combine data sources ---------------------------------------------------------------------
source_module("07_combine_data.R", "07")
check_objects(c("All.RR.sub.trans", "All.Resp.sub"), "07")

# --- 08: Effect sizes -----------------------------------------------------------------------------
source_module("08_effect_sizes.R", "08")
check_objects(c("SumStats", "SumStats.Final"), "08")

# --- 09: Meta-analysis ---------------------------------------------------------------------------
source_module("09_meta_analysis.R", "09")
check_objects(c("Table2", "meta_biomass", "meta_density"), "09")

# --- 10: Figures ----------------------------------------------------------------------------------
source_module("10_figures.R", "10")

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
  "Swath.join.sub"   = "Swath.join.sub",
  "KFM.join.ave"     = "KFM.join.ave",
  "LTER.join.ave"    = "LTER.join.ave",
  "All.RR.sub.trans" = "All.RR.sub.trans",
  "All.Resp.sub"     = "All.Resp.sub",
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
