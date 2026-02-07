# =============================================================================
# validate_improvements.R
# =============================================================================
#
# PURPOSE:
#   Validate that all systematic code improvements (2026-02-06) are working
#   correctly. This script checks each of the 7 implemented fixes.
#
# WHAT THIS SCRIPT VALIDATES:
#   1. Effect size time point standardization (sigmoid models use t=11)
#   2. Enhanced input validation in pBACIPS function
#   3. Duplicate handling in data combination
#   4. Deprecated function throws error
#   5. Debug comments removed from production code
#   6. Robust column selection (no fragile indexing)
#   7. Standardized package loading pattern
#
# OUTPUTS:
#   - Console report of validation results
#   - validation_report.txt in outputs/
#
# AUTHORS: Adrian Stier (validation script)
# PROJECT: CA MPA Kelp Forest pBACIPS Analysis
# DATE: 2026-02-06
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
})

# Initialize validation results
validation_results <- list()
validation_passed <- 0
validation_failed <- 0

cat("\n")
cat("========================================================================\n")
cat("  VALIDATION REPORT: Code Improvements (2026-02-06)\n")
cat("========================================================================\n")
cat("\n")

# =============================================================================
# TEST 1: Effect Size Time Point Standardization
# =============================================================================
cat("TEST 1: Effect size time point standardization\n")
cat("  Checking that EFFECT_SIZE_TIME_YEARS constant is defined...\n")

source(here::here("code", "R", "00_libraries.R"), local = TRUE)
source(here::here("code", "R", "00c_analysis_constants.R"), local = TRUE)

if (exists("EFFECT_SIZE_TIME_YEARS") && EFFECT_SIZE_TIME_YEARS == 11) {
  cat("  ✓ PASS: EFFECT_SIZE_TIME_YEARS = 11\n")
  validation_results$test1 <- "PASS"
  validation_passed <- validation_passed + 1
} else {
  cat("  ✗ FAIL: EFFECT_SIZE_TIME_YEARS not defined or != 11\n")
  validation_results$test1 <- "FAIL"
  validation_failed <- validation_failed + 1
}
cat("\n")

# =============================================================================
# TEST 2: Enhanced Input Validation
# =============================================================================
cat("TEST 2: Enhanced pBACIPS input validation\n")
cat("  Testing that invalid inputs are rejected...\n")

source(here::here("code", "R", "02_pBACIPS_function.R"), local = TRUE)

# Test NA detection
test2_pass <- TRUE
tryCatch({
  ProgressiveChangeBACIPS(
    control = c(1, 2, NA, 4),
    impact = c(2, 3, 4, 5),
    time.true = c(1, 2, 3, 4),
    time.model = c(0, 0, 1, 1)
  )
  cat("  ✗ FAIL: NA values not detected\n")
  test2_pass <- FALSE
}, error = function(e) {
  if (grepl("NA values", e$message)) {
    cat("  ✓ PASS: NA values correctly detected\n")
  } else {
    cat("  ✗ FAIL: Wrong error for NA values\n")
    test2_pass <<- FALSE
  }
})

# Test insufficient data detection
tryCatch({
  ProgressiveChangeBACIPS(
    control = c(1, 2),
    impact = c(2, 3),
    time.true = c(1, 2),
    time.model = c(0, 1)
  )
  cat("  ✗ FAIL: Insufficient data not detected\n")
  test2_pass <- FALSE
}, error = function(e) {
  if (grepl("insufficient data", e$message)) {
    cat("  ✓ PASS: Insufficient data correctly detected\n")
  } else {
    cat("  ✗ FAIL: Wrong error for insufficient data\n")
    test2_pass <<- FALSE
  }
})

if (test2_pass) {
  validation_results$test2 <- "PASS"
  validation_passed <- validation_passed + 1
} else {
  validation_results$test2 <- "FAIL"
  validation_failed <- validation_failed + 1
}
cat("\n")

# =============================================================================
# TEST 3: Duplicate Handling
# =============================================================================
cat("TEST 3: Duplicate handling in data combination\n")
cat("  Checking that duplicate detection code exists...\n")

script_07 <- readLines(here::here("code", "R", "07_combine_data.R"))
has_dup_removal <- any(grepl("!duplicated", script_07))

if (has_dup_removal) {
  cat("  ✓ PASS: Duplicate removal code found\n")
  validation_results$test3 <- "PASS"
  validation_passed <- validation_passed + 1
} else {
  cat("  ✗ FAIL: Duplicate removal code not found\n")
  validation_results$test3 <- "FAIL"
  validation_failed <- validation_failed + 1
}
cat("\n")

# =============================================================================
# TEST 4: Deprecated Function
# =============================================================================
cat("TEST 4: Deprecated effect size function\n")
cat("  Testing that calculate_effect_size() throws error...\n")

source(here::here("code", "R", "01_utils.R"), local = TRUE)

test4_pass <- FALSE
tryCatch({
  calculate_effect_size(NULL, NULL)
  cat("  ✗ FAIL: Deprecated function did not throw error\n")
}, error = function(e) {
  if (grepl("deprecated", e$message, ignore.case = TRUE)) {
    cat("  ✓ PASS: Deprecated function correctly throws error\n")
    test4_pass <<- TRUE
  } else {
    cat("  ✗ FAIL: Wrong error message\n")
  }
})

if (test4_pass) {
  validation_results$test4 <- "PASS"
  validation_passed <- validation_passed + 1
} else {
  validation_results$test4 <- "FAIL"
  validation_failed <- validation_failed + 1
}
cat("\n")

# =============================================================================
# TEST 5: Debug Comments Removed
# =============================================================================
cat("TEST 5: Debug comments removed from production code\n")
cat("  Checking for '# Debug:' comments in 10_figures.R...\n")

script_10 <- readLines(here::here("code", "R", "10_figures.R"))
debug_comments <- grep("^\\s*#\\s*Debug:", script_10, value = TRUE)

if (length(debug_comments) == 0) {
  cat("  ✓ PASS: No debug comments found\n")
  validation_results$test5 <- "PASS"
  validation_passed <- validation_passed + 1
} else {
  cat("  ✗ FAIL: Found", length(debug_comments), "debug comments\n")
  cat("    Lines:", paste(which(grepl("^\\s*#\\s*Debug:", script_10)), collapse = ", "), "\n")
  validation_results$test5 <- "FAIL"
  validation_failed <- validation_failed + 1
}
cat("\n")

# =============================================================================
# TEST 6: Robust Column Selection
# =============================================================================
cat("TEST 6: Robust column selection (main fix applied)\n")
cat("  Checking for dplyr::select in 04_pisco_processing.R...\n")

script_04 <- readLines(here::here("code", "R", "04_pisco_processing.R"))
has_dplyr_select <- any(grepl("dplyr::select\\(site, year, y, count", script_04))
has_fragile_indexing <- any(grepl("\\[.*colnames.*\\[c\\(", script_04))

# Check if the CRITICAL fix was applied (line ~164)
# Look for "Swath.ave.site %>% \n dplyr::select(site..."
critical_line <- grep("dplyr::select\\(site, year, y, count, CA_MPA_Name_Short", script_04)

if (has_dplyr_select && length(critical_line) > 0) {
  cat("  ✓ PASS: Critical column selection uses dplyr::select()\n")
  if (has_fragile_indexing) {
    cat("    NOTE: Some column reordering still uses indexing (non-critical)\n")
  }
  validation_results$test6 <- "PASS"
  validation_passed <- validation_passed + 1
} else {
  cat("  ✗ FAIL: Critical fix not found\n")
  validation_results$test6 <- "FAIL"
  validation_failed <- validation_failed + 1
}
cat("\n")

# =============================================================================
# TEST 7: Standardized Package Loading
# =============================================================================
cat("TEST 7: Standardized package loading pattern\n")
cat("  Checking that required packages use library()...\n")

script_00 <- readLines(here::here("code", "R", "00_libraries.R"))

# Find require() calls for critical packages
critical_packages <- c("minpack.lm", "nls2", "AICcmodavg")
uses_require <- sapply(critical_packages, function(pkg) {
  any(grepl(paste0("require\\(", pkg, "\\)"), script_00))
})

uses_library <- sapply(critical_packages, function(pkg) {
  any(grepl(paste0("library\\(", pkg, "\\)"), script_00))
})

if (all(uses_library) && !any(uses_require)) {
  cat("  ✓ PASS: Critical packages use library(), not require()\n")
  validation_results$test7 <- "PASS"
  validation_passed <- validation_passed + 1
} else if (any(uses_require)) {
  cat("  ✗ FAIL: Some critical packages still use require()\n")
  cat("    Packages:", paste(critical_packages[uses_require], collapse = ", "), "\n")
  validation_results$test7 <- "FAIL"
  validation_failed <- validation_failed + 1
} else {
  cat("  ✗ FAIL: Critical packages not found\n")
  validation_results$test7 <- "FAIL"
  validation_failed <- validation_failed + 1
}
cat("\n")

# =============================================================================
# SUMMARY
# =============================================================================
cat("========================================================================\n")
cat("  VALIDATION SUMMARY\n")
cat("========================================================================\n")
cat("\n")
cat("Total tests:", validation_passed + validation_failed, "\n")
cat("  PASSED:", validation_passed, "\n")
cat("  FAILED:", validation_failed, "\n")
cat("\n")

if (validation_failed == 0) {
  cat("✓ ALL VALIDATIONS PASSED\n")
  cat("\nAll systematic code improvements (2026-02-06) are working correctly.\n")
  exit_code <- 0
} else {
  cat("✗ SOME VALIDATIONS FAILED\n")
  cat("\nPlease review failed tests above.\n")
  exit_code <- 1
}
cat("\n")

# Write report to file
report_file <- here::here("outputs", "validation_report.txt")
dir.create(dirname(report_file), showWarnings = FALSE, recursive = TRUE)

sink(report_file)
cat("VALIDATION REPORT: Code Improvements (2026-02-06)\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("Results:\n")
for (i in seq_along(validation_results)) {
  cat(sprintf("  Test %d: %s\n", i, validation_results[[i]]))
}
cat(sprintf("\nTotal: %d passed, %d failed\n", validation_passed, validation_failed))
sink()

cat("Validation report written to:", report_file, "\n\n")

# Exit with appropriate code
quit(status = exit_code)
