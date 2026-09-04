# =============================================================================
# 00_download_data.R - Public data availability check
# =============================================================================
#
# PURPOSE:
#   This public analysis repository does not download or transform the raw
#   monitoring data. The analysis-ready harmonized CSVs needed by run_all.R are
#   tracked in data/harmonized/. Raw-data processing lives in the superseded
#   data-processing repository and will be replaced by a formal data archive when
#   the manuscript package is released.
#
# USAGE:
#   source(here::here("code", "R", "00_download_data.R"))
#
#   This verifies that the repo-controlled harmonized inputs are present. To run
#   the analysis itself, use:
#
#   source(here::here("code", "R", "run_all.R"))
#
# =============================================================================

library(here)

HARMONIZED_FILES <- c(
  "harmonized_response_ratios.csv",
  "harmonized_raw_responses.csv",
  "harmonized_landsat_rr.csv",
  "harmonized_site_metadata.csv"
)

check_harmonized_data <- function(dir = here::here("data", "harmonized")) {
  paths <- file.path(dir, HARMONIZED_FILES)
  existing <- file.exists(paths)
  list(
    all_present = all(existing),
    directory = dir,
    present = HARMONIZED_FILES[existing],
    missing = HARMONIZED_FILES[!existing],
    n_present = sum(existing),
    n_total = length(HARMONIZED_FILES)
  )
}

download_from_dryad <- function(...) {
  stop(
    "Dryad/raw-data download is not configured in this public analysis repo.\n",
    "The tracked data/harmonized/*.csv files are the reproducible inputs for run_all.R.\n",
    "Use the sibling data-processing repository only when regenerating the raw-data handoff.",
    call. = FALSE
  )
}

verify_checksums <- function(checksum_file = here::here("dryad_staging", "checksums.md5")) {
  if (!file.exists(checksum_file)) {
    message("Checksum file not found: ", checksum_file)
    message("The current public analysis handoff is verified by git-tracked harmonized CSVs.")
    return(invisible(FALSE))
  }

  lines <- readLines(checksum_file)
  lines <- lines[nzchar(trimws(lines))]

  n_ok <- 0L
  n_fail <- 0L
  for (line in lines) {
    parts <- strsplit(trimws(line), "\\s+", perl = TRUE)[[1]]
    expected_md5 <- parts[1]
    rel_path <- parts[2]
    full_path <- here::here(rel_path)

    if (!file.exists(full_path)) {
      message("MISSING: ", rel_path)
      n_fail <- n_fail + 1L
      next
    }

    actual_md5 <- unname(tools::md5sum(full_path))
    if (identical(actual_md5, expected_md5)) {
      n_ok <- n_ok + 1L
    } else {
      message("MISMATCH: ", rel_path)
      message("  Expected: ", expected_md5)
      message("  Actual:   ", actual_md5)
      n_fail <- n_fail + 1L
    }
  }

  message(sprintf("Checksum verification: %d OK, %d failed", n_ok, n_fail))
  invisible(n_fail == 0L)
}

status <- check_harmonized_data()
if (status$all_present) {
  message("Harmonized data OK: ", status$n_present, " of ", status$n_total,
          " files present in ", status$directory)
} else {
  stop(
    "Harmonized analysis inputs are missing from ", status$directory, ":\n  - ",
    paste(status$missing, collapse = "\n  - "),
    "\nRestore data/harmonized/ from git before running run_all.R.",
    call. = FALSE
  )
}
