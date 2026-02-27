# =============================================================================
# 00_download_data.R — Download data from Dryad and set up data/ directory
# =============================================================================
#
# This script downloads the data archived on Dryad, including:
#   - Raw monitoring data (PISCO, MBON, LTER, Landsat, MPA boundaries)
#   - Harmonized analysis-ready CSVs (response ratios, raw responses, etc.)
#
# Usage:
#   source(here::here("code", "R", "00_download_data.R"))
#
# You do NOT need this script if:
#   - You only want to run the analysis (harmonized CSVs are tracked in git)
#   - You only want to regenerate figures (use run_figures_only.R instead)
#   - You already have data/ populated via Google Drive symlinks
#
# =============================================================================

library(here)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# TODO: Replace with actual Dryad DOI after upload
DRYAD_DOI <- "doi:10.5061/dryad.XXXXXXX"

# Data directory
DATA_DIR <- here::here("data")

# Expected subdirectories
EXPECTED_DIRS <- c("PISCO", "MBON", "LTER", "LANDSAT", "MPA")

# ---------------------------------------------------------------------------
# Check if data already exists
# ---------------------------------------------------------------------------

check_data_exists <- function() {
  required_files <- c(
    file.path(DATA_DIR, "PISCO", "MLPA_kelpforest_swath_2024.csv"),
    file.path(DATA_DIR, "PISCO", "MLPA_kelpforest_fish_2024.csv"),
    file.path(DATA_DIR, "PISCO", "master_site_table_Emilyedit.csv"),
    file.path(DATA_DIR, "PISCO", "spp_attribute_table_downloaded_9-13-22_SHSPUL.csv"),
    file.path(DATA_DIR, "MBON", "SBCMBON_kelp_forest_integrated_quad_and_swath_20231022.csv"),
    file.path(DATA_DIR, "MBON", "SBCMBON_kelp_forest_integrated_fish_20231022.csv"),
    file.path(DATA_DIR, "MBON", "KFM_Macrocystis_RawData_1984-2023.csv"),
    file.path(DATA_DIR, "MBON", "SBCMBON_kelp_forest_site_geolocation_20210120_KFM_LTER.csv"),
    file.path(DATA_DIR, "LTER", "Annual_fish_comb_20240307.csv"),
    file.path(DATA_DIR, "LTER", "Annual_Kelp_All_Years_20240305.csv"),
    file.path(DATA_DIR, "LTER", "Lobster_Abundance_All_Years_20230831.csv"),
    file.path(DATA_DIR, "LANDSAT", "MPA_Runs_new.csv"),
    file.path(DATA_DIR, "MPA", "California_Marine_Protected_Areas_[ds582].shp"),
    file.path(DATA_DIR, "ALL_sizefreq_2024.csv")
  )

  existing <- file.exists(required_files)
  list(
    all_present = all(existing),
    missing = required_files[!existing],
    n_present = sum(existing),
    n_total = length(required_files)
  )
}

# ---------------------------------------------------------------------------
# Download from Dryad
# ---------------------------------------------------------------------------

download_from_dryad <- function(doi = DRYAD_DOI, dest = DATA_DIR) {

  cat("========================================================================\n")
  cat("  Downloading raw data from Dryad\n")
  cat("  DOI:", doi, "\n")
  cat("  Destination:", dest, "\n")
  cat("========================================================================\n\n")

  # Check if data already exists
  status <- check_data_exists()
  if (status$all_present) {
    cat("All", status$n_total, "required data files already present. Skipping download.\n")
    cat("To force re-download, delete the data subdirectories and re-run.\n")
    return(invisible(TRUE))
  }

  cat("Found", status$n_present, "of", status$n_total, "required files.\n")
  if (status$n_present > 0) {
    cat("Missing files:\n")
    for (f in status$missing) cat("  -", f, "\n")
  }
  cat("\n")

  # Create subdirectories
  for (d in EXPECTED_DIRS) {
    dir.create(file.path(dest, d), showWarnings = FALSE, recursive = TRUE)
  }

  # Construct Dryad download URL
  # Dryad API: https://datadryad.org/api/v2/docs/
  # For a dataset DOI like doi:10.5061/dryad.XXXXXXX, the download URL is:
  dryad_id <- sub("^doi:", "", doi)
  download_url <- paste0("https://datadryad.org/api/v2/datasets/", dryad_id, "/download")

  zip_path <- file.path(tempdir(), "dryad_data.zip")

  cat("Downloading from:", download_url, "\n")
  cat("This may take several minutes (~300 MB compressed)...\n\n")

  tryCatch({
    download.file(download_url, zip_path, mode = "wb", quiet = FALSE)
  }, error = function(e) {
    stop(
      "Download failed: ", e$message, "\n\n",
      "If the DOI has not been published yet, download manually from:\n",
      "  https://datadryad.org/stash/dataset/", dryad_id, "\n",
      "and extract into: ", dest, "\n",
      call. = FALSE
    )
  })

  cat("\nExtracting archive...\n")
  tryCatch(
    unzip(zip_path, exdir = dest, overwrite = TRUE),
    error = function(e) stop("Failed to extract archive: ", e$message, call. = FALSE)
  )
  unlink(zip_path)

  # Copy harmonized CSVs from the Dryad archive's harmonized/ subdirectory

  # into data/harmonized/ (which is tracked in git for the analysis repo)
  harmonized_src <- file.path(dest, "harmonized")
  harmonized_dest <- here::here("data", "harmonized")
  if (dir.exists(harmonized_src) && normalizePath(harmonized_src) != normalizePath(harmonized_dest, mustWork = FALSE)) {
    csv_files <- list.files(harmonized_src, pattern = "\\.csv$", full.names = TRUE)
    if (length(csv_files) > 0) {
      dir.create(harmonized_dest, showWarnings = FALSE, recursive = TRUE)
      file.copy(csv_files, harmonized_dest, overwrite = TRUE)
      cat("Extracted", length(csv_files), "harmonized CSVs to data/harmonized/\n")
    }
  }

  # Verify
  status <- check_data_exists()
  if (status$all_present) {
    cat("\nAll", status$n_total, "data files successfully downloaded and verified.\n")
  } else {
    warning("Download completed but some files are still missing:\n",
            paste("  -", status$missing, collapse = "\n"))
  }

  invisible(status$all_present)
}

# ---------------------------------------------------------------------------
# Verify checksums
# ---------------------------------------------------------------------------

verify_checksums <- function(checksum_file = here::here("dryad_staging", "checksums.md5")) {
  if (!file.exists(checksum_file)) {
    cat("Checksum file not found:", checksum_file, "\n")
    return(invisible(FALSE))
  }

  cat("Verifying data file checksums...\n")
  lines <- readLines(checksum_file)
  lines <- lines[nzchar(trimws(lines))]

  n_ok <- 0
  n_fail <- 0

  for (line in lines) {
    parts <- strsplit(trimws(line), "\\s+", perl = TRUE)[[1]]
    expected_md5 <- parts[1]
    rel_path <- parts[2]
    full_path <- file.path(DATA_DIR, rel_path)

    if (!file.exists(full_path)) {
      cat("  MISSING:", rel_path, "\n")
      n_fail <- n_fail + 1
      next
    }

    actual_md5 <- tools::md5sum(full_path)
    if (unname(actual_md5) == expected_md5) {
      n_ok <- n_ok + 1
    } else {
      cat("  MISMATCH:", rel_path, "\n")
      cat("    Expected:", expected_md5, "\n")
      cat("    Actual:  ", actual_md5, "\n")
      n_fail <- n_fail + 1
    }
  }

  cat(sprintf("Checksum verification: %d OK, %d failed (of %d total)\n",
              n_ok, n_fail, n_ok + n_fail))
  invisible(n_fail == 0)
}

# ---------------------------------------------------------------------------
# Main: run when sourced
# ---------------------------------------------------------------------------

status <- check_data_exists()

if (status$all_present) {
  cat("Data setup OK: all", status$n_total, "required files present.\n")
  cat("Run verify_checksums() to validate file integrity.\n")
} else {
  cat("\n")
  cat("========================================================================\n")
  cat("  Raw data not found (", status$n_present, "/", status$n_total, " files present)\n")
  cat("========================================================================\n\n")
  cat("Options:\n")
  cat("  1. Download from Dryad (raw + harmonized data):\n")
  cat("       download_from_dryad()\n\n")
  cat("  2. Symlink from Google Drive (collaborators with shared access):\n")
  cat("       See README.md 'Data Setup' section for symlink commands\n\n")
  cat("  3. Run the analysis pipeline (harmonized CSVs are tracked in git):\n")
  cat("       source(here::here('code', 'R', 'run_all.R'))\n")
  cat("       (No raw data needed — uses data/harmonized/ CSVs)\n\n")
  cat("  4. Regenerate figures only:\n")
  cat("       source(here::here('code', 'R', 'run_figures_only.R'))\n")
  cat("       (Uses cached intermediate results — no raw data needed)\n\n")

  if (status$n_present > 0) {
    cat("Missing files:\n")
    for (f in status$missing) cat("  -", basename(f), "\n")
  }
}
