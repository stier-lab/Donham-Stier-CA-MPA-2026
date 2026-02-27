# =============================================================================
# 03_load_harmonized_data.R
# =============================================================================
#
# PURPOSE:
#   Load harmonized CSV data produced by the data-processing repo.
#   Replaces scripts 03-07 (data import, PISCO/KFM/LTER/Landsat processing,
#   and data combination).
#
# DATA SOURCES:
#   Looks for harmonized CSVs in this order:
#   1. data/harmonized/ (tracked in git, preferred)
#   2. Sibling data repo output/harmonized/ (fallback for development)
#
# OUTPUTS:
#   Creates the same objects that scripts 03-07 previously produced:
#   - All.RR.sub.trans: Combined response ratios (3,368 rows x 14 cols)
#   - All.Resp.sub: Combined raw responses (6,244 rows x 9 cols)
#   - Landsat.RR: Satellite kelp response ratios (719 rows x 14 cols)
#   - Site: MPA site metadata (34 rows x 9 cols)
#
# AUTHORS: Emily Donham & Adrian Stier
# PROJECT: CA MPA Kelp Forest pBACIPS Analysis
# =============================================================================

####################################################################################################
## Locate harmonized data directory ################################################################
####################################################################################################

harmonized_dir <- here::here("data", "harmonized")

if (!dir.exists(harmonized_dir)) {
  # Fallback: check sibling data repo
  sibling <- file.path(dirname(here::here()), "Donham-Stier-CA-MPA-Data-2026", "output", "harmonized")
  if (dir.exists(sibling)) {
    harmonized_dir <- sibling
    cat("  Using harmonized data from sibling repo:", sibling, "\n")
  } else {
    stop("Harmonized data not found.\n",
         "  Expected: ", here::here("data", "harmonized"), "\n",
         "  Or sibling: ", sibling, "\n",
         "  See README.md 'Data Setup' for instructions.", call. = FALSE)
  }
}

####################################################################################################
## Load harmonized CSVs ############################################################################
####################################################################################################

cat("Loading harmonized data from:", harmonized_dir, "\n")

All.RR.sub.trans <- read.csv(file.path(harmonized_dir, "harmonized_response_ratios.csv"),
                              stringsAsFactors = FALSE)
All.Resp.sub <- read.csv(file.path(harmonized_dir, "harmonized_raw_responses.csv"),
                          stringsAsFactors = FALSE)
Site <- read.csv(file.path(harmonized_dir, "harmonized_site_metadata.csv"),
                  stringsAsFactors = FALSE)

# Landsat.RR is consumed directly by 08_effect_sizes.R (bypasses main RR pipeline)
landsat_path <- file.path(harmonized_dir, "harmonized_landsat_rr.csv")
if (!file.exists(landsat_path)) {
  stop("harmonized_landsat_rr.csv not found in ", harmonized_dir,
       "\n  This file is required by 08_effect_sizes.R for Landsat kelp analysis.",
       call. = FALSE)
}
Landsat.RR <- read.csv(landsat_path, stringsAsFactors = FALSE)
Landsat.RR$year <- as.integer(Landsat.RR$year)
Landsat.RR$time <- as.integer(Landsat.RR$time)

####################################################################################################
## Restore types that CSV round-trip may lose ######################################################
####################################################################################################

All.RR.sub.trans$year <- as.integer(All.RR.sub.trans$year)
All.RR.sub.trans$time <- as.integer(All.RR.sub.trans$time)

All.Resp.sub$year <- as.integer(All.Resp.sub$year)
All.Resp.sub$time <- as.integer(All.Resp.sub$time)

Site$MPA_Start <- as.integer(Site$MPA_Start)
Site$Quarter <- as.integer(Site$Quarter)

####################################################################################################
## Schema validation ###############################################################################
####################################################################################################

# Verify expected columns exist
stopifnot("All.RR.sub.trans schema" =
  all(c("CA_MPA_Name_Short", "year", "y", "lnDiff", "mpa", "reference", "Diff",
        "resp", "time", "type", "Location", "Hectares", "source", "BA") %in% names(All.RR.sub.trans)))

stopifnot("All.Resp.sub schema" =
  all(c("CA_MPA_Name_Short", "year", "taxon_name", "source", "status", "value",
        "resp", "BA", "time") %in% names(All.Resp.sub)))

stopifnot("Site schema" =
  all(c("CA_MPA_Name_Short", "MPA_Start", "Quarter", "Lat", "Lon", "type", "Location") %in% names(Site)))

stopifnot("Landsat.RR schema" =
  all(c("CA_MPA_Name_Short", "year", "lnDiff", "BA", "time") %in% names(Landsat.RR)))

####################################################################################################
## Summary #########################################################################################
####################################################################################################

cat("Loaded harmonized data:\n")
cat("  All.RR.sub.trans:", nrow(All.RR.sub.trans), "rows x", ncol(All.RR.sub.trans), "cols\n")
cat("  All.Resp.sub:    ", nrow(All.Resp.sub), "rows x", ncol(All.Resp.sub), "cols\n")
cat("  Site:            ", nrow(Site), "rows x", ncol(Site), "cols\n")
cat("  Landsat.RR:      ", nrow(Landsat.RR), "rows x", ncol(Landsat.RR), "cols\n")
