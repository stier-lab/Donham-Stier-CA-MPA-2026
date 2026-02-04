# =============================================================================
# 06b_landsat_processing.R
# =============================================================================
#
# PURPOSE:
#   Process satellite-derived kelp canopy biomass data from Landsat imagery
#   for the MPA pBACIPS analysis.
#
# WHAT THIS SCRIPT DOES:
#   1. Imports Landsat kelp canopy area/biomass data
#   2. Transforms data to long format for analysis
#   3. Calculates annual mean biomass by MPA and status (MPA vs reference)
#   4. Converts to proportions of time series maximum
#   5. Calculates log response ratios
#   6. Assigns time since MPA and Before/After labels
#
# DATA SOURCE:
#   Landsat-derived kelp canopy area from remote sensing analysis.
#   This provides an independent measure of kelp abundance that:
#   - Covers a broader spatial extent than diver surveys
#   - Has consistent methodology back to 1984
#   - Measures canopy area (proxy for biomass) from above
#
# ADVANTAGES OF SATELLITE DATA:
#   - No diver bias or sampling error
#   - Continuous time series without gaps
#   - Can detect large-scale patterns
#
# LIMITATIONS:
#   - Only measures surface canopy (misses understory kelp)
#   - Affected by water clarity and atmospheric conditions
#   - Lower taxonomic resolution (just "kelp", not species)
#
# INPUTS:
#   - data/LANDSAT/MPA_Runs_new.csv
#   - Site object (from 03_data_import.R)
#   - sites.short.edit object
#
# OUTPUTS:
#   - Landsat.RR: Kelp biomass response ratios with metadata
#   - Landsat.resp: Raw kelp biomass data (long format)
#
# DEPENDENCIES:
#   Requires 00-03 scripts to be sourced first
#
# AUTHORS: Emily Donham & Adrian Stier
# PROJECT: CA MPA Kelp Forest pBACIPS Analysis
# =============================================================================

# Source utility functions
source(here::here("code", "R", "01_utils.R"))

####################################################################################################
## 1. Import Landsat kelp biomass data #############################################################
####################################################################################################
# Landsat satellite imagery allows us to estimate kelp canopy area remotely.
# This provides an independent measure of kelp abundance that:
# - Is not affected by diver sampling bias
# - Has consistent methodology back to 1984
# - Covers a broader spatial extent than diver surveys

landsat_data_path <- here::here("data", "LANDSAT", "MPA_Runs_new.csv")

# Check if file exists before trying to read it
# file.exists() returns TRUE/FALSE; this makes the script robust to missing data
if (!file.exists(landsat_data_path)) {
  warning("Landsat data file not found: ", landsat_data_path,
          "\nLandsat.RR and Landsat.resp will not be created.")
} else {

  fn <- read.csv(landsat_data_path)

  ####################################################################################################
  ## 2. Transform to long format #####################################################################
  ####################################################################################################
  # The raw data is in "wide" format with one column per year.
  # We need "long" format (one row per observation) for analysis.
  # gather() from tidyr converts wide to long: creates Year and Biomass columns.
  # The minus signs mean "don't gather these columns" (keep them as-is).

  bio.summarize.long <- gather(fn, Year, Biomass, -MPA, -status, -rep, -lat, -lon)

  # Dynamically extract years from column names and calculate rows per year group
  year_cols <- grep("^X[0-9]{4}$", names(fn), value = TRUE)
  years <- as.numeric(gsub("X", "", year_cols))
  n_rows_per_year <- nrow(fn)  # Each row in fn represents one MPA x status x rep combination

  # Validate dimensions before assignment
  expected_total_rows <- length(years) * n_rows_per_year
  if (nrow(bio.summarize.long) != expected_total_rows) {
    warning("Landsat data dimension mismatch: expected ", expected_total_rows,
            " rows but got ", nrow(bio.summarize.long))
  }

  # Assign years dynamically based on actual data structure
  bio.summarize.long$year <- rep(years, each = n_rows_per_year)
  bio.summarize.long <- bio.summarize.long %>% dplyr::arrange(MPA)  # Sort by MPA name

  ####################################################################################################
  ## 3. Calculate annual mean biomass by MPA and status ##############################################
  ####################################################################################################

  bio.ave.mpa <- bio.summarize.long %>%
    dplyr::group_by(MPA, status, year) %>%
    dplyr::summarise_at(c("Biomass"), mean, na.rm = TRUE) %>%
    dplyr::ungroup()

  # Spread to wide format (one column per status) and remove incomplete cases
  bio.sum.max.short <- bio.ave.mpa %>% spread(status, Biomass)
  bio.sum.max.short <- bio.sum.max.short[complete.cases(bio.sum.max.short), ]

  ####################################################################################################
  ## 4. Convert to proportions of maximum ############################################################
  ####################################################################################################

  bio.sum.max.long <- gather(bio.sum.max.short, status, meanBio, -MPA, -year)

  # Add required columns for calculate_proportions() utility function
  bio.sum.max.long$CA_MPA_Name_Short <- bio.sum.max.long$MPA
  bio.sum.max.long$y <- "Macrocystis pyrifera"

  # Use the standardized calculate_proportions function from 01_utils.R
  # This ensures consistent proportion calculation and zero-correction across all data sources
  bio.sum.max.long <- calculate_proportions(
    bio.sum.max.long,
    value_col = "meanBio",
    correction_method = "adaptive"  # Use adaptive instead of fixed 0.01 for consistency
  )

  ####################################################################################################
  ## 5. Calculate log response ratios (MPA vs reference) #############################################
  ####################################################################################################

  bio.sum.max.long$taxon_name <- "Macrocystis pyrifera"
  # Select columns needed for spread: MPA, year, taxon_name, status, PropCorr
  bio.sum.max.long.sub <- bio.sum.max.long[, c("MPA", "year", "taxon_name", "status", "PropCorr")]
  bio.sum.max.LANDSAT.diff <- bio.sum.max.long.sub %>%
    spread(status, PropCorr)
  bio.sum.max.LANDSAT.diff <- calculate_log_response_ratio(bio.sum.max.LANDSAT.diff)

  ####################################################################################################
  ## 6. Assign Before/After and time since MPA implementation ########################################
  ####################################################################################################

  bio.sum.max.LANDSAT.diff <- assign_ba_from_site_table(
    data.frame(CA_MPA_Name_Short = bio.sum.max.LANDSAT.diff$MPA,
               year = bio.sum.max.LANDSAT.diff$year,
               lnDiff = bio.sum.max.LANDSAT.diff$lnDiff),
    Site
  )
  bio.sum.max.LANDSAT.diff <- assign_time_from_site_table(bio.sum.max.LANDSAT.diff, Site)

  ####################################################################################################
  ## 7. Build Landsat.RR output (matching format of LTER/PISCO/KFM RR datasets) #####################
  ####################################################################################################

  # Reconstruct full response ratio dataframe with columns matching other modules
  Landsat.RR <- bio.sum.max.LANDSAT.diff
  Landsat.RR$y <- "Macrocystis pyrifera"
  Landsat.RR$resp <- "Bio"
  Landsat.RR$source <- "Landsat"

  # Merge with site metadata for type, Location, Hectares
  Landsat.RR <- merge(Landsat.RR, sites.short.edit,
                       by.x = "CA_MPA_Name_Short",
                       by.y = "CA_MPA_Name_Short",
                       all.x = TRUE)

  ####################################################################################################
  ## 8. Build Landsat.resp output (raw response data in long format) #################################
  ####################################################################################################

  # Use the wide-format biomass (mpa vs reference columns) before log-ratio transformation
  Landsat.resp.wide <- bio.sum.max.short
  colnames(Landsat.resp.wide)[colnames(Landsat.resp.wide) == "MPA"] <- "CA_MPA_Name_Short"
  Landsat.resp.wide$taxon_name <- "Macrocystis pyrifera"
  Landsat.resp.wide$source <- "Landsat"

  Landsat.resp <- gather(Landsat.resp.wide, status, value, -CA_MPA_Name_Short, -year,
                         -taxon_name, -source)
  Landsat.resp$resp <- "Bio"

  ####################################################################################################
  ## Clean up intermediate objects ###################################################################
  ####################################################################################################

  rm(fn, bio.summarize.long, bio.ave.mpa, bio.sum.max.short,
     bio.sum.max.long, bio.sum.max.LANDSAT.diff, Landsat.resp.wide)

}

###########################################################################################
## Summary of outputs available for downstream scripts
###########################################################################################
# Landsat.RR   - Kelp biomass log response ratios with time, BA, site info
#                Columns: CA_MPA_Name_Short, year, lnDiff, BA, time, y, resp, source,
#                         type, Location, Hectares
# Landsat.resp - Raw kelp biomass response data (mpa vs reference, long format)
#                Columns: CA_MPA_Name_Short, year, taxon_name, source, status, value, resp
