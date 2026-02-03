# =============================================================================
# 03_data_import.R
# =============================================================================
#
# PURPOSE:
#   Import raw data files and prepare foundational datasets for the MPA kelp
#   forest pBACIPS analysis. This script loads and cleans the data that all
#   subsequent processing modules depend on.
#
# WHAT THIS SCRIPT DOES:
#   1. Imports the combined size frequency dataset (ALL_sizefreq_2024.csv)
#   2. Filters to target taxa and removes problematic sites
#   3. Standardizes species codes across programs (PISCO, KFM, LTER)
#   4. Imports MPA feature metadata (area, type, location)
#   5. Imports and merges site lists with MPA assignments
#   6. Imports MBON site geolocation data
#
# INPUTS (from data/ directory):
#   - ALL_sizefreq_2024.csv: Size frequency data from all programs
#   - MPAfeatures_subset.csv: MPA characteristics (area, type, etc.)
#   - Site_List_All.csv: Master site list with MPA assignments
#   - MBON/SBCMBON_kelp_forest_site_geolocation_*.csv: Site coordinates
#
# OUTPUTS (objects created in R environment):
#   - SizeFreq: Cleaned size frequency data for biomass estimation
#   - Site: Site list merged with MPA features
#   - Site.size: MPA features (area, type, etc.)
#   - sites.short.edit: Subset of Site with key MPA attributes
#   - Sites2: MBON site geolocation data (KFM + LTER)
#
# DEPENDENCIES:
#   Requires 00_libraries.R and 01_utils.R to be sourced first
#
# AUTHORS: Emily Donham & Adrian Stier
# PROJECT: CA MPA Kelp Forest pBACIPS Analysis
# =============================================================================

cat("Importing raw data files...\n")


# =============================================================================
# SECTION 1: SIZE FREQUENCY DATA
# =============================================================================
# Size frequency data records the size distribution of organisms at each site.
# This is essential for converting counts to biomass via bootstrap resampling.
# The data combines records from multiple monitoring programs:
#   - PISCO (UCSB and VRG campuses)
#   - KFM (NPS Kelp Forest Monitoring)
#   - LTER (Santa Barbara Coastal Long-Term Ecological Research)

# Import the combined size frequency file
SizeFreq <- read.csv(here::here("data", "ALL_sizefreq_2024.csv"))

cat("  Loaded size frequency data:", nrow(SizeFreq), "rows\n")

# -----------------------------------------------------------------------------
# Filter to target taxa across all programs
# -----------------------------------------------------------------------------
# Each program uses different species codes:
#   PISCO (VRG/UCSB): MESFRA, STRPUR, PANINT
#   KFM: Full scientific names (Strongylocentrotus franciscanus, etc.)
#   LTER: SFL, SFS, SPS, SPL (Size class codes for urchins)
#
# We subset to records that match our target taxa in each program's format

SizeFreq.south <- subset(SizeFreq,
  # Red urchin (Mesocentrotus franciscanus)
  (campus == "VRG"  & classcode == "MESFRA") |
  (campus == "UCSB" & classcode == "MESFRA") |
  (campus == "KFM"  & classcode == "Strongylocentrotus franciscanus") |

  # Purple urchin (Strongylocentrotus purpuratus)
  (campus == "VRG"  & classcode == "STRPUR") |
  (campus == "UCSB" & classcode == "STRPUR") |
  (campus == "KFM"  & classcode == "Strongylocentrotus purpuratus") |

  # Spiny lobster (Panulirus interruptus)
  (campus == "UCSB" & classcode == "PANINT") |
  (campus == "VRG"  & classcode == "PANINT") |

  # LTER urchin size classes
  # SFL/SFS = Strongylocentrotus franciscanus Large/Small
  # SPL/SPS = Strongylocentrotus purpuratus Large/Small
  (campus == "LTER" & classcode == "SFL") |
  (campus == "LTER" & classcode == "SFS") |
  (campus == "LTER" & classcode == "SPS") |
  (campus == "LTER" & classcode == "SPL")
)

cat("  Filtered to target taxa:", nrow(SizeFreq.south), "rows\n")

# -----------------------------------------------------------------------------
# Remove problematic KFM sites
# -----------------------------------------------------------------------------
# These sites are excluded because they:
#   - Are not proper control/impact pairs
#   - Have short time series where longer ones exist at the same location
# The exclusion list is defined in 01_utils.R (EXCLUDED_KFM_SITES)

SizeFreq.south <- SizeFreq.south[!(SizeFreq.south$site %in% EXCLUDED_KFM_SITES), ]

# -----------------------------------------------------------------------------
# Remove poor reference sites
# -----------------------------------------------------------------------------
# Based on consultation with Scott Hamilton (PISCO), these reference sites
# don't provide good comparisons due to habitat or sampling issues
# The exclusion list is defined in 01_utils.R (EXCLUDED_REFERENCE_SITES)

SizeFreq.south <- SizeFreq.south[!(SizeFreq.south$site %in% EXCLUDED_REFERENCE_SITES), ]

cat("  After site exclusions:", nrow(SizeFreq.south), "rows\n")

# -----------------------------------------------------------------------------
# Standardize species codes
# -----------------------------------------------------------------------------
# Convert all codes to PISCO short format for consistency:
#   STRPUR = Strongylocentrotus purpuratus (purple urchin)
#   MESFRA = Mesocentrotus franciscanus (red urchin)
#   PANINT = Panulirus interruptus (spiny lobster)

SizeFreq <- mutate(SizeFreq.south, classcode = case_when(
  # KFM full names -> PISCO codes
  classcode == "Strongylocentrotus purpuratus"    ~ "STRPUR",
  classcode == "Strongylocentrotus franciscanus"  ~ "MESFRA",

  # LTER size class codes -> PISCO codes
  # Note: KFM incorrectly has "Strongylocentrotus franciscanus" for red urchin
  # The species was reclassified to Mesocentrotus franciscanus
  classcode == "SFL" ~ "MESFRA",  # franciscanus Large
  classcode == "SFS" ~ "MESFRA",  # franciscanus Small
  classcode == "SPS" ~ "STRPUR",  # purpuratus Small
  classcode == "SPL" ~ "STRPUR",  # purpuratus Large

  # Keep all other codes as-is (PISCO codes already correct)
  TRUE ~ classcode
))

# Clean up intermediate object
rm(SizeFreq.south)

cat("  Species codes standardized.\n")


# =============================================================================
# SECTION 2: SITE METADATA
# =============================================================================
# Load MPA characteristics and site assignments

# -----------------------------------------------------------------------------
# MPA features (area, type, protection level)
# -----------------------------------------------------------------------------

Site.size <- read.csv(here::here("data", "MPAfeatures_subset.csv"))

cat("  Loaded MPA features:", nrow(Site.size), "MPAs\n")

# -----------------------------------------------------------------------------
# Master site list with MPA assignments
# -----------------------------------------------------------------------------
# This file links monitoring sites to their associated MPAs

Site <- read.csv(here::here("data", "Site_List_All.csv"))

cat("  Loaded site list:", nrow(Site), "sites\n")

# -----------------------------------------------------------------------------
# Merge site list with MPA features
# -----------------------------------------------------------------------------
# This gives us MPA characteristics (area, type) for each monitoring site

Site <- merge(Site, Site.size,
  by.x = "Site",        # Column in Site
  by.y = "NAME",        # Column in Site.size
  all = TRUE            # Keep all rows from both tables
)

# Remove rows that don't have valid location data
# (these are MPAs without monitoring sites or vice versa)
Site <- Site[!is.na(Site$Lat), ]

cat("  Merged site-MPA data:", nrow(Site), "rows\n")

# -----------------------------------------------------------------------------
# Create reduced site table for downstream joins
# -----------------------------------------------------------------------------
# This subset contains just the columns needed for merging with processed data

sites.short.edit <- subset(Site, select = c(
  CA_MPA_Name_Short,  # MPA name (short version)
  type,               # MPA type (SMR, SMCA, etc.)
  Location,           # Geographic region
  Hectares            # MPA area
))


# =============================================================================
# SECTION 3: MBON SITE GEOLOCATION
# =============================================================================
# The Marine Biodiversity Observation Network (MBON) synthesis includes
# harmonized site location data for KFM and LTER monitoring programs

Sites2 <- read.csv(here::here("data", "MBON",
                               "SBCMBON_kelp_forest_site_geolocation_20210120_KFM_LTER.csv"))

# Convert site_id to factor for consistent handling
Sites2$site_id <- as.factor(Sites2$site_id)

cat("  Loaded MBON site geolocation:", nrow(Sites2), "sites\n")


# =============================================================================
# SUMMARY
# =============================================================================

cat("\nData import complete. Objects created:\n")
cat("  - SizeFreq:", nrow(SizeFreq), "rows (size frequency for biomass estimation)\n")
cat("  - Site:", nrow(Site), "rows (site list with MPA features)\n")
cat("  - Site.size:", nrow(Site.size), "rows (MPA features)\n")
cat("  - sites.short.edit:", nrow(sites.short.edit), "rows (key MPA attributes)\n")
cat("  - Sites2:", nrow(Sites2), "rows (MBON site coordinates)\n")
