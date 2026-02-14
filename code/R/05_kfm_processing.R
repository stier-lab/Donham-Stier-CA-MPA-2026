# =============================================================================
# 05_kfm_processing.R
# =============================================================================
#
# PURPOSE:
#   Process KFM (Kelp Forest Monitoring) data from the National Park Service
#   Channel Islands monitoring program for the MPA pBACIPS analysis.
#
# WHAT THIS SCRIPT DOES:
#   1. Imports MBON-integrated quad and swath data for KFM sites
#   2. Processes urchin data (S. purpuratus, M. franciscanus)
#   3. Processes Macrocystis kelp stipe data
#   4. Processes sheephead fish data (visual fish surveys)
#   5. Calculates biomass via bootstrap resampling
#   6. Calculates log response ratios (MPA vs reference sites)
#   7. Assigns time since MPA and Before/After labels
#
# DATA SOURCES:
#   KFM = NPS Kelp Forest Monitoring, Channel Islands National Park
#   Data comes via the MBON (Marine Biodiversity Observation Network) synthesis
#   which integrates multiple long-term monitoring programs.
#
# KEY DIFFERENCES FROM PISCO:
#   - KFM uses different transect sizes (quads are 2m2 vs PISCO 60m2)
#   - KFM has longer time series (since 1984 for some sites)
#   - All KFM sites are in the Channel Islands
#   - Sheephead surveyed via visual fish counts and roving diver fish counts
#
# INPUTS:
#   - data/MBON/SBCMBON_kelp_forest_integrated_quad_and_swath_20231022.csv
#   - data/MBON/SBCMBON_kelp_forest_integrated_fish_20231022.csv
#   - data/MBON/KFM_Macrocystis_RawData_1984-2023.csv
#   - data/MBON/SBCMBON_kelp_forest_site_geolocation_*.csv
#   - SizeFreq.Urch.OG (from 04_pisco_processing.R - unfiltered urchin sizes)
#
# OUTPUTS:
#   - KFM.join.ave: Response ratios for urchins and kelp
#   - KFM.resp: Raw density/biomass data (MPA vs reference)
#   - kfm.fish: Sheephead response ratios
#   - KFM.fish.den.long: Raw fish density data
#   - lter: LTER subset of MBON (passed to 06_lter_processing.R)
#
# DEPENDENCIES:
#   Requires 00-04 scripts to be sourced first
#
# AUTHORS: Emily Donham & Adrian Stier
# PROJECT: CA MPA Kelp Forest pBACIPS Analysis
# =============================================================================

cat("Processing KFM/MBON data...\n")

###########################################################################################
## 1. Import MBON quad/swath data and site geolocation
###########################################################################################
# MBON (Marine Biodiversity Observation Network) integrates data from multiple programs.
# This file contains both KFM (Kelp Forest Monitoring - Channel Islands NPS) and
# LTER (Long-Term Ecological Research - mainland) data.
# We'll split them apart and process KFM here, LTER in 06_lter_processing.R.

MBON <- read.csv(here::here("data", "MBON", "SBCMBON_kelp_forest_integrated_quad_and_swath_20231022.csv"))

# Sites2 contains site coordinates and MPA information
# This table links site IDs to MPA names and geographic metadata
Sites2 <- read.csv(here::here("data", "MBON", "SBCMBON_kelp_forest_site_geolocation_20210120_KFM_LTER.csv"))
Sites2$site_id <- as.factor(Sites2$site_id)  # Convert to factor for proper merging

# Treat FORCE_BOOTSTRAP consistently across scripts: only force recompute when TRUE.
force_boot <- exists("FORCE_BOOTSTRAP", envir = .GlobalEnv) && isTRUE(get("FORCE_BOOTSTRAP", envir = .GlobalEnv))

###########################################################################################
## 2. Process KFM subset: type conversions, date formatting, missing values
###########################################################################################
# Data imported from CSV often needs type conversions.
# R may read text as "character" strings when we want "factor" (categorical variables)
# or numeric. This section standardizes column types.

# Factor conversions - factors are R's way of representing categorical variables
# They're more memory-efficient and have defined levels for grouping/plotting
factor_cols <- c("data_source", "sample_method", "site_id", "subsite_id",
                 "transect_id", "replicate_id", "proj_taxon_id",
                 "auth_taxon_id", "auth_name", "taxon_name",
                 "site_name", "subsite_name")

# Loop through columns and convert to factor if not already
# [[col]] accesses a column by name (like $col but works with variables)
for (col in factor_cols) {
  if (!is.factor(MBON[[col]])) MBON[[col]] <- as.factor(MBON[[col]])
}

# Numeric conversions
numeric_cols <- c("area", "count", "latitude", "longitude")
for (col in numeric_cols) {
  if (is.factor(MBON[[col]])) {
    MBON[[col]] <- as.numeric(levels(MBON[[col]]))[as.integer(MBON[[col]])]
  } else if (is.character(MBON[[col]])) {
    MBON[[col]] <- as.numeric(MBON[[col]])
  }
}

# Date conversion
tmpDateFormat <- "%Y-%m-%d"
tmp1date <- as.Date(MBON$date, format = tmpDateFormat)
if (length(tmp1date) == length(tmp1date[!is.na(tmp1date)])) {
  MBON$date <- tmp1date
} else {
  print("Date conversion failed for MBON$date. Please inspect the data and do the date conversion yourself.")
}
rm(tmpDateFormat, tmp1date)

# Convert missing value placeholders to NA
MBON$count <- ifelse(trimws(as.character(MBON$count)) == trimws("."), NA, MBON$count)
suppressWarnings(MBON$count <- ifelse(!is.na(as.numeric(".")) & (trimws(as.character(MBON$count)) == as.character(as.numeric("."))), NA, MBON$count))
MBON$auth_taxon_id <- as.factor(ifelse(trimws(as.character(MBON$auth_taxon_id)) == trimws("."), NA, as.character(MBON$auth_taxon_id)))
MBON$auth_name <- as.factor(ifelse(trimws(as.character(MBON$auth_name)) == trimws("."), NA, as.character(MBON$auth_name)))
MBON$subsite_name <- as.factor(ifelse(trimws(as.character(MBON$subsite_name)) == trimws("."), NA, as.character(MBON$subsite_name)))
MBON$latitude <- ifelse(trimws(as.character(MBON$latitude)) == trimws("."), NA, MBON$latitude)
suppressWarnings(MBON$latitude <- ifelse(!is.na(as.numeric(".")) & (trimws(as.character(MBON$latitude)) == as.character(as.numeric("."))), NA, MBON$latitude))
MBON$longitude <- ifelse(trimws(as.character(MBON$longitude)) == trimws("."), NA, MBON$longitude)
suppressWarnings(MBON$longitude <- ifelse(!is.na(as.numeric(".")) & (trimws(as.character(MBON$longitude)) == as.character(as.numeric("."))), NA, MBON$longitude))

# Split into KFM and LTER subsets based on data_source column
# KFM data is processed here; LTER data is passed to 06_lter_processing.R
kfm <- subset(MBON, data_source == "kfm")
lter <- subset(MBON, data_source == "lter")  # Used by 06_lter_processing.R

###########################################################################################
## 3. Remove excluded KFM sites and join to site table
###########################################################################################
# merge() joins the survey data with site metadata using site_id as the key
# This adds columns like CA_MPA_Name_Short, status (mpa/reference), etc.
kfm.site <- merge(kfm, Sites2, by = c("site_id"))

# Extract year from date column using lubridate's year() function
kfm.site$year <- year(kfm.site$date)

# Remove empty MPA names (sites not associated with any MPA)
kfm.site <- subset(kfm.site, CA_MPA_Name_Short != "")

# Remove excluded KFM sites using constant from 01_utils.R
# EXCLUDED_KFM_SITES contains site IDs that don't meet data quality standards
# The %in% operator checks membership; ! negates to keep sites NOT in the list
kfm.site <- kfm.site[!(kfm.site$site_id %in% EXCLUDED_KFM_SITES), ]

###########################################################################################
## 4. Assign time since MPA
###########################################################################################

kfm.site <- assign_time_from_site_table(kfm.site, Sites2)

###########################################################################################
## 5. Subset target taxa
###########################################################################################

kfm.site <- subset(kfm.site, taxon_name == "Strongylocentrotus purpuratus" |
                     taxon_name == "Mesocentrotus franciscanus" |
                     taxon_name == "Macrocystis pyrifera" |
                     taxon_name == "Panulirus interruptus")

# Remove juvenile giant kelp (not done by other programs, not comparable)
kfm.site <- subset(kfm.site, taxon_name != "Macrocystis pyrifera" | proj_taxon_id != "t-k-002")

# Sum across different Macro names (considering all together)
kfm.sum <- kfm.site %>%
  dplyr::group_by(site_id, taxon_name, status, CA_MPA_Name_Short, area, ChannelIsland, replicate_id, MPA_Start, time, year) %>%
  dplyr::summarise_at(c("count"), sum) %>%
  dplyr::ungroup()

# Transform to density since there are different survey methods
# Use safe_divide to handle any zero area values
kfm.sum$den <- safe_divide(kfm.sum$count, kfm.sum$area, context = "KFM density calculation")

# Average across sites and replicates
kfm.ave <- kfm.sum %>%
  dplyr::group_by(site_id, taxon_name, status, CA_MPA_Name_Short, area, ChannelIsland, MPA_Start, time, year) %>%
  dplyr::summarise_at(c("den"), mean) %>%
  dplyr::ungroup()

###########################################################################################
## 6. Calculate urchin biomass via bootstrap resampling
###########################################################################################

URCHINS <- subset(kfm.site, taxon_name == "Strongylocentrotus purpuratus" |
                    taxon_name == "Mesocentrotus franciscanus")
URCHINS <- dplyr::mutate(URCHINS, taxon_name = dplyr::case_when(
  taxon_name == "Strongylocentrotus purpuratus" ~ "STRPURAD",
  taxon_name == "Mesocentrotus franciscanus" ~ "MESFRAAD",
  TRUE ~ taxon_name))

u <- unique(URCHINS[, c("year", "site_id", "CA_MPA_Name_Short", "status", "area", "taxon_name")])

# Reset size frequency data to original (before PISCO 25mm filter)
SizeFreq.Urch <- SizeFreq.Urch.OG

# Cache bootstrap results to avoid re-running the slow loop
.cache_kfm_urch <- here::here("data", "cache", "kfm_urchin_bootstrap.rds")
if (file.exists(.cache_kfm_urch) && !force_boot) {
  cat("    Loading cached KFM urchin bootstrap results...\n")
  Urchin.site <- safe_read_rds(.cache_kfm_urch)
  if (is.null(Urchin.site)) {
    cat("    Cache invalid, will recompute...\n")
  } else {
    # Cache schema can evolve; invalidate cache if required columns are missing.
    required_urch_cols <- c(
      "site", "CA_MPA_Name_Short", "site_status",
      "year", "area", "transect",
      "count", "biomass", "y"
    )
    missing_urch_cols <- setdiff(required_urch_cols, names(Urchin.site))
    if (length(missing_urch_cols) > 0) {
      cat("    Cached KFM urchin bootstrap missing columns (",
          paste(missing_urch_cols, collapse = ", "),
          "); will recompute...\n", sep = "")
      Urchin.site <- NULL
    }
  }
}

if (!exists("Urchin.site") || is.null(Urchin.site)) {
  # Pre-allocate list for results (PERFORMANCE: avoids O(n²) rbind copies)
  results_list <- vector("list", nrow(u))

  # Set random seed for reproducibility of bootstrap resampling
  set.seed(12347)
  for (i in seq_len(nrow(u))) {
    t <- which(URCHINS$site_id == u$site_id[i] &
                 URCHINS$CA_MPA_Name_Short == u$CA_MPA_Name_Short[i] &
                 URCHINS$year == u$year[i] &
                 URCHINS$taxon_name == u$taxon_name[i])
    n <- sum(URCHINS$count[t], na.rm = TRUE)
    t2 <- which(SizeFreq.Urch$CA_MPA_Name_Short == u$CA_MPA_Name_Short[i] &
                  SizeFreq.Urch$site_status == u$status[i] &
                  SizeFreq.Urch$year == u$year[i] &
                  SizeFreq.Urch$classcode == u$taxon_name[i])

    # Select biomass function based on species
    bio_fun <- if (u$taxon_name[i] == "MESFRAAD") bio_redurch else bio_purpurch

    # Wrap bootstrap in tryCatch for robustness
    result <- tryCatch(
      bootstrap_biomass(
        count = n,
        size_freq_indices = t2,
        size_freq_table = SizeFreq.Urch,
        biomass_fun = bio_fun
      ),
      error = function(e) {
        warning("Bootstrap failed for site ", u$site_id[i], " year ", u$year[i],
                " species ", u$taxon_name[i], ": ", e$message)
        list(biomass = NA, se = NA, count = n)
      }
    )

    results_list[[i]] <- data.frame(
      site = u$site_id[i],
      CA_MPA_Name_Short = u$CA_MPA_Name_Short[i],
      site_status = u$status[i],
      year = u$year[i],
      area = u$area[i],
      transect = length(t),
      count = result$count,
      biomass = result$biomass,
      y = u$taxon_name[i]
    )
  }

  # Combine all results at once (efficient)
  Urchin.site <- do.call(rbind, results_list)

  # Save cache
  dir.create(dirname(.cache_kfm_urch), recursive = TRUE, showWarnings = FALSE)
  saveRDS(Urchin.site, .cache_kfm_urch)
  cat("    Cached KFM urchin bootstrap results.\n")
}

# Urchin.site already contains KFM/MBON site_status and MPA assignment.
# Do NOT merge with PISCO's `sites.short` (which may exist in the session and has different site IDs).
KFM.Urchin.site.merge <- Urchin.site

# Remove Old Anacapa MPA and the SMCA since they don't meet criteria for inclusion
KFM.Urchin.site.merge <- subset(
  KFM.Urchin.site.merge,
  !(CA_MPA_Name_Short %in% c("Anacapa Island SMCA", "Anacapa Island SMR 1978"))
)

# Find sites where urchin was found on transect but no size was recorded
KFM.Urchin.null <- subset(KFM.Urchin.site.merge, count > 0 & (is.na(biomass) | biomass == 0))
KFM.Urchin.site.merge.sub <- KFM.Urchin.site.merge

# Remove sites where urchins were seen but not measured (can't get accurate biomass)
for (i in seq_len(nrow(KFM.Urchin.null))) {
  KFM.Urchin.site.merge.sub <- subset(KFM.Urchin.site.merge.sub,
                                        site != KFM.Urchin.null$site[i] |
                                          year != KFM.Urchin.null$year[i])
}

# Account for different number of transects at sites
KFM.Urchin.site.merge.sub$count.ave <- safe_divide(
  KFM.Urchin.site.merge.sub$count,
  KFM.Urchin.site.merge.sub$transect,
  context = "KFM urchin count/transect"
)
KFM.Urchin.site.merge.sub$biomass.ave <- safe_divide(
  KFM.Urchin.site.merge.sub$biomass,
  KFM.Urchin.site.merge.sub$transect,
  context = "KFM urchin biomass/transect"
)

KFM.Urchin.site.all <- KFM.Urchin.site.merge.sub

# Quads are 2m2, so divide by 2 to get per m2
KFM.Urchin.site.all$Density <- KFM.Urchin.site.all$count.ave / 2
KFM.Urchin.site.all$bio.m2 <- KFM.Urchin.site.all$biomass.ave / 2

# Restore full species names
KFM.Urchin.site.all <- dplyr::mutate(KFM.Urchin.site.all, y = dplyr::case_when(
  y == "STRPURAD" ~ "Strongylocentrotus purpuratus",
  y == "MESFRAAD" ~ "Mesocentrotus franciscanus",
  TRUE ~ y))

###########################################################################################
## 7. Process KFM Macrocystis stipe data from raw stipe counts file
###########################################################################################

kfm.stipe <- read.csv(here::here("data", "MBON", "KFM_Macrocystis_RawData_1984-2023.csv"))
kfm.stipe <- kfm.stipe[!duplicated(kfm.stipe), ] # Data issues caused duplicates

# Create site_id matching site table format
kfm.stipe$site_id <- ifelse(kfm.stipe$SiteNumber > 9,
                              paste0("a-k-", kfm.stipe$SiteNumber),
                              paste0("a-k-0", kfm.stipe$SiteNumber))
kfm.stipe.site <- merge(kfm.stipe, Sites2, by = c("site_id"))
kfm.stipe.site <- kfm.stipe.site %>%
  dplyr::mutate(date = as.Date(SurveyDate, format = "%m/%d/%Y"))
kfm.stipe.site$year <- year(kfm.stipe.site$date)

# Remove excluded KFM sites
kfm.stipe.site <- kfm.stipe.site[!(kfm.stipe.site$site_id %in% EXCLUDED_KFM_SITES), ]

# Calculate average number of stipes per site
kfm.stipe.mean <- kfm.stipe.site %>%
  dplyr::group_by(site_id, status, CA_MPA_Name_Short, ChannelIsland, MPA_Start, year) %>%
  dplyr::summarise_at(c("Stipe_Count"), mean)

# Calculate average number of stipes per MPA/Reference
kfm.stipe.max <- kfm.stipe.mean %>%
  dplyr::group_by(status, CA_MPA_Name_Short, ChannelIsland, MPA_Start, year) %>%
  dplyr::summarise_at(c("Stipe_Count"), mean)

kfm.stipe.ave <- subset(kfm.stipe.max, CA_MPA_Name_Short != "")

# Assign time since MPA using utility function
kfm.stipe.ave <- assign_time_from_site_table(kfm.stipe.ave, Sites2)

###########################################################################################
## 8. Calculate macrocystis biomass via bootstrap resampling
###########################################################################################

MACRO <- subset(kfm.sum, taxon_name == "Macrocystis pyrifera")
# Remove Old Anacapa MPA and the SMCA
MACRO <- subset(MACRO, CA_MPA_Name_Short != "Anacapa Island SMCA" &
                  CA_MPA_Name_Short != "Anacapa Island SMR 1978")

u <- unique(MACRO[, c("year", "site_id", "CA_MPA_Name_Short", "status", "area", "taxon_name")])

kfm.stipe.OG <- kfm.stipe.site

# Cache bootstrap results to avoid re-running the slow loop
.cache_kfm_macro <- here::here("data", "cache", "kfm_macro_bootstrap.rds")
if (file.exists(.cache_kfm_macro) && !force_boot) {
  cat("    Loading cached KFM Macrocystis bootstrap results...\n")
  Macro.site <- safe_read_rds(.cache_kfm_macro)
  if (is.null(Macro.site)) {
    cat("    Cache invalid, will recompute...\n")
  } else {
    required_macro_cols <- c(
      "site", "CA_MPA_Name_Short", "site_status",
      "year", "transect", "area",
      "y", "count", "biomass", "ind"
    )
    missing_macro_cols <- setdiff(required_macro_cols, names(Macro.site))
    if (length(missing_macro_cols) > 0) {
      cat("    Cached KFM Macrocystis bootstrap missing columns (",
          paste(missing_macro_cols, collapse = ", "),
          "); will recompute...\n", sep = "")
      Macro.site <- NULL
    }
  }
}

if (!exists("Macro.site") || is.null(Macro.site)) {
  # Pre-allocate list for results (PERFORMANCE: avoids O(n²) rbind copies)
  results_list <- vector("list", nrow(u))

  # Set random seed for reproducibility of bootstrap resampling
  set.seed(12348)
  for (i in seq_len(nrow(u))) {
    t <- which(MACRO$site_id == u$site_id[i] &
                 MACRO$year == u$year[i] &
                 MACRO$area == u$area[i] &
                 MACRO$taxon_name == u$taxon_name[i])
    n <- sum(MACRO$count[t], na.rm = TRUE)
    t2 <- which(kfm.stipe.site$site_id == u$site_id[i] &
                  kfm.stipe.site$SurveyYear == u$year[i] &
                  kfm.stipe.site$ScientificName == u$taxon_name[i])

    # Wrap bootstrap in tryCatch for robustness
    result_row <- tryCatch({
      if (n != 0 & length(t2) != 0) {
        # Kelp found on transect and in stipe data
        a <- kfm.stipe.site[t2, ]
        s <- data.frame(matrix(NA, nrow = 1000, ncol = n))
        for (k in 1:1000) {
          s[k, ] <- sample(a$Stipe_Count, n, replace = TRUE)
        }
        s_bio <- s %>%
          dplyr::mutate(dplyr::across(starts_with("X"), ~ bio_macro(.)))
        s_bio <- s_bio %>%
          dplyr::mutate(sum = rowSums(.))
        s <- s %>%
          dplyr::mutate(sum = rowSums(.))
        aveB <- mean(s_bio$sum)
        aveD <- mean(s$sum)

        data.frame(site = u$site_id[i], CA_MPA_Name_Short = u$CA_MPA_Name_Short[i],
                   site_status = u$status[i],
                   year = u$year[i], transect = length(t), area = u$area[i],
                   y = u$taxon_name[i], count = aveD, biomass = aveB, ind = n)

      } else if (n == 0 & length(t2) == 0) {
        # No kelp on transect or in stipe data
        data.frame(site = u$site_id[i], CA_MPA_Name_Short = u$CA_MPA_Name_Short[i],
                   site_status = u$status[i],
                   year = u$year[i], transect = length(t), area = u$area[i],
                   y = u$taxon_name[i], count = 0, biomass = 0, ind = 0)

      } else if (n != 0 & length(t2) == 0) {
        # Kelp on transect but no stipe data - mark as NA
        data.frame(site = u$site_id[i], CA_MPA_Name_Short = u$CA_MPA_Name_Short[i],
                   site_status = u$status[i],
                   year = u$year[i], transect = length(t), area = u$area[i],
                   y = u$taxon_name[i], count = NA, biomass = NA, ind = n)

      } else {
        data.frame(site = u$site_id[i], CA_MPA_Name_Short = u$CA_MPA_Name_Short[i],
                   site_status = u$status[i],
                   year = u$year[i], transect = length(t), area = u$area[i],
                   y = u$taxon_name[i], count = 0, biomass = 0, ind = n)
      }
    }, error = function(e) {
      warning("KFM Macro bootstrap failed for site ", u$site_id[i], " year ", u$year[i], ": ", e$message)
      data.frame(site = u$site_id[i], CA_MPA_Name_Short = u$CA_MPA_Name_Short[i],
                 site_status = u$status[i],
                 year = u$year[i], transect = NA, area = u$area[i],
                 y = u$taxon_name[i], count = NA, biomass = NA, ind = NA)
    })

    results_list[[i]] <- result_row
  }

  # Combine all results at once (efficient)
  Macro.site <- do.call(rbind, results_list)

  # Save cache
  dir.create(dirname(.cache_kfm_macro), recursive = TRUE, showWarnings = FALSE)
  saveRDS(Macro.site, .cache_kfm_macro)
  cat("    Cached KFM Macrocystis bootstrap results.\n")
}

# Macro.site already contains KFM/MBON site_status and MPA assignment.
# Do NOT merge with PISCO's `sites.short` (which may exist in the session and has different site IDs).
KFM.Macro.site.merge <- Macro.site

# Convert to per m2 by dividing by transects and area
# Use safe_divide to handle any zero area/transect values
area_transect <- KFM.Macro.site.merge$area * KFM.Macro.site.merge$transect
KFM.Macro.site.merge$biomassCorr <- safe_divide(KFM.Macro.site.merge$biomass, area_transect, context = "KFM macro biomass/m2")
KFM.Macro.site.merge$densityCorr <- safe_divide(KFM.Macro.site.merge$count, area_transect, context = "KFM macro density/m2")
KFM.Macro.site.merge$densityIndCorr <- safe_divide(KFM.Macro.site.merge$ind, area_transect, context = "KFM macro ind density/m2")

###########################################################################################
## 9. Combine urchin and macro data, calculate density per m2
###########################################################################################

# Wrangle macro columns
KFM.Macro.site.all.sub <- KFM.Macro.site.merge[, c(
  "site", "site_status", "year", "y", "CA_MPA_Name_Short", "area", "transect",
  "count", "biomass", "ind",
  "biomassCorr", "densityCorr", "densityIndCorr"
)]
names(KFM.Macro.site.all.sub)[names(KFM.Macro.site.all.sub) == "biomass"] <- "biomassRaw"
names(KFM.Macro.site.all.sub)[names(KFM.Macro.site.all.sub) == "biomassCorr"] <- "biomass"
names(KFM.Macro.site.all.sub)[names(KFM.Macro.site.all.sub) == "densityCorr"] <- "density"

# Wrangle urchin columns
KFM.Urchin.site.all.sub <- KFM.Urchin.site.all[, c(
  "site", "CA_MPA_Name_Short", "site_status", "year", "area", "transect",
  "count", "biomass", "y",
  "Density", "bio.m2"
)]
names(KFM.Urchin.site.all.sub)[names(KFM.Urchin.site.all.sub) == "Density"] <- "density"
names(KFM.Urchin.site.all.sub)[names(KFM.Urchin.site.all.sub) == "biomass"] <- "biomassRaw"
names(KFM.Urchin.site.all.sub)[names(KFM.Urchin.site.all.sub) == "bio.m2"] <- "biomass"
KFM.Urchin.site.all.sub$ind <- NA_real_
KFM.Urchin.site.all.sub$densityIndCorr <- NA

# Ensure identical column sets + ordering before rbind (prevents 'numbers of columns' errors).
common_cols <- c(
  "site", "site_status", "year", "y", "CA_MPA_Name_Short", "area", "transect",
  "count", "biomassRaw", "ind", "biomass", "density", "densityIndCorr"
)
KFM.Macro.site.all.sub <- KFM.Macro.site.all.sub[, common_cols]
KFM.Urchin.site.all.sub <- KFM.Urchin.site.all.sub[, common_cols]

# Combine
KFM.bio.site.all <- rbind(KFM.Macro.site.all.sub, KFM.Urchin.site.all.sub)

# Merge with average density data
# After merging: den = density of individuals; for kelp, density = stipe density
KFM.allbio <- merge(kfm.ave, KFM.bio.site.all,
                     by.x = c("site_id", "status", "year", "taxon_name", "CA_MPA_Name_Short", "area"),
                     by.y = c("site", "site_status", "year", "y", "CA_MPA_Name_Short", "area"),
                     all.x = TRUE)

Swath.KFM <- KFM.allbio
names(Swath.KFM)[names(Swath.KFM) == "taxon_name"] <- "y"
names(Swath.KFM)[names(Swath.KFM) == "count"] <- "countRaw"
names(Swath.KFM)[names(Swath.KFM) == "density"] <- "den.y"
names(Swath.KFM)[names(Swath.KFM) == "site_id"] <- "site"

Swath.KFM.Corr <- Swath.KFM
# Make density the stipe density (not density of ind) for macro - for cross-program comparison
Swath.KFM.Corr$den <- ifelse(Swath.KFM.Corr$y == "Macrocystis pyrifera",
                               Swath.KFM.Corr$den.y, Swath.KFM.Corr$den)

# Calculate the mean across sites
kfm.ave.ave <- Swath.KFM.Corr %>%
  dplyr::group_by(y, status, CA_MPA_Name_Short, ChannelIsland, MPA_Start, time, year, area) %>%
  dplyr::summarise_at(c("den", "biomass"), mean)

# Remove MPAs that don't meet criteria for inclusion
kfm.ave.ave <- subset(kfm.ave.ave,
                        CA_MPA_Name_Short != "Anacapa Island SMCA" &
                          CA_MPA_Name_Short != "Anacapa Island SMR 1978" &
                          CA_MPA_Name_Short != "Carrington Point SMR" &
                          CA_MPA_Name_Short != "Painted Cave SMCA")

###########################################################################################
## 10. Calculate proportions and log response ratios
###########################################################################################
# This section converts raw density/biomass values to LOG RESPONSE RATIOS (lnRR).
# lnRR compares MPA sites to reference sites, standardized across different MPAs.
#
# Steps:
# 1. calculate_proportions(): Convert to proportion of maximum (0-1 scale + 0.01)
# 2. spread(): Reshape from long (one row per status) to wide (mpa/reference columns)
# 3. calculate_log_response_ratio(): Calculate lnRR = ln(mpa / reference)

# -- KFM Density proportions and response ratios --
All.den <- kfm.ave.ave
All.den <- subset(All.den, den != "N/A")  # Remove invalid values

# calculate_proportions() from 01_utils.R normalizes by max value within each MPA x taxa
All.den <- calculate_proportions(All.den, "den")

# Select columns needed for reshaping and calculate response ratio
All.den.sub <- All.den[, c("CA_MPA_Name_Short", "year", "y", "status", "area", "PropCorr")]

# spread() from tidyr: converts long format to wide format
# status column values (mpa, reference) become new columns with PropCorr values
Short.den <- All.den.sub %>%
  spread(status, PropCorr)

# calculate_log_response_ratio() from 01_utils.R: lnRR = ln(mpa / reference)
Short.den.diff <- calculate_log_response_ratio(Short.den)
Short.den.diff$resp <- "Den"  # Label as density response

# -- KFM Biomass proportions and response ratios --
# Same process as density, but for biomass values
All.bio <- kfm.ave.ave
All.bio <- subset(All.bio, biomass != "N/A")
All.bio <- calculate_proportions(All.bio, "biomass")

All.bio.sub <- All.bio[, c("CA_MPA_Name_Short", "year", "y", "status", "area", "PropCorr")]
Short.bio <- All.bio.sub %>%
  spread(status, PropCorr)
Short.bio.diff <- calculate_log_response_ratio(Short.bio)
Short.bio.diff$resp <- "Bio"  # Label as biomass response

# Combine density and biomass response ratios into one dataframe
# rbind() stacks dataframes vertically (row bind)
KFM.join.ave <- rbind(Short.den.diff, Short.bio.diff)

###########################################################################################
## 11. Create KFM.resp dataframe with raw density/biomass
###########################################################################################

# Only want 10m area for Macro
kfm.edit.den <- subset(kfm.ave.ave, y != "Macrocystis pyrifera" | area != 2)
kfm.edit.den <- subset(kfm.edit.den, y != "Macrocystis pyrifera" | area != 1)

kfm.edit.den <- kfm.edit.den[, c("status", "CA_MPA_Name_Short", "MPA_Start", "year", "den", "y")]
kfm.edit.den <- kfm.edit.den[complete.cases(kfm.edit.den), ]
KFM.den <- kfm.edit.den %>%
  spread(status, den)
KFM.den$source <- "KFM"
KFM.den <- KFM.den[complete.cases(KFM.den), ]
# FIXED: Replaced positional column indices 5:6 with explicit column names
KFM.den.long <- gather(KFM.den, status, value, mpa, reference)
KFM.den.long$resp <- "Den"

kfm.edit.bio <- subset(kfm.ave.ave, y != "Macrocystis pyrifera" | area != 2)
kfm.edit.bio <- subset(kfm.edit.bio, y != "Macrocystis pyrifera" | area != 1)

kfm.edit.bio <- kfm.edit.bio[, c("status", "CA_MPA_Name_Short", "MPA_Start", "year", "biomass", "y")]
kfm.edit.bio <- kfm.edit.bio[complete.cases(kfm.edit.bio), ]
KFM.bio <- kfm.edit.bio %>%
  spread(status, biomass)
KFM.bio$source <- "KFM"
KFM.bio <- KFM.bio[complete.cases(KFM.bio), ]
# FIXED: Replaced positional column indices 5:6 with explicit column names
KFM.bio.long <- gather(KFM.bio, status, value, mpa, reference)
KFM.bio.long$resp <- "Bio"

KFM.resp <- rbind(KFM.den.long, KFM.bio.long)

###########################################################################################
## 12. Assign time and Before/After labels to KFM.join.ave
###########################################################################################

# Assign time since MPA
KFM.join.ave <- assign_time_from_site_table(KFM.join.ave, Sites2)

# Merge with site metadata (type, location, hectares)
KFM.join.ave <- merge(KFM.join.ave, sites.short.edit,
                        by.x = c("CA_MPA_Name_Short"),
                        by.y = c("CA_MPA_Name_Short"),
                        all = TRUE)
KFM.join.ave <- KFM.join.ave[complete.cases(KFM.join.ave$year), ]
KFM.join.ave$source <- "KFM"

# Assign Before/After labels using utility function
KFM.join.ave <- assign_ba_from_site_table(KFM.join.ave, Sites2)

###########################################################################################
## 13. Process KFM fish data (sheephead from MBON fish file)
###########################################################################################

mbon.fish <- read.csv(here::here("data", "MBON", "SBCMBON_kelp_forest_integrated_fish_20231022.csv"))
mbon.fish <- subset(mbon.fish, data_source == "kfm")

# Subset sheephead (Semicossyphus pulcher)
mbon.fish.spul <- subset(mbon.fish, proj_taxon_id == "t-k-127" | proj_taxon_id == "t-k-233")

# Date conversion
tmpDateFormat <- "%Y-%m-%d"
tmp1date <- as.Date(mbon.fish.spul$date, format = tmpDateFormat)
if (length(tmp1date) == length(tmp1date[!is.na(tmp1date)])) {
  mbon.fish.spul$date <- tmp1date
} else {
  print("Date conversion failed for mbon.fish.spul$date.")
}
rm(tmpDateFormat, tmp1date)
mbon.fish.spul$year <- year(mbon.fish.spul$date)

# Subset relevant survey methods
mbon.fish.spul.sub <- subset(mbon.fish.spul,
                               (sample_method == "rdfc" & year >= 2003) |
                                 (sample_method == "visualfish"))

# Join to site table
mbon.spul.site <- merge(mbon.fish.spul.sub, Sites2, by = c("site_id"))
mbon.spul.site$count <- as.numeric(mbon.spul.site$count)
mbon.spul.site <- subset(mbon.spul.site, CA_MPA_Name_Short != "")
mbon.spul.site$den <- mbon.spul.site$count / mbon.spul.site$area

# Calculate sheephead density per MPA/Reference - averaging hierarchy
mbon.spul.sum <- mbon.spul.site %>%
  dplyr::group_by(site_id, sample_method, transect_id, status, CA_MPA_Name_Short, ChannelIsland, MPA_Start, year) %>%
  dplyr::summarise_at(c("den"), mean)

mbon.spul.ave <- mbon.spul.sum %>%
  dplyr::group_by(site_id, sample_method, status, CA_MPA_Name_Short, ChannelIsland, MPA_Start, year) %>%
  dplyr::summarise_at(c("den"), mean)

mbon.spul.max <- mbon.spul.ave %>%
  dplyr::group_by(status, sample_method, CA_MPA_Name_Short, ChannelIsland, MPA_Start, year) %>%
  dplyr::summarise_at(c("den"), mean)

mbon.spul.max <- subset(mbon.spul.max, CA_MPA_Name_Short != "Anacapa Island SMR 1978")
colnames(mbon.spul.max)[colnames(mbon.spul.max) == "den"] <- "count"

mbon.spul.max$Prop <- 0
mbon.spul.max$PropCorr <- 0
mbon.spul.max$taxon_name <- "SPUL"

mbon.spul.max.rd <- subset(mbon.spul.max, sample_method == "rdfc")
mbon.spul.max.vf <- subset(mbon.spul.max, sample_method == "visualfish")

# -- KFM fish raw density dataframe (visual fish method) --
mbon.fish.edit <- mbon.spul.max.vf[, c("status", "CA_MPA_Name_Short", "MPA_Start", "year", "count", "taxon_name")]
KFM.fish.den <- mbon.fish.edit %>%
  spread(status, count)
KFM.fish.den$source <- "KFM"
KFM.fish.den.sub <- KFM.fish.den[complete.cases(KFM.fish.den), ]
# FIXED: Replaced positional column indices 5:6 with explicit column names
KFM.fish.den.long <- gather(KFM.fish.den.sub, status, value, mpa, reference)
KFM.fish.den.long$resp <- "Den"

# -- Proportions for RDFC method --
# Use standardized calculate_proportions() function from 01_utils.R
# Add 'y' column required by the function (taxon_name already set to "SPUL")
# NOTE: Using adaptive correction (min/2) for consistency across all taxa/programs.
# This is statistically preferred over fixed +0.01 (Aitchison 1986) as it scales
# appropriately with the data and avoids inflating effect sizes for rare species.
mbon.spul.max.rd$y <- mbon.spul.max.rd$taxon_name
mbon.spul.max.rd <- calculate_proportions(mbon.spul.max.rd, "count")

All.den.spul.kfm.sub.rd <- mbon.spul.max.rd[, c("CA_MPA_Name_Short", "sample_method", "year", "taxon_name", "status", "PropCorr")]
Short.den.spul.kfm.sub.rd <- All.den.spul.kfm.sub.rd %>%
  spread(status, PropCorr)
Short.den.spul.kfm.sub.rd.diff <- calculate_log_response_ratio(Short.den.spul.kfm.sub.rd)
names(Short.den.spul.kfm.sub.rd.diff)[names(Short.den.spul.kfm.sub.rd.diff) == "taxon_name"] <- "y"
Short.den.spul.kfm.sub.rd.diff$resp <- "RD"

# -- Proportions for visual fish method --
# Use standardized calculate_proportions() function from 01_utils.R
# NOTE: Using adaptive correction (min/2) for consistency - see RDFC method comment above.
mbon.spul.max.vf$y <- mbon.spul.max.vf$taxon_name
mbon.spul.max.vf <- calculate_proportions(mbon.spul.max.vf, "count")

All.den.spul.kfm.sub.vf <- mbon.spul.max.vf[, c("CA_MPA_Name_Short", "sample_method", "year", "taxon_name", "status", "PropCorr")]
Short.den.spul.kfm.sub.vf <- All.den.spul.kfm.sub.vf %>%
  spread(status, PropCorr)
Short.den.spul.kfm.sub.vf.diff <- calculate_log_response_ratio(Short.den.spul.kfm.sub.vf)
names(Short.den.spul.kfm.sub.vf.diff)[names(Short.den.spul.kfm.sub.vf.diff) == "taxon_name"] <- "y"
Short.den.spul.kfm.sub.vf.diff$resp <- "Den"

# Combine RDFC and visual fish response ratios
Short.den.spul.kfm.sub.diff <- rbind(Short.den.spul.kfm.sub.vf.diff, Short.den.spul.kfm.sub.rd.diff)

# Assign time since MPA
Short.den.spul.kfm.sub.diff <- assign_time_from_site_table(Short.den.spul.kfm.sub.diff, Sites2)

# Merge with site metadata
kfm.fish <- merge(Short.den.spul.kfm.sub.diff, sites.short.edit,
                    by.x = c("CA_MPA_Name_Short"),
                    by.y = c("CA_MPA_Name_Short"),
                    all = TRUE)
kfm.fish <- kfm.fish[complete.cases(kfm.fish$year), ]
kfm.fish$source <- "KFM"

# Assign Before/After labels
kfm.fish <- assign_ba_from_site_table(kfm.fish, Sites2)

###########################################################################################
## Final outputs:
##   KFM.join.ave  - Response ratio data for urchins, macro (density + biomass)
##   KFM.resp      - Raw density/biomass data inside/outside MPAs
##   kfm.fish      - Sheephead fish response ratio data
##   KFM.fish.den.long - Raw fish density data
##   Sites2        - Site geolocation table (used by downstream scripts)
##   lter          - LTER subset of MBON (used by 06_lter_processing.R)
##   kfm.stipe.site, kfm.stipe.OG - Stipe data (used downstream)
###########################################################################################

cat("KFM/MBON processing complete.\n")
