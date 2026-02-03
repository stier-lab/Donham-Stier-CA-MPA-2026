# =============================================================================
# 06_lter_processing.R
# =============================================================================
#
# PURPOSE:
#   Process SBC LTER (Santa Barbara Coastal Long-Term Ecological Research)
#   data for the MPA pBACIPS analysis.
#
# WHAT THIS SCRIPT DOES:
#   1. Processes LTER urchin data (S. purpuratus, M. franciscanus)
#   2. Processes LTER lobster abundance and biomass data
#   3. Processes LTER Macrocystis kelp frond data
#   4. Processes LTER sheephead fish data
#   5. Calculates biomass via bootstrap resampling (urchins) or
#      allometric equations (lobster, fish)
#   6. Calculates log response ratios (MPA vs reference sites)
#   7. Assigns time since MPA and Before/After labels
#
# DATA SOURCES:
#   SBC LTER (https://sbclter.msi.ucsb.edu/) has monitored kelp forests
#   since 2000. Key datasets used:
#   - Annual_Kelp_All_Years: Frond counts (1984-present)
#   - Annual_fish_comb: Fish abundance by size (1984-present)
#   - Lobster_Abundance_All_Years: Lobster size-abundance (2012-present)
#
# KEY FEATURES OF LTER DATA:
#   - Only covers Campus Point SMCA and Naples SMCA (mainland MPAs)
#   - Has complementary datasets for same sites (fish, kelp, inverts)
#   - Standardized protocols since program inception
#   - Lobster data has individual sizes (no bootstrap needed)
#
# INPUTS:
#   - lter object from 05_kfm_processing.R (MBON subset)
#   - data/LTER/Annual_Kelp_All_Years_20240305.csv
#   - data/LTER/Annual_fish_comb_20240307.csv
#   - data/LTER/Lobster_Abundance_All_Years_20230831.csv
#   - SizeFreq.Urch.OG (from 04_pisco_processing.R)
#
# OUTPUTS:
#   - LTER.join.ave: Urchin density & biomass response ratios
#   - LTER.lob: Lobster density & biomass response ratios
#   - Short.lter.macro: Kelp frond density & biomass response ratios
#   - LTER.fish: Sheephead density & biomass response ratios
#   - LTER.resp, LTER.lob.resp, LTER.macro.resp, LTER.fish.resp: Raw data
#
# DEPENDENCIES:
#   Requires 00-05 scripts to be sourced first
#
# AUTHORS: Emily Donham & Adrian Stier
# PROJECT: CA MPA Kelp Forest pBACIPS Analysis
# =============================================================================

# Source utility functions
source(here::here("code", "R", "01_utils.R"))

###########################################################################################
## 1. Process LTER urchin/macro data from MBON integrated dataset
###########################################################################################
# The 'lter' object was created in 05_kfm_processing.R by subsetting MBON data.
# LTER = Santa Barbara Coastal Long-Term Ecological Research program.
# They monitor kelp forest sites along the mainland coast (not Channel Islands).

# Merge LTER survey data with site metadata
# merge() joins two tables like a database JOIN or Excel VLOOKUP
lter.site <- merge(lter, Sites2, by = c("site_id"))

# Remove sites with no MPA association and Arroyo Honda (data quality issues per Bob Miller)
lter.site <- subset(lter.site, CA_MPA_Name_Short != "" & site_id != "a-l-02")

# Extract year from date column using lubridate's year() function
lter.site$year <- year(lter.site$date)

# Assign time since MPA implementation using utility function from 01_utils.R
# This calculates how many years after MPA start each observation occurred
# (replaces hundreds of if-else lines with a single function call)
lter.site <- assign_time_from_site_table(lter.site, Sites2)

# Subset to only the two LTER MPAs we're analyzing
# These are the mainland MPAs monitored by LTER with sufficient before/after data
lter.site <- subset(lter.site, CA_MPA_Name_Short == "Campus Point SMCA" |
                      CA_MPA_Name_Short == "Naples SMCA")

# Subset target taxa: urchins from quads, macrocystis from swath
# Different taxa are surveyed using different methods:
# - Urchins: counted in quadrats (small fixed areas)
# - Kelp: counted along swath transects (longer strips)
# The | operator means "OR" - keep any row that matches any condition
lter.sub <- subset(lter.site,
                   (taxon_name == "Strongylocentrotus purpuratus" & sample_method == "quad") |
                     (taxon_name == "Mesocentrotus franciscanus" & sample_method == "quad") |
                     (taxon_name == "Macrocystis pyrifera" & sample_method == "swath"))

# Calculate density: count divided by survey area gives individuals per square meter
lter.sub$den <- lter.sub$count / lter.sub$area

# Calculate site means using two-step aggregation (nested sampling design)
# Step 1: Average across replicates within each transect
# group_by() defines the grouping structure; summarise_at() calculates mean for 'den'
lter.ave <- lter.sub %>%
  group_by(site_id, transect_id, taxon_name, status, CA_MPA_Name_Short,
           ChannelIsland, MPA_Start, time, year) %>%
  summarise_at(c("den"), mean)

# Step 2: Average across transects to get site-level annual means
lter.ave <- lter.ave %>%
  group_by(site_id, taxon_name, status, CA_MPA_Name_Short,
           ChannelIsland, MPA_Start, time, year) %>%
  summarise_at(c("den"), mean)

###########################################################################################
## 2. Bootstrap urchin biomass from size frequency distributions
###########################################################################################

URCHINS <- subset(lter.site,
                  taxon_name == "Strongylocentrotus purpuratus" |
                    taxon_name == "Mesocentrotus franciscanus")
URCHINS <- mutate(URCHINS, taxon_name = case_when(
  taxon_name == "Strongylocentrotus purpuratus" ~ "STRPURAD",
  taxon_name == "Mesocentrotus franciscanus"    ~ "MESFRAAD",
  TRUE ~ taxon_name))
URCHINS$replicate_id <- as.numeric(URCHINS$replicate_id)

# Get unique site-year-transect-taxon combinations
u <- unique(URCHINS[, c('year', 'site_id', 'CA_MPA_Name_Short', 'status',
                        'area', 'taxon_name', 'transect_id')])

# Prepare size frequency data
SizeFreq.Urch <- SizeFreq.Urch.OG
SizeFreq.Urch <- SizeFreq.Urch[complete.cases(SizeFreq.Urch$size), ]

# Bootstrap urchin biomass using utility function
Urchin.site <- data.frame(site = NA, CA_MPA_Name_Short = NA, site_status = NA,
                          year = NA, transect = NA, quad = NA,
                          biomass = NA, count = NA, y = NA, area = NA)
Urchin.SF <- data.frame(site = NA, CA_MPA_Name_Short = NA, site_status = NA,
                        year = NA, transect = NA, quad = NA,
                        biomass = NA, count = NA, y = NA, area = NA)

# Cache bootstrap results to avoid re-running the slow loop
.cache_lter_urch <- here::here("data", "cache", "lter_urchin_bootstrap.rds")
if (file.exists(.cache_lter_urch) && !exists("FORCE_BOOTSTRAP")) {
  cat("    Loading cached LTER urchin bootstrap results...\n")
  Urchin.site <- readRDS(.cache_lter_urch)
} else {
  for (i in 1:nrow(u)) {
    t <- which(URCHINS$site_id == u$site_id[i] &
                 URCHINS$year == u$year[i] &
                 URCHINS$transect_id == u$transect_id[i] &
                 URCHINS$taxon_name == u$taxon_name[i])
    n <- sum(URCHINS$count[t])
    t2 <- which(SizeFreq.Urch$CA_MPA_Name_Short == u$CA_MPA_Name_Short[i] &
                  SizeFreq.Urch$site_status == u$status[i] &
                  SizeFreq.Urch$year == u$year[i] &
                  SizeFreq.Urch$classcode == u$taxon_name[i])

    # Select appropriate biomass function based on species
    bio_fun <- if (u$taxon_name[i] == "MESFRAAD") bio_redurch else bio_purpurch

    # Use bootstrap_biomass utility function
    result <- bootstrap_biomass(count = n,
                                size_freq_indices = t2,
                                size_freq_table = SizeFreq.Urch,
                                biomass_fun = bio_fun,
                                n_resamples = 1000)

    Urchin.SF$site           <- u$site_id[i]
    Urchin.SF$CA_MPA_Name_Short <- u$CA_MPA_Name_Short[i]
    Urchin.SF$year           <- u$year[i]
    Urchin.SF$transect       <- u$transect_id[i]
    Urchin.SF$quad           <- length(t)
    Urchin.SF$site_status    <- u$status[i]
    Urchin.SF$area           <- u$area[i]
    Urchin.SF$y              <- u$taxon_name[i]
    Urchin.SF$biomass        <- result$biomass
    Urchin.SF$count          <- result$count
    Urchin.site <- rbind(Urchin.site, Urchin.SF)
  }

  # Save cache
  dir.create(dirname(.cache_lter_urch), recursive = TRUE, showWarnings = FALSE)
  saveRDS(Urchin.site, .cache_lter_urch)
  cat("    Cached LTER urchin bootstrap results.\n")
}

# Merge bootstrapped biomass with site info
LTER.Urchin.site.merge <- merge(Urchin.site, sites.short,
                                by.x = c("site", "site_status", "CA_MPA_Name_Short"),
                                by.y = c("site", "site_status", "CA_MPA_Name_Short"),
                                all.x = TRUE)
LTER.Urchin.site.merge <- LTER.Urchin.site.merge[, colnames(LTER.Urchin.site.merge)[c(3, 1, 2, 4, 17, 20, 5:10)]]
# Zero-fill only numeric columns; preserve NAs in character/factor columns
num_cols <- sapply(LTER.Urchin.site.merge, is.numeric)
LTER.Urchin.site.merge[num_cols][is.na(LTER.Urchin.site.merge[num_cols])] <- 0

# Remove site-years where urchins were counted but not sized (can't estimate biomass)
LTER.Urchin.null <- subset(LTER.Urchin.site.merge, count > 0 & biomass == 0)
LTER.Urchin.site.merge.sub <- LTER.Urchin.site.merge
for (i in 1:nrow(LTER.Urchin.null)) {
  LTER.Urchin.site.merge.sub <- subset(LTER.Urchin.site.merge.sub,
                                       site != LTER.Urchin.null$site[i] |
                                         year != LTER.Urchin.null$year[i])
}

# Normalize by number of quads per transect
LTER.Urchin.site.transectSum <- LTER.Urchin.site.merge.sub
LTER.Urchin.site.transectSum$count   <- LTER.Urchin.site.transectSum$count / LTER.Urchin.site.transectSum$quad
LTER.Urchin.site.transectSum$biomass <- LTER.Urchin.site.transectSum$biomass / LTER.Urchin.site.transectSum$quad

# Average across transects
LTER.Urchin.site.all <- LTER.Urchin.site.transectSum %>%
  group_by(year, site, CA_MPA_Name_Short, site_designation, site_status,
           BaselineRegion, y, area) %>%
  summarise_at(c("count", "biomass"), mean)
LTER.Urchin.site.all$Density <- LTER.Urchin.site.all$count / LTER.Urchin.site.all$area

# Recode taxa in lter.ave to match
lter.ave <- mutate(lter.ave, taxon_name = case_when(
  taxon_name == "Strongylocentrotus purpuratus" ~ "STRPURAD",
  taxon_name == "Mesocentrotus franciscanus"    ~ "MESFRAAD",
  TRUE ~ taxon_name))

# Merge density and biomass data
LTER.allbio <- merge(lter.ave, LTER.Urchin.site.all,
                     by.x = c("site_id", "status", "year", "taxon_name", "CA_MPA_Name_Short"),
                     by.y = c("site", "site_status", "year", "y", "CA_MPA_Name_Short"),
                     all = TRUE)

Swath.LTER <- LTER.allbio
Swath.LTER <- Swath.LTER[, colnames(Swath.LTER)[c(1:5, 12, 6:10, 13:15)]]
names(Swath.LTER)[names(Swath.LTER) == "Density"]   <- "count.y"
names(Swath.LTER)[names(Swath.LTER) == "status"]    <- "site_status"
names(Swath.LTER)[names(Swath.LTER) == "den"]       <- "Density"
names(Swath.LTER)[names(Swath.LTER) == "site_id"]   <- "site"
names(Swath.LTER)[names(Swath.LTER) == "taxon_name"] <- "y"

###########################################################################################
## 3. Calculate proportions and log response ratios for urchins
###########################################################################################

# Average across sites
lter.ave.ave <- Swath.LTER %>%
  group_by(y, site_status, CA_MPA_Name_Short, ChannelIsland, MPA_Start, time, year) %>%
  summarise_at(vars("Density", "biomass"), mean)

# Restore full species names
lter.ave.ave <- mutate(lter.ave.ave, y = case_when(
  y == "STRPURAD" ~ "Strongylocentrotus purpuratus",
  y == "MESFRAAD" ~ "Mesocentrotus franciscanus",
  TRUE ~ y))
names(lter.ave.ave)[names(lter.ave.ave) == "Density"] <- "den"

## --- LTER density proportions and log response ratios ---
All.den <- lter.ave.ave
All.den <- All.den[complete.cases(All.den$den), ]
All.den <- calculate_proportions(All.den, "den")

All.den.sub <- All.den[, colnames(All.den)[c(3, 7, 1, 2, 11)]]
Short.den <- All.den.sub %>% spread(site_status, PropCorr)
Short.den.diff <- calculate_log_response_ratio(Short.den)
Short.den.diff$resp <- "Den"

## --- LTER biomass proportions and log response ratios ---
All.bio <- lter.ave.ave
All.bio <- All.bio[complete.cases(All.bio$biomass), ]
All.bio <- calculate_proportions(All.bio, "biomass")

All.bio.sub <- All.bio[, colnames(All.bio)[c(3, 7, 1, 2, 11)]]
Short.bio <- All.bio.sub %>% spread(site_status, PropCorr)
Short.bio.diff <- calculate_log_response_ratio(Short.bio)
Short.bio.diff$resp <- "Bio"

# Combine density and biomass response ratios
LTER.join.ave <- rbind(Short.den.diff, Short.bio.diff)

## --- Create LTER urchin response dataframes (density/biomass, long format) ---
# These dataframes store the raw density/biomass values for MPA vs reference sites.
# They're used for plotting raw data time series and as input to other analyses.

# Select relevant columns for density response dataframe
lter.edit.den <- lter.ave.ave[, colnames(lter.ave.ave)[c(2, 3, 5, 7, 8, 1)]]
lter.edit.den <- lter.edit.den[complete.cases(lter.edit.den), ]  # Remove rows with NAs

# spread() converts from long to wide format: site_status values become columns
LTER.den <- lter.edit.den %>% spread(site_status, den)
LTER.den <- subset(LTER.den, y != "Macrocystis pyrifera")  # Kelp processed separately
LTER.den$source <- "LTER"  # Add source label for tracking
LTER.den <- LTER.den[complete.cases(LTER.den), ]

# gather() converts from wide back to long format for storage
# Columns 5:6 (mpa, reference) are gathered into key-value pairs
# key = "site_status" (mpa or reference), value = "value" (the density)
LTER.den.long <- gather(LTER.den, site_status, value, 5:6)
LTER.den.long$resp <- "Den"  # Label as density response

# Same process for biomass
lter.edit.bio <- lter.ave.ave[, colnames(lter.ave.ave)[c(2, 3, 5, 7, 9, 1)]]
lter.edit.bio <- lter.edit.bio[complete.cases(lter.edit.bio), ]
LTER.bio <- lter.edit.bio %>% spread(site_status, biomass)
LTER.bio$source <- "LTER"
LTER.bio <- LTER.bio[complete.cases(LTER.bio), ]
LTER.bio.long <- gather(LTER.bio, site_status, value, 5:6)
LTER.bio.long$resp <- "Bio"

# Combine density and biomass into one response dataframe
# rbind() stacks dataframes vertically (row bind)
LTER.resp <- rbind(LTER.bio.long, LTER.den.long)

## --- Assign time and BA labels to LTER.join.ave ---
LTER.join.ave <- assign_time_from_site_table(LTER.join.ave, Sites2)
LTER.join.ave <- merge(LTER.join.ave, sites.short.edit,
                       by.x = c("CA_MPA_Name_Short"),
                       by.y = c("CA_MPA_Name_Short"),
                       all = TRUE)
LTER.join.ave <- LTER.join.ave[complete.cases(LTER.join.ave$year), ]
LTER.join.ave$source <- "LTER"
LTER.join.ave <- assign_ba_from_site_table(LTER.join.ave, Sites2)

###########################################################################################
## 4. Import and process LTER lobster data
###########################################################################################
# LTER has a separate lobster monitoring dataset with individual size measurements.
# Unlike PISCO/KFM, we don't need bootstrap resampling because LTER measures
# each individual lobster's carapace length directly.

lter.lob <- read.csv(here::here("data", "LTER", "Lobster_Abundance_All_Years_20230831.csv"))

# Merge with site metadata using SITE column as the key
lter.lob.site <- merge(lter.lob, Sites2, by.x = c("SITE"), by.y = c("site_id"))
lter.lob.site <- subset(lter.lob.site, CA_MPA_Name_Short != "")

# Calculate biomass from size using allometric equation
# bio_lobster() from 01_utils.R converts carapace length (mm) to weight (g)
# Total biomass = count * weight per individual
lter.lob.site$biomass <- lter.lob.site$COUNT * bio_lobster(lter.lob.site$SIZE)
lter.lob.site$biomass[is.na(lter.lob.site$biomass)] <- 0  # NA sizes -> 0 biomass

# Summarise: replicate -> transect -> site -> MPA level
lter.lob.site.sum <- lter.lob.site %>%
  group_by(YEAR, SITE, TRANSECT, REPLICATE, status, CA_MPA_Name_Short,
           ChannelIsland, MPA_Start) %>%
  summarise_at(c("COUNT", "biomass"), sum)

lter.lob.site.ave <- lter.lob.site.sum %>%
  group_by(YEAR, SITE, TRANSECT, status, CA_MPA_Name_Short,
           ChannelIsland, MPA_Start) %>%
  summarise_at(c("COUNT", "biomass"), mean)

lter.lob.site.ave <- lter.lob.site.ave %>%
  group_by(YEAR, SITE, status, CA_MPA_Name_Short, ChannelIsland, MPA_Start) %>%
  summarise_at(c("COUNT", "biomass"), mean)

lter.lob.site.max <- lter.lob.site.ave %>%
  group_by(status, CA_MPA_Name_Short, MPA_Start, YEAR) %>%
  summarise_at(c("COUNT", "biomass"), mean)

## --- Lobster biomass proportions and log response ratios ---
lter.lob.site.max$Prop <- 0
lter.lob.site.max$PropCorr <- 0
lter.lob.site.max$taxon_name <- "Panulirus interruptus"

# Biomass proportions
x_lob <- unique(lter.lob.site.max$CA_MPA_Name_Short)
for (i in seq_along(x_lob)) {
  idx <- which(lter.lob.site.max$CA_MPA_Name_Short == x_lob[i])
  m <- max(lter.lob.site.max$biomass[idx])
  if (m > 0) {
    lter.lob.site.max$Prop[idx] <- lter.lob.site.max$biomass[idx] / m
    lter.lob.site.max$PropCorr[idx] <- lter.lob.site.max$Prop[idx] + 0.01
  }
}

All.bio.panLTERall.sub <- lter.lob.site.max[, colnames(lter.lob.site.max)[c(2, 4, 9, 1, 8)]]
Short.bio.panLTERall <- All.bio.panLTERall.sub %>% spread(status, PropCorr)
Short.bio.panLTERall.diff <- calculate_log_response_ratio(Short.bio.panLTERall)
Short.bio.panLTERall.diff$resp <- "Bio"

# Count/density proportions
for (i in seq_along(x_lob)) {
  idx <- which(lter.lob.site.max$CA_MPA_Name_Short == x_lob[i])
  m <- max(lter.lob.site.max$COUNT[idx])
  if (m > 0) {
    lter.lob.site.max$Prop[idx] <- lter.lob.site.max$COUNT[idx] / m
    lter.lob.site.max$PropCorr[idx] <- lter.lob.site.max$Prop[idx] + 0.01
  }
}

All.den.panLTERall.sub <- lter.lob.site.max[, colnames(lter.lob.site.max)[c(2, 4, 9, 1, 8)]]
Short.den.panLTERall <- All.den.panLTERall.sub %>% spread(status, PropCorr)
Short.den.panLTERall.diff <- calculate_log_response_ratio(Short.den.panLTERall)
Short.den.panLTERall.diff$resp <- "Den"

# Combine lobster results
LTER.lob <- rbind(Short.bio.panLTERall.diff, Short.den.panLTERall.diff)
colnames(LTER.lob)[2] <- "year"
colnames(LTER.lob)[3] <- "y"

# Assign time and BA labels using utility functions
LTER.lob <- assign_time_from_site_table(LTER.lob, Sites2)
LTER.lob <- merge(LTER.lob, sites.short.edit,
                  by.x = c("CA_MPA_Name_Short"),
                  by.y = c("CA_MPA_Name_Short"),
                  all = TRUE)
LTER.lob <- LTER.lob[complete.cases(LTER.lob$year), ]
LTER.lob$source <- "LTER lob surveys"
LTER.lob$BA <- "N/A"
LTER.lob <- assign_ba_from_site_table(LTER.lob, Sites2)

## --- Create LTER lobster response dataframes ---
lter.lob.site.max$Den <- lter.lob.site.max$COUNT / 300  # plots are 300 m2
lter.lob.site.max$Bio <- lter.lob.site.max$biomass / 300

lter.lob.edit.den <- lter.lob.site.max[, colnames(lter.lob.site.max)[c(1:4, 10, 9)]]
LTER.lob.den <- lter.lob.edit.den %>% spread(status, Den)
LTER.lob.den$source <- "LTER"
LTER.lob.den.sub <- LTER.lob.den[complete.cases(LTER.lob.den), ]
LTER.lob.den.long <- gather(LTER.lob.den.sub, status, value, 5:6)
LTER.lob.den.long$resp <- "Den"

lter.lob.edit.bio <- lter.lob.site.max[, colnames(lter.lob.site.max)[c(1:4, 11, 9)]]
LTER.lob.bio <- lter.lob.edit.bio %>% spread(status, Bio)
LTER.lob.bio$source <- "LTER"
LTER.lob.bio.sub <- LTER.lob.bio[complete.cases(LTER.lob.bio), ]
LTER.lob.bio.long <- gather(LTER.lob.bio.sub, status, value, 5:6)
LTER.lob.bio.long$resp <- "Bio"

LTER.lob.resp <- rbind(LTER.lob.den.long, LTER.lob.bio.long)

###########################################################################################
## 5. Import and process LTER kelp stipe data
###########################################################################################

lter.macro <- read.csv(here::here("data", "LTER", "Annual_Kelp_All_Years_20240305.csv"))
lter.macro.site <- merge(lter.macro, Sites2, by.x = c("SITE"), by.y = c("site_id"))

# Filter: remove empty MPA names, Arroyo Honda, invalid frond counts, non-target MPAs
lter.macro.site <- subset(lter.macro.site, CA_MPA_Name_Short != "" & SITE != "AHND")
lter.macro.site <- subset(lter.macro.site, FRONDS != -99999)
lter.macro.site <- subset(lter.macro.site, CA_MPA_Name_Short == "Campus Point SMCA" |
                            CA_MPA_Name_Short == "Naples SMCA")

# Calculate frond density and biomass using bio_macro() utility function
lter.macro.site$frondDen <- lter.macro.site$FRONDS / lter.macro.site$AREA
lter.macro.site$biomass  <- bio_macro(lter.macro.site$frondDen)

# Summarise: quad -> transect -> site -> MPA level
lter.macro.site.sum <- lter.macro.site %>%
  group_by(YEAR, SITE, TRANSECT, QUAD, SIDE, AREA, status,
           CA_MPA_Name_Short, ChannelIsland, MPA_Start) %>%
  summarise_at(c("FRONDS", "biomass"), sum)
lter.macro.site.sum$den <- lter.macro.site.sum$FRONDS / lter.macro.site.sum$AREA

lter.macro.site.ave <- lter.macro.site.sum %>%
  group_by(YEAR, TRANSECT, SITE, status, CA_MPA_Name_Short,
           ChannelIsland, MPA_Start) %>%
  summarise_at(c("den", "biomass"), mean)

lter.macro.site.ave <- lter.macro.site.ave %>%
  group_by(YEAR, SITE, status, CA_MPA_Name_Short, ChannelIsland, MPA_Start) %>%
  summarise_at(c("den", "biomass"), mean)

lter.macro.site.max <- lter.macro.site.ave %>%
  group_by(status, CA_MPA_Name_Short, MPA_Start, YEAR) %>%
  summarise_at(c("den", "biomass"), mean)

## --- Macrocystis frond density proportions and log response ratios ---
lter.macro.site.max$Prop <- 0
lter.macro.site.max$PropCorr <- 0
lter.macro.site.max$taxon_name <- "Macrocystis pyrifera"

x_macro <- unique(lter.macro.site.max$CA_MPA_Name_Short)

# Frond density proportions
for (i in seq_along(x_macro)) {
  idx <- which(lter.macro.site.max$CA_MPA_Name_Short == x_macro[i])
  m <- max(lter.macro.site.max$den[idx])
  if (m > 0) {
    lter.macro.site.max$Prop[idx] <- lter.macro.site.max$den[idx] / m
    lter.macro.site.max$PropCorr[idx] <- lter.macro.site.max$Prop[idx] + 0.01
  }
}

All.lter.macro.frond.sub <- lter.macro.site.max[, colnames(lter.macro.site.max)[c(2, 4, 9, 1, 8)]]
Short.lter.macro.frond <- All.lter.macro.frond.sub %>% spread(status, PropCorr)
Short.lter.macro.frond.diff <- calculate_log_response_ratio(Short.lter.macro.frond)
Short.lter.macro.frond.diff$resp <- "Den"
colnames(Short.lter.macro.frond.diff)[2] <- "year"
colnames(Short.lter.macro.frond.diff)[3] <- "y"

# Biomass proportions
for (i in seq_along(x_macro)) {
  idx <- which(lter.macro.site.max$CA_MPA_Name_Short == x_macro[i])
  m <- max(lter.macro.site.max$biomass[idx])
  if (m > 0) {
    lter.macro.site.max$Prop[idx] <- lter.macro.site.max$biomass[idx] / m
    lter.macro.site.max$PropCorr[idx] <- lter.macro.site.max$Prop[idx] + 0.01
  }
}

All.lter.macro.bio.sub <- lter.macro.site.max[, colnames(lter.macro.site.max)[c(2, 4, 9, 1, 8)]]
Short.lter.macro.bio <- All.lter.macro.bio.sub %>% spread(status, PropCorr)
Short.lter.macro.bio.diff <- calculate_log_response_ratio(Short.lter.macro.bio)
Short.lter.macro.bio.diff$resp <- "Bio"
colnames(Short.lter.macro.bio.diff)[2] <- "year"
colnames(Short.lter.macro.bio.diff)[3] <- "y"

# Combine macrocystis results
Short.lter.macro <- rbind(Short.lter.macro.frond.diff, Short.lter.macro.bio.diff)

# Assign time and BA labels using utility functions
Short.lter.macro <- assign_time_from_site_table(Short.lter.macro, Sites2)
Short.lter.macro <- merge(Short.lter.macro, sites.short.edit,
                          by.x = c("CA_MPA_Name_Short"),
                          by.y = c("CA_MPA_Name_Short"),
                          all = TRUE)
Short.lter.macro <- Short.lter.macro[complete.cases(Short.lter.macro$year), ]
Short.lter.macro$source <- "LTER macro surveys"
Short.lter.macro <- assign_ba_from_site_table(Short.lter.macro, Sites2)

## --- Create LTER macrocystis response dataframes ---
lter.macro.edit.den <- lter.macro.site.max[, colnames(lter.macro.site.max)[c(1:4, 5, 9)]]
LTER.macro.den <- lter.macro.edit.den %>% spread(status, den)
LTER.macro.den$source <- "LTER"
LTER.macro.den.sub <- LTER.macro.den[complete.cases(LTER.macro.den), ]
LTER.macro.den.long <- gather(LTER.macro.den.sub, status, value, 5:6)
LTER.macro.den.long$resp <- "Den"

lter.macro.edit.bio <- lter.macro.site.max[, colnames(lter.macro.site.max)[c(1:4, 6, 9)]]
LTER.macro.bio <- lter.macro.edit.bio %>% spread(status, biomass)
LTER.macro.bio$source <- "LTER"
LTER.macro.bio.sub <- LTER.macro.bio[complete.cases(LTER.macro.bio), ]
LTER.macro.bio.long <- gather(LTER.macro.bio.sub, status, value, 5:6)
LTER.macro.bio.long$resp <- "Bio"

LTER.macro.resp <- rbind(LTER.macro.den.long, LTER.macro.bio.long)

###########################################################################################
## 6. Import and process LTER fish data (sheephead only)
###########################################################################################
# LTER conducts visual fish surveys along with their benthic monitoring.
# We focus on California sheephead (Semicossyphus pulcher), a key urchin predator.
# Sheephead are heavily fished outside MPAs, making them a good indicator species.

lter.fish <- read.csv(here::here("data", "LTER", "Annual_fish_comb_20240307.csv"))

# Subset to California sheephead only (Semicossyphus pulcher)
lter.fish.sub <- subset(lter.fish, SCIENTIFIC_NAME == "Semicossyphus pulcher")

# Merge with site metadata and filter to our study MPAs
lter.fish.sub.site <- merge(lter.fish.sub, Sites2, by.x = c("SITE"), by.y = c("site_id"))
lter.fish.sub.site <- subset(lter.fish.sub.site, CA_MPA_Name_Short != "" & SITE != "AHND")
lter.fish.sub.site <- subset(lter.fish.sub.site, CA_MPA_Name_Short == "Campus Point SMCA" |
                                CA_MPA_Name_Short == "Naples SMCA")

# Remove invalid count values (COUNT = -99999 is a placeholder for missing data)
lter.fish.sub.site <- subset(lter.fish.sub.site, COUNT != -99999)

# Calculate sheephead biomass using length-weight relationship
# This is an allometric equation: W = a * L^b where L = total length (cm)
# The coefficients 0.0144 and 3.04 are from published literature for this species
lter.fish.sub.site$biomass <- 0.0144 * (lter.fish.sub.site$SIZE ^ 3.04)
lter.fish.sub.site$biomass[is.na(lter.fish.sub.site$biomass)] <- 0

# Summarise: transect -> site -> MPA level
lter.fish.sum <- lter.fish.sub.site %>%
  group_by(YEAR, SITE, TRANSECT, status, CA_MPA_Name_Short,
           ChannelIsland, MPA_Start) %>%
  summarise_at(c("COUNT", "biomass"), sum)

lter.fish.ave <- lter.fish.sum %>%
  group_by(YEAR, SITE, status, CA_MPA_Name_Short, ChannelIsland, MPA_Start) %>%
  summarise_at(c("COUNT", "biomass"), mean)

lter.fish.max <- lter.fish.ave %>%
  group_by(status, CA_MPA_Name_Short, MPA_Start, YEAR) %>%
  summarise_at(c("COUNT", "biomass"), mean)

## --- Fish count proportions and log response ratios ---
lter.fish.max$Prop <- 0
lter.fish.max$PropCorr <- 0
lter.fish.max$taxon_name <- "SPUL"

x_fish <- unique(lter.fish.max$CA_MPA_Name_Short)

# Count proportions
for (i in seq_along(x_fish)) {
  idx <- which(lter.fish.max$CA_MPA_Name_Short == x_fish[i])
  m <- max(lter.fish.max$COUNT[idx])
  if (m > 0) {
    lter.fish.max$Prop[idx] <- lter.fish.max$COUNT[idx] / m
    lter.fish.max$PropCorr[idx] <- lter.fish.max$Prop[idx] + 0.01
  }
}

All.lter.fish.count.sub <- lter.fish.max[, colnames(lter.fish.max)[c(2, 4, 9, 1, 8)]]
Short.lter.fish.count <- All.lter.fish.count.sub %>% spread(status, PropCorr)
Short.lter.fish.count.diff <- calculate_log_response_ratio(Short.lter.fish.count)
Short.lter.fish.count.diff$resp <- "Den"

# Biomass proportions
for (i in seq_along(x_fish)) {
  idx <- which(lter.fish.max$CA_MPA_Name_Short == x_fish[i])
  m <- max(lter.fish.max$biomass[idx])
  if (m > 0) {
    lter.fish.max$Prop[idx] <- lter.fish.max$biomass[idx] / m
    lter.fish.max$PropCorr[idx] <- lter.fish.max$Prop[idx] + 0.01
  }
}

All.lter.spul.bio.sub <- lter.fish.max[, colnames(lter.fish.max)[c(2, 4, 9, 1, 8)]]
Short.lter.spul.bio <- All.lter.spul.bio.sub %>% spread(status, PropCorr)
Short.lter.spul.bio.diff <- calculate_log_response_ratio(Short.lter.spul.bio)
Short.lter.spul.bio.diff$resp <- "Bio"

# Combine fish results
LTER.fish <- rbind(Short.lter.fish.count.diff, Short.lter.spul.bio.diff)
colnames(LTER.fish)[2] <- "year"
colnames(LTER.fish)[3] <- "y"

# Assign time and BA labels using utility functions
LTER.fish <- assign_time_from_site_table(LTER.fish, Sites2)
LTER.fish <- merge(LTER.fish, sites.short.edit,
                   by.x = c("CA_MPA_Name_Short"),
                   by.y = c("CA_MPA_Name_Short"),
                   all = TRUE)
LTER.fish <- LTER.fish[complete.cases(LTER.fish$year), ]
LTER.fish$source <- "LTER"
LTER.fish <- assign_ba_from_site_table(LTER.fish, Sites2)

## --- Create LTER fish response dataframes ---
lter.fish.max$Bio <- lter.fish.max$biomass / 80
lter.fish.max$Den <- lter.fish.max$COUNT / 80

lter.fish.edit.den <- lter.fish.max[, colnames(lter.fish.max)[c(1:4, 11, 9)]]
LTER.fish.den <- lter.fish.edit.den %>% spread(status, Den)
LTER.fish.den$source <- "LTER"
LTER.fish.den.sub <- LTER.fish.den[complete.cases(LTER.fish.den), ]
LTER.fish.den.long <- gather(LTER.fish.den.sub, status, value, 5:6)
LTER.fish.den.long$resp <- "Den"

lter.fish.edit.bio <- lter.fish.max[, colnames(lter.fish.max)[c(1:4, 10, 9)]]
LTER.fish.bio <- lter.fish.edit.bio %>% spread(status, Bio)
LTER.fish.bio$source <- "LTER"
LTER.fish.bio.sub <- LTER.fish.bio[complete.cases(LTER.fish.bio), ]
LTER.fish.bio.long <- gather(LTER.fish.bio.sub, status, value, 5:6)
LTER.fish.bio.long$resp <- "Bio"

LTER.fish.resp <- rbind(LTER.fish.den.long, LTER.fish.bio.long)

###########################################################################################
## Summary of outputs available for downstream scripts
###########################################################################################
# LTER.join.ave    - Urchin density & biomass log response ratios with time, BA, site info
# LTER.lob         - Lobster density & biomass log response ratios with time, BA, site info
# Short.lter.macro - Macrocystis frond density & biomass log response ratios with time, BA
# LTER.fish        - Sheephead density & biomass log response ratios with time, BA
# LTER.resp        - Urchin response (den & bio, mpa vs reference, long format)
# LTER.lob.resp    - Lobster response (den & bio, mpa vs reference, long format)
# LTER.macro.resp  - Macrocystis response (den & bio, mpa vs reference, long format)
# LTER.fish.resp   - Fish response (den & bio, mpa vs reference, long format)
