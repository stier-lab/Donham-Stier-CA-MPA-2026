# =============================================================================
# 07_combine_data.R
# =============================================================================
#
# PURPOSE:
#   Combine response ratio and raw response data from all monitoring programs
#   (LTER, PISCO, KFM) into unified datasets for the meta-analysis.
#
# WHAT THIS SCRIPT DOES:
#   1. Combines response ratio data from all programs into All.RR
#   2. Applies MPA exclusion rules (removes sites with insufficient data)
#   3. Handles special cases for sheephead-only MPAs
#   4. Standardizes species names to full scientific names
#   5. Standardizes source names (e.g., "LTER macro surveys" -> "LTER")
#   6. Combines raw response data (density/biomass) into All.Resp
#   7. Assigns Before/After labels and time since MPA to all data
#   8. Exports summary statistics
#
# EXCLUSION LOGIC:
#   Some MPAs only have data for certain taxa. We:
#   - Exclude MPAs from SHEEPHEAD_ONLY_MPAS list from main analysis
#   - Re-include those MPAs specifically for sheephead analysis
#   - This prevents sheephead from being over-represented in cross-taxa summaries
#
# INPUTS:
#   Response ratio objects:
#   - LTER.join.ave, Swath.join.sub (PISCO), KFM.join.ave
#   - LTER.lob, Short.lter.macro, LTER.fish, kfm.fish
#
#   Raw response objects:
#   - LTER.resp, KFM.resp, PISCO.resp
#   - KFM.fish.den.long, LTER.fish.resp, LTER.macro.resp, LTER.lob.resp
#
#   Site metadata: Site, Sites2
#
# OUTPUTS:
#   - All.RR.sub.trans: Combined response ratios with standardized names
#   - All.Resp.sub: Combined raw data with time and BA columns
#   - average_responses.csv: Summary statistics (written to project root)
#
# DEPENDENCIES:
#   Requires 00-06b scripts to be sourced first
#
# AUTHORS: Emily Donham & Adrian Stier
# PROJECT: CA MPA Kelp Forest pBACIPS Analysis
# =============================================================================

####################################################################################################
## Combine all response ratio data #################################################################
####################################################################################################
# At this point we have response ratio dataframes from each monitoring program:
# - LTER.join.ave (LTER urchins), Swath.join.sub (PISCO), KFM.join.ave (KFM urchins/kelp)
# - LTER.lob (LTER lobsters), Short.lter.macro (LTER kelp), LTER.fish (LTER sheephead)
# - kfm.fish (KFM sheephead)
# Now we combine them into one master dataset for the meta-analysis.

# Ensure consistent column structure before binding
# Different sources may have columns in different orders; we standardize them here
#
# TECHNICAL DEBT: Column indices are fragile - if upstream processing changes column
# order, these selections may break silently. Future improvement: convert to named
# column selection using dplyr::select(). The target columns for rbind are:
#   CA_MPA_Name_Short, year, y, lnDiff, mpa, reference, Diff, resp, time, type,
#   Location, Hectares, source, BA
#
# Current indices drop column 4 from each dataframe (typically a duplicate or unused column)
Swath.join.sub <- Swath.join.sub[, colnames(Swath.join.sub)[c(1:3, 5:15)]]
kfm.fish <- kfm.fish[, colnames(kfm.fish)[c(1, 3:15)]]
KFM.join.ave <- KFM.join.ave[, colnames(KFM.join.ave)[c(1:3, 5:15)]]

# Combine all response ratio datasets using rbind() (row bind)
# This stacks all dataframes vertically into one master dataframe
All.RR <- rbind(LTER.join.ave, Swath.join.sub, KFM.join.ave, LTER.lob,
                Short.lter.macro, LTER.fish, kfm.fish)

####################################################################################################
## Filter excluded MPAs with special handling for sheephead-only sites #############################
####################################################################################################
# IMPORTANT: Some MPAs only have data for certain taxa (usually just sheephead from PISCO fish surveys).
# If we include these MPAs in cross-taxa summaries, sheephead would be over-represented.
# SOLUTION: Remove these "sheephead-only" MPAs from the main analysis, then selectively
# add them back ONLY for sheephead-specific analyses.
#
# SHEEPHEAD_ONLY_MPAS is a constant defined in 01_utils.R containing these MPA names.

# Remove sheephead-only MPAs from the main dataset
# The ! negates the %in% check, so we KEEP rows that are NOT in the exclusion list
All.RR.sub <- subset(All.RR, !(CA_MPA_Name_Short %in% SHEEPHEAD_ONLY_MPAS))

# Re-include sheephead-only MPAs for SPUL (sheephead) taxa only
# SPUL = species code for Semicossyphus pulcher (California sheephead)
# This selective re-inclusion lets us analyze sheephead at these MPAs without
# biasing the multi-species summaries
All.spul <- subset(All.RR,
  (CA_MPA_Name_Short == "Blue Cavern Onshore SMCA" & y == "SPUL") |
  (CA_MPA_Name_Short == "Dana Point SMCA" & y == "SPUL") |
  (CA_MPA_Name_Short == "Farnsworth Onshore SMCA" & y == "SPUL") |
  (CA_MPA_Name_Short == "Swamis SMCA" & y == "SPUL") |
  (CA_MPA_Name_Short == "Cat Harbor SMCA" & y == "SPUL") |
  (CA_MPA_Name_Short == "Long Point SMR" & y == "SPUL") |
  (CA_MPA_Name_Short == "Santa Barbara Island SMR" & y == "SPUL" & source == "PISCO") |
  (CA_MPA_Name_Short == "Santa Barbara Island SMR" & source == "KFM")
)

# Add the sheephead-specific rows back to the main dataset
All.RR.sub <- rbind(All.RR.sub, All.spul)

####################################################################################################
## Duplicate detection for response ratio data #####################################################
####################################################################################################
# Check for duplicates that might arise from overlapping data sources or processing errors.
# A duplicate is defined as same MPA + year + taxa + response type + source combination.

dup_key_cols <- c("CA_MPA_Name_Short", "year", "y", "resp", "source")
dup_check <- duplicated(All.RR.sub[, dup_key_cols])
n_dups <- sum(dup_check)

if (n_dups > 0) {
  warning("Found ", n_dups, " duplicate rows in All.RR.sub (by MPA/year/taxa/resp/source). ",
          "Review data processing pipeline for overlapping records.")
  cat("\nDuplicate rows detected in All.RR.sub:", n_dups, "\n")
  cat("First 5 duplicates:\n")
  print(head(All.RR.sub[dup_check, dup_key_cols], 5))
} else {
  cat("Duplicate check passed: No duplicate MPA/year/taxa/resp/source combinations in All.RR.sub.\n")
}

####################################################################################################
## Standardize species and source names ############################################################
####################################################################################################
# Different monitoring programs use different species codes:
# - PISCO: STRPURAD, MESFRAAD, MACPYRAD, PANINT, SPUL
# - KFM/LTER: Full names like "Strongylocentrotus purpuratus"
# For consistent analysis and plotting, we standardize to full scientific names.

All.RR.sub.trans <- All.RR.sub
All.RR.sub.trans$y <- as.character(All.RR.sub.trans$y)

# standardize_species_names() from 01_utils.R converts codes to full names:
# STRPURAD -> Strongylocentrotus purpuratus (purple urchin)
# MESFRAAD -> Mesocentrotus franciscanus (red urchin)
# MACPYRAD -> Macrocystis pyrifera (giant kelp)
# PANINT -> Panulirus interruptus (California spiny lobster)
# SPUL -> Semicossyphus pulcher (California sheephead)
All.RR.sub.trans$y <- standardize_species_names(All.RR.sub.trans$y)

# Standardize source names to simple program names
# Some datasets have more descriptive source names that we simplify
All.RR.sub.trans$source[All.RR.sub.trans$source == "LTER macro surveys"] <- "LTER"
All.RR.sub.trans$source[All.RR.sub.trans$source == "LTER lob surveys"] <- "LTER"

####################################################################################################
## Combine all raw response data ###################################################################
####################################################################################################

# Standardize species names in individual response datasets before combining
PISCO.resp$y <- standardize_species_names(as.character(PISCO.resp$y))
colnames(PISCO.resp)[3] <- "taxon_name"

KFM.fish.den.long <- KFM.fish.den.long[, colnames(KFM.fish.den.long)[c(1, 3:8)]]
KFM.fish.den.long$taxon_name[KFM.fish.den.long$taxon_name == "SPUL"] <- "Semicossyphus pulcher"

KFM.resp <- KFM.resp[, colnames(KFM.resp)[c(1, 3:8)]]
colnames(KFM.resp)[3] <- "taxon_name"
KFM.resp$taxon_name <- standardize_species_names(as.character(KFM.resp$taxon_name))

LTER.fish.resp <- LTER.fish.resp[, colnames(LTER.fish.resp)[c(1, 3:8)]]
LTER.fish.resp$taxon_name[LTER.fish.resp$taxon_name == "SPUL"] <- "Semicossyphus pulcher"
colnames(LTER.fish.resp)[2] <- "year"

LTER.macro.resp <- LTER.macro.resp[, colnames(LTER.macro.resp)[c(1, 3:8)]]
colnames(LTER.macro.resp)[2] <- "year"

LTER.lob.resp <- LTER.lob.resp[, colnames(LTER.lob.resp)[c(1, 3:8)]]
colnames(LTER.lob.resp)[2] <- "year"

LTER.resp <- LTER.resp[, colnames(LTER.resp)[c(1, 3:8)]]

# Combine all raw response datasets
All.Resp <- rbind(KFM.fish.den.long, LTER.fish.resp, LTER.macro.resp,
                  LTER.lob.resp, LTER.resp, KFM.resp, PISCO.resp)

####################################################################################################
## Apply same exclusion logic to raw response data #################################################
####################################################################################################

All.Resp.sub <- subset(All.Resp, !(CA_MPA_Name_Short %in% SHEEPHEAD_ONLY_MPAS))

# Re-include sheephead-only MPAs for Semicossyphus pulcher
All.Resp.spul <- subset(All.Resp,
  (CA_MPA_Name_Short == "Blue Cavern Onshore SMCA" & taxon_name == "Semicossyphus pulcher") |
  (CA_MPA_Name_Short == "Dana Point SMCA" & taxon_name == "Semicossyphus pulcher") |
  (CA_MPA_Name_Short == "Farnsworth Onshore SMCA" & taxon_name == "Semicossyphus pulcher") |
  (CA_MPA_Name_Short == "Swamis SMCA" & taxon_name == "Semicossyphus pulcher") |
  (CA_MPA_Name_Short == "Cat Harbor SMCA" & taxon_name == "Semicossyphus pulcher") |
  (CA_MPA_Name_Short == "Long Point SMR" & taxon_name == "Semicossyphus pulcher") |
  (CA_MPA_Name_Short == "Santa Barbara Island SMR" & taxon_name == "Semicossyphus pulcher" & source == "PISCO") |
  (CA_MPA_Name_Short == "Santa Barbara Island SMR" & source == "KFM")
)

All.Resp.sub <- rbind(All.Resp.sub, All.Resp.spul)

####################################################################################################
## Duplicate detection for raw response data #######################################################
####################################################################################################
# Check for duplicates in raw response data.

dup_key_resp <- c("CA_MPA_Name_Short", "year", "taxon_name", "resp", "source", "status")
dup_check_resp <- duplicated(All.Resp.sub[, dup_key_resp])
n_dups_resp <- sum(dup_check_resp)

if (n_dups_resp > 0) {
  warning("Found ", n_dups_resp, " duplicate rows in All.Resp.sub. ",
          "Review data processing pipeline for overlapping records.")
  cat("\nDuplicate rows detected in All.Resp.sub:", n_dups_resp, "\n")
} else {
  cat("Duplicate check passed: No duplicate rows in All.Resp.sub.\n")
}

####################################################################################################
## Assign Before/After and time to raw response data ###############################################
####################################################################################################

# Use vectorized utility functions instead of year-by-year if-else blocks
All.Resp.sub <- assign_ba_from_site_table(All.Resp.sub, Site)
All.Resp.sub <- assign_time_from_site_table(All.Resp.sub, Site)

# Keep only essential columns
All.Resp.sub <- All.Resp.sub[, colnames(All.Resp.sub)[c(1:7, 10:11)]]

####################################################################################################
## Calculate average responses and write summary ###################################################
####################################################################################################

AveResponse <- summarySE(data = All.Resp.sub, measurevar = "value",
                         groupvars = c("taxon_name", "source", "resp"))
write_csv(AveResponse, here::here("data", "average_responses.csv"))

# Clean up intermediate objects
rm(All.RR, All.spul, All.Resp, All.Resp.spul)
