# =============================================================================
# 01_utils.R
# =============================================================================
#
# PURPOSE:
#   Define utility functions and constants used across all data processing
#   modules in the MPA kelp forest pBACIPS analysis pipeline.
#
# WHAT THIS SCRIPT DOES:
#   1. Defines functions to calculate time since MPA implementation
#   2. Defines functions to assign Before/After periods
#   3. Defines species-specific biomass conversion functions
#   4. Defines bootstrap resampling for biomass estimation
#   5. Defines functions for calculating response ratios and effect sizes
#   6. Defines constants for site/MPA exclusions
#
# USAGE:
#   source(here::here("code", "R", "01_utils.R"))
#   # Then use functions like:
#   df <- assign_time_from_site_table(df, site_table)
#   biomass <- bio_lobster(carapace_length_mm)
#
# AUTHORS: Emily Donham & Adrian Stier
# PROJECT: CA MPA Kelp Forest pBACIPS Analysis
# =============================================================================


# =============================================================================
# SECTION 1: TIME AND PERIOD ASSIGNMENT FUNCTIONS
# =============================================================================
# These functions handle the temporal aspects of the BACI design:
# - Calculating years since MPA implementation
# - Labeling observations as "Before" or "After" MPA establishment

#' Calculate time since MPA implementation
#'
#' For the pBACIPS design, we need to know how many years have passed since
#' each MPA was established. Years before implementation get time = 0.
#'
#' @param years Numeric vector of survey years (e.g., 2005, 2010, 2015)
#' @param mpa_start Integer, the year the MPA was legally established
#' @return Numeric vector: 0 for years <= mpa_start, otherwise (year - mpa_start)
#'
#' @examples
#' # MPA established in 2012
#' calculate_time_since_mpa(c(2010, 2012, 2014, 2016), mpa_start = 2012)
#' # Returns: c(0, 0, 2, 4)
calculate_time_since_mpa <- function(years, mpa_start) {
  ifelse(years <= mpa_start, 0, years - mpa_start)
}


#' Assign Before/After labels based on MPA implementation year
#'
#' Labels each observation as occurring "Before" or "After" the MPA was
#' established. The implementation year itself is considered "Before".
#'
#' @param years Numeric vector of survey years
#' @param mpa_start Integer, the year the MPA was established
#' @return Character vector: "Before" for years <= mpa_start, "After" otherwise
#'
#' @examples
#' assign_ba_period(c(2010, 2012, 2014), mpa_start = 2012)
#' # Returns: c("Before", "Before", "After")
assign_ba_period <- function(years, mpa_start) {
  ifelse(years <= mpa_start, "Before", "After")
}


#' Apply time since MPA to a dataframe using a site lookup table
#'
#' This function joins a dataframe with a site lookup table to get MPA start
#' years, then calculates time since implementation for each observation.
#'
#' @param df Dataframe with columns 'CA_MPA_Name_Short' and 'year'
#' @param site_table Dataframe with 'MPA_Start' column and a column matching mpa_col
#' @param mpa_col Name of column in site_table that matches CA_MPA_Name_Short
#' @return Modified dataframe with 'time' column added
assign_time_from_site_table <- function(df, site_table, mpa_col = "CA_MPA_Name_Short") {
  # Initialize time column
  df$time <- 0

  # Get unique MPAs in the data
  mpas <- unique(df$CA_MPA_Name_Short)

  # Loop through each MPA and calculate time
  for (mpa in mpas) {
    # Find rows in df for this MPA
    idx <- which(df$CA_MPA_Name_Short == mpa)

    # Find matching row in site table
    j <- which(site_table[[mpa_col]] == mpa)

    if (length(j) > 0) {
      # Get the MPA start year from the site table
      mpa_start <- site_table$MPA_Start[j[1]]

      # Calculate time since MPA for all rows of this MPA
      df$time[idx] <- calculate_time_since_mpa(df$year[idx], mpa_start)
    }
  }

  df
}


#' Apply Before/After labels to a dataframe using a site lookup table
#'
#' Similar to assign_time_from_site_table, but assigns "Before"/"After" labels.
#'
#' @param df Dataframe with columns CA_MPA_Name_Short and year
#' @param site_table Dataframe with MPA_Start column
#' @param mpa_col Name of column in site_table matching CA_MPA_Name_Short
#' @return Modified dataframe with 'BA' column added
assign_ba_from_site_table <- function(df, site_table, mpa_col = "CA_MPA_Name_Short") {
  df$BA <- NA

  mpas <- unique(df$CA_MPA_Name_Short)

  for (mpa in mpas) {
    idx <- which(df$CA_MPA_Name_Short == mpa)
    j <- which(site_table[[mpa_col]] == mpa)

    if (length(j) > 0) {
      mpa_start <- site_table$MPA_Start[j[1]]
      df$BA[idx] <- assign_ba_period(df$year[idx], mpa_start)
    }
  }

  df
}


# =============================================================================
# SECTION 2: BIOMASS CONVERSION FUNCTIONS
# =============================================================================
# Species-specific allometric relationships to convert size measurements to
# biomass (grams). These relationships come from published literature and
# LTER long-term research.
#
# General form: Biomass = a * (Size)^b
# where 'a' and 'b' are species-specific constants

#' Lobster carapace length to biomass (grams)
#'
#' Converts California spiny lobster (Panulirus interruptus) carapace length
#' to wet weight biomass using the allometric equation from LTER research.
#'
#' @param size Carapace length in MILLIMETERS
#' @return Biomass in grams
#'
#' @details
#' Equation: Biomass = 0.001352821 * CL^2.913963113
#' Source: SBC LTER long-term lobster monitoring
#'
#' @examples
#' bio_lobster(80)  # 80mm lobster = ~538 grams
bio_lobster <- function(size) {
  0.001352821 * (size) ^ 2.913963113
}


#' Red sea urchin test diameter to biomass (grams)
#'
#' Converts red urchin (Mesocentrotus franciscanus) test diameter to wet weight.
#'
#' @param size Test diameter in MILLIMETERS
#' @return Biomass in grams
#'
#' @details
#' Equation: Biomass = 0.00059 * TD^2.917
#' Source: SBC LTER urchin size-weight relationships
bio_redurch <- function(size) {
  0.00059 * (size) ^ 2.917
}


#' Purple sea urchin test diameter to biomass (grams)
#'
#' Converts purple urchin (Strongylocentrotus purpuratus) test diameter to wet weight.
#'
#' @param size Test diameter in MILLIMETERS
#' @return Biomass in grams
#'
#' @details
#' Equation: Biomass = 0.00059 * TD^2.870
#' Note: Slightly lower exponent than red urchin (smaller species)
bio_purpurch <- function(size) {
  0.00059 * (size) ^ 2.870
}


#' Giant kelp stipe count to biomass (grams per m2)
#'
#' Converts Macrocystis pyrifera frond/stipe density to biomass using the
#' average seasonal slope from LTER monitoring data.
#'
#' @param stipe Number of stipes (or stipe density per area)
#' @return Biomass in grams (or grams per m2 if input is density)
#'
#' @details
#' The conversion uses the average slope from May-October LTER data.
#' This seasonal average accounts for variation in frond weight through
#' the growing season.

# Average slope from May-October LTER data (kg dry weight per frond)
MACRO_AVE_SLOPE <- mean(c(0.10386, 0.10103, 0.09267, 0.09204, 0.08054, 0.08505))

bio_macro <- function(stipe) {
  # Multiply by 1000 to convert kg to grams
  stipe * MACRO_AVE_SLOPE * 1000
}


# =============================================================================
# SECTION 3: BOOTSTRAP BIOMASS ESTIMATION
# =============================================================================
# When we have count data but not individual sizes for every observation,
# we use bootstrap resampling from size frequency distributions to estimate
# biomass. This approach:
# 1. Finds matching size frequency records for the same MPA/year/site_status
# 2. Creates a "population" of sizes from the size frequency data
# 3. Resamples N individuals from this population (where N = observed count)
# 4. Calculates biomass for each resample
# 5. Returns the mean biomass across all resamples

#' Bootstrap biomass estimation from size frequency data
#'
#' Estimates total biomass by resampling from a size frequency distribution.
#' This is used when we know how many individuals were observed but don't have
#' individual size measurements for that specific transect.
#'
#' @param count Number of individuals observed on the transect
#' @param size_freq_indices Row indices in size_freq_table for matching records
#'        (same MPA, year, site status, species)
#' @param size_freq_table The size frequency dataframe (needs 'size' and 'count' columns)
#' @param biomass_fun Function to convert size to biomass (e.g., bio_lobster)
#' @param n_resamples Number of bootstrap iterations (default 1000)
#'
#' @return List with components:
#'   - biomass: Mean estimated biomass across bootstrap resamples
#'   - se: Standard error of the bootstrap distribution (NEW - for uncertainty)
#'   - count: Total count (echoed back for convenience)
#'
#' @details
#' STATISTICAL FIX (2026-02-03): Now returns bootstrap SE for uncertainty
#' quantification. The SE represents variability due to unknown individual
#' sizes when only aggregate counts are available.
#'
#' @examples
#' # See usage in 04_pisco_processing.R for VRG lobster bootstrap
bootstrap_biomass <- function(count, size_freq_indices, size_freq_table,
                              biomass_fun, n_resamples = 1000) {
  n <- count
  t2 <- size_freq_indices

  # Case 1: We have both counts and size frequency data
  if (n != 0 & length(t2) != 0) {
    # Build a "population" of sizes from the size frequency records
    # Each size is repeated by its frequency count
    a <- NULL
    for (j in seq_along(t2)) {
      reps <- rep(size_freq_table$size[t2[j]], size_freq_table$count[t2[j]])
      a <- c(a, reps)
    }
    a <- as.numeric(na.omit(a))

    # If no valid sizes found, return NA
    if (length(a) == 0) {
      return(list(biomass = NA, se = NA, count = n))
    }

    # Bootstrap resample: draw n individuals from the size population
    # Repeat n_resamples times
    s <- matrix(NA, nrow = n_resamples, ncol = n)
    for (k in 1:n_resamples) {
      s[k, ] <- sample(a, n, replace = TRUE)
    }

    # Convert all sizes to biomass
    s_bio <- apply(s, c(1, 2), biomass_fun)

    # Sum biomass for each resample
    bootstrap_sums <- rowSums(s_bio)

    # Calculate mean AND standard error from bootstrap distribution
    # SE quantifies uncertainty due to unknown individual sizes
    ave_biomass <- mean(bootstrap_sums)
    se_biomass <- sd(bootstrap_sums)

    return(list(biomass = ave_biomass, se = se_biomass, count = n))

  # Case 2: No individuals observed (count = 0)
  } else if (n == 0) {
    return(list(biomass = 0, se = 0, count = 0))

  # Case 3: Individuals observed but no size data available
  } else if (n != 0 & length(t2) == 0) {
    return(list(biomass = NA, se = NA, count = n))

  # Case 4: Fallback
  } else {
    return(list(biomass = 0, se = 0, count = 0))
  }
}


# =============================================================================
# SECTION 4: RESPONSE RATIO CALCULATIONS
# =============================================================================
# The pBACIPS analysis uses log response ratios to compare MPA vs reference
# sites. These functions standardize data and calculate the ratios.

#' Convert data to proportions of time series maximum
#'
#' Before calculating response ratios, we standardize each time series to
#' range from 0 to 1 (proportion of maximum). This makes different MPAs
#' and taxa comparable regardless of absolute abundance.
#'
#' @param df Dataframe with columns CA_MPA_Name_Short, y (taxa name), and value_col
#' @param value_col Name of the column containing the values to convert
#' @param correction_method Method for zero-correction. Options:
#'   - "adaptive" (default): Uses half of the minimum non-zero proportion
#'   - "fixed": Uses a fixed value of 0.01 (legacy behavior)
#'   - "none": No correction (PropCorr = Prop, may cause log(0) issues)
#'
#' @return Modified dataframe with two new columns:
#'   - Prop: Raw proportion (value / max)
#'   - PropCorr: Corrected proportion to avoid log(0) issues
#'
#' @details
#' STATISTICAL FIX (2026-02-03): The original fixed correction of +0.01 is
#' arbitrary and can bias results, especially for rare species. The new
#' "adaptive" method uses half the minimum observed non-zero proportion,
#' which scales appropriately with the data.
#'
#' Rationale for half-minimum: This is a standard approach in compositional
#' data analysis (Aitchison, 1986) that maintains the relative structure of
#' the data while avoiding log(0).
calculate_proportions <- function(df, value_col, correction_method = "adaptive") {
  df$Prop <- NA_real_
  df$PropCorr <- NA_real_

  mpas <- unique(df$CA_MPA_Name_Short)
  taxa <- unique(df$y)

  # Calculate proportions separately for each MPA x taxa combination
  for (mpa in mpas) {
    for (taxon in taxa) {
      idx <- which(df$CA_MPA_Name_Short == mpa & df$y == taxon)

      if (length(idx) > 0) {
        max_val <- max(df[[value_col]][idx], na.rm = TRUE)

        if (is.finite(max_val) && max_val > 0) {
          df$Prop[idx] <- df[[value_col]][idx] / max_val

          # Apply zero-correction based on method
          if (correction_method == "adaptive") {
            # Use half of the minimum non-zero proportion
            # This scales with the data and is less arbitrary
            non_zero_props <- df$Prop[idx][df$Prop[idx] > 0]
            if (length(non_zero_props) > 0) {
              correction <- min(non_zero_props, na.rm = TRUE) / 2
            } else {
              correction <- 0.01  # Fallback if all zeros
            }
            df$PropCorr[idx] <- df$Prop[idx] + correction
          } else if (correction_method == "fixed") {
            # Legacy behavior: fixed 0.01 correction
            df$PropCorr[idx] <- df$Prop[idx] + 0.01
          } else if (correction_method == "none") {
            # No correction (use with caution)
            df$PropCorr[idx] <- df$Prop[idx]
          } else {
            warning("Unknown correction_method '", correction_method,
                    "'. Using 'adaptive'.")
            non_zero_props <- df$Prop[idx][df$Prop[idx] > 0]
            correction <- ifelse(length(non_zero_props) > 0,
                                  min(non_zero_props, na.rm = TRUE) / 2, 0.01)
            df$PropCorr[idx] <- df$Prop[idx] + correction
          }
        } else {
          # If max is 0, NA, or -Inf, mark as NA
          df$Prop[idx] <- NA_real_
          df$PropCorr[idx] <- NA_real_
        }
      }
    }
  }

  df
}


#' Calculate log response ratio from wide-format MPA vs reference data
#'
#' The log response ratio (lnRR) is: ln(MPA / Reference)
#' Positive values indicate higher abundance inside the MPA.
#'
#' @param df Dataframe in wide format with 'mpa' and 'reference' columns
#'        (these should be the proportions, not raw values)
#'
#' @return Dataframe with two new columns:
#'   - Diff: Raw ratio (mpa / reference)
#'   - lnDiff: Natural log of the ratio
calculate_log_response_ratio <- function(df) {
  # Remove rows with any missing values
  df <- df[complete.cases(df), ]

  # Calculate ratio and log ratio
  df$Diff <- df$mpa / df$reference
  df$lnDiff <- log(df$Diff)

  df
}


# =============================================================================
# SECTION 5: EFFECT SIZE CALCULATION
# =============================================================================
# For the meta-analysis, we need effect sizes (mean differences) and their
# standard errors from the emmeans (estimated marginal means) output.
#
# STATISTICAL NOTE (2026-02-03):
# When before and after estimates come from the same model, they share error
# variance and are correlated. The correct SE calculation should use:
#   Var(A - B) = Var(A) + Var(B) - 2*Cov(A,B)
#
# The function below uses the independence assumption (omits covariance) which
# is CONSERVATIVE (overestimates SE) when covariance is positive, but may
# UNDERESTIMATE SE when covariance is negative. For predictions from the same
# regression at different x-values, covariance is typically positive, so the
# independence assumption gives slightly larger SEs than the true value.
#
# For more precise SE calculation, use calculate_effect_size_from_contrast()
# which uses emmeans::pairs() to properly handle the covariance.

#' Calculate effect size from emmeans estimates (independence assumption)
#'
#' Returns the mean difference on the log-ratio scale between two time points
#' (typically time=0 for Before and time=T for After).
#'
#' @param before_emmeans Tidy emmeans object for the "before" period (time=0)
#' @param after_emmeans Tidy emmeans object for the "after" period (time=T)
#'
#' @return List with components:
#'   - mean: Effect size (after - before)
#'   - SE: Standard error of the difference (assumes independence - see note)
#'   - SD: Pooled standard deviation (for legacy compatibility)
#'   - CI: 95% confidence interval half-width (SE * 1.96)
#'
#' @note This function assumes independence between before and after estimates.
#'   This is appropriate when estimates come from different models or datasets.
#'   When estimates come from the same model, consider using
#'   calculate_effect_size_from_contrast() instead.
calculate_effect_size <- function(before_emmeans, after_emmeans) {
  # Effect size is simply the difference in estimates
  mean_es <- after_emmeans$estimate[1] - before_emmeans$estimate[1]

  # SE of the difference assuming independence of the two predictions
  # Var(A - B) = Var(A) + Var(B) when independent
  # NOTE: This ignores covariance. See function documentation.
  se_diff <- sqrt(before_emmeans$std.error[1]^2 + after_emmeans$std.error[1]^2)

  # 95% CI half-width
  ci_diff <- se_diff * 1.96

  # Back-calculate SD for legacy compatibility (not used in meta-analysis)
  n_before <- before_emmeans$df[1] + 1
  n_after <- after_emmeans$df[1] + 1
  sd_before <- before_emmeans$std.error[1] * sqrt(n_before)
  sd_after <- after_emmeans$std.error[1] * sqrt(n_after)

  # Pooled standard deviation
  pSD <- sqrt(((n_before - 1) * sd_before^2 + (n_after - 1) * sd_after^2) /
                (n_before + n_after - 2))

  list(mean = mean_es, SE = se_diff, SD = pSD, CI = ci_diff)
}


#' Calculate effect size using proper contrast (handles covariance)
#'
#' This function calculates the effect size using emmeans::pairs() which
#' properly accounts for the covariance between estimates from the same model.
#'
#' @param model A fitted model object (lm, nls, etc.)
#' @param time_var Character name of the time variable in the model
#' @param time_before Numeric value for "before" time point (typically 0)
#' @param time_after Numeric value for "after" time point
#'
#' @return List with components:
#'   - mean: Effect size (after - before)
#'   - SE: Standard error of the contrast (properly accounts for covariance)
#'   - CI: 95% confidence interval half-width
#'   - df: Degrees of freedom for the contrast
#'
#' @details
#' STATISTICAL FIX (2026-02-03): This function uses emmeans::pairs() to compute
#' the contrast, which extracts the variance-covariance matrix from the model
#' and correctly calculates Var(A - B) = Var(A) + Var(B) - 2*Cov(A,B).
#'
#' @examples
#' # model <- lm(lnDiff ~ time, data = dat)
#' # es <- calculate_effect_size_from_contrast(model, "time", 0, 10)
calculate_effect_size_from_contrast <- function(model, time_var, time_before = 0, time_after) {
  # Create emmeans grid at before and after time points
  at_list <- list()
  at_list[[time_var]] <- c(time_before, time_after)

  em <- emmeans::emmeans(model, as.formula(paste("~", time_var)), at = at_list)

  # Use pairs() to compute the contrast with proper covariance handling
  contrast_result <- emmeans::pairs(em, reverse = TRUE)  # after - before
  contrast_summary <- summary(contrast_result)

  # Extract results
  mean_es <- contrast_summary$estimate[1]
  se_es <- contrast_summary$SE[1]
  df_es <- contrast_summary$df[1]
  ci_es <- se_es * qt(0.975, df_es)  # Use t-distribution for CI

 list(mean = mean_es, SE = se_es, CI = ci_es, df = df_es)
}


# =============================================================================
# SECTION 6: SPECIES NAME STANDARDIZATION
# =============================================================================
# Different data sources use different naming conventions. This function
# converts all names to full scientific names.

#' Standardize species names across programs
#'
#' Converts PISCO short codes and other abbreviations to full scientific names.
#'
#' @param names Character vector of species codes/names
#' @return Character vector with standardized scientific names
standardize_species_names <- function(names) {
  # PISCO short codes
  names[names == "PANINT"] <- "Panulirus interruptus"
  names[names == "SPUL"] <- "Semicossyphus pulcher"
  names[names == "MACPYRAD"] <- "Macrocystis pyrifera"
  names[names == "MESFRAAD"] <- "Mesocentrotus franciscanus"
  names[names == "STRPURAD"] <- "Strongylocentrotus purpuratus"

  names
}


# =============================================================================
# SECTION 7: SITE AND MPA EXCLUSION CONSTANTS
# =============================================================================
# These constants define which sites and MPAs to exclude from analysis.
# Exclusions are based on data quality, sampling issues, or MPA characteristics
# that make them unsuitable for the pBACIPS comparison.

#' Reference sites to exclude from comparisons
#'
#' These reference sites are excluded based on consultation with Scott Hamilton.
#' Reasons include: poor habitat match, unusual location, incomplete time series.
EXCLUDED_REFERENCE_SITES <- c(
  "ANACAPA_BLACK_SEA_BASS",
  "SMI_BAY_POINT",
  "SCAI_SHIP_ROCK",
  "SMI_CROOK_POINT_E",
  "SMI_CROOK_POINT_W",
  "SBI_SUTIL",
  "SCI_YELLOWBANKS_CEN",
  "SCI_YELLOWBANKS_W",
  "SRI_FORD_POINT",
  "SMI_PRINCE_ISLAND_CEN",
  "SMI_PRINCE_ISLAND_N",
  "SMI_HARE_ROCK",
  "SMI_TYLER_BIGHT_E",
  "SMI_TYLER_BIGHT_W",
  "SBI_GRAVEYARD_CANYON",
  "SBI_GRAVEYARD_CANYON_N"
)


#' MPAs to exclude from the main analysis
#'
#' Excluded for various reasons:
#' - Some allow certain types of take (not fully protected)
#' - Some have insufficient before-period data
#' - Some have weird geography or don't fit the study design
EXCLUDED_MPAS <- c(
  "Carrington Point SMR",
  "N/A",
  "Arrow Point to Lion Head Point SMCA",
  "Crystal Cove SMCA",
  "Laguna Beach SMR",
  "South La Jolla SMR",
  "Vandenberg SMR",
  "Point Conception SMR",
  "Anacapa Island SMR 1978",  # Old MPA, different protection regime

"Painted Cave SMCA",
  "Anacapa Island SMCA"
)


#' KFM sites to exclude
#'
#' These NPS Kelp Forest Monitoring sites are excluded because they:
#' - Are not proper control/impact sites
#' - Have short time series where longer series exist at the same location
EXCLUDED_KFM_SITES <- c(
  "a-k-21", "a-k-05", "a-k-07", "a-k-26", "a-k-27", "a-k-28",
  "a-k-29", "a-k-30", "a-k-36", "a-k-37", "a-k-12", "a-k-13",
  "a-k-31", "a-k-35"
)


#' MPAs that only have sheephead data
#'
#' These MPAs only have sufficient data for sheephead (S. pulcher).
#' They are excluded from the main analysis but included when analyzing
#' sheephead specifically.
SHEEPHEAD_ONLY_MPAS <- c(
  "Blue Cavern Onshore SMCA",
  "Painted Cave SMCA",
  "Dana Point SMCA",
  "Farnsworth Onshore SMCA",
  "Point Dume SMR",
  "Cat Harbor SMCA",
  "Swamis SMCA",
  "Anacapa Island SMCA",
  "Long Point SMR",
  "Point Dume SMCA",
  "Santa Barbara Island SMR"
)


# =============================================================================
# Confirmation message
# =============================================================================
cat("Utility functions and constants loaded.\n")
