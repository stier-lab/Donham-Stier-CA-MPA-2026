# =============================================================================
# 00c_analysis_constants.R
# =============================================================================
#
# PURPOSE:
#   Define named constants for all hardcoded values used across the analysis
#   pipeline. Centralizing these values improves maintainability and makes
#   assumptions explicit.
#
# USAGE:
#   source(here::here("code", "R", "00c_analysis_constants.R"))
#
# AUTHORS: Emily Donham & Adrian Stier
# PROJECT: CA MPA Kelp Forest pBACIPS Analysis
# =============================================================================


# =============================================================================
# SECTION 1: SURVEY AREA NORMALIZATION FACTORS
# =============================================================================
# Each monitoring program uses different sampling units. These constants
# convert counts to per-square-meter densities.

#' PISCO swath transect area (2m wide × 30m long = 60 m²)
PISCO_SWATH_AREA_M2 <- 60

#' KFM quadrat area (default 2 m²; data-repo code uses row-level `area`)
#' KFM urchin quads were 1 m² in 1982-1984, then 2 m² from 1985 onward.
#' Macro quads were 1, 2, or 10 m² (analysis uses 10 m² subset only).
#' See data-repo 05_kfm_processing.R for density normalization that uses
#' the row-level `area` column rather than this constant.
KFM_QUAD_AREA_M2 <- 2

#' LTER lobster plot area (300 m²)
LTER_LOBSTER_PLOT_AREA_M2 <- 300

#' LTER fish survey transect area (80 m²)
LTER_FISH_SURVEY_AREA_M2 <- 80


# =============================================================================
# SECTION 2: SIZE THRESHOLDS
# =============================================================================
# Minimum sizes for including organisms in analyses, based on field protocols.

#' PISCO minimum urchin size (mm test diameter)
#' RATIONALE: PISCO field protocols only record urchins >= 25mm because
#' smaller urchins are difficult to see and identify consistently.
PISCO_URCHIN_MIN_SIZE_MM <- 25


# =============================================================================
# SECTION 3: TEMPORAL PARAMETERS
# =============================================================================
# Year ranges and time-based filters for analysis.

#' Survey season months (May-October)
#' RATIONALE: Only summer months for consistency across survey methods
SURVEY_MONTH_START <- 5   # May
SURVEY_MONTH_END <- 10    # October

#' Year lobster size measurements began by campus
#' RATIONALE: Carapace length measurements started at different times
UCSB_LOBSTER_SIZE_START_YEAR <- 2010
VRG_LOBSTER_SIZE_START_YEAR <- 2011

#' KFM RDFC (Roving Diver Fish Count) survey methodology start year.
#' Per Davis et al. 1997 (NPS Kelp Forest Monitoring Handbook Vol 1, p.48),
#' RDFC was implemented in 1996. Earlier value of 2003 was incorrect;
#' corrected 2026-05-04. See data-repo 00c_analysis_constants.R for full note.
KFM_RDFC_SURVEY_START_YEAR <- 1996

#' Minimum years of data required for site inclusion
#' RATIONALE: Sites with fewer years provide unreliable trend estimates
MINIMUM_YEARS_OF_DATA <- 5

#' Standardized time point for effect size extraction (years post-MPA implementation)
#' RATIONALE: Using a fixed time point (t=11) rather than the maximum observed time
#' allows for comparable effect sizes across MPAs with different establishment dates.
#' We chose 11 years because it corresponds to the age of our youngest MPA in 2023
#' (MLPA South Coast MPAs implemented in 2012). This controls for differences in
#' effect size that could arise simply from longer protection duration.
#' See: Thiault et al. (2017) Methods in Ecology and Evolution for pBACIPS methodology.
EFFECT_SIZE_TIME_YEARS <- 11


# =============================================================================
# SECTION 4: BOOTSTRAP PARAMETERS
# =============================================================================
# Parameters for bootstrap resampling of biomass estimates.

#' Number of bootstrap resamples
N_BOOTSTRAP_RESAMPLES <- 1000

#' Random seeds for reproducibility
#' Sequential numbering (12345-12349) for easy identification
SEED_PISCO_LOBSTER <- 12345
SEED_PISCO_URCHIN <- 12346
SEED_KFM_URCHIN <- 12347
SEED_KFM_MACRO <- 12348
SEED_LTER_URCHIN <- 12349


# =============================================================================
# SECTION 5: STATISTICAL THRESHOLDS
# =============================================================================
# Significance levels and diagnostic thresholds.

#' Standard significance level (α = 0.05)
SIGNIFICANCE_ALPHA <- 0.05

#' Z-critical value for 95% confidence intervals
Z_CRITICAL_95_CI <- 1.96

#' Shapiro-Wilk test sample size bounds
MIN_OBS_SHAPIRO_TEST <- 3
MAX_OBS_SHAPIRO_TEST <- 5000

#' DHARMa simulation count for residual diagnostics
DHARMA_N_SIMULATIONS <- 250

#' Heteroscedasticity correlation thresholds
#' LM (linear models): more strict
#' NLS (non-linear models): more lenient due to inherent complexity
MAX_HETEROSCEDASTICITY_CORRELATION_LM <- 0.5
MAX_HETEROSCEDASTICITY_CORRELATION_NLS <- 0.6

#' Maximum outliers allowed in NLS diagnostics
MAX_OUTLIERS_NLS <- 1

#' NLS normality threshold (more lenient than default 0.05)
NLS_NORMALITY_THRESHOLD <- 0.01

#' Cook's distance threshold numerator for meta-analysis outlier detection
#' Standard rule: outliers have Cook's D > 4/n
COOKS_DISTANCE_NUMERATOR <- 4


# =============================================================================
# SECTION 6: ALLOMETRIC COEFFICIENTS
# =============================================================================
# Length-weight relationships from published literature.

#' California sheephead (Semicossyphus pulcher) length-weight relationship
#' biomass = a * length^b (where length is in cm and biomass in grams)
#' Source: Published literature
SHEEPHEAD_BIOMASS_ALLOMETRIC_A <- 0.0144
SHEEPHEAD_BIOMASS_ALLOMETRIC_B <- 3.04


# =============================================================================
# SECTION 7: TROPHIC LEVEL ASSIGNMENTS
# =============================================================================
# Maps species names (both full and abbreviated forms) to trophic levels.
# Used by temporal analysis (10), figures (11), and other scripts.

trophic_assignment <- c(
  "Panulirus interruptus" = "Predators",
  "Semicossyphus pulcher" = "Predators",
  "P. interruptus" = "Predators",
  "S. pulcher" = "Predators",
  "Strongylocentrotus purpuratus" = "Urchins",
  "Mesocentrotus franciscanus" = "Urchins",
  "S. purpuratus" = "Urchins",
  "M. franciscanus" = "Urchins",
  "Macrocystis pyrifera" = "Kelp",
  "M. pyrifera" = "Kelp"
)


# =============================================================================
# SECTION 8: RESILIENCE / HEATWAVE MODULE PARAMETERS (scripts 14-23)
# =============================================================================
# Single source of truth for the constants shared across the resilience suite
# (heatwave, resistance/recovery, replication, memory, stability). Previously each
# of scripts 14/19/21/22/23 re-typed these; centralizing prevents drift.

#' Marine heatwave event windows for the Santa Barbara Coastal region
#' (Hobday et al. 2016 heatwaveR catalog on daily NOAA OISST; see 14_heatwave_analysis.R).
MHW1_YEARS <- 2014:2016   # first heatwave ("the Blob" + 2015-16 El Nino); most extreme on record
MHW2_YEARS <- 2018:2020   # second heatwave ("Blob 2.0")

#' State-variable resistance/recovery windows, relative to a pre-heatwave baseline (script 22).
RESILIENCE_BASELINE_YEARS <- 2010:2013   # pre-heatwave baseline
RESILIENCE_AFTER_YEARS    <- 2017:2019   # immediate post-MHW1 recovery
RESILIENCE_RECENT_YEARS   <- 2020:2023   # recent / recovery-completeness window

#' Southern California Bight scope. Resilience analyses are restricted to reefs
#' south of Point Conception; this is the maximum latitude (deg N), asserted in
#' scripts 14 and 19 via assert_socal_scope() (01_utils.R).
SOCAL_MAX_LAT <- 34.45

#' Canonical resilience taxa (top predator -> primary producer), with abbreviated
#' display names, trophic role, and the response metric used (biomass for giant
#' kelp, density for the animals). Consumed by scripts 14/19/21/22/23.
RESILIENCE_TAXA <- c("Panulirus interruptus", "Semicossyphus pulcher",
                     "Strongylocentrotus purpuratus", "Mesocentrotus franciscanus",
                     "Macrocystis pyrifera")
RESILIENCE_TAXA_SHORT <- c(
  "Panulirus interruptus" = "P. interruptus", "Semicossyphus pulcher" = "S. pulcher",
  "Strongylocentrotus purpuratus" = "S. purpuratus",
  "Mesocentrotus franciscanus" = "M. franciscanus", "Macrocystis pyrifera" = "M. pyrifera")
RESILIENCE_TAXA_ROLE <- c(
  "Panulirus interruptus" = "Predator", "Semicossyphus pulcher" = "Predator",
  "Strongylocentrotus purpuratus" = "Herbivore", "Mesocentrotus franciscanus" = "Herbivore",
  "Macrocystis pyrifera" = "Producer")
RESILIENCE_RESP_OF <- c(
  "Panulirus interruptus" = "Den", "Semicossyphus pulcher" = "Den",
  "Strongylocentrotus purpuratus" = "Den", "Mesocentrotus franciscanus" = "Den",
  "Macrocystis pyrifera" = "Bio")


# =============================================================================
# Confirmation message
# =============================================================================
cat("Analysis constants loaded.\n")
