# =============================================================================
# 08_effect_sizes.R
# =============================================================================
#
# PURPOSE:
#   Run pBACIPS analyses and calculate effect sizes for all combinations of
#   data source, taxa, response type, and MPA.
#
# WHAT THIS SCRIPT DOES:
#   1. For each data subset (taxa x source x response x MPA):
#      a. Tests for linear trends in the "Before" period (violates BACI assumptions)
#      b. Runs pBACIPS model selection (step, linear, asymptotic, sigmoid)
#      c. Calculates effect size (mean difference on log-ratio scale)
#      d. Records the best-fitting model and effect size with SE
#   2. Handles three types of analyses:
#      - pBACIPS: Full before-after comparison with model selection
#      - BACI: Simple step-change model (for sites with linear before period)
#      - CI: Control-Impact only (for sites without before period data)
#   3. Compiles all effect sizes into SumStats dataframe
#   4. Filters to final publication-ready subset (SumStats.Final)
#
# EFFECT SIZE INTERPRETATION:
#   - Effect size = change in log response ratio after MPA implementation
#   - Positive values: MPA site increased relative to reference
#   - Negative values: MPA site decreased relative to reference
#   - SE allows calculation of confidence intervals and meta-analysis weights
#
# STATISTICAL NOTES:
#   - Linear before period: If there's a significant trend before MPA,
#     we cannot attribute post-MPA changes to protection (confounded)
#   - pBACIPS vs BACI: pBACIPS allows progressive change; BACI assumes step
#   - CI analysis: Used when no before-period data exists
#
# INPUTS:
#   - All.RR.sub.trans: Combined response ratio data
#   - Site, Sites2: Site metadata with MPA start years
#   - ProgressiveChangeBACIPS function from 02_pBACIPS_function.R
#
# OUTPUTS:
#   - SumStats: All effect size estimates with metadata
#   - SumStats.Final: Filtered subset for publication
#     (excludes MPAs with data quality issues)
#
# DEPENDENCIES:
#   Requires 00-07 scripts to be sourced first
#
# AUTHORS: Emily Donham & Adrian Stier
# PROJECT: CA MPA Kelp Forest pBACIPS Analysis
# =============================================================================

####################################################################################################
## Source pBACIPS function and set up ##############################################################
####################################################################################################

source(here::here("code", "R", "02_pBACIPS_function.R"))

# Color palette for plots
cols <- c("Sheephead" = "#377eb8", "Kelp" = "#4daf4a", "Lobs" = "#ff7f00",
          "Purps" = "#984ea3", "Reds" = "#e41a1c")

# Rename Site columns for MPA_implement reference
MPA_implement <- Site
names(MPA_implement)[names(MPA_implement) == "Site"] <- "Site_ID"
names(MPA_implement)[names(MPA_implement) == "CA_MPA_Name_Short"] <- "Site"

####################################################################################################
## Initialize SumStats dataframe ###################################################################
####################################################################################################

SumStats <- data.frame(
  Taxa = NA, MPA = NA, Mean = NA, SE = NA, SD = NA, CI = NA,
  Model = NA, Source = NA, Resp = NA, BA = NA, Primary = NA,
  Type = NA, LinearBefore = NA
)

####################################################################################################
## Helper functions for repeated analysis patterns #################################################
####################################################################################################

#' Add time.model and time.true columns to a subset of RR data
#'
#' time.model: 0 for Before period, sequential starting at 0 for After
#' time.true: sequential index across all time points per MPA
#'
#' @param dat Dataframe with CA_MPA_Name_Short, BA, year columns
#' @return Modified dataframe with time.model and time.true columns
add_time_columns <- function(dat) {
  dat$time.model <- NA
  dat$time.true <- NA
  dat <- dat[order(as.numeric(dat$year), decreasing = FALSE), ]
  mpas <- unique(dat$CA_MPA_Name_Short)
  for (i in seq_along(mpas)) {
    idx <- which(dat$CA_MPA_Name_Short == mpas[i])
    n_before <- length(which(dat$BA[idx] == "Before"))
    n_after <- length(which(dat$BA[idx] == "After"))
    dat$time.model[idx] <- c(rep(0, n_before), seq(0, n_after - 1))
    dat$time.true[idx] <- seq(1, length(idx))
  }
  dat
}

#' Test for linear trend in the Before period for each MPA
#'
#' METHODOLOGICAL NOTE:
#' The BACI/pBACIPS design assumes that MPA and reference sites have parallel
#' trajectories before intervention. A significant pre-existing trend in the
#' log response ratio violates this assumption and can confound estimates.
#'
#' DECISION RATIONALE:
#' Sites with significant (p < 0.05) linear trends in the before period are
#' excluded from pBACIPS analysis and analyzed separately as CI (Control-Impact)
#' sites, which don't require the parallel trend assumption.
#'
#' LIMITATIONS AND CAVEATS:
#' 1. Low statistical power: With only 3-5 before-period points, power to detect
#'    trends is low. Sites may pass this test simply due to noise.
#' 2. Selection bias: This approach may favor inclusion of noisy sites over sites
#'    with detectable patterns, potentially biasing effect estimates.
#' 3. Arbitrary threshold: The p < 0.05 cutoff is conventional but arbitrary.
#'
#' POTENTIAL SENSITIVITY ANALYSES:
#' - Compare results using different p-value thresholds (0.01, 0.10)
#' - Include all sites and use detrending approaches
#' - Report results both with and without exclusions
#'
#' @param dat Dataframe with CA_MPA_Name_Short, lnDiff, time.true, BA columns
#' @return Dataframe with pInt (intercept p-value), pSlope (slope p-value),
#'         and site columns. Sites with pSlope < 0.05 have significant before trends.
test_linear_before <- function(dat) {
  mpas <- unique(dat$CA_MPA_Name_Short)
  LinearBefore <- data.frame(pInt = NA, pSlope = NA)
  for (j in mpas) {
    temp <- subset(dat, CA_MPA_Name_Short == j &  BA == "Before",
                   select = c(CA_MPA_Name_Short, mpa, reference, lnDiff, time.true, time.model, BA))
    lmBefore <- lm(lnDiff ~ time.true, data = temp)
    pl <- summary(lmBefore)$coefficients[, 4]
    LinearBefore <- rbind(LinearBefore, pl)
  }
  LinearBefore <- na.omit(LinearBefore)
  LinearBefore$site <- mpas
  LinearBefore
}

#' Run ProgressiveChangeBACIPS across MPAs and collect model weights/p-values
#'
#' @param dat Dataframe with lnDiff, lnmpa, lnreference, time.true, time.model
#' @return List with likelihood.poly, pvals.poly, Norm.poly dataframes
run_pbacips_loop <- function(dat) {
  mpas <- unique(dat$CA_MPA_Name_Short)
  likelihood.poly <- data.frame(matrix(ncol = 4, nrow = 0))
  pvals.poly <- data.frame(matrix(ncol = 9, nrow = 0))
  Norm.poly <- data.frame(matrix(ncol = 9, nrow = 0))

  dat$lnmpa <- log(dat$mpa)
  dat$lnreference <- log(dat$reference)

  for (j in mpas) {
    temp <- subset(dat, CA_MPA_Name_Short == j,
                   select = c(CA_MPA_Name_Short, lnmpa, lnreference, lnDiff, time.true, time.model))
    result <- ProgressiveChangeBACIPS(
      control = temp$lnreference,
      impact = temp$lnmpa,
      time.true = temp$time.true,
      time.model = temp$time.model
    )
    ps <- summary(result$step)[[1]][["Pr(>F)"]][1]
    pl <- summary(result$linear)$coefficients[, 4]
    psi <- summary(result$sigmoid)$coefficients[, 4]
    pa <- summary(result$asymptotic)$coefficients[, 4]
    pval <- c(ps, pl, psi, pa)
    likelihood.poly <- rbind(likelihood.poly, result$weights)
    pvals.poly <- rbind(pvals.poly, pval)
    t_test <- shapiro.test(temp$lnDiff)
    Norm.poly <- rbind(Norm.poly, t_test$p.value)
  }

  colnames(likelihood.poly) <- c("Step", "Linear", "Asymptotic", "Sigmoid")
  colnames(pvals.poly) <- c("Step", "Linear", "Asymptotic M", "Asymptotic B", "Asymptotic L",
                             "Sigmoid M", "Sigmoid B", "Sigmoid L", "Sigmoid K")
  likelihood.poly <- cbind(likelihood.poly, mpas)

  list(likelihood = likelihood.poly, pvals = pvals.poly, normality = Norm.poly)
}

#' Run CI analysis: calculate both linear and mean effect sizes
#'
#' For each MPA in the data, fits lnDiff ~ time, calculates effect size at time 0 vs 11,
#' and also calculates the simple mean effect size. Returns two rows per MPA.
#'
#' @param dat Dataframe subset with lnDiff, time columns
#' @param taxa_name Character, species name for SumStats
#' @param source_name Character, data source for SumStats
#' @param resp_name Character, response type ("Den" or "Bio")
#' @param time_var Character, name of the time variable to use (default "time")
#' @return Dataframe with two rows per MPA (linear and mean effect sizes)
run_ci_analysis <- function(dat, taxa_name, source_name, resp_name, time_var = "time") {
  mpas <- unique(dat$CA_MPA_Name_Short)
  rows <- list()
  for (i in seq_along(mpas)) {
    idx <- which(dat$CA_MPA_Name_Short == mpas[i])
    sub_dat <- dat[idx, ]

    # Fit linear model
    formula_str <- as.formula(paste("lnDiff ~", time_var))
    Lm.Ab <- lm(formula_str, data = sub_dat)
    p <- summary(Lm.Ab)$coefficients[2, 4]

    # Calculate linear effect size using emmeans at time 0 and 11
    before <- tidy(emmeans(Lm.Ab, as.formula(paste("~", time_var)),
                           at = setNames(list(0), time_var)))
    after <- tidy(emmeans(Lm.Ab, as.formula(paste("~", time_var)),
                          at = setNames(list(11), time_var)))

    es <- calculate_effect_size(before, after)
    mean_val <- es$mean

    # Calculate simple mean CI
    CP.mean <- summarySE(sub_dat, measurevar = "lnDiff")

    # Determine Primary flag based on p-value (linear sig -> linear is primary, else mean is primary)
    if (p <= 0.05) {
      rows[[length(rows) + 1]] <- data.frame(
        Taxa = taxa_name, MPA = mpas[i], Mean = mean_val, SE = es$SE, SD = es$SD, CI = es$CI,
        Model = "Linear", Source = source_name, Resp = resp_name, BA = "N", Primary = "Y",
        Type = "CI", LinearBefore = "NA", stringsAsFactors = FALSE
      )
    } else {
      rows[[length(rows) + 1]] <- data.frame(
        Taxa = taxa_name, MPA = mpas[i], Mean = mean_val, SE = es$SE, SD = es$SD, CI = es$CI,
        Model = "Linear", Source = source_name, Resp = resp_name, BA = "N", Primary = "N",
        Type = "CI", LinearBefore = "NA", stringsAsFactors = FALSE
      )
    }

    if (p <= 0.05) {
      rows[[length(rows) + 1]] <- data.frame(
        Taxa = taxa_name, MPA = mpas[i], Mean = CP.mean$lnDiff, SE = CP.mean$se,
        SD = CP.mean$sd, CI = CP.mean$ci, Model = "Mean", Source = source_name, Resp = resp_name,
        BA = "N", Primary = "N", Type = "CI", LinearBefore = "NA", stringsAsFactors = FALSE
      )
    } else {
      rows[[length(rows) + 1]] <- data.frame(
        Taxa = taxa_name, MPA = mpas[i], Mean = CP.mean$lnDiff, SE = CP.mean$se,
        SD = CP.mean$sd, CI = CP.mean$ci, Model = "Mean", Source = source_name, Resp = resp_name,
        BA = "N", Primary = "Y", Type = "CI", LinearBefore = "NA", stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

#' Add a single step-model (BACI) effect size row
#'
#' @param dat Dataframe for a single MPA with lnDiff and BA columns
#' @param taxa_name Character species name
#' @param mpa_name Character MPA name
#' @param source_name Character data source
#' @param resp_name Character response type
#' @return One-row dataframe with the same columns as SumStats
add_step_effect_size <- function(dat, taxa_name, mpa_name, source_name, resp_name) {
  mod1 <- lm(data = dat, lnDiff ~ BA)
  ba_emmeans <- tidy(emmeans(mod1, ~ BA))
  ba_emmeans$stdev <- ba_emmeans$std.error * sqrt(ba_emmeans$df + 2)

  pSD <- sqrt((((ba_emmeans$df[1] + 1) * (ba_emmeans$stdev[1]^2)) +
                 ((ba_emmeans$df[2] + 1) * (ba_emmeans$stdev[2]^2))) /
                (ba_emmeans$df[1] + 2))
  pSE <- pSD * sqrt((1 / (ba_emmeans$df[1] + 2)) + (1 / (ba_emmeans$df[2] + 2)))
  pCI <- pSE * 1.96
  mean_es <- ba_emmeans$estimate[1] - ba_emmeans$estimate[2]

  data.frame(
    Taxa = taxa_name, MPA = mpa_name, Mean = mean_es, SE = pSE, SD = pSD, CI = pCI,
    Model = "Step", Source = source_name, Resp = resp_name, BA = "Y", Primary = "Y",
    Type = "BACI", LinearBefore = "N", stringsAsFactors = FALSE
  )
}

#' Add a single linear-model (pBACIPS) effect size row
#'
#' @param dat Dataframe for a single MPA with lnDiff and time.model columns
#' @param taxa_name Character species name
#' @param mpa_name Character MPA name
#' @param source_name Character data source
#' @param resp_name Character response type
#' @param time_var Character, time variable name (default "time.model")
#' @param time_after Numeric, the time point for "after" prediction (default NULL = use max observed)
#' @return One-row dataframe with the same columns as SumStats
add_linear_effect_size <- function(dat, taxa_name, mpa_name, source_name, resp_name,
                                    time_var = "time.model", time_after = NULL) {
  # Use the maximum observed time value if not specified
  if (is.null(time_after)) {
    time_after <- max(dat[[time_var]], na.rm = TRUE)
  }
  formula_str <- as.formula(paste("lnDiff ~", time_var))
  mod1 <- lm(formula_str, data = dat)

  before <- tidy(emmeans(mod1, as.formula(paste("~", time_var)),
                         at = setNames(list(0), time_var)))
  after <- tidy(emmeans(mod1, as.formula(paste("~", time_var)),
                        at = setNames(list(time_after), time_var)))

  es <- calculate_effect_size(before, after)

  data.frame(
    Taxa = taxa_name, MPA = mpa_name, Mean = es$mean, SE = es$SE, SD = es$SD, CI = es$CI,
    Model = "Linear", Source = source_name, Resp = resp_name, BA = "Y", Primary = "Y",
    Type = "pBACIPS", LinearBefore = "N", stringsAsFactors = FALSE
  )
}

#' Add a mean-only effect size row (for cases without enough data for pBACIPS)
#'
#' @param dat Dataframe for a single MPA
#' @param taxa_name Character species name
#' @param mpa_name Character MPA name
#' @param source_name Character data source
#' @param resp_name Character response type
#' @return One-row dataframe with the same columns as SumStats
add_mean_effect_size <- function(dat, taxa_name, mpa_name, source_name, resp_name) {
  CP.mean <- summarySE(dat, measurevar = "lnDiff")
  data.frame(
    Taxa = taxa_name, MPA = mpa_name, Mean = CP.mean$lnDiff, SE = CP.mean$se,
    SD = CP.mean$sd, CI = CP.mean$ci, Model = "Mean", Source = source_name, Resp = resp_name,
    BA = "N", Primary = "Y", Type = "CI", LinearBefore = "NA", stringsAsFactors = FALSE
  )
}

####################################################################################################
## LTER ANALYSES ###################################################################################
####################################################################################################

#----- S. purpuratus Density (LTER) -----#
LTER.purps <- subset(All.RR.sub.trans, source == "LTER" &
                       y == "Strongylocentrotus purpuratus" & resp == "Den")
LTER.purps <- add_time_columns(LTER.purps)
LinearBefore.purps <- test_linear_before(LTER.purps)
# Note: Linear in before period -- excluded from pBACIPS (conservative: exclude den & bio)

#----- M. franciscanus Density (LTER) -----#
LTER.reds <- subset(All.RR.sub.trans, source == "LTER" &
                      y == "Mesocentrotus franciscanus" & resp == "Den")
LTER.reds <- add_time_columns(LTER.reds)
LinearBefore.reds <- test_linear_before(LTER.reds)
pbacips.reds.den <- run_pbacips_loop(LTER.reds)

# Only Naples (Campus Point is linear in before); step was best fit model
LTER.redsNap <- subset(LTER.reds, CA_MPA_Name_Short == "Naples SMCA")
SumStats <- rbind(SumStats, add_step_effect_size(LTER.redsNap, "M. franciscanus", "Naples SMCA", "LTER", "Den"))

#----- M. franciscanus Biomass (LTER) -----#
LTER.reds.bio <- subset(All.RR.sub.trans, source == "LTER" &
                           y == "Mesocentrotus franciscanus" & resp == "Bio")
LTER.reds.bio <- add_time_columns(LTER.reds.bio)

# Naples only since Campus Point was linear in density in the before
LTER.nap.red <- filter(LTER.reds.bio, CA_MPA_Name_Short == "Naples SMCA" & time > 0)
SumStats <- rbind(SumStats, add_mean_effect_size(LTER.nap.red, "M. franciscanus", "Naples SMCA", "LTER", "Bio"))

#----- M. pyrifera Biomass (LTER) -----#
# Note: only biomass since density is the exact same proportional data (bio is a multiplier of density)
LTER.macro <- subset(All.RR.sub.trans, source == "LTER" &
                       y == "Macrocystis pyrifera" & resp == "Bio")
LTER.macro <- add_time_columns(LTER.macro)
LinearBefore.macro <- test_linear_before(LTER.macro)
pbacips.macro <- run_pbacips_loop(LTER.macro)

# Naples
LTER.nap.macro <- filter(LTER.macro, CA_MPA_Name_Short == "Naples SMCA")
SumStats <- rbind(SumStats, add_linear_effect_size(LTER.nap.macro, "M. pyrifera", "Naples SMCA", "LTER", "Bio"))

# Campus Point
LTER.cp.macro <- filter(LTER.macro, CA_MPA_Name_Short == "Campus Point SMCA")
SumStats <- rbind(SumStats, add_linear_effect_size(LTER.cp.macro, "M. pyrifera", "Campus Point SMCA", "LTER", "Bio"))

#----- S. pulcher Density (LTER) -----#
LTER.SPUL.den <- subset(All.RR.sub.trans, source == "LTER" &
                           y == "Semicossyphus pulcher" & resp == "Den")
LTER.SPUL.den <- add_time_columns(LTER.SPUL.den)
LinearBefore.spul.den <- test_linear_before(LTER.SPUL.den)
pbacips.spul.den <- run_pbacips_loop(LTER.SPUL.den)

# Naples was linear in the before; Campus Point only
LTER.cp.SPUL.den <- filter(LTER.SPUL.den, CA_MPA_Name_Short == "Campus Point SMCA")
SumStats <- rbind(SumStats, add_linear_effect_size(LTER.cp.SPUL.den, "S. pulcher", "Campus Point SMCA", "LTER", "Den"))

#----- S. pulcher Biomass (LTER) -----#
LTER.SPUL.bio <- subset(All.RR.sub.trans, source == "LTER" &
                           y == "Semicossyphus pulcher" & resp == "Bio")
LTER.SPUL.bio <- add_time_columns(LTER.SPUL.bio)
LinearBefore.spul.bio <- test_linear_before(LTER.SPUL.bio)
pbacips.spul.bio <- run_pbacips_loop(LTER.SPUL.bio)

# Campus Point (Naples biomass was linear in the before)
LTER.cp.SPUL.bio <- filter(LTER.SPUL.bio, CA_MPA_Name_Short == "Campus Point SMCA")
SumStats <- rbind(SumStats, add_linear_effect_size(LTER.cp.SPUL.bio, "S. pulcher", "Campus Point SMCA", "LTER", "Bio"))

#----- P. interruptus Density (LTER) -- CI data only -----#
LTER.PANINT.den <- subset(All.RR.sub.trans, source == "LTER" &
                             y == "Panulirus interruptus" & resp == "Den")

# Campus Point
LTER.cp.LOB <- filter(LTER.PANINT.den, CA_MPA_Name_Short == "Campus Point SMCA")
SumStats <- rbind(SumStats, add_mean_effect_size(LTER.cp.LOB, "P. interruptus", "Campus Point SMCA", "LTER", "Den"))

# Naples
LTER.nap.LOB <- filter(LTER.PANINT.den, CA_MPA_Name_Short == "Naples SMCA")
SumStats <- rbind(SumStats, add_mean_effect_size(LTER.nap.LOB, "P. interruptus", "Naples SMCA", "LTER", "Den"))

#----- P. interruptus Biomass (LTER) -----#
LTER.PANINT.bio <- subset(All.RR.sub.trans, source == "LTER" &
                             y == "Panulirus interruptus" & resp == "Bio")

# Campus Point (linear)
LTER.cp.LOB.bio <- filter(LTER.PANINT.bio, CA_MPA_Name_Short == "Campus Point SMCA")
SumStats <- rbind(SumStats, add_linear_effect_size(LTER.cp.LOB.bio, "P. interruptus", "Campus Point SMCA", "LTER", "Bio",
                       time_var = "time"))

# Naples (mean only)
LTER.nap.LOB.bio <- filter(LTER.PANINT.bio, CA_MPA_Name_Short == "Naples SMCA")
SumStats <- rbind(SumStats, add_mean_effect_size(LTER.nap.LOB.bio, "P. interruptus", "Naples SMCA", "LTER", "Bio"))


####################################################################################################
## PISCO ANALYSES ##################################################################################
####################################################################################################

# The PISCO loop pattern: for each taxa x response, loop over MPAs and calculate both
# linear and mean effect sizes using run_ci_analysis()

#----- M. franciscanus Density (PISCO) -----#
Mes.den <- subset(All.RR.sub.trans, y == "Mesocentrotus franciscanus" &
                    source == "PISCO" & resp == "Den")
SumStats <- rbind(SumStats, run_ci_analysis(Mes.den, "M. franciscanus", "PISCO", "Den"))

#----- M. franciscanus Biomass (PISCO) -----#
Mes.bio <- subset(All.RR.sub.trans, y == "Mesocentrotus franciscanus" &
                    source == "PISCO" & resp == "Bio")
Mes.bio <- subset(Mes.bio, CA_MPA_Name_Short != "Matlahuayl SMR")
SumStats <- rbind(SumStats, run_ci_analysis(Mes.bio, "M. franciscanus", "PISCO", "Bio"))

#----- S. purpuratus Density (PISCO) -----#
Str.den <- subset(All.RR.sub.trans, y == "Strongylocentrotus purpuratus" &
                    source == "PISCO" & resp == "Den")
SumStats <- rbind(SumStats, run_ci_analysis(Str.den, "S. purpuratus", "PISCO", "Den"))

#----- S. purpuratus Biomass (PISCO) -----#
Str.bio <- subset(All.RR.sub.trans, y == "Strongylocentrotus purpuratus" &
                    source == "PISCO" & resp == "Bio")
Str.bio <- subset(Str.bio, CA_MPA_Name_Short != "Matlahuayl SMR")
SumStats <- rbind(SumStats, run_ci_analysis(Str.bio, "S. purpuratus", "PISCO", "Bio"))

#----- M. pyrifera Biomass (PISCO) -----#
Mac <- subset(All.RR.sub.trans, y == "Macrocystis pyrifera" &
                source == "PISCO" & resp == "Bio")
SumStats <- rbind(SumStats, run_ci_analysis(Mac, "M. pyrifera", "PISCO", "Bio"))

#----- P. interruptus Density (PISCO) -----#
Pan.den <- subset(All.RR.sub.trans, y == "Panulirus interruptus" &
                    source == "PISCO" & resp == "Den")
SumStats <- rbind(SumStats, run_ci_analysis(Pan.den, "P. interruptus", "PISCO", "Den"))

#----- P. interruptus Biomass (PISCO) -----#
Pan.bio <- subset(All.RR.sub.trans, y == "Panulirus interruptus" &
                    resp == "Bio" & source == "PISCO")
# Exclude MPAs with insufficient lobster biomass data
Pan.bio <- subset(Pan.bio,
  CA_MPA_Name_Short != "Gull Island SMR" &
  CA_MPA_Name_Short != "Harris Point SMR" &
  CA_MPA_Name_Short != "South Point SMR" &
  CA_MPA_Name_Short != "Scorpion SMR" &
  CA_MPA_Name_Short != "Anacapa Island SMR 2003" &
  CA_MPA_Name_Short != "Santa Barbara Island SMR" &
  CA_MPA_Name_Short != "Point Dume SMCA"
)
SumStats <- rbind(SumStats, run_ci_analysis(Pan.bio, "P. interruptus", "PISCO", "Bio"))

#----- S. pulcher Density (PISCO) -----#
Sheep.den <- subset(All.RR.sub.trans, y == "Semicossyphus pulcher" &
                      resp == "Den" & source == "PISCO")
SumStats <- rbind(SumStats, run_ci_analysis(Sheep.den, "S. pulcher", "PISCO", "Den"))

#----- S. pulcher Biomass (PISCO) -----#
Sheep.bio <- subset(All.RR.sub.trans, y == "Semicossyphus pulcher" &
                      resp == "Bio" & source == "PISCO")
SumStats <- rbind(SumStats, run_ci_analysis(Sheep.bio, "S. pulcher", "PISCO", "Bio"))


####################################################################################################
## KFM ANALYSES ####################################################################################
####################################################################################################

# KFM sites are split into:
# BA sites (Scorpion, SBI, Harris Point, Gull Island) -> pBACIPS with before/after data
# CI sites (Anacapa 2003, South Point) -> CI-only linear regression

# KFM BA site lists
KFM_BA_SITES <- c("Scorpion SMR", "Santa Barbara Island SMR",
                   "Harris Point SMR", "Gull Island SMR")
KFM_CI_SITES <- c("Anacapa Island SMR 2003", "South Point SMR")

#----- S. purpuratus Density (KFM) -----#
KFM.purps.den <- subset(All.RR.sub, source == "KFM" &
                           y == "Strongylocentrotus purpuratus" & resp == "Den")
KFM.purps.den <- add_time_columns(KFM.purps.den)

# BA sites: test for linear before
KFM.purps.BA <- subset(KFM.purps.den, CA_MPA_Name_Short %in% KFM_BA_SITES)
KFM.purps.CI <- subset(KFM.purps.den, CA_MPA_Name_Short %in% KFM_CI_SITES)
LinearBefore.kfm.purps <- test_linear_before(KFM.purps.BA)

# Santa Barbara Island SMR is only site not linear in before period
KFM.SB.purp <- filter(KFM.purps.den, CA_MPA_Name_Short == "Santa Barbara Island SMR")
SumStats <- rbind(SumStats, add_step_effect_size(KFM.SB.purp, "S. purpuratus", "Santa Barbara Island SMR", "KFM", "Den"))

# CI sites
SumStats <- rbind(SumStats, run_ci_analysis(KFM.purps.CI, "S. purpuratus", "KFM", "Den"))

#----- S. purpuratus Biomass (KFM) -----#
KFM.purps.bio <- subset(All.RR.sub, source == "KFM" &
                           y == "Strongylocentrotus purpuratus" & resp == "Bio")
KFM.purps.bio <- add_time_columns(KFM.purps.bio)

# All BA sites with before data were linear in before for density, so removed
# CI sites only
KFM.purps.bio.CI <- subset(KFM.purps.bio, CA_MPA_Name_Short %in% KFM_CI_SITES)
SumStats <- rbind(SumStats, run_ci_analysis(KFM.purps.bio.CI, "S. purpuratus", "KFM", "Bio"))

#----- M. franciscanus Density (KFM) -----#
KFM.reds.den <- subset(All.RR.sub, source == "KFM" &
                          y == "Mesocentrotus franciscanus" & resp == "Den")
KFM.reds.den <- add_time_columns(KFM.reds.den)

KFM.reds.den.BA <- subset(KFM.reds.den, CA_MPA_Name_Short %in% KFM_BA_SITES)
KFM.reds.den.CI <- subset(KFM.reds.den, CA_MPA_Name_Short %in% KFM_CI_SITES)
LinearBefore.kfm.reds <- test_linear_before(KFM.reds.den.BA)

# Scorpion SMR - step model
SumStats <- rbind(SumStats, add_step_effect_size(
  filter(KFM.reds.den, CA_MPA_Name_Short == "Scorpion SMR"),
  "M. franciscanus", "Scorpion SMR", "KFM", "Den"
))

# Gull Island SMR - step model
SumStats <- rbind(SumStats, add_step_effect_size(
  filter(KFM.reds.den, CA_MPA_Name_Short == "Gull Island SMR"),
  "M. franciscanus", "Gull Island SMR", "KFM", "Den"
))

# Harris Point SMR - linear model
SumStats <- rbind(SumStats, add_linear_effect_size(
  subset(KFM.reds.den, CA_MPA_Name_Short == "Harris Point SMR"),
  "M. franciscanus", "Harris Point SMR", "KFM", "Den"
))

# CI sites
SumStats <- rbind(SumStats, run_ci_analysis(KFM.reds.den.CI, "M. franciscanus", "KFM", "Den"))

#----- M. franciscanus Biomass (KFM) -----#
KFM.reds.bio <- subset(All.RR.sub, source == "KFM" &
                          y == "Mesocentrotus franciscanus" & resp == "Bio")
KFM.reds.bio <- add_time_columns(KFM.reds.bio)

KFM.reds.bio.BA <- subset(KFM.reds.bio, CA_MPA_Name_Short %in% KFM_BA_SITES)
KFM.reds.bio.CI <- subset(KFM.reds.bio, CA_MPA_Name_Short %in% KFM_CI_SITES)

# Run pBACIPS on BA sites
pbacips.kfm.reds.bio <- run_pbacips_loop(KFM.reds.bio.BA)

# Scorpion SMR - step model
SumStats <- rbind(SumStats, add_step_effect_size(
  filter(KFM.reds.bio.BA, CA_MPA_Name_Short == "Scorpion SMR"),
  "M. franciscanus", "Scorpion SMR", "KFM", "Bio"
))

# Gull Island SMR - linear model
SumStats <- rbind(SumStats, add_linear_effect_size(
  filter(KFM.reds.bio.BA, CA_MPA_Name_Short == "Gull Island SMR"),
  "M. franciscanus", "Gull Island SMR", "KFM", "Bio"
))

# Harris Point SMR - sigmoid model
KFM.mfran.harris <- subset(KFM.reds.bio.BA, CA_MPA_Name_Short == "Harris Point SMR")
time.model.of.impact <- max(which(KFM.mfran.harris$time.model == 0))
time.true <- KFM.mfran.harris$time.true
sigmoid.Model <- mySIGfun_standalone(
  delta = KFM.mfran.harris$lnDiff,
  time.model = KFM.mfran.harris$time.model,
  time.model.of.impact = time.model.of.impact,
  time.true = time.true
)
# Use confidence interval (not prediction interval) for effect sizes
# Confidence intervals quantify uncertainty in the mean, which is what effect sizes estimate
# Prediction intervals include residual variance and would overestimate uncertainty
interval <- data.frame(predFit(sigmoid.Model, newdata = KFM.mfran.harris, interval = "confidence"))
interval$se <- abs((interval$lwr - interval$fit) / 1.96)
n_obs <- nrow(interval)
n_params <- length(coef(sigmoid.Model))
n_df <- n_obs - n_params
interval$stdev <- interval$se * sqrt(n_obs)
pSD <- sqrt((((n_df) * (interval$stdev[1]^2)) + ((n_df) * (interval$stdev[n_obs]^2))) / (n_df + 1))
pSE <- pSD * sqrt((1 / n_obs) + (1 / n_obs))
pCI <- pSE * 1.96
mean_es <- interval$fit[n_obs] - interval$fit[1]
SumStats[nrow(SumStats) + 1, ] <- c("M. franciscanus", "Harris Point SMR", mean_es, pSE, pSD, pCI,
                                      "Sigmoid", "KFM", "Bio", "Y", "Y", "pBACIPS", "N")

# CI sites
SumStats <- rbind(SumStats, run_ci_analysis(KFM.reds.bio.CI, "M. franciscanus", "KFM", "Bio"))

#----- M. pyrifera Biomass (KFM) -----#
KFM.macro <- subset(All.RR.sub.trans, source == "KFM" &
                      y == "Macrocystis pyrifera" & resp == "Bio")
KFM.macro <- add_time_columns(KFM.macro)

KFM.macro.BA <- subset(KFM.macro, CA_MPA_Name_Short == "Scorpion SMR" |
                          CA_MPA_Name_Short == "Harris Point SMR" |
                          CA_MPA_Name_Short == "Gull Island SMR")
KFM.macro.CI <- subset(KFM.macro, CA_MPA_Name_Short %in%
                          c("Anacapa Island SMR 2003", "South Point SMR", "Santa Barbara Island SMR"))

LinearBefore.kfm.macro <- test_linear_before(KFM.macro.BA)
# Harris Point is linear in before
pbacips.kfm.macro <- run_pbacips_loop(KFM.macro.BA)

# Gull Island SMR - linear
SumStats <- rbind(SumStats, add_linear_effect_size(
  filter(KFM.macro.BA, CA_MPA_Name_Short == "Gull Island SMR"),
  "M. pyrifera", "Gull Island SMR", "KFM", "Bio"
))

# Scorpion SMR - sigmoid
KFM.macro.stipe <- filter(KFM.macro.BA, CA_MPA_Name_Short == "Scorpion SMR")
time.model.of.impact <- max(which(KFM.macro.stipe$time.model == 0))
time.true <- KFM.macro.stipe$time.true
sigmoid.Model <- mySIGfun_standalone(
  delta = KFM.macro.stipe$lnDiff,
  time.model = KFM.macro.stipe$time.model,
  time.model.of.impact = time.model.of.impact,
  time.true = time.true
)
# Use confidence interval for effect size estimation (uncertainty in mean, not prediction)
interval <- data.frame(predFit(sigmoid.Model, newdata = KFM.macro.stipe, interval = "confidence"))
interval$se <- abs((interval$lwr - interval$fit) / 1.96)
n_obs <- nrow(interval)
n_params <- length(coef(sigmoid.Model))
n_df <- n_obs - n_params
interval$stdev <- interval$se * sqrt(n_obs)
pSD <- sqrt((((n_df) * (interval$stdev[1]^2)) + ((n_df) * (interval$stdev[n_obs]^2))) / (n_df + 1))
pSE <- pSD * sqrt((1 / n_obs) + (1 / n_obs))
pCI <- pSE * 1.96
mean_es <- interval$fit[n_obs] - interval$fit[1]
SumStats[nrow(SumStats) + 1, ] <- c("M. pyrifera", "Scorpion SMR", mean_es, pSE, pSD, pCI,
                                      "Sigmoid", "KFM", "Bio", "Y", "Y", "pBACIPS", "N")

# CI sites
SumStats <- rbind(SumStats, run_ci_analysis(KFM.macro.CI, "M. pyrifera", "KFM", "Bio"))

#----- P. interruptus Density (KFM) -----#
KFM.lob <- subset(All.RR.sub.trans, source == "KFM" &
                     y == "Panulirus interruptus")
KFM.lob$CA_MPA_Name_Short <- as.character(KFM.lob$CA_MPA_Name_Short)
KFM.lob <- add_time_columns(KFM.lob)

# Dropped Harris, no lobsters in before suggests not good sampling
KFM.lob.BA <- subset(KFM.lob, CA_MPA_Name_Short %in%
                        c("Scorpion SMR", "Santa Barbara Island SMR", "Gull Island SMR"))
KFM.lob.CI <- subset(KFM.lob, CA_MPA_Name_Short %in% KFM_CI_SITES)

LinearBefore.kfm.lob <- test_linear_before(KFM.lob.BA)

# Run pBACIPS on subset (exclude Harris and SBI from pBACIPS loop)
KFM.lob.sub <- subset(KFM.lob.BA, CA_MPA_Name_Short != "Harris Point SMR" &
                         CA_MPA_Name_Short != "Santa Barbara Island SMR")
pbacips.kfm.lob <- run_pbacips_loop(KFM.lob.sub)

# Gull Island SMR - sigmoid
KFM.lob.gull <- filter(KFM.lob.sub, CA_MPA_Name_Short == "Gull Island SMR")
time.model.of.impact <- max(which(KFM.lob.gull$time == 0))
time.true <- KFM.lob.gull$time
sigmoid.Model <- mySIGfun_standalone(
  delta = KFM.lob.gull$lnDiff,
  time.model = KFM.lob.gull$time.model,
  time.model.of.impact = time.model.of.impact,
  time.true = time.true
)
# Use confidence interval for effect size estimation (uncertainty in mean, not prediction)
interval <- data.frame(predFit(sigmoid.Model, newdata = KFM.lob.gull, interval = "confidence"))
interval$se <- abs((interval$lwr - interval$fit) / 1.96)
n_obs <- nrow(interval)
n_params <- length(coef(sigmoid.Model))
n_df <- n_obs - n_params
interval$stdev <- interval$se * sqrt(n_obs)
pSD <- sqrt((((n_df) * (interval$stdev[1]^2)) +
               ((n_df) * (interval$stdev[n_obs]^2))) / (n_df + 1))
pSE <- pSD * sqrt((1 / n_obs) + (1 / n_obs))
pCI <- pSE * 1.96
mean_es <- interval$fit[n_obs] - interval$fit[1]
SumStats[nrow(SumStats) + 1, ] <- c("P. interruptus", "Gull Island SMR", mean_es, pSE, pSD, pCI,
                                      "Sigmoid", "KFM", "Den", "Y", "Y", "pBACIPS", "N")

# Scorpion SMR - linear
SumStats <- rbind(SumStats, add_linear_effect_size(
  filter(KFM.lob.sub, CA_MPA_Name_Short == "Scorpion SMR"),
  "P. interruptus", "Scorpion SMR", "KFM", "Den"
))

# Santa Barbara Island SMR - linear
KFM.lobs.sbi <- subset(KFM.lob.BA, CA_MPA_Name_Short == "Santa Barbara Island SMR")
SumStats <- rbind(SumStats, add_linear_effect_size(KFM.lobs.sbi, "P. interruptus", "Santa Barbara Island SMR", "KFM", "Den"))

# CI sites
SumStats <- rbind(SumStats, run_ci_analysis(KFM.lob.CI, "P. interruptus", "KFM", "Den"))

#----- S. pulcher Density (KFM) -----#
Sheep.den.kfm <- subset(All.RR.sub.trans, y == "Semicossyphus pulcher" &
                           resp == "Den" & source == "KFM")
Sheep.den.kfm <- add_time_columns(Sheep.den.kfm)

KFM.sheep.BA <- subset(Sheep.den.kfm, CA_MPA_Name_Short %in% KFM_BA_SITES)
KFM.sheep.CI <- subset(Sheep.den.kfm, CA_MPA_Name_Short %in% KFM_CI_SITES)

LinearBefore.kfm.sheep <- test_linear_before(KFM.sheep.BA)

# Harris Point SMR - linear
SumStats <- rbind(SumStats, add_linear_effect_size(
  filter(KFM.sheep.BA, CA_MPA_Name_Short == "Harris Point SMR"),
  "S. pulcher", "Harris Point SMR", "KFM", "Den"
))

# Scorpion SMR - step
SumStats <- rbind(SumStats, add_step_effect_size(
  filter(KFM.sheep.BA, CA_MPA_Name_Short == "Scorpion SMR"),
  "S. pulcher", "Scorpion SMR", "KFM", "Den"
))

# CI sites
SumStats <- rbind(SumStats, run_ci_analysis(KFM.sheep.CI, "S. pulcher", "KFM", "Den"))


####################################################################################################
## LANDSAT Kelp Biomass Effect Sizes ###############################################################
####################################################################################################

# Landsat data processing is now handled by 06b_landsat_processing.R which produces Landsat.RR.
# Here we only run the effect size analysis loop over the pre-processed response ratios.

if (exists("Landsat.RR") && nrow(Landsat.RR) > 0) {
  mpas_landsat <- unique(Landsat.RR$CA_MPA_Name_Short)
  for (mpa in mpas_landsat) {
    mpa_dat <- filter(Landsat.RR, CA_MPA_Name_Short == mpa)
    if (nrow(mpa_dat) > 3) {
      SumStats <- rbind(SumStats, add_linear_effect_size(mpa_dat, "M. pyrifera", mpa, "Landsat", "Bio",
                             time_var = "time"))
    }
  }
}

####################################################################################################
## Finalize SumStats ###############################################################################
####################################################################################################

# Convert numeric columns
SumStats$Mean <- as.numeric(SumStats$Mean)
SumStats$SE <- as.numeric(SumStats$SE)
SumStats$CI <- as.numeric(SumStats$CI)
SumStats$SD <- as.numeric(SumStats$SD)

# Remove incomplete rows
SumStats.sub <- SumStats[complete.cases(SumStats), ]

# Set factor levels for plotting
SumStats.sub$Taxa <- factor(SumStats.sub$Taxa,
  levels = c("S. purpuratus", "M. franciscanus", "M. pyrifera", "P. interruptus", "S. pulcher"),
  ordered = TRUE
)
SumStats.sub$Resp <- factor(SumStats.sub$Resp, levels = c("Bio", "Den"), ordered = TRUE)
SumStats.sub$Source <- factor(SumStats.sub$Source,
  levels = c("KFM", "LTER", "PISCO", "Landsat"), ordered = TRUE
)
SumStats.sub$Type <- factor(SumStats.sub$Type,
  levels = c("pBACIPS", "BACI", "CI"), ordered = TRUE
)

# Merge with site metadata
SumStats.sub <- merge(SumStats.sub, Site, by.x = "MPA", by.y = "CA_MPA_Name_Short")

####################################################################################################
## Filter SumStats.Final for publication ###########################################################
####################################################################################################

# Include pBACIPS results, plus CI results where before period was not linear and is primary
SumStats.Final <- subset(SumStats.sub,
  Type.x == "pBACIPS" | (Type.x == "CI" & LinearBefore == "NA" & Primary == "Y")
)

# Remove problematic MPAs:
# Painted Cave SMCA (lobster take allowed), San Miguel Island SC (weird location),
# Arrow Point (finfish fishing allowed), Judith Rk SMR (overlaps San Miguel Island SC),
# Point Conception SMR
SumStats.Final <- subset(SumStats.Final,
  MPA != "Painted Cave SMCA" &
  MPA != "San Miguel Island SC" &
  MPA != "Arrow Point to Lion Head Point SMCA" &
  MPA != "Judith Rk SMR" &
  MPA != "Point Conception SMR"
)

cat("Effect size calculation complete.\n")
cat("SumStats rows:", nrow(SumStats), "\n")
cat("SumStats.Final rows:", nrow(SumStats.Final), "\n")
