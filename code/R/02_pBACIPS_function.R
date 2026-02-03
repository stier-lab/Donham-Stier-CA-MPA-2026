# =============================================================================
# 02_pBACIPS_function.R
# =============================================================================
#
# PURPOSE:
#   Implement the Progressive-Change BACIPS (Before-After-Control-Impact-Pairs
#   with Symmetry) analysis method. This is the core statistical approach for
#   detecting MPA effects that may change progressively over time.
#
# WHAT THIS SCRIPT DOES:
#   1. Defines the main ProgressiveChangeBACIPS() function
#   2. Fits four candidate temporal models to the MPA vs reference difference:
#      - Step model: Abrupt change after MPA implementation
#      - Linear model: Gradual linear change over time
#      - Asymptotic model: Change that levels off (Michaelis-Menten form)
#      - Sigmoid model: S-shaped change (logistic form)
#   3. Selects the best model using AICc (corrected Akaike Information Criterion)
#   4. Provides standalone model fitting functions for effect size calculation
#
# BACKGROUND:
#   Traditional BACI designs assume an instantaneous "step" change after
#   intervention. But ecological responses to protection often unfold gradually:
#   - Fish populations may take years to recover
#   - Trophic cascades propagate slowly through the food web
#   - Habitat structure changes incrementally
#
#   The pBACIPS approach tests whether the data support a step change or one
#   of several progressive change models, providing more realistic estimates
#   of MPA effects.
#
# REFERENCE:
#   Thiault, L., Kernaléguen, L., Osenberg, C.W. & Claudet, J. (2017)
#   "Progressive-Change BACIPS: a flexible approach for environmental impact
#   assessment." Methods in Ecology and Evolution, 8(3), 288-296.
#
# AUTHORS: L. Thiault, L. Kernaléguen, C.W. Osenberg & J. Claudet
#          (adapted by Emily Donham & Adrian Stier)
# PROJECT: CA MPA Kelp Forest pBACIPS Analysis
# =============================================================================


# =============================================================================
# REQUIRED PACKAGES
# =============================================================================
# These should already be loaded by 00_libraries.R, but we ensure they're
# available for this specific functionality

require(minpack.lm)    # Levenberg-Marquardt nonlinear least squares
require(nls2)          # Enhanced NLS with grid search
require(AICcmodavg)    # AICc calculation


# =============================================================================
# MAIN FUNCTION: ProgressiveChangeBACIPS
# =============================================================================

#' Progressive-Change BACIPS Analysis
#'
#' Fits four candidate models to BACI delta (Impact - Control) data and
#' selects the best model using AICc weights.
#'
#' @param control Numeric vector of response values at Control (reference) site
#'        at each sampling time. Should be on the same scale as impact.
#' @param impact Numeric vector of response values at Impact (MPA) site at each
#'        sampling time. Same length as control.
#' @param time.true Sequential time index for ALL sampling dates (1, 2, 3, ...)
#'        This represents the full time series regardless of before/after.
#' @param time.model Time since intervention: 0 for Before period, sequential
#'        values (0, 1, 2, ...) for After period. This is the "treatment time".
#'
#' @return A named list with components:
#'   - weights: AICc weights for each model (% likelihood)
#'   - step: The fitted step (ANOVA) model object
#'   - linear: The fitted linear model object
#'   - asymptotic: The fitted asymptotic (NLS) model object
#'   - sigmoid: The fitted sigmoid (NLS) model object
#'   - best: Integer index of the best model (1=step, 2=linear, 3=asymp, 4=sig)
#'
#' @details
#' The function proceeds in 4 steps:
#'
#' STEP 1: Calculate delta at each time point
#'   delta = impact - control
#'   This is the difference between MPA and reference at each sampling time.
#'
#' STEP 2: Define the four candidate models
#'   All models predict delta as a function of time.model:
#'
#'   Step Model: delta ~ period (categorical Before/After)
#'     - Tests for a simple mean shift between periods
#'     - Fit using ANOVA (aov)
#'
#'   Linear Model: delta ~ time.model
#'     - Tests for a gradual linear trend starting at MPA implementation
#'     - Fit using ordinary least squares (lm)
#'
#'   Asymptotic Model: delta ~ (M * time.model) / (L + time.model) + B
#'     - Michaelis-Menten form: rapid initial change that levels off
#'     - M = maximum change, L = half-saturation time, B = baseline
#'     - Fit using nonlinear least squares (nls.lm + nls2)
#'
#'   Sigmoid Model: delta ~ (M * (time.model/L)^K) / (1 + (time.model/L)^K) + B
#'     - Hill function form: S-shaped change with flexible steepness
#'     - M = maximum change, L = midpoint, K = steepness, B = baseline
#'     - Fit using nonlinear least squares (nls.lm + nls2)
#'
#' STEP 3: Compare models using AICc
#'   - AICc penalizes for model complexity more than AIC
#'   - Appropriate for small sample sizes (typical in BACI studies)
#'   - Weights sum to 100% and represent relative likelihood
#'
#' STEP 4: Report best model
#'   - Prints AICc weights and summary of best-fitting model
#'   - Returns all model objects for further analysis
#'
#' @examples
#' # See 08_effect_sizes.R for full usage examples
#' # Basic usage:
#' # result <- ProgressiveChangeBACIPS(control, impact, time.true, time.model)
#' # result$weights  # Check which model is best
#' # result$best     # Index of best model
ProgressiveChangeBACIPS <- function(control, impact, time.true, time.model) {

  # ---------------------------------------------------------------------------
  # STEP 1: Calculate delta at each sampling date
  # ---------------------------------------------------------------------------
  # Delta is the difference between impact (MPA) and control (reference)
  # Positive delta = higher values inside MPA than outside

  delta <- impact - control

  # Find the last time point in the "Before" period
  # This is needed for setting initial parameter values
  time.model.of.impact <- max(which(time.model == 0))

  # ---------------------------------------------------------------------------
  # STEP 2: Create period factor for step model
  # ---------------------------------------------------------------------------
  # Convert time.model to categorical Before/After

  period <- ifelse(time.model == 0, "Before", "After")

  # ---------------------------------------------------------------------------
  # STEP 3: Fit the four candidate models
  # ---------------------------------------------------------------------------

  ## --- Step Model (ANOVA) ---
  # H0: mean delta is the same Before and After
  # This is the traditional BACI analysis approach

  step.Model <- aov(delta ~ period)

  ## --- Linear Model (OLS) ---
  # H0: no linear trend in delta over time since MPA
  # Tests for gradual, constant-rate change

  linear.Model <- lm(delta ~ time.model)

  ## --- Asymptotic Model (Nonlinear) ---
  # Michaelis-Menten form: rapid initial change that saturates
  # Equation: delta = (M * time) / (L + time) + B
  # M = maximum effect, L = time to reach half of maximum, B = baseline

  myASYfun <- function(delta, time.model) {
    # Define the model formula
    foAsy <- delta ~ (M * time.model) / (L + time.model) + B

    # Define prediction function for nls.lm
    funAsy <- function(parS, time.model) {
      (parS$M * time.model) / (parS$L + time.model) + parS$B
    }

    # Residual function: observed - predicted
    residFun <- function(p, observed, time.model) {
      observed - funAsy(p, time.model)
    }

    # Set starting values based on data
    # M: mean of After period values (expected maximum change)
    # B: mean of Before period values (baseline)
    # L: start with 1 (half-saturation at 1 year)
    parStart <- list(
      M = mean(delta[time.model.of.impact:length(time.true)]),
      B = mean(delta[1:time.model.of.impact]),
      L = 1
    )

    # First pass: use Levenberg-Marquardt algorithm (more robust)
    nls_ASY_out <- nls.lm(
      par = parStart,
      fn = residFun,
      observed = delta,
      time.model = time.model,
      control = nls.lm.control(maxfev = integer(), maxiter = 1000)
    )

    # Refine with nls2 grid search around the LM solution
    startPar <- as.list(coef(nls_ASY_out))
    gridPar <- expand.grid(
      M = startPar$M * c(0.5, 1.0, 1.5),
      B = startPar$B * c(0.5, 1.0, 1.5),
      L = pmax(0.01, startPar$L * c(0.5, 1.0, 1.5))  # L must be positive
    )

    # Try grid search, fall back to single point if it fails
    nls2_fit <- tryCatch(
      nls2(foAsy, start = gridPar, algorithm = "brute-force"),
      error = function(e) {
        nls2(foAsy, start = data.frame(t(unlist(startPar))),
             algorithm = "brute-force")
      }
    )

    # Final refinement with default NLS algorithm
    tryCatch(
      nls(foAsy, start = coef(nls2_fit)),
      error = function(e) nls2_fit
    )
  }

  asymptotic.Model <- myASYfun(delta = delta, time.model = time.model)

  ## --- Sigmoid Model (Nonlinear) ---
  # Hill function form: S-shaped change with adjustable steepness
  # Equation: delta = (M * (time/L)^K) / (1 + (time/L)^K) + B
  # M = maximum effect, L = midpoint (time to half-max), K = steepness, B = baseline

  mySIGfun <- function(delta, time.model) {
    # Define the model formula
    foSIG <- delta ~ (M * (time.model / L)^K) / (1 + (time.model / L)^K) + B

    # Define prediction function for nls.lm
    funSIG <- function(parS, time.model) {
      (parS$M * (time.model / parS$L)^parS$K) /
        (1 + (time.model / parS$L)^parS$K) + parS$B
    }

    # Residual function: observed - predicted
    residFun <- function(p, observed, time.model) {
      observed - funSIG(p, time.model)
    }

    # Set starting values
    # K = 5 gives a moderately steep sigmoid
    parStart <- list(
      M = mean(delta[time.model.of.impact:length(time.true)]),
      B = mean(delta[1:time.model.of.impact]),
      L = mean(time.model),  # Midpoint at average time
      K = 5
    )

    # First pass: Levenberg-Marquardt
    nls_SIG_out <- nls.lm(
      par = parStart,
      fn = residFun,
      observed = delta,
      time.model = time.model,
      control = nls.lm.control(maxfev = integer(), maxiter = 1000)
    )

    # Refine with nls2 grid search
    startPar <- as.list(coef(nls_SIG_out))
    gridPar <- expand.grid(
      M = startPar$M * c(0.5, 1.0, 1.5),
      B = startPar$B * c(0.5, 1.0, 1.5),
      L = pmax(0.01, startPar$L * c(0.5, 1.0, 1.5)),
      K = pmax(0.1, startPar$K * c(0.5, 1.0, 1.5))
    )

    nls2_fit <- tryCatch(
      nls2(foSIG, start = gridPar, algorithm = "brute-force"),
      error = function(e) {
        nls2(foSIG, start = data.frame(t(unlist(startPar))),
             algorithm = "brute-force")
      }
    )

    # Final refinement
    tryCatch(
      nls(foSIG, start = coef(nls2_fit)),
      error = function(e) nls2_fit
    )
  }

  sigmoid.Model <- mySIGfun(delta = delta, time.model = time.model)

  # ---------------------------------------------------------------------------
  # STEP 4: Model comparison using AICc
  # ---------------------------------------------------------------------------
  # AICc is preferred over AIC for small samples (n/K < 40)
  # where n = sample size and K = number of parameters

  # Get AIC from all models
  AIC.test <- AIC(step.Model, linear.Model, asymptotic.Model, sigmoid.Model)

  # Calculate AICc for each model
  AICc.test <- as.data.frame(cbind(
    AIC.test[, 1],  # df
    c(AICc(step.Model), AICc(linear.Model),
      AICc(asymptotic.Model), AICc(sigmoid.Model))  # AICc
  ))
  rownames(AICc.test) <- rownames(AIC.test)
  names(AICc.test) <- names(AIC.test)

  # Calculate AICc weights (Akaike weights)
  # Lower AICc = better fit; weights sum to 1 (or 100%)
  for (i in 1:dim(AICc.test)[1]) {
    AICc.test$diff[i] <- AICc.test$AIC[i] - min(AICc.test$AIC)
  }

  # Relative likelihood: exp(-0.5 * delta_AICc)
  AICc.test$RL <- exp(-0.5 * AICc.test$diff)

  # Normalize to get weights (as percentages)
  RL_sum <- sum(AICc.test$RL)
  AICc.test$aicWeights <- (AICc.test$RL / RL_sum) * 100

  # Extract weights vector
  w <- AICc.test$aicWeights
  names(w) <- rownames(AICc.test)

  # ---------------------------------------------------------------------------
  # STEP 5: Display results
  # ---------------------------------------------------------------------------

  # Print AICc weights
  print(w)

  # Identify best model
  best.Model <- which(w == max(w))

  # Print summary of best model
  if (best.Model == 1) {
    writeLines(paste("\n\nSTEP MODEL SELECTED - Likelihood = ",
                     round(w[1], 1), "%\n\n", sep = ""))
    print(summary(step.Model))
  }
  if (best.Model == 2) {
    writeLines(paste("\n\nLINEAR MODEL SELECTED - Likelihood = ",
                     round(w[2], 1), "%\n\n", sep = ""))
    print(summary(linear.Model))
  }
  if (best.Model == 3) {
    writeLines(paste("\n\nASYMPTOTIC MODEL SELECTED - Likelihood = ",
                     round(w[3], 1), "%\n\n", sep = ""))
    print(asymptotic.Model)
  }
  if (best.Model == 4) {
    writeLines(paste("\n\nSIGMOID MODEL SELECTED - Likelihood = ",
                     round(w[4], 1), "%\n\n", sep = ""))
    print(sigmoid.Model)
  }

  # ---------------------------------------------------------------------------
  # Return all models and weights
  # ---------------------------------------------------------------------------

  list(
    weights = w,
    step = step.Model,
    linear = linear.Model,
    asymptotic = asymptotic.Model,
    sigmoid = sigmoid.Model,
    best = best.Model
  )
}


# =============================================================================
# STANDALONE MODEL FITTING FUNCTIONS
# =============================================================================
# These functions are used in 08_effect_sizes.R to fit specific models
# independently of the full pBACIPS comparison. Useful when you already
# know which model form is appropriate.

#' Fit asymptotic model standalone
#'
#' Fits the Michaelis-Menten form asymptotic model outside of the main
#' ProgressiveChangeBACIPS function. Used for effect size calculation.
#'
#' @param delta Numeric vector of Impact - Control differences
#' @param time.model Time since intervention (0 for Before, sequential for After)
#' @param time.model.of.impact Index of last Before-period observation
#' @param time.true Full time index vector
#' @return Fitted NLS model object
myASYfun_standalone <- function(delta, time.model, time.model.of.impact, time.true) {
  foAsy <- delta ~ (M * time.model) / (L + time.model) + B

  funAsy <- function(parS, time.model) {
    (parS$M * time.model) / (parS$L + time.model) + parS$B
  }

  residFun <- function(p, observed, time.model) {
    observed - funAsy(p, time.model)
  }

  parStart <- list(
    M = mean(delta[time.model.of.impact:length(time.true)]),
    B = mean(delta[1:time.model.of.impact]),
    L = 1
  )

  nls_ASY_out <- nls.lm(
    par = parStart,
    fn = residFun,
    observed = delta,
    time.model = time.model,
    control = nls.lm.control(maxfev = integer(), maxiter = 1000)
  )

  startPar <- as.list(coef(nls_ASY_out))
  gridPar <- expand.grid(
    M = startPar$M * c(0.5, 1.0, 1.5),
    B = startPar$B * c(0.5, 1.0, 1.5),
    L = pmax(0.01, startPar$L * c(0.5, 1.0, 1.5))
  )

  nls2_fit <- tryCatch(
    nls2(foAsy, start = gridPar, algorithm = "brute-force"),
    error = function(e) {
      nls2(foAsy, start = data.frame(t(unlist(startPar))),
           algorithm = "brute-force")
    }
  )

  tryCatch(
    nls(foAsy, start = coef(nls2_fit)),
    error = function(e) nls2_fit
  )
}


#' Fit sigmoid model standalone
#'
#' Fits the Hill function form sigmoid model outside of the main
#' ProgressiveChangeBACIPS function. Used for effect size calculation.
#'
#' @param delta Numeric vector of Impact - Control differences
#' @param time.model Time since intervention (0 for Before, sequential for After)
#' @param time.model.of.impact Index of last Before-period observation
#' @param time.true Full time index vector
#' @return Fitted NLS model object
mySIGfun_standalone <- function(delta, time.model, time.model.of.impact, time.true) {
  foSIG <- delta ~ (M * (time.model / L)^K) / (1 + (time.model / L)^K) + B

  funSIG <- function(parS, time.model) {
    (parS$M * (time.model / parS$L)^parS$K) /
      (1 + (time.model / parS$L)^parS$K) + parS$B
  }

  residFun <- function(p, observed, time.model) {
    observed - funSIG(p, time.model)
  }

  parStart <- list(
    M = mean(delta[time.model.of.impact:length(time.true)]),
    B = mean(delta[1:time.model.of.impact]),
    L = mean(time.model),
    K = 5
  )

  nls_SIG_out <- nls.lm(
    par = parStart,
    fn = residFun,
    observed = delta,
    time.model = time.model,
    control = nls.lm.control(maxfev = integer(), maxiter = 1000)
  )

  startPar <- as.list(coef(nls_SIG_out))
  gridPar <- expand.grid(
    M = startPar$M * c(0.5, 1.0, 1.5),
    B = startPar$B * c(0.5, 1.0, 1.5),
    L = pmax(0.01, startPar$L * c(0.5, 1.0, 1.5)),
    K = pmax(0.1, startPar$K * c(0.5, 1.0, 1.5))
  )

  nls2_fit <- tryCatch(
    nls2(foSIG, start = gridPar, algorithm = "brute-force"),
    error = function(e) {
      nls2(foSIG, start = data.frame(t(unlist(startPar))),
           algorithm = "brute-force")
    }
  )

  tryCatch(
    nls(foSIG, start = coef(nls2_fit)),
    error = function(e) nls2_fit
  )
}


# =============================================================================
# Confirmation message
# =============================================================================
cat("pBACIPS functions loaded.\n")
cat("Main function: ProgressiveChangeBACIPS(control, impact, time.true, time.model)\n")
