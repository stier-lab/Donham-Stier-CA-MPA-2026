# =============================================================================
# 01_utils.R
# =============================================================================
#
# PURPOSE:
#   Shared utility functions and constants used across the MPA kelp forest
#   analysis pipeline, including core pBACIPS scripts and resilience modules.
#
# CONTENTS:
#   Section 4: EXCLUDED_MPAS, EXCLUDED_KFM_SITES, SHEEPHEAD_ONLY_MPAS
#              -> Which MPAs are dropped from analysis and why
#   Section 8: calculate_effect_size_from_contrast()
#              -> Calculates effect sizes using emmeans (covariance-aware)
#   Section 9: should_render(), save_fig(), scale_y_rr(), scale_x_rr()
#              -> Figure rendering, RR axis scales, and display helpers
#   Section 10: MPA name abbreviations for figure labels
#
# HOW EFFECT SIZES WORK:
#   For a linear model (lnRR ~ time), the effect size is:
#     predicted(t=11) - predicted(t=0)
#   This is the estimated change in log response ratio over 11 years of
#   MPA protection. emmeans::pairs() computes this with proper covariance
#   handling (predictions at t=0 and t=11 share model parameters).
#
# NOTE: Data-processing functions (biomass conversion, bootstrap, response
#   ratios, species standardization) are in the sibling data-processing repo
#   (Donham-Stier-CA-MPA-Data-2026).
#
# USAGE:
#   source(here::here("code", "R", "01_utils.R"))
#
# AUTHORS: Emily Donham & Adrian Stier
# PROJECT: CA MPA Kelp Forest pBACIPS Analysis
# =============================================================================


# =============================================================================
# SECTION 4: SITE/MPA EXCLUSION CONSTANTS
# =============================================================================
# These constants define which MPAs and sites are excluded from analysis.
# Used by 08_effect_sizes.R and 12_results_summary.R.

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

EXCLUDED_MPAS <- c(
  # NOTE: names match the spelling in THIS repo's inputs, including the Landsat
  # product (e.g. "Carrington Pt SMR", "San Miguel Island SC", "Judith Rk SMR").
  # The data-processing repo excludes the in-water spelling ("Carrington Point
  # SMR"); both refer to the same excluded MPAs (verified excluded from final results).
  "Carrington Pt SMR",
  "N/A",
  "Arrow Point to Lion Head Point SMCA",
  "Crystal Cove SMCA",
  "Laguna Beach SMR",
  "South La Jolla SMR",
  "Vandenberg SMR",
  "Point Conception SMR",
  "Anacapa Island SMR 1978",
  "Painted Cave SMCA",
  "Anacapa Island SMCA",
  # Additional MPAs excluded in 08_effect_sizes.R (lines 1583-1588):
  # San Miguel Island SC (unusual location), Judith Rk SMR (overlaps San Miguel Island SC)
  "San Miguel Island SC",
  "Judith Rk SMR"
)

EXCLUDED_KFM_SITES <- c(
  "a-k-21", "a-k-05", "a-k-07", "a-k-26", "a-k-27", "a-k-28",
  "a-k-29", "a-k-30", "a-k-36", "a-k-37", "a-k-12", "a-k-13",
  "a-k-31", "a-k-35"
)

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
# The function below uses emmeans::pairs() to properly handle covariance
# between before and after predictions from the same model.

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
  # Note: pairs() works on emmeans objects (method dispatch), can't be called as emmeans::pairs()
  contrast_result <- pairs(em, reverse = TRUE)  # after - before
  contrast_summary <- summary(contrast_result)

  # Extract results
  mean_es <- contrast_summary$estimate[1]
  se_es <- contrast_summary$SE[1]
  df_es <- contrast_summary$df[1]
  ci_es <- se_es * qt(0.975, df_es)  # Use t-distribution for CI

 list(mean = mean_es, SE = se_es, CI = ci_es, df = df_es)
}


# =============================================================================
# SECTION 9: FIGURE RENDERING AND DISPLAY UTILITIES
# =============================================================================

# --- Selective rendering control ---
# Set RENDER_FIGURES before sourcing to render a subset (e.g., c("fig03"))
# Default: render all figures
should_render <- function(fig_name) {
  if (!exists("RENDER_FIGURES", envir = .GlobalEnv)) return(TRUE)
  rf <- get("RENDER_FIGURES", envir = .GlobalEnv)
  identical(rf, "all") || fig_name %in% rf
}

# --- save_fig(): save a ggplot/patchwork figure to PDF + PNG ---
save_fig <- function(plot, name, w, h, dpi = 600) {
  pdf_path <- here::here("plots", paste0(name, ".pdf"))
  png_path <- here::here("plots", paste0(name, ".png"))

  # Detect if this is a patchwork plot (complex multi-panel)
  is_patchwork <- inherits(plot, "patchwork") || (inherits(plot, "gg") && !is.null(plot$patches))

  # Remove old PDF to prevent stale files if a save strategy silently fails
  if (file.exists(pdf_path)) file.remove(pdf_path)

  # Save PDF with error handling and fallback strategies
  pdf_success <- FALSE

  # Strategy 1: For patchwork plots, use pdf() device (most reliable)
  if (is_patchwork && !pdf_success) {
    pdf_success <- tryCatch({
      pdf(pdf_path, width = w / 2.54, height = h / 2.54, bg = "white")
      print(plot)
      dev.off()
      file.exists(pdf_path) && file.size(pdf_path) > 0
    }, error = function(e) {
      if (dev.cur() > 1) dev.off()
      FALSE
    })
  }

  # Strategy 2: Try cairo_pdf via ggsave (catch warnings too. cairo can
  # fail with a warning rather than an error, leaving the file unwritten)
  if (!pdf_success && capabilities("cairo")) {
    pdf_success <- tryCatch(
      withCallingHandlers({
        ggsave(pdf_path, plot, width = w, height = h, units = "cm",
               device = cairo_pdf, bg = "white", limitsize = FALSE)
        file.exists(pdf_path) && file.size(pdf_path) > 0
      },
      warning = function(w) {
        if (grepl("cairo", conditionMessage(w), ignore.case = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }),
      error = function(e) { FALSE }
    )
  }

  # Strategy 3: Try standard pdf device via ggsave
  if (!pdf_success) {
    pdf_success <- tryCatch({
      ggsave(pdf_path, plot, width = w, height = h, units = "cm",
             device = "pdf", bg = "white", limitsize = FALSE)
      file.exists(pdf_path) && file.size(pdf_path) > 0
    }, error = function(e) { FALSE })
  }

  # Strategy 4: Last resort - use pdf() device directly
  if (!pdf_success) {
    pdf_success <- tryCatch({
      pdf(pdf_path, width = w / 2.54, height = h / 2.54, bg = "white")
      print(plot)
      dev.off()
      file.exists(pdf_path) && file.size(pdf_path) > 0
    }, error = function(e) {
      if (dev.cur() > 1) dev.off()
      warning("All PDF save strategies failed for ", name, ": ", e$message)
      FALSE
    })
  }

  # Save PNG with error handling
  png_success <- tryCatch({
    ggsave(png_path, plot, width = w, height = h, units = "cm",
           dpi = dpi, bg = "white", limitsize = FALSE)
    TRUE
  }, error = function(e) {
    warning("Failed to save PNG for ", name, ": ", e$message)
    FALSE
  })

  # Verify files were created successfully
  pdf_exists <- file.exists(pdf_path) && file.size(pdf_path) > 0
  png_exists <- file.exists(png_path) && file.size(png_path) > 0

  # Report results
  if (pdf_exists && png_exists) {
    pdf_size <- format(file.size(pdf_path) / 1024, digits = 1, nsmall = 1)
    png_size <- format(file.size(png_path) / 1024, digits = 1, nsmall = 1)
    cat(sprintf("  Saved: %s (PDF: %s KB, PNG: %s KB @ %d DPI)\n",
                name, pdf_size, png_size, dpi))
  } else if (!png_exists) {
    stop("CRITICAL: Failed to create PNG for ", name,
         " at ", png_path, "\n",
         "  Check disk space and write permissions.")
  } else if (!pdf_exists) {
    png_size <- format(file.size(png_path) / 1024, digits = 1, nsmall = 1)
    cat(sprintf("  Saved: %s (PNG: %s KB @ %d DPI) - PDF generation failed\n",
                name, png_size, dpi))
    warning("PDF generation failed for ", name, " but PNG was created. ",
            "This is a known issue with some complex ggplot2 figures. ",
            "PNG can be converted to PDF externally if needed.")
  }

  invisible(list(pdf = pdf_path, png = png_path))
}

# --- shorten_mpa_name(): abbreviate MPA names for display ---
shorten_mpa_name <- function(mpa_name) {
  replacements <- c(
    " SMCA| SMR| SC| 2003" = "",
    "Anacapa Island"       = "Anacapa Is.",
    "Santa Barbara Island" = "Santa Barbara Is.",
    "San Miguel Island"    = "San Miguel Is.",
    "Campus Point"         = "Campus Pt.",
    "Point Vicente"        = "Pt. Vicente",
    "Harris Point"         = "Harris Pt.",
    "South Point"          = "South Pt.",
    "Carrington Pt"        = "Carrington Pt.",
    "Skunk Pt"             = "Skunk Pt.",
    "Gull Island"          = "Gull Is."
  )
  stringr::str_replace_all(mpa_name, replacements)
}


# --- scale_y_rr(): back-transformed Response Ratio y-axis ---
# Data stays on lnRR scale; tick labels show RR values at log-spaced positions.
# Provides a consistent presentation across all figures that display effect sizes.
#
# @param y_lo  Lower bound of visible lnRR range
# @param y_hi  Upper bound of visible lnRR range
# @param name  Axis title (default: "MPA / Reference")
# @return      A scale_y_continuous() layer
scale_y_rr <- function(y_lo, y_hi,
                       name = "MPA / Reference",
                       sparse = FALSE) {
  rr_pool <- if (sparse) {
    c(0.001, 0.01, 0.1, 1, 8, 100, 1000)
  } else {
    c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8, 20, 100, 500, 1000)
  }
  rr_vis <- rr_pool[log(rr_pool) >= y_lo & log(rr_pool) <= y_hi]
  if (length(rr_vis) == 0) rr_vis <- 1
  # Ensure outermost ticks bracket the full data range: add nearest pool
  # values just outside (or at) the axis limits
  lo_candidates <- rr_pool[log(rr_pool) <= y_lo]
  hi_candidates <- rr_pool[log(rr_pool) >= y_hi]
  if (length(lo_candidates) > 0) rr_vis <- union(tail(lo_candidates, 1), rr_vis)
  if (length(hi_candidates) > 0) rr_vis <- union(rr_vis, hi_candidates[1])
  rr_vis <- sort(rr_vis)
  # Cap at 5 ticks max to prevent overlap in compact panels
  if (length(rr_vis) > 5) {
    idx <- round(seq(1, length(rr_vis), length.out = 5))
    rr_vis <- rr_vis[idx]
  }
  ggplot2::scale_y_continuous(
    breaks = log(rr_vis),
    labels = as.character(rr_vis),
    name   = name
  )
}

# Identical helper for x-axis (used in bivariate scatter plots)
scale_x_rr <- function(x_lo, x_hi,
                       name = "MPA / Reference",
                       sparse = FALSE) {
  rr_pool <- if (sparse) {
    c(0.01, 0.1, 0.5, 1, 2, 8, 100)
  } else {
    c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8, 20, 100, 500)
  }
  rr_vis <- rr_pool[log(rr_pool) >= x_lo & log(rr_pool) <= x_hi]
  if (length(rr_vis) == 0) rr_vis <- 1
  # Ensure outermost ticks bracket the full data range: add nearest pool
  # values just outside (or at) the axis limits
  lo_candidates <- rr_pool[log(rr_pool) <= x_lo]
  hi_candidates <- rr_pool[log(rr_pool) >= x_hi]
  if (length(lo_candidates) > 0) rr_vis <- union(tail(lo_candidates, 1), rr_vis)
  if (length(hi_candidates) > 0) rr_vis <- union(rr_vis, hi_candidates[1])
  rr_vis <- sort(rr_vis)
  # Cap at 5 ticks max to prevent overlap in compact panels
  if (length(rr_vis) > 5) {
    idx <- round(seq(1, length(rr_vis), length.out = 5))
    rr_vis <- rr_vis[idx]
  }
  ggplot2::scale_x_continuous(
    breaks = log(rr_vis),
    labels = as.character(rr_vis),
    name   = name
  )
}


# --- scale_y_rr_free(): RR-labelled y-axis for facets with scales = "free" ---
# Uses a breaks function so ggplot2 auto-selects RR ticks per panel.
# Thins to max_breaks when the range is wide to prevent label overlap.
scale_y_rr_free <- function(name = "MPA / Reference",
                            max_breaks = 6) {
  rr_pool <- c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8, 20, 100, 500)
  lnrr_pool <- log(rr_pool)
  ggplot2::scale_y_continuous(
    breaks = function(limits) {
      sel <- lnrr_pool[lnrr_pool >= limits[1] & lnrr_pool <= limits[2]]
      if (length(sel) == 0) return(0)
      # Thin evenly if too many breaks
      iter <- 0
      while (length(sel) > max_breaks && iter < 20) {
        iter <- iter + 1
        # Always keep first, last, and 0 (RR=1); drop every other
        keep <- rep(c(TRUE, FALSE), length.out = length(sel))
        keep[1] <- TRUE
        keep[length(sel)] <- TRUE
        # Always keep RR = 1 (lnRR = 0) if present
        zero_idx <- which(abs(sel) < 1e-8)
        if (length(zero_idx) > 0) keep[zero_idx] <- TRUE
        sel <- sel[keep]
      }
      sel
    },
    labels = function(breaks) {
      rr_vals <- exp(breaks)
      ifelse(rr_vals >= 1, as.character(round(rr_vals)), as.character(rr_vals))
    },
    name = name
  )
}

# =============================================================================
# SECTION 11: RESILIENCE MODULE HELPERS (scripts 14-23)
# =============================================================================
# Shared loaders and helpers for the resilience suite. Constants they rely on
# (RESILIENCE_TAXA, SOCAL_MAX_LAT, MHW*_YEARS, ...) live in 00c_analysis_constants.R,
# which resilience scripts source before this file.

#' Load the harmonized raw abundance table (density/biomass per MPA x year x taxon
#' x status x source). Single source of truth for the file path.
load_harmonized_raw <- function() {
  read.csv(here::here("data", "harmonized", "harmonized_raw_responses.csv"),
           stringsAsFactors = FALSE)
}

#' Load the harmonized response-ratio table (lnDiff = ln[MPA / Reference]).
load_harmonized_rr <- function() {
  read.csv(here::here("data", "harmonized", "harmonized_response_ratios.csv"),
           stringsAsFactors = FALSE)
}

#' Response metric for a taxon in the resilience analyses: biomass for giant kelp,
#' density for the animals. Mirrors the RESILIENCE_RESP_OF named vector.
resilience_resp <- function(tx) if (tx == "Macrocystis pyrifera") "Bio" else "Den"

#' Mean of `values` over rows whose `years` fall in `target_years` (NA if none).
#' Replaces the per-script win()/winb() window helpers.
window_mean <- function(values, years, target_years) {
  v <- values[years %in% target_years]
  if (length(v)) mean(v) else NA_real_
}

#' Assert that a site-metadata frame lies within the Southern California Bight
#' scope (<= SOCAL_MAX_LAT). Used by scripts 14 and 19.
assert_socal_scope <- function(meta, lat_col = "Lat") {
  stopifnot("Non-Southern-California sites (lat > SOCAL_MAX_LAT)" =
              all(meta[[lat_col]] <= SOCAL_MAX_LAT, na.rm = TRUE))
  invisible(TRUE)
}

# =============================================================================
# Confirmation message
# =============================================================================
cat("Utility functions and constants loaded.\n")
