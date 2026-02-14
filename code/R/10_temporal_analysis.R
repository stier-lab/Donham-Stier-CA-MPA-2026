# =============================================================================
# 10_temporal_analysis.R
# =============================================================================
# PURPOSE:
#   Supplemental temporal dynamics analysis for the Conservation Letters manuscript.
#   Explores how MPA effects develop over time, complementing the main analysis
#   which standardizes effect sizes at t=11 years.
#
# ANALYSES:
#   1. Species-level GAM recovery curves with MPA-level spaghetti
#   2. Temporal meta-regression (lmer with random slopes)
#   2b. AR1 temporal autocorrelation sensitivity analysis (nlme::lme)
#   3. Trophic cascade phase portrait
#   4. Triptych space-time heatmap (extends Fig S4 to all trophic levels)
#   5. Rate-of-change comparison and cascade consistency scores
#
# INPUTS:
#   - All.RR.sub.trans: Annual log response ratios (from 07_combine_data.R)
#   - SumStats.Final: Effect size estimates (from 08_effect_sizes.R)
#   - Site: MPA metadata (from 03_data_import.R)
#   - trophic_assignment: Species-to-trophic-level map (from 00c_analysis_constants.R)
#   - Color palette objects from 00b_color_palette.R
#
# OUTPUTS:
#   Figures: fig_s03_recovery_curves, fig_s04_cascade_phase,
#            fig_s05_triptych_heatmap, fig_s06_slope_comparison
#   Tables:  table_s_temporal_meta_regression.csv, table_s_cascade_consistency.csv,
#            table_s_ar1_sensitivity.csv, table_s_temporal_meta_regression_ar1.csv
#
# AUTHORS: Emily Donham & Adrian Stier
# =============================================================================


# =============================================================================
# Section A: Setup
# =============================================================================

dir.create(here::here("plots"), showWarnings = FALSE)

# --- Selective rendering control ---
# Set RENDER_FIGURES before sourcing to render a subset (e.g., c("fig_s03"))
# Default: render all figures
if (!exists("RENDER_FIGURES", envir = .GlobalEnv)) {
  RENDER_FIGURES <- "all"
}
# --- Figure dimension constants (Conservation Letters supplemental) ---
FIG_S03_DIMS <- c(w = 17.8, h = 22)    # Recovery curves: 5 species panels
FIG_S04_DIMS <- c(w = 17, h = 15)      # Phase portrait: 4-panel (2x2)
FIG_S05_DIMS <- c(w = 17.8, h = 28)    # Triptych heatmap: 5 species panels
FIG_S06_DIMS <- c(w = 17, h = 12)      # Slope comparison

# --- Validate required data objects ---
required_objects <- c("All.RR.sub.trans", "SumStats.Final", "Site",
                      "trophic_assignment", "col_taxa", "theme_mpa")
missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), envir = globalenv())
]
if (length(missing_objects) > 0) {
  stop("Missing required objects: ", paste(missing_objects, collapse = ", "),
       "\nPlease run scripts 00-09 first.")
}

# --- Detect taxa column name in All.RR.sub.trans ---
if ("y" %in% names(All.RR.sub.trans)) {
  taxa_col <- "y"
} else if ("Taxa" %in% names(All.RR.sub.trans)) {
  taxa_col <- "Taxa"
} else {
  taxa_col <- names(All.RR.sub.trans)[1]
  warning("Could not find 'y' or 'Taxa' column; falling back to '",
          taxa_col, "'")
}

cat("  taxa_col detected as: '", taxa_col, "'\n", sep = "")

# --- Full-name to abbreviated-name mapping (used across all analyses) ---
full_to_abbrev <- c(
  "Macrocystis pyrifera"          = "M. pyrifera",
  "Mesocentrotus franciscanus"    = "M. franciscanus",
  "Strongylocentrotus purpuratus" = "S. purpuratus",
  "Panulirus interruptus"         = "P. interruptus",
  "Semicossyphus pulcher"         = "S. pulcher"
)

# Species-level colors keyed by full name
species_colors_full <- setNames(
  col_taxa[full_to_abbrev],
  names(full_to_abbrev)
)

# Ordered species list: predators → urchins → kelp
species_order_full <- c(
  "Panulirus interruptus", "Semicossyphus pulcher",
  "Strongylocentrotus purpuratus", "Mesocentrotus franciscanus",
  "Macrocystis pyrifera"
)
species_order_abbrev <- full_to_abbrev[species_order_full]

# Trophic level colors (kept for any legacy use)
trophic_colors <- c(
  "Predators" = "#5C8A70",
  "Urchins"   = "#956079",
  "Kelp"      = unname(col_taxa["M. pyrifera"])
)

# Expected direction for cascade consistency (per species)
expected_direction <- c(
  "Panulirus interruptus"         = "positive",
  "Semicossyphus pulcher"         = "positive",
  "Strongylocentrotus purpuratus" = "negative",
  "Mesocentrotus franciscanus"    = "negative",
  "Macrocystis pyrifera"          = "positive"
)

# --- Build temporal working dataset (After period only) ---
# NOTE: This dataset includes both biomass (Bio) and density (Den) response
# types. The main meta-analysis (09_meta_analysis.R) runs separate models for
# Bio and Den because they are non-independent measures of the same populations.
# Here we include both but account for the non-independence by including resp
# as a random effect in the lmer model (Analysis 2). For descriptive figures
# (GAMs, phase portraits), we note that both response types are pooled.
temporal_after <- All.RR.sub.trans %>%
  dplyr::filter(BA == "After", time >= 0) %>%
  dplyr::mutate(
    Trophic_Level = trophic_assignment[.data[[taxa_col]]],
    Species = .data[[taxa_col]],
    Species_abbrev = full_to_abbrev[.data[[taxa_col]]],
    time = as.numeric(time)
  ) %>%
  dplyr::filter(!is.na(Trophic_Level))

cat(sprintf("  Temporal dataset: %d observations, %d MPAs, %d species\n",
            nrow(temporal_after),
            length(unique(temporal_after$CA_MPA_Name_Short)),
            length(unique(temporal_after[[taxa_col]]))))

cat("=== Temporal analysis setup complete ===\n")


# =============================================================================
# Section B: Analysis 1 -- Species-Level GAM Recovery Curves
# =============================================================================
# Figure S3: Individual species recovery trajectories over time since MPA
# implementation. Faint spaghetti lines show MPA-level lnRR trajectories;
# bold GAM smooth reveals the population-level trend.

if (should_render("fig_s03")) {
cat("\n--- Figure S3: Recovery Curves ---\n")

# Require mgcv for GAM fitting
if (!requireNamespace("mgcv", quietly = TRUE)) {
  stop("Package 'mgcv' is required for GAM recovery curves. ",
       "Install with: install.packages('mgcv')")
}

# Species present in the temporal data
species_in_data <- unique(temporal_after[[taxa_col]])
cat(sprintf("  Species in temporal data: %s\n",
            paste(species_in_data, collapse = ", ")))

# Build color mapping from full species names
species_colors <- species_colors_full[species_in_data]
# Fill in any unmatched species with a neutral grey
missing_color <- is.na(species_colors)
if (any(missing_color)) {
  species_colors[missing_color] <- "grey40"
  warning("No col_taxa entry for: ",
          paste(species_in_data[missing_color], collapse = ", "))
}

# Construct the faceted figure
# Use linetype aesthetic to create a legend differentiating individual MPAs
# from the population-level GAM smooth
fig_s03 <- ggplot(temporal_after,
                  aes(x = time, y = lnDiff)) +
  # Individual MPA trajectories (spaghetti)
  geom_line(aes(group = CA_MPA_Name_Short, color = .data[[taxa_col]],
                linetype = "Individual MPA"),
            alpha = 0.25, linewidth = 0.4, show.legend = TRUE) +
  # Reference line at zero effect
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
             linewidth = 0.4) +
  # Population-level GAM smooth with 95% CI ribbon
  geom_smooth(aes(color = .data[[taxa_col]], fill = .data[[taxa_col]],
                  linetype = "GAM smooth (95% CI)"),
              method = "gam",
              formula = y ~ s(x, k = 5),
              linewidth = 1.1, alpha = 0.2,
              se = TRUE, level = 0.95) +
  # Facet by species
  facet_wrap(as.formula(paste("~", taxa_col)),
             ncol = 2, scales = "free_y") +
  # Use species-specific colors (suppress their legends; color = facet label)
  scale_color_manual(values = species_colors, guide = "none") +
  scale_fill_manual(values = species_colors, guide = "none") +
  # Linetype legend for the two line types
  scale_linetype_manual(
    name = NULL,
    values = c("Individual MPA" = "solid", "GAM smooth (95% CI)" = "solid"),
    guide = guide_legend(
      override.aes = list(
        linewidth = c(0.4, 1.1),
        alpha = c(0.35, 1.0),
        color = c("grey50", "grey20")
      )
    )
  ) +
  scale_x_continuous(breaks = seq(0, 20, by = 5)) +
  labs(
    x = "Years since MPA implementation",
    y = "Log response ratio (lnRR)"
  ) +
  theme_mpa(base_size = 9) +
  theme(
    strip.text = element_text(face = "italic", size = 9),
    panel.grid.major = element_blank(),
    plot.margin = margin(6, 8, 6, 6, "pt"),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.margin = margin(t = 0)
  )

save_fig(fig_s03, "fig_s03_recovery_curves",
         w = FIG_S03_DIMS["w"], h = FIG_S03_DIMS["h"])

# --- GAM non-linearity diagnostics ---
cat("\n  GAM EDF diagnostics (non-linearity test):\n")
for (sp in species_in_data) {
  sp_data <- dplyr::filter(temporal_after, .data[[taxa_col]] == sp)
  if (nrow(sp_data) >= 10) {
    gam_fit <- tryCatch(
      mgcv::gam(lnDiff ~ s(time, k = 5), data = sp_data),
      error = function(e) NULL
    )
    if (!is.null(gam_fit)) {
      s_table <- summary(gam_fit)$s.table
      cat(sprintf("  %s: EDF = %.2f (p = %.4f)\n",
                  sp, s_table[1, "edf"], s_table[1, "p-value"]))
    } else {
      cat(sprintf("  %s: GAM fitting failed\n", sp))
    }
  } else {
    cat(sprintf("  %s: Too few observations (n = %d, need >= 10)\n",
                sp, nrow(sp_data)))
  }
}

} # end fig_s03


# =============================================================================
# Section C: Analysis 2 -- Temporal Meta-Regression (Species-Level)
# =============================================================================
# Fits a multilevel model to test whether individual species differ in their
# rate of change (slope) over time since MPA implementation.
# Uses species (5 levels) rather than trophic level to preserve the distinct
# recovery trajectories of each predator and urchin species.

cat("\n--- Temporal Meta-Regression (Species-Level) ---\n")

if (!requireNamespace("lme4", quietly = TRUE)) {
  stop("Package 'lme4' is required. Install with: install.packages('lme4')")
}

# Factor species for interpretable contrasts (reference = first predator)
temporal_after$Species <- factor(
  temporal_after$Species,
  levels = species_order_full
)

# Attempt full random-slopes model first
cat("  Fitting random-slopes model (Species * time)...\n")
temporal_lmer <- tryCatch({
  lme4::lmer(
    lnDiff ~ Species * time +
      (1 + time | CA_MPA_Name_Short) + (1 | source) + (1 | resp),
    data = temporal_after,
    control = lme4::lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 20000)
    )
  )
}, warning = function(w) {
  cat("  WARNING (random slopes): ", conditionMessage(w), "\n")
  tryCatch(
    suppressWarnings(
      lme4::lmer(
        lnDiff ~ Species * time +
          (1 + time | CA_MPA_Name_Short) + (1 | source) + (1 | resp),
        data = temporal_after,
        control = lme4::lmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = 20000)
        )
      )
    ),
    error = function(e) NULL
  )
}, error = function(e) {
  cat("  ERROR (random slopes): ", conditionMessage(e), "\n")
  NULL
})

# Check convergence; fall back to random-intercepts if needed
use_fallback <- FALSE
if (is.null(temporal_lmer)) {
  use_fallback <- TRUE
} else if (length(temporal_lmer@optinfo$conv$lme4) > 0) {
  cat("  Convergence issues detected; falling back to random-intercepts model\n")
  use_fallback <- TRUE
}

if (use_fallback) {
  cat("  Fitting random-intercepts fallback model...\n")
  temporal_lmer <- tryCatch(
    lme4::lmer(
      lnDiff ~ Species * time +
        (1 | CA_MPA_Name_Short) + (1 | source) + (1 | resp),
      data = temporal_after
    ),
    error = function(e) {
      cat("  ERROR (random intercepts): ", conditionMessage(e), "\n")
      NULL
    }
  )
}

# Check for singular fit (near-zero variance components)
if (!is.null(temporal_lmer) && lme4::isSingular(temporal_lmer)) {
  cat("  NOTE: Model has singular fit (near-zero variance for one or more random effects).\n")
  cat("  Model variant:", ifelse(use_fallback, "random-intercepts", "random-slopes"), "\n")
}

if (!is.null(temporal_lmer)) {

  # Extract fixed effects using lmerTest for p-values (Satterthwaite)
  if (requireNamespace("lmerTest", quietly = TRUE)) {
    lmer_summary <- summary(lmerTest::as_lmerModLmerTest(temporal_lmer))
  } else {
    lmer_summary <- summary(temporal_lmer)
  }

  fixed_ef <- as.data.frame(coef(lmer_summary))
  fixed_ef$Term <- rownames(fixed_ef)
  rownames(fixed_ef) <- NULL

  if ("Pr(>|t|)" %in% names(fixed_ef)) {
    names(fixed_ef)[names(fixed_ef) == "Pr(>|t|)"] <- "p_value"
  }
  names(fixed_ef)[names(fixed_ef) == "Estimate"] <- "Estimate"
  names(fixed_ef)[names(fixed_ef) == "Std. Error"] <- "SE"
  names(fixed_ef)[names(fixed_ef) == "t value"] <- "t_value"

  col_order <- c("Term", "Estimate", "SE", "t_value")
  if ("df" %in% names(fixed_ef)) col_order <- c(col_order, "df")
  if ("p_value" %in% names(fixed_ef)) col_order <- c(col_order, "p_value")
  remaining <- setdiff(names(fixed_ef), col_order)
  fixed_ef <- fixed_ef[, c(col_order, remaining)]

  # Marginal and conditional R-squared
  r2_text <- ""
  if (requireNamespace("MuMIn", quietly = TRUE)) {
    r2 <- tryCatch(MuMIn::r.squaredGLMM(temporal_lmer), error = function(e) NULL)
    if (!is.null(r2)) {
      r2_text <- sprintf("    Marginal R2: %.4f, Conditional R2: %.4f\n",
                         r2[1, "R2m"], r2[1, "R2c"])
    }
  }

  # Save fixed effects table
  write.csv(fixed_ef,
            here::here("data", "table_s_temporal_meta_regression.csv"),
            row.names = FALSE)
  cat("  Saved: data/table_s_temporal_meta_regression.csv\n")

  # --- Print species-level slopes ---
  cat("\n  Temporal Meta-Regression Results (Species-Level):\n")
  if (use_fallback) {
    cat("    (Random-intercepts model; random slopes did not converge)\n")
  } else {
    cat("    (Random-slopes model converged)\n")
  }

  # Reference species slope (P. interruptus = first level)
  ref_slope <- fixed_ef$Estimate[fixed_ef$Term == "time"]
  ref_se    <- fixed_ef$SE[fixed_ef$Term == "time"]

  # Extract vcov matrix for proper SE of summed coefficients
  vcov_mat <- as.matrix(vcov(temporal_lmer))
  term_names <- rownames(vcov_mat)

  for (sp in species_order_full) {
    sp_abbrev <- full_to_abbrev[sp]
    if (sp == species_order_full[1]) {
      # Reference level
      cat(sprintf("    %s: slope = %.4f +/- %.4f\n", sp_abbrev, ref_slope, ref_se))
    } else {
      # Clean the species name for matching term names (R replaces spaces/dots)
      sp_clean <- gsub(" ", ".", sp)
      interaction_term <- grep(paste0("Species", sp_clean, ":time|time:Species", sp_clean),
                               term_names, value = TRUE)
      if (length(interaction_term) > 0) {
        interaction_term <- interaction_term[1]
        sp_slope <- ref_slope + fixed_ef$Estimate[fixed_ef$Term == interaction_term]
        # Correct SE: Var(a+b) = Var(a) + Var(b) + 2*Cov(a,b)
        ref_term <- "time"
        var_sum <- vcov_mat[ref_term, ref_term] +
                   vcov_mat[interaction_term, interaction_term] +
                   2 * vcov_mat[ref_term, interaction_term]
        sp_se <- sqrt(max(var_sum, 0))  # guard against tiny negative from numerics
        cat(sprintf("    %s: slope = %.4f +/- %.4f\n", sp_abbrev, sp_slope, sp_se))
      }
    }
  }
  if (nchar(r2_text) > 0) cat(r2_text)

} else {
  cat("  WARNING: Temporal meta-regression model failed. Skipping.\n")
}


# =============================================================================
# Section C2: Temporal Autocorrelation Sensitivity Analysis (AR1 correction)
# =============================================================================
# MOTIVATION: The lmer model above assumes independent residuals within each
# MPA time series. However, annual ecological observations spanning 10-20 years
# likely exhibit positive temporal autocorrelation (from environmental regime
# shifts, multi-year recruitment pulses, persistent demographic structure).
# Positive autocorrelation inflates effective sample size, deflating standard
# errors and potentially overstating statistical significance.
#
# APPROACH: We fit a parallel model using nlme::lme() with an AR1 correlation
# structure (corAR1) that estimates and accounts for first-order autocorrelation
# within each MPA's time series. We compare the AR1-corrected SEs and slopes
# against the original lmer to assess sensitivity of conclusions.
#
# NOTE: nlme is a base R recommended package (always available). The lme()
# formula syntax differs from lmer(): random effects use a list specification
# and nested random effects require explicit nesting.
# =============================================================================

cat("\n--- AR1 Temporal Autocorrelation Sensitivity Analysis ---\n")

if (!is.null(temporal_lmer)) {

  # --- Step 1: Durbin-Watson diagnostic on lmer residuals ---
  # Provides a quick test for autocorrelation in the pooled residuals.
  # DW ~ 2 indicates no autocorrelation; DW < 2 indicates positive autocorrelation.
  cat("  Durbin-Watson test on lmer residuals:\n")
  dw_result <- tryCatch({
    resids <- residuals(temporal_lmer)
    # Order residuals by MPA then time for meaningful DW test
    resid_df <- data.frame(
      resid = resids,
      MPA = temporal_after$CA_MPA_Name_Short,
      time = temporal_after$time
    )
    resid_df <- resid_df[order(resid_df$MPA, resid_df$time), ]
    # Compute DW statistic on ordered residuals
    n <- nrow(resid_df)
    dw_stat <- sum(diff(resid_df$resid)^2) / sum(resid_df$resid^2)
    cat(sprintf("    DW statistic (ordered by MPA, time) = %.4f\n", dw_stat))
    cat(sprintf("    Interpretation: %s\n",
                ifelse(dw_stat < 1.5, "Substantial positive autocorrelation (DW << 2)",
                       ifelse(dw_stat < 1.8, "Moderate positive autocorrelation (DW < 2)",
                              ifelse(dw_stat > 2.5, "Negative autocorrelation (DW > 2)",
                                     "Weak or no autocorrelation (DW ~ 2)")))))
    dw_stat
  }, error = function(e) {
    cat("    DW test failed:", conditionMessage(e), "\n")
    NA_real_
  })

  # --- Step 2: Fit nlme::lme with AR1 correlation structure ---
  cat("\n  Fitting nlme::lme model with AR1 correction...\n")

  # Ensure required columns are not missing
  temporal_after$source_f <- factor(temporal_after$source)
  temporal_after$resp_f   <- factor(temporal_after$resp)
  temporal_after$MPA_f    <- factor(temporal_after$CA_MPA_Name_Short)

  # Remove any rows with NA in key columns (nlme is less tolerant than lmer)
  ar1_data <- temporal_after %>%
    dplyr::filter(!is.na(lnDiff), !is.na(Species), !is.na(time),
                  !is.na(MPA_f), !is.na(source_f), !is.na(resp_f))

  # Try progressively simpler random effects structures until convergence
  temporal_lme_ar1 <- NULL
  ar1_model_description <- ""

  # Attempt 1: Full model matching lmer — random slopes for MPA + random intercepts
  # for source and resp. nlme requires nested or crossed random effects specified

  # via pdBlocked or groupedData; for crossed designs we nest source and resp
  # within a dummy grouping variable. However, nlme only supports nested (not
  # fully crossed) random effects. We use the most important grouping (MPA with
  # random slopes) plus source as a nested or additional intercept.
  #
  # Strategy: Use MPA as the primary grouping with random slopes, and add
  # source and resp as fixed covariates rather than random effects if the full
  # random structure won't converge. This is a pragmatic compromise — the key
  # goal is estimating the AR1 coefficient (rho).

  temporal_lme_ar1 <- tryCatch({
    cat("    Attempt 1: Random slopes (MPA) + random intercepts (source, resp)...\n")
    nlme::lme(
      fixed = lnDiff ~ Species * time,
      random = list(
        source_f = nlme::pdDiag(~ 1),
        MPA_f = ~ 1 + time,
        resp_f = nlme::pdDiag(~ 1)
      ),
      correlation = nlme::corAR1(form = ~ time | MPA_f),
      data = ar1_data,
      method = "REML",
      control = nlme::lmeControl(
        maxIter = 200, msMaxIter = 200,
        opt = "optim", msVerbose = FALSE,
        returnObject = TRUE
      )
    )
  }, error = function(e) {
    cat("    Attempt 1 failed:", conditionMessage(e), "\n")
    NULL
  })
  if (!is.null(temporal_lme_ar1)) {
    ar1_model_description <- "random slopes (MPA) + random intercepts (source, resp)"
  }

  # Attempt 2: Random slopes for MPA only (drop source and resp random effects)
  if (is.null(temporal_lme_ar1)) {
    temporal_lme_ar1 <- tryCatch({
      cat("    Attempt 2: Random slopes (MPA) only, AR1 within MPA...\n")
      nlme::lme(
        fixed = lnDiff ~ Species * time,
        random = ~ 1 + time | MPA_f,
        correlation = nlme::corAR1(form = ~ time | MPA_f),
        data = ar1_data,
        method = "REML",
        control = nlme::lmeControl(
          maxIter = 200, msMaxIter = 200,
          opt = "optim", msVerbose = FALSE,
          returnObject = TRUE
        )
      )
    }, error = function(e) {
      cat("    Attempt 2 failed:", conditionMessage(e), "\n")
      NULL
    })
    if (!is.null(temporal_lme_ar1)) {
      ar1_model_description <- "random slopes (MPA) only"
    }
  }

  # Attempt 3: Random intercepts only for MPA (simplest structure that allows AR1)
  if (is.null(temporal_lme_ar1)) {
    temporal_lme_ar1 <- tryCatch({
      cat("    Attempt 3: Random intercepts (MPA) only, AR1 within MPA...\n")
      nlme::lme(
        fixed = lnDiff ~ Species * time,
        random = ~ 1 | MPA_f,
        correlation = nlme::corAR1(form = ~ time | MPA_f),
        data = ar1_data,
        method = "REML",
        control = nlme::lmeControl(
          maxIter = 200, msMaxIter = 200,
          opt = "optim", msVerbose = FALSE,
          returnObject = TRUE
        )
      )
    }, error = function(e) {
      cat("    Attempt 3 failed:", conditionMessage(e), "\n")
      NULL
    })
    if (!is.null(temporal_lme_ar1)) {
      ar1_model_description <- "random intercepts (MPA) only"
    }
  }

  # --- Step 3: Compare lmer vs AR1 model ---
  if (!is.null(temporal_lme_ar1)) {
    cat(sprintf("    AR1 model converged (%s)\n", ar1_model_description))

    # Extract estimated AR1 coefficient (rho)
    ar1_rho <- tryCatch({
      cs <- coef(temporal_lme_ar1$modelStruct$corStruct,
                 unconstrained = FALSE)
      cs[1]
    }, error = function(e) NA_real_)
    cat(sprintf("    Estimated AR1 coefficient (rho) = %.4f\n",
                ifelse(is.na(ar1_rho), NA_real_, ar1_rho)))
    if (!is.na(ar1_rho)) {
      cat(sprintf("    Interpretation: %s\n",
                  ifelse(abs(ar1_rho) < 0.1,
                         "Negligible autocorrelation; lmer results are robust",
                         ifelse(ar1_rho > 0.3,
                                "Substantial positive autocorrelation; lmer SEs likely underestimated",
                                ifelse(ar1_rho > 0.1,
                                       "Moderate positive autocorrelation; some SE inflation possible",
                                       "Negative autocorrelation (unusual in ecology)")))))
    }

    # Extract AR1 fixed effects for comparison
    ar1_summary <- summary(temporal_lme_ar1)
    ar1_fixed <- as.data.frame(ar1_summary$tTable)
    ar1_fixed$Term <- rownames(ar1_fixed)
    rownames(ar1_fixed) <- NULL
    names(ar1_fixed)[names(ar1_fixed) == "Value"] <- "Estimate"
    names(ar1_fixed)[names(ar1_fixed) == "Std.Error"] <- "SE"
    names(ar1_fixed)[names(ar1_fixed) == "t-value"] <- "t_value"
    names(ar1_fixed)[names(ar1_fixed) == "p-value"] <- "p_value"

    # --- Compare slopes side-by-side ---
    cat("\n  Comparison of species slopes: lmer vs AR1-corrected lme\n")
    cat("  ---------------------------------------------------------------\n")
    cat(sprintf("  %-20s %10s %8s | %10s %8s | %s\n",
                "Species", "lmer_slope", "lmer_SE",
                "AR1_slope", "AR1_SE", "SE_ratio"))
    cat("  ---------------------------------------------------------------\n")

    # Build comparison table for export
    comparison_rows <- list()

    for (sp in species_order_full) {
      sp_abbrev <- full_to_abbrev[sp]

      # lmer slope and SE
      if (sp == species_order_full[1]) {
        lmer_sl <- fixed_ef$Estimate[fixed_ef$Term == "time"]
        lmer_se <- fixed_ef$SE[fixed_ef$Term == "time"]
      } else {
        sp_clean <- gsub(" ", ".", sp)
        int_term <- grep(paste0("Species", sp_clean, ":time|time:Species", sp_clean),
                         fixed_ef$Term, value = TRUE)
        if (length(int_term) > 0) {
          lmer_sl <- fixed_ef$Estimate[fixed_ef$Term == "time"] +
                     fixed_ef$Estimate[fixed_ef$Term == int_term[1]]
          # Proper SE from vcov
          vm <- as.matrix(vcov(temporal_lmer))
          v <- vm["time", "time"] + vm[int_term[1], int_term[1]] +
               2 * vm["time", int_term[1]]
          lmer_se <- sqrt(max(v, 0))
        } else {
          lmer_sl <- NA_real_; lmer_se <- NA_real_
        }
      }

      # AR1 slope and SE
      if (sp == species_order_full[1]) {
        ar1_sl <- ar1_fixed$Estimate[ar1_fixed$Term == "time"]
        ar1_se <- ar1_fixed$SE[ar1_fixed$Term == "time"]
      } else {
        sp_clean <- gsub(" ", ".", sp)
        int_term_ar1 <- grep(paste0("Species", sp_clean, ":time|time:Species", sp_clean),
                             ar1_fixed$Term, value = TRUE)
        if (length(int_term_ar1) > 0) {
          ar1_sl <- ar1_fixed$Estimate[ar1_fixed$Term == "time"] +
                    ar1_fixed$Estimate[ar1_fixed$Term == int_term_ar1[1]]
          # SE from nlme vcov
          vm_ar1 <- vcov(temporal_lme_ar1)
          v_ar1 <- vm_ar1["time", "time"] + vm_ar1[int_term_ar1[1], int_term_ar1[1]] +
                   2 * vm_ar1["time", int_term_ar1[1]]
          ar1_se <- sqrt(max(v_ar1, 0))
        } else {
          ar1_sl <- NA_real_; ar1_se <- NA_real_
        }
      }

      se_ratio <- ifelse(!is.na(lmer_se) && !is.na(ar1_se) && lmer_se > 0,
                         ar1_se / lmer_se, NA_real_)
      cat(sprintf("  %-20s %10.4f %8.4f | %10.4f %8.4f | %.2f\n",
                  sp_abbrev,
                  ifelse(is.na(lmer_sl), NA_real_, lmer_sl),
                  ifelse(is.na(lmer_se), NA_real_, lmer_se),
                  ifelse(is.na(ar1_sl), NA_real_, ar1_sl),
                  ifelse(is.na(ar1_se), NA_real_, ar1_se),
                  ifelse(is.na(se_ratio), NA_real_, se_ratio)))

      comparison_rows[[sp_abbrev]] <- data.frame(
        Species = sp_abbrev,
        lmer_slope = lmer_sl, lmer_SE = lmer_se,
        AR1_slope = ar1_sl, AR1_SE = ar1_se,
        SE_ratio = se_ratio,
        stringsAsFactors = FALSE
      )
    }
    cat("  ---------------------------------------------------------------\n")
    cat("  SE_ratio > 1 means AR1 model has larger (more conservative) SEs\n")

    comparison_df <- dplyr::bind_rows(comparison_rows)

    # Assess whether conclusions change
    # If the AR1 SE_ratio is substantially > 1 for most species, flag it
    median_se_ratio <- median(comparison_df$SE_ratio, na.rm = TRUE)
    cat(sprintf("\n  Median SE ratio (AR1 / lmer) = %.2f\n", median_se_ratio))
    if (median_se_ratio > 1.3) {
      cat("  CAUTION: AR1 correction substantially inflates SEs (>30%).\n")
      cat("  Temporal autocorrelation is non-trivial; consider using AR1 model.\n")
    } else if (median_se_ratio > 1.1) {
      cat("  NOTE: AR1 correction moderately inflates SEs (10-30%).\n")
      cat("  Results are qualitatively robust but SEs should be interpreted with care.\n")
    } else {
      cat("  GOOD: AR1 correction has minimal effect on SEs (<10%).\n")
      cat("  Temporal autocorrelation does not materially affect conclusions.\n")
    }

    # --- Save AR1 comparison table ---
    ar1_comparison_out <- comparison_df
    ar1_comparison_out$AR1_rho <- ar1_rho
    ar1_comparison_out$DW_statistic <- ifelse(!is.na(dw_result), dw_result, NA_real_)
    ar1_comparison_out$AR1_model <- ar1_model_description
    write.csv(ar1_comparison_out,
              here::here("data", "table_s_ar1_sensitivity.csv"),
              row.names = FALSE)
    cat("  Saved: data/table_s_ar1_sensitivity.csv\n")

    # --- Save AR1 fixed effects table (parallel to lmer table) ---
    ar1_fixed_out <- ar1_fixed[, intersect(c("Term", "Estimate", "SE",
                                              "t_value", "DF", "p_value"),
                                            names(ar1_fixed))]
    write.csv(ar1_fixed_out,
              here::here("data", "table_s_temporal_meta_regression_ar1.csv"),
              row.names = FALSE)
    cat("  Saved: data/table_s_temporal_meta_regression_ar1.csv\n")

  } else {
    cat("    WARNING: All AR1 model attempts failed to converge.\n")
    cat("    Reporting Durbin-Watson diagnostic only.\n")
    if (!is.na(dw_result)) {
      cat(sprintf("    DW = %.4f — %s\n", dw_result,
                  ifelse(dw_result < 1.5,
                         "suggests positive autocorrelation (lmer SEs may be underestimated)",
                         "suggests autocorrelation is not severe")))
    }
    # Save a minimal diagnostic CSV
    dw_out <- data.frame(
      diagnostic = "Durbin-Watson",
      statistic = ifelse(!is.na(dw_result), dw_result, NA_real_),
      note = ifelse(!is.na(dw_result) && dw_result < 1.5,
                    "Positive autocorrelation detected; AR1 model did not converge",
                    "Weak autocorrelation; lmer results likely robust"),
      stringsAsFactors = FALSE
    )
    write.csv(dw_out,
              here::here("data", "table_s_ar1_sensitivity.csv"),
              row.names = FALSE)
    cat("  Saved: data/table_s_ar1_sensitivity.csv (DW diagnostic only)\n")
  }

} else {
  cat("  Skipping AR1 analysis — lmer model did not fit.\n")
}

cat("=== AR1 sensitivity analysis complete ===\n")


# =============================================================================
# Section D: Analysis 3 -- Trophic Cascade Phase Portrait (Species-Level)
# =============================================================================
# Figure S4: Four-panel species-level phase portrait showing how individual
# species effects co-evolve over time. Each panel pairs one species against
# M. pyrifera (kelp) or against an urchin species, revealing whether
# cascade sequencing operates at the species level.

if (should_render("fig_s04")) {
cat("\n--- Figure S4: Cascade Phase Portrait (Species-Level) ---\n")

library(patchwork)

# Aggregate yearly means by species using two-stage approach:
# Stage 1: compute MPA-level means (avoids data-rich programs dominating)
# Stage 2: average across MPAs
phase_yearly <- temporal_after %>%
  dplyr::group_by(Species, time, CA_MPA_Name_Short) %>%
  dplyr::summarise(mpa_mean = mean(lnDiff, na.rm = TRUE), .groups = "drop") %>%
  dplyr::group_by(Species, time) %>%
  dplyr::summarise(
    mean_lnRR = mean(mpa_mean, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::filter(n >= 2) %>%
  dplyr::arrange(time)

# Pivot wide: one column per species
phase_wide <- phase_yearly %>%
  dplyr::mutate(Species_abbrev = full_to_abbrev[Species]) %>%
  dplyr::select(Species_abbrev, time, mean_lnRR) %>%
  tidyr::pivot_wider(names_from = Species_abbrev, values_from = mean_lnRR)

# Helper: build one phase portrait panel
build_phase_panel <- function(data, x_col, y_col, x_lab, y_lab,
                              q_top_right, q_bot_left) {
  d <- data %>%
    dplyr::filter(!is.na(.data[[x_col]]), !is.na(.data[[y_col]]))
  if (nrow(d) < 3) return(NULL)

  ggplot(d, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70",
               linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey70",
               linewidth = 0.4) +
    annotate("text", x = Inf, y = Inf, label = q_top_right,
             hjust = 1.1, vjust = 1.5, size = 2.1, color = "grey60") +
    annotate("text", x = -Inf, y = -Inf, label = q_bot_left,
             hjust = -0.1, vjust = -0.5, size = 2.1, color = "grey60") +
    geom_path(color = "grey40", linewidth = 0.6,
              arrow = arrow(length = unit(1.5, "mm"), type = "closed",
                            ends = "last"),
              lineend = "round") +
    geom_point(aes(color = time), size = 2.5) +
    scale_color_viridis_c(name = "Years since\nMPA", option = "D") +
    labs(x = x_lab, y = y_lab) +
    theme_mpa(base_size = 8) +
    theme(
      legend.position = "bottom",
      legend.key.width = unit(10, "mm"),
      legend.key.height = unit(2.5, "mm"),
      panel.grid.major = element_blank()
    )
}

# Invert urchin columns (positive = urchin decline = cascade direction)
phase_wide_inv <- phase_wide %>%
  dplyr::mutate(
    `S. purpuratus_inv`  = -1 * `S. purpuratus`,
    `M. franciscanus_inv` = -1 * `M. franciscanus`
  )

# Panel (a): S. purpuratus decline vs M. pyrifera increase
pa <- build_phase_panel(
  phase_wide_inv, "S. purpuratus_inv", "M. pyrifera",
  x_lab = expression(italic("S. purpuratus") ~ "lnRR \u00d7 (-1)"),
  y_lab = expression(italic("M. pyrifera") ~ "lnRR"),
  q_top_right = "Purples \u2193  Kelp \u2191",
  q_bot_left  = "Purples \u2191  Kelp \u2193"
)

# Panel (b): M. franciscanus decline vs M. pyrifera increase
pb <- build_phase_panel(
  phase_wide_inv, "M. franciscanus_inv", "M. pyrifera",
  x_lab = expression(italic("M. franciscanus") ~ "lnRR \u00d7 (-1)"),
  y_lab = expression(italic("M. pyrifera") ~ "lnRR"),
  q_top_right = "Reds \u2193  Kelp \u2191",
  q_bot_left  = "Reds \u2191  Kelp \u2193"
)

# Panel (c): P. interruptus increase vs S. purpuratus decline
pc <- build_phase_panel(
  phase_wide_inv, "P. interruptus", "S. purpuratus_inv",
  x_lab = expression(italic("P. interruptus") ~ "lnRR"),
  y_lab = expression(italic("S. purpuratus") ~ "lnRR \u00d7 (-1)"),
  q_top_right = "Lobster \u2191  Purples \u2193",
  q_bot_left  = "Lobster \u2193  Purples \u2191"
)

# Panel (d): S. pulcher increase vs M. franciscanus decline
pd <- build_phase_panel(
  phase_wide_inv, "S. pulcher", "M. franciscanus_inv",
  x_lab = expression(italic("S. pulcher") ~ "lnRR"),
  y_lab = expression(italic("M. franciscanus") ~ "lnRR \u00d7 (-1)"),
  q_top_right = "Sheephead \u2191  Reds \u2193",
  q_bot_left  = "Sheephead \u2193  Reds \u2191"
)

# Assemble available panels
panels <- list(pa, pb, pc, pd)
panels <- panels[!vapply(panels, is.null, logical(1))]

if (length(panels) >= 2) {
  fig_s04 <- wrap_plots(panels, ncol = 2) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(
      legend.position = "bottom",
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.line.x.bottom = element_line(colour = "black", linewidth = 0.3),
      axis.line.y.left = element_line(colour = "black", linewidth = 0.3)
    )

  save_fig(fig_s04, "fig_s04_cascade_phase",
           w = FIG_S04_DIMS["w"], h = FIG_S04_DIMS["h"])
} else {
  cat("  WARNING: Too few species-level pairings for phase portrait.\n")
}

} # end fig_s04


# =============================================================================
# Analysis 4: Species-Level Space-Time Heatmap (Figure S5)
# =============================================================================
# Shows all five species as stacked heatmap panels (replaces old urchin-only version)
# stacked vertically, preserving distinct recovery patterns for each
# predator and urchin species.

if (should_render("fig_s05")) {
cat("\n--- Figure S5: Species-Level Heatmap ---\n")

# (FIG_S05_DIMS defined in header constants section)

# ---------------------------------------------------------------------------
# 4a. Prepare heatmap data by species
# ---------------------------------------------------------------------------

heatmap_data <- temporal_after %>%
  dplyr::group_by(Species, Species_abbrev, CA_MPA_Name_Short, time) %>%
  dplyr::summarise(
    mean_lnRR = mean(lnDiff, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::filter(n >= 1, time >= 0, time <= 15)

# Get MPA implementation year for ordering
mpa_years <- Site %>%
  dplyr::select(CA_MPA_Name_Short, MPA_Start) %>%
  dplyr::distinct() %>%
  dplyr::filter(!is.na(MPA_Start)) %>%
  dplyr::arrange(MPA_Start)

heatmap_data <- heatmap_data %>%
  dplyr::left_join(mpa_years, by = "CA_MPA_Name_Short") %>%
  dplyr::filter(!is.na(MPA_Start))

# Shorten MPA names
heatmap_data$CA_MPA_Name_Short <- shorten_mpa_name(heatmap_data$CA_MPA_Name_Short)
mpa_years$CA_MPA_Name_Short    <- shorten_mpa_name(mpa_years$CA_MPA_Name_Short)

# Unified MPA ordering (oldest at top)
mpa_order <- heatmap_data %>%
  dplyr::left_join(
    mpa_years %>% dplyr::distinct(CA_MPA_Name_Short, MPA_Start),
    by = "CA_MPA_Name_Short"
  ) %>%
  dplyr::distinct(CA_MPA_Name_Short,
                  MPA_Start = dplyr::coalesce(MPA_Start.x, MPA_Start.y)) %>%
  dplyr::arrange(dplyr::desc(MPA_Start)) %>%
  dplyr::pull(CA_MPA_Name_Short)

heatmap_data <- heatmap_data %>%
  dplyr::mutate(CA_MPA_Name_Short = factor(CA_MPA_Name_Short, levels = mpa_order))

cat("  Heatmap data: ", nrow(heatmap_data), " cells across ",
    length(mpa_order), " MPAs\n", sep = "")

# ---------------------------------------------------------------------------
# 4b. Build individual heatmap panels per species
# ---------------------------------------------------------------------------

color_lim <- c(-3, 3)

build_heatmap_panel <- function(sp_full, show_x_axis = FALSE,
                                show_legend = FALSE) {
  sp_abbrev <- full_to_abbrev[sp_full]
  d <- dplyr::filter(heatmap_data, Species == sp_full)

  if (nrow(d) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5,
                          label = paste0("No data for ", sp_abbrev),
                          size = 3.5) +
        theme_mpa()
    )
  }

  all_mpas  <- levels(heatmap_data$CA_MPA_Name_Short)
  all_times <- seq(min(d$time), max(d$time))
  bg_grid   <- tidyr::expand_grid(
    CA_MPA_Name_Short = factor(all_mpas, levels = all_mpas),
    time = all_times
  )

  p <- ggplot2::ggplot(d, ggplot2::aes(x = time, y = CA_MPA_Name_Short)) +
    ggplot2::geom_tile(data = bg_grid, fill = "grey85",
                       color = "white", linewidth = 0.2) +
    ggplot2::geom_tile(ggplot2::aes(fill = mean_lnRR),
                       color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(
      low = "#B2182B", mid = "white", high = "#2166AC",
      midpoint = 0, limits = color_lim, oob = scales::squish,
      na.value = "grey85", name = "lnRR"
    ) +
    ggplot2::scale_x_continuous(breaks = seq(0, 15, by = 3),
                                expand = c(0, 0)) +
    ggplot2::labs(title = sp_abbrev, y = NULL) +
    theme_mpa(base_size = 9) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 7),
      panel.grid  = ggplot2::element_blank(),
      plot.title  = ggplot2::element_text(size = 9, face = "bold.italic",
                                          hjust = 0)
    )

  if (!show_x_axis) {
    p <- p + ggplot2::labs(x = NULL) +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank())
  } else {
    p <- p + ggplot2::labs(x = "Years since MPA implementation")
  }

  if (!show_legend) {
    p <- p + ggplot2::theme(legend.position = "none")
  } else {
    p <- p + ggplot2::theme(legend.position = "right")
  }

  p
}

# ---------------------------------------------------------------------------
# 4c. Assemble 5-panel figure (predators → urchins → kelp)
# ---------------------------------------------------------------------------

n_sp <- length(species_order_full)
heatmap_panels <- lapply(seq_along(species_order_full), function(i) {
  build_heatmap_panel(
    species_order_full[i],
    show_x_axis = (i == n_sp),
    show_legend = (i == n_sp)
  )
})

fig_s05 <- Reduce(`/`, heatmap_panels) +
  patchwork::plot_layout(guides = "collect",
                         heights = c(rep(1, n_sp - 1), 1.15)) +
  patchwork::plot_annotation(
    tag_levels = "a",
    caption = "Grey cells = no data available"
  ) &
  ggplot2::theme(
    plot.tag     = ggplot2::element_text(size = 9, face = "bold"),
    plot.caption = ggplot2::element_text(size = 7, color = "grey40", hjust = 0)
  )

save_fig(fig_s05, "fig_s05_triptych_heatmap", FIG_S05_DIMS["w"], FIG_S05_DIMS["h"])
cat("  Figure S5 saved.\n")

} # end fig_s05


# =============================================================================
# Analysis 5: Rate-of-Change Comparison & Cascade Consistency (Figure S06)
# =============================================================================
# Two complementary views at the species level:
#   (a) Per-MPA slopes of lnRR ~ time, grouped by species (5 rows).
#   (b) Cascade consistency: does each MPA show the expected direction
#       (predators +, urchins -, kelp +) for each of the 5 species?

if (should_render("fig_s06")) {
cat("\n--- Figure S06: Slope Comparison & Cascade Consistency (Species-Level) ---\n")

FIG_S06_DIMS <- c(w = 17, h = 16)

# ---------------------------------------------------------------------------
# 5a. Calculate per-MPA slopes by species
# ---------------------------------------------------------------------------

# NOTE: Per-MPA slopes use OLS (lm). Standard errors may be underestimated
# if positive temporal autocorrelation is present. These are descriptive
# summaries; the formal test is the lmer temporal meta-regression above.
slope_data <- temporal_after %>%
  dplyr::group_by(CA_MPA_Name_Short, Species, Species_abbrev) %>%
  dplyr::filter(dplyr::n() >= 5) %>%
  dplyr::summarise(
    slope = tryCatch(
      coef(lm(lnDiff ~ time))[2],
      error = function(e) NA_real_
    ),
    slope_se = tryCatch(
      summary(lm(lnDiff ~ time))$coefficients[2, 2],
      error = function(e) NA_real_
    ),
    n_years = dplyr::n_distinct(time),
    .groups = "drop"
  ) %>%
  dplyr::filter(!is.na(slope))

cat("  Slopes calculated for ", nrow(slope_data), " MPA x species combos\n",
    sep = "")

# Species-level summaries
slope_summary <- slope_data %>%
  dplyr::group_by(Species, Species_abbrev) %>%
  dplyr::summarise(
    mean_slope = mean(slope, na.rm = TRUE),
    se_slope   = sd(slope, na.rm = TRUE) / sqrt(dplyr::n()),
    n_mpas     = dplyr::n(),
    .groups    = "drop"
  )

cat("  Slope summaries:\n")
for (i in seq_len(nrow(slope_summary))) {
  cat(sprintf("    %s: mean = %.4f +/- %.4f (n = %d MPAs)\n",
              slope_summary$Species_abbrev[i],
              slope_summary$mean_slope[i],
              slope_summary$se_slope[i],
              slope_summary$n_mpas[i]))
}

# ---------------------------------------------------------------------------
# 5b. Panel (a): Slopes by species — jitter + species summary
# ---------------------------------------------------------------------------

# Factor species in cascade order
slope_data$Species_abbrev    <- factor(slope_data$Species_abbrev,
                                       levels = rev(species_order_abbrev))
slope_summary$Species_abbrev <- factor(slope_summary$Species_abbrev,
                                       levels = rev(species_order_abbrev))

# Build color map keyed by abbreviated name
sp_color_abbrev <- setNames(
  unname(species_colors_full[species_order_full]),
  species_order_abbrev
)

p_slopes <- ggplot2::ggplot(slope_data,
                            ggplot2::aes(x = Species_abbrev, y = slope)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
                      linewidth = 0.4) +
  ggplot2::geom_jitter(
    ggplot2::aes(color = Species_abbrev),
    width = 0.15, size = 1.8, alpha = 0.6, shape = 16
  ) +
  ggplot2::geom_pointrange(
    data = slope_summary,
    ggplot2::aes(x = Species_abbrev, y = mean_slope,
                 ymin = mean_slope - se_slope,
                 ymax = mean_slope + se_slope),
    size = 0.6, linewidth = 0.8, color = "black",
    position = ggplot2::position_nudge(x = 0.25)
  ) +
  ggplot2::scale_color_manual(values = sp_color_abbrev, guide = "none") +
  ggplot2::coord_flip() +
  ggplot2::labs(
    x = NULL,
    y = expression("Annual rate of change (slope of lnRR" ~ "~" ~ "time)")
  ) +
  theme_mpa(base_size = 9) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 8, face = "italic")
  )

# ---------------------------------------------------------------------------
# 5c. Panel (b): Cascade consistency tile chart (species-level)
# ---------------------------------------------------------------------------

# Pivot to wide: one column per species
cascade_wide <- slope_data %>%
  dplyr::select(CA_MPA_Name_Short, Species_abbrev, slope) %>%
  tidyr::pivot_wider(names_from = Species_abbrev, values_from = slope)

# Determine which species each MPA has data for
available_species <- intersect(species_order_abbrev, names(cascade_wide))

# For cascade consistency, check expected direction per species
# expected_direction from setup: "positive" for predators/kelp, "negative" for urchins
# NOTE: Consistency is defined by the SIGN of the slope only, regardless of
# statistical significance. This is appropriate for a descriptive summary
# of MPA-level replication. The formal statistical test is the lmer model.
cascade_long_list <- list()
for (sp_full in species_order_full) {
  sp_abbrev <- full_to_abbrev[sp_full]
  if (!(sp_abbrev %in% names(cascade_wide))) next

  direction <- expected_direction[sp_full]
  cascade_long_list[[sp_abbrev]] <- cascade_wide %>%
    dplyr::filter(!is.na(.data[[sp_abbrev]])) %>%
    dplyr::mutate(
      Species_label = sp_abbrev,
      consistent = if (direction == "positive") {
        .data[[sp_abbrev]] > 0
      } else {
        .data[[sp_abbrev]] < 0
      }
    ) %>%
    dplyr::select(CA_MPA_Name_Short, Species_label, consistent)
}
cascade_long <- dplyr::bind_rows(cascade_long_list)

# Build direction labels for x-axis
direction_labels <- setNames(
  ifelse(expected_direction == "positive", "\u2191", "\u2193"),
  full_to_abbrev[names(expected_direction)]
)
cascade_long <- cascade_long %>%
  dplyr::mutate(
    Species_dir = paste0(Species_label, " ", direction_labels[Species_label]),
    Species_dir = factor(Species_dir,
                         levels = paste0(species_order_abbrev, " ",
                                         direction_labels[species_order_abbrev]))
  )

# Shorten MPA names
cascade_long$MPA_short <- shorten_mpa_name(cascade_long$CA_MPA_Name_Short)

# Order MPAs by cascade consistency score (most consistent at top)
cascade_scores <- cascade_long %>%
  dplyr::group_by(CA_MPA_Name_Short, MPA_short) %>%
  dplyr::summarise(
    cascade_score = sum(consistent, na.rm = TRUE),
    n_species = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(cascade_score, CA_MPA_Name_Short)

mpa_cascade_order_short <- cascade_scores$MPA_short
cascade_long$MPA_short <- factor(cascade_long$MPA_short,
                                 levels = mpa_cascade_order_short)

fill_cascade <- c("TRUE" = "#5C8A70", "FALSE" = "#C47272")

p_cascade <- ggplot2::ggplot(
  cascade_long,
  ggplot2::aes(x = Species_dir, y = MPA_short,
               fill = as.character(consistent))
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.5) +
  ggplot2::scale_fill_manual(
    values = fill_cascade,
    labels = c("FALSE" = "No", "TRUE" = "Yes"),
    name   = "Expected direction"
  ) +
  ggplot2::labs(x = NULL, y = NULL) +
  theme_mpa(base_size = 9) +
  ggplot2::theme(
    axis.text.y     = ggplot2::element_text(size = 7),
    axis.text.x     = ggplot2::element_text(size = 6.5, face = "italic",
                                            angle = 45, hjust = 1),
    legend.position = "bottom",
    panel.grid      = ggplot2::element_blank()
  )

# ---------------------------------------------------------------------------
# 5d. Assemble 2-panel figure
# ---------------------------------------------------------------------------

fig_s06 <- p_slopes + p_cascade +
  patchwork::plot_layout(widths = c(1, 1.3)) +
  patchwork::plot_annotation(tag_levels = "a") &
  ggplot2::theme(
    plot.tag = ggplot2::element_text(size = 9, face = "bold")
  )

save_fig(fig_s06, "fig_s06_slope_comparison", FIG_S06_DIMS["w"], FIG_S06_DIMS["h"])
cat("  Figure S06 saved.\n")

# ---------------------------------------------------------------------------
# 5e. Save cascade consistency table (species-level)
# ---------------------------------------------------------------------------

cascade_out <- cascade_long %>%
  tidyr::pivot_wider(
    id_cols = c(CA_MPA_Name_Short, MPA_short),
    names_from = Species_label,
    values_from = consistent
  ) %>%
  dplyr::left_join(cascade_scores %>% dplyr::select(CA_MPA_Name_Short,
                                                     cascade_score, n_species),
                   by = "CA_MPA_Name_Short") %>%
  dplyr::arrange(dplyr::desc(cascade_score), CA_MPA_Name_Short)

write.csv(cascade_out,
          here::here("data", "table_s_cascade_consistency.csv"),
          row.names = FALSE)
cat("  Table saved: data/table_s_cascade_consistency.csv\n")

} # end fig_s06


# =============================================================================
# Completion
# =============================================================================
cat("\n=== 10_temporal_analysis.R complete ===\n")
