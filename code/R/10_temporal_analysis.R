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
#   4. [DROPPED] Triptych space-time heatmap (disabled; redundant with S3+S6)
#   5. Rate-of-change comparison and cascade consistency scores
#
# INPUTS:
#   - All.RR.sub.trans: Annual log response ratios (from 03_load_harmonized_data.R)
#   - SumStats.Final: Effect size estimates (from 08_effect_sizes.R)
#   - Site: MPA metadata (from 03_load_harmonized_data.R)
#   - trophic_assignment: Species-to-trophic-level map (from 00c_analysis_constants.R)
#   - Color palette objects from 00b_color_palette.R
#
# OUTPUTS:
#   Figures: fig_s04_recovery_curves, fig_s05_cascade_phase,
#            fig_s06_slope_comparison (heatmap dropped)
#   Tables:  table_s_temporal_meta_regression.csv      — pooled lmer fixed effects
#            table_s_temporal_meta_regression_bio.csv   — biomass lmer fixed effects
#            table_s_temporal_meta_regression_den.csv   — density lmer fixed effects
#            table_s_lmer_model_fit.csv                 — model-level stats (AIC, BIC, R2, var components)
#            table_s_species_slopes.csv                 — per-species derived slopes with SE & 95% CI
#            table_s_gam_linearity_diagnostics.csv      — GAM EDF and non-linearity classification
#            table_s_ar1_sensitivity.csv                — AR1 vs lmer slope/SE comparison (or DW only)
#            table_s_temporal_meta_regression_ar1.csv   — AR1 fixed effects (when AR1 converges)
#            table_s_cascade_consistency.csv            — per-MPA cascade direction consistency
#            table_s_per_mpa_slopes.csv                 — per-MPA OLS slopes by species
#            table_s_slope_summary_by_species.csv       — species-level mean slopes across MPAs
#            table_s_phase_portrait_yearly_means.csv    — yearly mean lnRR by species (for phase plots)
#            table_s_phase_portrait_correlations.csv    — pairwise species correlations from phase portraits
#
# AUTHORS: Emily Donham & Adrian Stier
# =============================================================================


# =============================================================================
# Section A: Setup
# =============================================================================

dir.create(here::here("plots"), showWarnings = FALSE)

# --- Selective rendering control ---
# Set RENDER_FIGURES before sourcing to render a subset (e.g., c("fig_s04"))
# Default: render all figures
if (!exists("RENDER_FIGURES", envir = .GlobalEnv)) {
  RENDER_FIGURES <- "all"
}
# --- Figure dimension constants (Conservation Letters supplemental) ---
FIG_S04_DIMS <- c(w = 17.8, h = 25)    # Recovery curves: 5 species × 2 response types
FIG_S05_DIMS <- c(w = 17, h = 16)      # Phase portrait: 4-panel (2x2)
# FIG_S06_DIMS <- c(w = 17.8, h = 28)    # Triptych heatmap: 5 species panels  [DROPPED]
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
# NOTE: The pooled model does not include resp as either a fixed or random
# effect (k=2 levels is too few for a random effect, and the separate
# response-type models in Section C1b provide properly stratified analysis).
# Bio/Den non-independence is a known limitation documented in 09_meta_analysis.R.
# For descriptive figures (GAMs, phase portraits), we note that both response
# types are pooled.
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
# Figure S4: Individual species recovery trajectories over time since MPA
# implementation. Faint spaghetti lines show MPA-level lnRR trajectories;
# bold GAM smooth reveals the population-level trend.

if (should_render("fig_s04")) {
cat("\n--- Figure S4: Recovery Curves ---\n")

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

# Factor taxa column in canonical trophic order for consistent facet ordering
temporal_after[[taxa_col]] <- factor(
  temporal_after[[taxa_col]],
  levels = species_order_full
)

# Create labelled copy for SI Fig S3 (disk: fig_s04) facet strips
# Cap at t<=15 to match Figure 4 (main text) time range
# Keep resp column for Bio/Den faceting
s03_data <- temporal_after %>%
  dplyr::filter(time <= 15, !is.na(resp), resp %in% c("Bio", "Den")) %>%
  dplyr::mutate(
    resp_label = dplyr::case_when(
      resp == "Bio" ~ "Biomass",
      resp == "Den" ~ "Density",
      TRUE ~ resp
    )
  )

# Remove M. pyrifera Density rows (no density data for kelp)
s03_data <- s03_data %>%
  dplyr::filter(!(Species == "Macrocystis pyrifera" & resp == "Den"))

# Create species factor with canonical ordering
s03_data$Taxa_Full_factor <- factor(
  s03_data[[taxa_col]], levels = species_order_full
)
s03_data$resp_label <- factor(s03_data$resp_label, levels = c("Biomass", "Density"))

cat(sprintf("  SI Fig S3 (fig_s04) data: %d obs, species x resp combinations:\n", nrow(s03_data)))
cat("  ", paste(unique(paste(s03_data$Species, s03_data$resp_label)), collapse = ", "), "\n")

# --- Load response-specific lmer slopes for linear overlay ---
# Helper: extract per-species intercepts and slopes from a response-specific CSV
build_lmer_preds <- function(csv_path, sp_order, resp_val) {
  if (!file.exists(csv_path)) {
    cat(sprintf("  NOTE: %s not found; skipping lmer overlay for %s\n",
                basename(csv_path), resp_val))
    return(NULL)
  }
  coefs <- read.csv(csv_path, stringsAsFactors = FALSE)
  ref_int   <- coefs$Estimate[coefs$Term == "(Intercept)"]
  ref_slope <- coefs$Estimate[coefs$Term == "time"]
  if (length(ref_int) == 0 || length(ref_slope) == 0) return(NULL)

  sp_ints   <- setNames(numeric(length(sp_order)), sp_order)
  sp_slopes <- setNames(numeric(length(sp_order)), sp_order)
  sp_ints[sp_order[1]]   <- ref_int
  sp_slopes[sp_order[1]] <- ref_slope
  for (sp in sp_order[-1]) {
    sp_clean <- gsub(" ", ".", sp)
    sp_main  <- paste0("Species", sp_clean)
    sp_iterm <- paste0("Species", sp_clean, ":time")
    main_row <- coefs[coefs$Term == sp_main, ]
    int_row  <- coefs[coefs$Term == sp_iterm, ]
    sp_ints[sp]   <- ref_int   + ifelse(nrow(main_row) > 0, main_row$Estimate, 0)
    sp_slopes[sp] <- ref_slope + ifelse(nrow(int_row) > 0, int_row$Estimate, 0)
  }

  # Build prediction data.frame (time 0-15)
  pred_df <- do.call(rbind, lapply(sp_order, function(sp) {
    data.frame(
      time = seq(0, 15, length.out = 50),
      lnDiff_pred = sp_ints[sp] + sp_slopes[sp] * seq(0, 15, length.out = 50),
      taxa_label = sp,
      resp_label = resp_val,
      stringsAsFactors = FALSE
    )
  }))
  pred_df
}

# Build predictions from response-specific CSVs
# NOTE: On first pipeline run, these CSVs don't exist yet (written in Section C1b below).
# build_lmer_preds() returns NULL gracefully; lmer overlay will appear on subsequent runs.
bio_csv <- here::here("tables", "table_s_temporal_meta_regression_bio.csv")
den_csv <- here::here("tables", "table_s_temporal_meta_regression_den.csv")
density_species <- species_order_full[species_order_full != "Macrocystis pyrifera"]

lmer_pred_bio <- build_lmer_preds(bio_csv, species_order_full, "Biomass")
lmer_pred_den <- build_lmer_preds(den_csv, density_species, "Density")

# Combine into a single prediction data.frame
lmer_pred_df <- dplyr::bind_rows(lmer_pred_bio, lmer_pred_den)
if (!is.null(lmer_pred_df) && nrow(lmer_pred_df) > 0) {
  names(lmer_pred_df)[names(lmer_pred_df) == "taxa_label"] <- taxa_col
  lmer_pred_df[[taxa_col]] <- factor(lmer_pred_df[[taxa_col]], levels = species_order_full)
  lmer_pred_df$Taxa_Full_factor <- factor(lmer_pred_df[[taxa_col]], levels = species_order_full)
  lmer_pred_df$resp_label <- factor(lmer_pred_df$resp_label, levels = c("Biomass", "Density"))
  cat("  Loaded response-specific lmer slopes for linear overlay on Fig S4\n")
} else {
  lmer_pred_df <- NULL
  cat("  NOTE: No response-specific lmer CSVs found; skipping linear overlay\n")
}

# Create combined facet label: "Species — Response" for facet_wrap
s03_data <- s03_data %>%
  dplyr::mutate(
    facet_label = paste0(.data[[taxa_col]], " \u2014 ", resp_label)
  )
# Order: species trophic order × biomass-first within species
facet_order <- as.vector(t(outer(
  species_order_full,
  c("Biomass", "Density"),
  function(sp, r) paste0(sp, " \u2014 ", r)
)))
# Drop M. pyrifera density (no data)
facet_order <- facet_order[facet_order != "Macrocystis pyrifera \u2014 Density"]
s03_data$facet_label <- factor(s03_data$facet_label, levels = facet_order)

# Also add facet labels to lmer predictions
if (!is.null(lmer_pred_df)) {
  lmer_pred_df <- lmer_pred_df %>%
    dplyr::mutate(
      facet_label = factor(
        paste0(.data[[taxa_col]], " \u2014 ", resp_label),
        levels = facet_order
      )
    )
}

# Construct the faceted figure: 2-column facet_wrap (no empty panel)
fig_s03 <- ggplot(s03_data,
                  aes(x = time, y = lnDiff)) +
  # Reference line at RR = 1 (no MPA effect)
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30",
             linewidth = 0.4) +
  # Individual MPA trajectories (spaghetti) — species-colored
  geom_line(aes(group = interaction(CA_MPA_Name_Short, resp),
                color = .data[[taxa_col]]),
            alpha = 0.28, linewidth = 0.5) +
  # Population-level GAM smooth with 95% CI ribbon
  geom_smooth(aes(color = .data[[taxa_col]], fill = .data[[taxa_col]]),
              method = "gam",
              formula = y ~ s(x, k = 5),
              linewidth = 1.2, alpha = 0.20,
              se = TRUE, level = 0.95) +
  # lmer linear prediction overlay (dashed black line)
  {if (!is.null(lmer_pred_df))
    geom_line(data = lmer_pred_df,
              aes(x = time, y = lnDiff_pred),
              color = "grey15", linewidth = 0.7, alpha = 0.8,
              linetype = "dashed",
              inherit.aes = FALSE)
  } +
  facet_wrap(~ facet_label, ncol = 2, scales = "free_y") +
  scale_color_manual(values = species_colors, guide = "none") +
  scale_fill_manual(values = species_colors, guide = "none") +
  scale_x_continuous(breaks = seq(0, 15, by = 5), limits = c(0, 16),
                     expand = c(0.02, 0)) +
  scale_y_rr_free() +
  labs(
    x = "Years since MPA implementation",
    y = "MPA / Reference"
  ) +
  theme_mpa(base_size = 9) +
  theme(
    strip.text = element_text(face = "italic", size = 9,
                              margin = margin(b = 3, t = 3)),
    axis.text  = element_text(size = 7),
    axis.title = element_text(size = 8),
    panel.grid.major = element_blank(),
    plot.margin = margin(6, 6, 6, 6, "pt"),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.margin = margin(t = 0)
  )

# Manual legend via cowplot
legend_grob <- cowplot::ggdraw() +
  cowplot::draw_line(x = c(0.06, 0.14), y = c(0.5, 0.5),
                     color = "grey30", linewidth = 1.2) +
  cowplot::draw_label("GAM smooth (95% CI)", x = 0.22, y = 0.5,
                      size = 8, fontfamily = "Helvetica", hjust = 0) +
  cowplot::draw_line(x = c(0.42, 0.50), y = c(0.5, 0.5),
                     color = "grey55", linewidth = 0.4) +
  cowplot::draw_label("Individual MPA", x = 0.52, y = 0.5,
                      size = 8, fontfamily = "Helvetica", hjust = 0) +
  cowplot::draw_line(x = c(0.68, 0.76), y = c(0.5, 0.5),
                     color = "grey15", linewidth = 0.7, linetype = "dashed") +
  cowplot::draw_label("lmer prediction", x = 0.78, y = 0.5,
                      size = 8, fontfamily = "Helvetica", hjust = 0)

# Placeholder for empty bottom-right cell (M. pyrifera has no density data)
note_grob <- cowplot::ggdraw() +
  cowplot::draw_label("No density metric\nfor M. pyrifera",
                      x = 0.5, y = 0.5, size = 8, color = "grey50",
                      fontface = "italic", fontfamily = "Helvetica")

# Combine main plot with note in bottom-right via inset
fig_s03_annotated <- cowplot::ggdraw(fig_s03) +
  cowplot::draw_plot(note_grob, x = 0.55, y = 0.035, width = 0.4, height = 0.1)

fig_s03_final <- cowplot::plot_grid(
  fig_s03_annotated, legend_grob, ncol = 1, rel_heights = c(1, 0.03)
)

save_fig(fig_s03_final, "fig_s04_recovery_curves",
         w = FIG_S04_DIMS["w"], h = FIG_S04_DIMS["h"])

} # end fig_s04

# --- GAM non-linearity diagnostics (per species per response type) ---
# NOTE: This table (Table S7) is generated unconditionally — it should always
# be written regardless of whether fig_s04 is rendered, since it is a
# standalone supplemental table used in the manuscript.
cat("\n  GAM EDF diagnostics (non-linearity test, per species x response):\n")

if (!requireNamespace("mgcv", quietly = TRUE)) {
  cat("  Skipping GAM EDF diagnostics: package 'mgcv' not available.\n")
} else {

# Prepare data for GAM diagnostics (self-contained, not dependent on fig_s04 rendering)
gam_diag_species <- unique(temporal_after[[taxa_col]])
gam_diag_data <- temporal_after %>%
  dplyr::filter(time <= 15, !is.na(resp), resp %in% c("Bio", "Den")) %>%
  dplyr::filter(!(Species == "Macrocystis pyrifera" & resp == "Den"))

edf_results <- list()
resp_types <- c("Bio", "Den")
resp_labels_map <- c(Bio = "Biomass", Den = "Density")

for (sp in gam_diag_species) {
  for (rt in resp_types) {
    # Skip M. pyrifera density (no data)
    if (sp == "Macrocystis pyrifera" && rt == "Den") next

    sp_rt_data <- dplyr::filter(gam_diag_data, .data[[taxa_col]] == sp, resp == rt)
    result_key <- paste(sp, rt, sep = "_")

    if (nrow(sp_rt_data) >= 10) {
      gam_fit <- tryCatch(
        mgcv::gam(lnDiff ~ s(time, k = 5), data = sp_rt_data),
        error = function(e) NULL
      )
      if (!is.null(gam_fit)) {
        s_table <- summary(gam_fit)$s.table
        cat(sprintf("  %s [%s]: EDF = %.2f (p = %.4f)\n",
                    sp, resp_labels_map[rt], s_table[1, "edf"], s_table[1, "p-value"]))
        edf_results[[result_key]] <- data.frame(
          Species = sp,
          Resp = rt,
          EDF = round(s_table[1, "edf"], 3),
          Ref_df = round(s_table[1, "Ref.df"], 3),
          F_stat = round(s_table[1, "F"], 3),
          p_value = s_table[1, "p-value"],
          n_obs = nrow(sp_rt_data),
          linear = ifelse(s_table[1, "edf"] < 1.5, "yes", "no"),
          stringsAsFactors = FALSE
        )
      } else {
        cat(sprintf("  %s [%s]: GAM fitting failed\n", sp, resp_labels_map[rt]))
      }
    } else {
      cat(sprintf("  %s [%s]: Too few observations (n = %d, need >= 10)\n",
                  sp, resp_labels_map[rt], nrow(sp_rt_data)))
    }
  }
}

# Save EDF diagnostics to CSV
if (length(edf_results) > 0) {
  edf_table <- do.call(rbind, edf_results)
  edf_csv_path <- here::here("tables", "table_s_gam_linearity_diagnostics.csv")
  write.csv(edf_table, edf_csv_path, row.names = FALSE)
  cat(sprintf("  GAM EDF diagnostics saved to: %s\n", edf_csv_path))
}

} # end mgcv check


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

# ---------------------------------------------------------------------------
# Helper: extract model-level fit statistics and per-species derived slopes
# from a fitted lmer model. Returns a list with two data.frames:
#   $fit_stats  — one-row summary (AIC, BIC, logLik, R2, variance, N, etc.)
#   $species_slopes — per-species total slopes with proper SEs and 95% CIs
# ---------------------------------------------------------------------------
extract_lmer_diagnostics <- function(lmer_obj, species_levels, full_to_abbrev,
                                      model_label = "pooled",
                                      is_fallback = FALSE) {
  # --- Model fit statistics ---
  aic_val  <- AIC(lmer_obj)
  bic_val  <- BIC(lmer_obj)
  ll_val   <- as.numeric(logLik(lmer_obj))
  n_obs    <- nobs(lmer_obj)
  mpa_col <- if ("CA_MPA_Name_Short" %in% names(lmer_obj@frame)) {
    lmer_obj@frame$CA_MPA_Name_Short
  } else {
    lmer_obj@frame$MPA_f
  }
  n_mpa <- length(unique(mpa_col))

  # R-squared (marginal and conditional)
  r2m <- r2c <- NA_real_
  if (requireNamespace("MuMIn", quietly = TRUE)) {
    r2 <- tryCatch(MuMIn::r.squaredGLMM(lmer_obj), error = function(e) NULL)
    if (!is.null(r2)) {
      r2m <- r2[1, "R2m"]
      r2c <- r2[1, "R2c"]
    }
  }

  # Random effects variance components
  vc <- as.data.frame(lme4::VarCorr(lmer_obj))
  vc_summary <- paste(
    sprintf("%s|%s: %.4f", vc$grp, vc$var1, vc$vcov),
    collapse = "; "
  )

  # Singular fit flag
  is_singular <- lme4::isSingular(lmer_obj)

  fit_stats <- data.frame(
    model = model_label,
    model_type = ifelse(is_fallback, "random-intercepts", "random-slopes"),
    n_obs = n_obs,
    n_mpa = n_mpa,
    AIC = round(aic_val, 2),
    BIC = round(bic_val, 2),
    logLik = round(ll_val, 2),
    R2_marginal = round(r2m, 4),
    R2_conditional = round(r2c, 4),
    variance_components = vc_summary,
    singular_fit = is_singular,
    stringsAsFactors = FALSE
  )

  # --- Per-species derived slopes with proper SEs ---
  # Get fixed effects with Satterthwaite p-values if available
  if (requireNamespace("lmerTest", quietly = TRUE)) {
    sm <- summary(lmerTest::as_lmerModLmerTest(lmer_obj))
  } else {
    sm <- summary(lmer_obj)
  }
  fe <- as.data.frame(coef(sm))
  fe$Term <- rownames(fe)
  if ("Pr(>|t|)" %in% names(fe)) names(fe)[names(fe) == "Pr(>|t|)"] <- "p_value"
  names(fe)[names(fe) == "Std. Error"] <- "SE"
  names(fe)[names(fe) == "t value"] <- "t_value"

  vcov_mat <- as.matrix(vcov(lmer_obj))
  term_names <- rownames(vcov_mat)

  ref_slope <- fe$Estimate[fe$Term == "time"]
  ref_se    <- fe$SE[fe$Term == "time"]
  ref_df    <- if ("df" %in% names(fe)) fe$df[fe$Term == "time"] else NA_real_
  ref_p     <- if ("p_value" %in% names(fe)) fe$p_value[fe$Term == "time"] else NA_real_

  slope_rows <- list()
  for (sp in species_levels) {
    sp_abbrev <- full_to_abbrev[sp]
    if (sp == species_levels[1]) {
      sp_slope <- ref_slope
      sp_se    <- ref_se
      sp_df    <- ref_df
      sp_t     <- sp_slope / sp_se
      sp_p     <- ref_p
    } else {
      sp_clean <- gsub(" ", ".", sp)
      int_term <- grep(paste0("Species", sp_clean, ":time|time:Species", sp_clean),
                       term_names, value = TRUE)
      if (length(int_term) == 0) next
      int_term <- int_term[1]
      sp_slope <- ref_slope + fe$Estimate[fe$Term == int_term]
      v_sum    <- vcov_mat["time", "time"] + vcov_mat[int_term, int_term] +
                  2 * vcov_mat["time", int_term]
      sp_se    <- sqrt(max(v_sum, 0))
      sp_t     <- sp_slope / sp_se
      sp_df    <- if ("df" %in% names(fe)) fe$df[fe$Term == int_term] else NA_real_
      # p-value from t-distribution if df available
      sp_p     <- if (!is.na(sp_df)) 2 * pt(abs(sp_t), df = sp_df, lower.tail = FALSE) else NA_real_
    }
    slope_rows[[sp_abbrev]] <- data.frame(
      model = model_label,
      Species = sp,
      Species_abbrev = sp_abbrev,
      slope = round(sp_slope, 6),
      SE = round(sp_se, 6),
      t_value = round(sp_t, 4),
      df = round(sp_df, 2),
      p_value = sp_p,
      CI_lower = round(sp_slope - qt(0.975, if (is.na(sp_df)) Inf else sp_df) * sp_se, 6),
      CI_upper = round(sp_slope + qt(0.975, if (is.na(sp_df)) Inf else sp_df) * sp_se, 6),
      stringsAsFactors = FALSE
    )
  }
  species_slopes <- dplyr::bind_rows(slope_rows)

  list(fit_stats = fit_stats, species_slopes = species_slopes)
}

# Factor species for interpretable contrasts (reference = first predator)
temporal_after$Species <- factor(
  temporal_after$Species,
  levels = species_order_full
)

# Attempt full random-slopes model first
cat("  Fitting random-slopes model (Species * time)...\n")
temporal_lmer <- tryCatch({
  withCallingHandlers(
    lme4::lmer(
      # Note: (1|resp) removed — with only k=2 levels (Bio/Den), this RE is unidentifiable
      lnDiff ~ Species * time +
        (1 + time | CA_MPA_Name_Short) + (1 | source),
      data = temporal_after,
      control = lme4::lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 20000)
      )
    ),
    warning = function(w) {
      cat("  WARNING (random slopes): ", conditionMessage(w), "\n")
      invokeRestart("muffleWarning")
    }
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
      # Note: (1|resp) removed — with only k=2 levels (Bio/Den), this RE is unidentifiable
      lnDiff ~ Species * time +
        (1 | CA_MPA_Name_Short) + (1 | source),
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
            here::here("tables", "table_s_temporal_meta_regression.csv"),
            row.names = FALSE)
  cat("  Saved: tables/table_s_temporal_meta_regression.csv\n")

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

  # --- Export model-level diagnostics and per-species derived slopes ---
  pooled_diag <- extract_lmer_diagnostics(
    temporal_lmer, species_order_full, full_to_abbrev,
    model_label = "pooled", is_fallback = use_fallback
  )
  # Accumulate across all models (pooled, bio, den) for a combined export
  all_fit_stats     <- list(pooled_diag$fit_stats)
  all_species_slopes <- list(pooled_diag$species_slopes)

} else {
  cat("  WARNING: Temporal meta-regression model failed. Skipping.\n")
  all_fit_stats     <- list()
  all_species_slopes <- list()
}


# =============================================================================
# Section C1b: Response-Specific Temporal Meta-Regression (Biomass & Density)
# =============================================================================
# Fits separate lmer models for biomass-only and density-only data.
# These complement the pooled model above by allowing response-type-specific
# slopes to be overlaid on Figure S4's biomass/density facets.
# NOTE: M. pyrifera has no density data, so the density model includes only
# 4 species (P. interruptus, S. pulcher, S. purpuratus, M. franciscanus).

cat("\n--- Response-Specific Temporal Meta-Regression (Bio & Den) ---\n")

# Helper: fit a response-specific lmer, extract coefficients, save CSV
fit_resp_lmer <- function(resp_type, resp_label, species_levels) {
  cat(sprintf("\n  Fitting %s-only model...\n", resp_label))

  resp_data <- temporal_after %>%
    dplyr::filter(resp == resp_type)

  cat(sprintf("    %s observations: %d\n", resp_label, nrow(resp_data)))

  if (nrow(resp_data) < 20) {
    cat(sprintf("    Too few observations for %s model. Skipping.\n", resp_label))
    return(NULL)
  }

  # Set species factor levels (subset to species present in this response type)
  sp_present <- intersect(species_levels, unique(as.character(resp_data$Species)))
  resp_data$Species <- factor(resp_data$Species, levels = sp_present)
  cat(sprintf("    Species in %s data: %s\n", resp_label,
              paste(sp_present, collapse = ", ")))

  # Attempt random-slopes model
  resp_lmer <- tryCatch({
    withCallingHandlers(
      lme4::lmer(
        lnDiff ~ Species * time +
          (1 + time | CA_MPA_Name_Short) + (1 | source),
        data = resp_data,
        control = lme4::lmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = 20000)
        )
      ),
      warning = function(w) {
        cat(sprintf("    WARNING (random slopes, %s): %s\n", resp_label, conditionMessage(w)))
        invokeRestart("muffleWarning")
      }
    )
  }, error = function(e) {
    cat(sprintf("    ERROR (random slopes, %s): %s\n", resp_label, conditionMessage(e)))
    NULL
  })

  # Check convergence; fall back to random-intercepts if needed
  resp_fallback <- FALSE
  if (is.null(resp_lmer)) {
    resp_fallback <- TRUE
  } else if (length(resp_lmer@optinfo$conv$lme4) > 0) {
    cat(sprintf("    Convergence issues in %s; falling back to random-intercepts\n", resp_label))
    resp_fallback <- TRUE
  }

  if (resp_fallback) {
    cat(sprintf("    Fitting random-intercepts fallback (%s)...\n", resp_label))
    resp_lmer <- tryCatch(
      lme4::lmer(
        lnDiff ~ Species * time +
          (1 | CA_MPA_Name_Short) + (1 | source),
        data = resp_data
      ),
      error = function(e) {
        cat(sprintf("    ERROR (random intercepts, %s): %s\n", resp_label, conditionMessage(e)))
        NULL
      }
    )
  }

  if (!is.null(resp_lmer) && lme4::isSingular(resp_lmer)) {
    cat(sprintf("    NOTE: %s model has singular fit.\n", resp_label))
  }

  if (is.null(resp_lmer)) {
    cat(sprintf("    WARNING: %s meta-regression model failed. Skipping.\n", resp_label))
    return(NULL)
  }

  # Extract fixed effects (same pattern as pooled model)
  if (requireNamespace("lmerTest", quietly = TRUE)) {
    resp_summary <- summary(lmerTest::as_lmerModLmerTest(resp_lmer))
  } else {
    resp_summary <- summary(resp_lmer)
  }

  resp_fixed <- as.data.frame(coef(resp_summary))
  resp_fixed$Term <- rownames(resp_fixed)
  rownames(resp_fixed) <- NULL

  if ("Pr(>|t|)" %in% names(resp_fixed)) {
    names(resp_fixed)[names(resp_fixed) == "Pr(>|t|)"] <- "p_value"
  }
  names(resp_fixed)[names(resp_fixed) == "Estimate"] <- "Estimate"
  names(resp_fixed)[names(resp_fixed) == "Std. Error"] <- "SE"
  names(resp_fixed)[names(resp_fixed) == "t value"] <- "t_value"

  col_order <- c("Term", "Estimate", "SE", "t_value")
  if ("df" %in% names(resp_fixed)) col_order <- c(col_order, "df")
  if ("p_value" %in% names(resp_fixed)) col_order <- c(col_order, "p_value")
  remaining <- setdiff(names(resp_fixed), col_order)
  resp_fixed <- resp_fixed[, c(col_order, remaining)]

  # Save CSV
  csv_path <- here::here("tables", paste0("table_s_temporal_meta_regression_",
                                          tolower(resp_type), ".csv"))
  write.csv(resp_fixed, csv_path, row.names = FALSE)
  cat(sprintf("    Saved: %s\n", csv_path))

  # Print species-level slopes
  cat(sprintf("\n    %s Meta-Regression Slopes:\n", resp_label))
  if (resp_fallback) {
    cat("      (Random-intercepts model)\n")
  } else {
    cat("      (Random-slopes model converged)\n")
  }

  ref_slope_r <- resp_fixed$Estimate[resp_fixed$Term == "time"]
  ref_se_r    <- resp_fixed$SE[resp_fixed$Term == "time"]
  vcov_r      <- as.matrix(vcov(resp_lmer))
  term_names_r <- rownames(vcov_r)

  for (sp in sp_present) {
    sp_abbrev <- full_to_abbrev[sp]
    if (sp == sp_present[1]) {
      cat(sprintf("      %s: slope = %.4f +/- %.4f\n", sp_abbrev, ref_slope_r, ref_se_r))
    } else {
      sp_clean <- gsub(" ", ".", sp)
      int_term <- grep(paste0("Species", sp_clean, ":time|time:Species", sp_clean),
                       term_names_r, value = TRUE)
      if (length(int_term) > 0) {
        int_term <- int_term[1]
        sp_slope_r <- ref_slope_r + resp_fixed$Estimate[resp_fixed$Term == int_term]
        var_sum_r <- vcov_r["time", "time"] + vcov_r[int_term, int_term] +
                     2 * vcov_r["time", int_term]
        sp_se_r <- sqrt(max(var_sum_r, 0))
        # Check significance
        p_int <- resp_fixed$p_value[resp_fixed$Term == int_term]
        sig_marker <- if (!is.null(p_int) && length(p_int) > 0 && !is.na(p_int) && p_int < 0.05) " *" else ""
        cat(sprintf("      %s: slope = %.4f +/- %.4f%s\n", sp_abbrev, sp_slope_r, sp_se_r, sig_marker))
      }
    }
  }

  # Extract diagnostics using shared helper
  resp_diag <- extract_lmer_diagnostics(
    resp_lmer, sp_present, full_to_abbrev,
    model_label = tolower(resp_type), is_fallback = resp_fallback
  )

  return(list(model = resp_lmer, diagnostics = resp_diag))
}

# Fit biomass-only model (all 5 species have biomass data)
bio_result <- fit_resp_lmer("Bio", "Biomass", species_order_full)
temporal_lmer_bio <- if (!is.null(bio_result)) bio_result$model else NULL

# Fit density-only model (4 species — M. pyrifera has no density data)
# species_order_full without M. pyrifera
density_species <- species_order_full[species_order_full != "Macrocystis pyrifera"]
den_result <- fit_resp_lmer("Den", "Density", density_species)
temporal_lmer_den <- if (!is.null(den_result)) den_result$model else NULL

# --- Accumulate response-specific diagnostics ---
if (!is.null(bio_result)) {
  all_fit_stats[[length(all_fit_stats) + 1]]      <- bio_result$diagnostics$fit_stats
  all_species_slopes[[length(all_species_slopes) + 1]] <- bio_result$diagnostics$species_slopes
}
if (!is.null(den_result)) {
  all_fit_stats[[length(all_fit_stats) + 1]]      <- den_result$diagnostics$fit_stats
  all_species_slopes[[length(all_species_slopes) + 1]] <- den_result$diagnostics$species_slopes
}

# --- Export vcov-derived prediction parameters for Figure 4/S10 CI ribbons ---
# For each response-specific lmer model, extract per-species intercept, slope,
# their SEs, and the covariance between intercept and slope. These are needed
# to compute prediction intervals: Var(y_hat) = Var(int) + t^2*Var(slope) + 2*t*Cov(int,slope).

extract_prediction_params <- function(lmer_obj, species_levels, resp_label) {
  if (is.null(lmer_obj)) return(NULL)

  # Get fixed effects with Satterthwaite df
  if (requireNamespace("lmerTest", quietly = TRUE)) {
    sm <- summary(lmerTest::as_lmerModLmerTest(lmer_obj))
  } else {
    sm <- summary(lmer_obj)
  }
  fe <- as.data.frame(coef(sm))
  fe$Term <- rownames(fe)
  if ("Pr(>|t|)" %in% names(fe)) names(fe)[names(fe) == "Pr(>|t|)"] <- "p_value"
  names(fe)[names(fe) == "Std. Error"] <- "SE"

  vcov_mat <- as.matrix(vcov(lmer_obj))
  term_names <- rownames(vcov_mat)

  # Reference species values (first level in species_levels)
  ref_int   <- fe$Estimate[fe$Term == "(Intercept)"]
  ref_slope <- fe$Estimate[fe$Term == "time"]

  pred_rows <- list()
  for (sp in species_levels) {
    if (sp == species_levels[1]) {
      # Reference species: parameters come directly from model
      sp_int    <- ref_int
      sp_slope  <- ref_slope
      sp_int_se <- sqrt(vcov_mat["(Intercept)", "(Intercept)"])
      sp_slope_se <- sqrt(vcov_mat["time", "time"])
      sp_cov    <- vcov_mat["(Intercept)", "time"]
      # Satterthwaite df for the slope term
      sp_df_slope <- if ("df" %in% names(fe)) fe$df[fe$Term == "time"] else NA_real_
    } else {
      # Non-reference species: add species main effect and interaction
      sp_clean <- gsub(" ", ".", sp)
      sp_main_term <- grep(paste0("^Species", sp_clean, "$"), term_names, value = TRUE)
      sp_int_term  <- grep(paste0("Species", sp_clean, ":time|time:Species", sp_clean),
                           term_names, value = TRUE)

      if (length(sp_main_term) == 0 || length(sp_int_term) == 0) next
      sp_main_term <- sp_main_term[1]
      sp_int_term  <- sp_int_term[1]

      # Combined intercept and slope
      sp_int   <- ref_int   + fe$Estimate[fe$Term == sp_main_term]
      sp_slope <- ref_slope + fe$Estimate[fe$Term == sp_int_term]

      # Var(intercept) = Var(Intercept) + Var(SpeciesX) + 2*Cov(Intercept, SpeciesX)
      sp_int_var <- vcov_mat["(Intercept)", "(Intercept)"] +
                    vcov_mat[sp_main_term, sp_main_term] +
                    2 * vcov_mat["(Intercept)", sp_main_term]
      sp_int_se <- sqrt(max(sp_int_var, 0))

      # Var(slope) = Var(time) + Var(SpeciesX:time) + 2*Cov(time, SpeciesX:time)
      sp_slope_var <- vcov_mat["time", "time"] +
                      vcov_mat[sp_int_term, sp_int_term] +
                      2 * vcov_mat["time", sp_int_term]
      sp_slope_se <- sqrt(max(sp_slope_var, 0))

      # Cov(combined intercept, combined slope) =
      #   Cov(Intercept, time) + Cov(Intercept, SpeciesX:time) +
      #   Cov(SpeciesX, time) + Cov(SpeciesX, SpeciesX:time)
      sp_cov <- vcov_mat["(Intercept)", "time"] +
                vcov_mat["(Intercept)", sp_int_term] +
                vcov_mat[sp_main_term, "time"] +
                vcov_mat[sp_main_term, sp_int_term]

      # Satterthwaite df for the interaction slope term
      sp_df_slope <- if ("df" %in% names(fe)) fe$df[fe$Term == sp_int_term] else NA_real_
    }

    pred_rows[[sp]] <- data.frame(
      resp = resp_label,
      Species = sp,
      intercept = sp_int,
      intercept_se = sp_int_se,
      slope = sp_slope,
      slope_se = sp_slope_se,
      cov_int_slope = sp_cov,
      df_slope = sp_df_slope,
      stringsAsFactors = FALSE
    )
  }

  dplyr::bind_rows(pred_rows)
}

pred_params_list <- list()
if (!is.null(bio_result)) {
  bio_sp <- intersect(species_order_full, unique(as.character(bio_result$model@frame$Species)))
  pred_params_list[["bio"]] <- extract_prediction_params(
    bio_result$model, bio_sp, "Bio"
  )
}
if (!is.null(den_result)) {
  den_sp <- intersect(density_species, unique(as.character(den_result$model@frame$Species)))
  pred_params_list[["den"]] <- extract_prediction_params(
    den_result$model, den_sp, "Den"
  )
}

if (length(pred_params_list) > 0) {
  pred_params_df <- dplyr::bind_rows(pred_params_list)
  write.csv(pred_params_df,
            here::here("tables", "table_s_lmer_prediction_params.csv"),
            row.names = FALSE)
  cat("  Saved: tables/table_s_lmer_prediction_params.csv\n")
}

# --- Export combined model diagnostics across pooled/bio/den ---
if (length(all_fit_stats) > 0) {
  combined_fit <- dplyr::bind_rows(all_fit_stats)
  write.csv(combined_fit,
            here::here("tables", "table_s_lmer_model_fit.csv"),
            row.names = FALSE)
  cat("  Saved: tables/table_s_lmer_model_fit.csv\n")
}
if (length(all_species_slopes) > 0) {
  combined_slopes <- dplyr::bind_rows(all_species_slopes)
  write.csv(combined_slopes,
            here::here("tables", "table_s_species_slopes.csv"),
            row.names = FALSE)
  cat("  Saved: tables/table_s_species_slopes.csv\n")
}

cat("\n=== Response-specific meta-regression complete ===\n")


# ---------------------------------------------------------------------------
# Export lmer residual diagnostics for downstream plotting (fig_s13)
# ---------------------------------------------------------------------------
# Save residuals, fitted values, and random effects from all three lmer models
# (pooled, bio, den) so that 11_figures.R can render diagnostics without
# re-fitting models.

extract_resid_diag <- function(lmer_obj, model_label) {
  if (is.null(lmer_obj)) return(NULL)
  md <- lmer_obj@frame
  required_cols <- c("Species", "time", "CA_MPA_Name_Short")
  missing_cols <- setdiff(required_cols, names(md))
  if (length(missing_cols) > 0) {
    warning("lmer model frame missing columns: ", paste(missing_cols, collapse = ", "),
            " — skipping ", model_label)
    return(NULL)
  }
  n <- nrow(md)
  resid_raw <- residuals(lmer_obj)
  # source and resp may not be in model frame for response-specific models
  src_col <- if ("source" %in% names(md)) as.character(md$source) else rep(NA_character_, n)
  resp_col <- if ("resp" %in% names(md)) as.character(md$resp) else rep(model_label, n)
  data.frame(
    model         = model_label,
    residual      = resid_raw,
    fitted        = fitted(lmer_obj),
    std_residual  = resid_raw / sigma(lmer_obj),
    species       = as.character(md$Species),
    time          = md$time,
    MPA           = as.character(md$CA_MPA_Name_Short),
    source        = src_col,
    resp          = resp_col,
    stringsAsFactors = FALSE
  )
}

resid_diag_list <- list(
  extract_resid_diag(temporal_lmer, "pooled"),
  extract_resid_diag(temporal_lmer_bio, "bio"),
  extract_resid_diag(temporal_lmer_den, "den")
)
resid_diag_all <- dplyr::bind_rows(resid_diag_list)
if (nrow(resid_diag_all) == 0) {
  warning("All lmer models are NULL — lmer_residual_diagnostics.csv will be empty. ",
          "Check model fitting above.")
}
write.csv(resid_diag_all,
          here::here("outputs", "lmer_residual_diagnostics.csv"),
          row.names = FALSE)
cat("  Saved: outputs/lmer_residual_diagnostics.csv\n")

# Random effects (BLUPs) for all lmer models (pooled, bio, den)
extract_ranef_diag <- function(lmer_obj, model_label) {
  if (is.null(lmer_obj)) return(NULL)
  re_list <- lme4::ranef(lmer_obj)
  if (!"CA_MPA_Name_Short" %in% names(re_list)) {
    warning("Random effect 'CA_MPA_Name_Short' not found in ", model_label,
            " model. Available groups: ", paste(names(re_list), collapse = ", "))
    return(NULL)
  }
  re <- re_list$CA_MPA_Name_Short
  re$MPA <- rownames(re)
  re$model <- model_label
  rownames(re) <- NULL
  re
}

ranef_diag_list <- list(
  extract_ranef_diag(temporal_lmer, "pooled"),
  extract_ranef_diag(temporal_lmer_bio, "bio"),
  extract_ranef_diag(temporal_lmer_den, "den")
)
ranef_diag_all <- dplyr::bind_rows(ranef_diag_list)
if (nrow(ranef_diag_all) == 0) {
  warning("All lmer models are NULL — lmer_ranef_diagnostics.csv will be empty.")
}
write.csv(ranef_diag_all,
          here::here("outputs", "lmer_ranef_diagnostics.csv"),
          row.names = FALSE)
cat("  Saved: outputs/lmer_ranef_diagnostics.csv\n")


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
  # NOTE: This DW test is computed on pooled cross-MPA residuals ordered by
  # MPA then time. The diff() operation at MPA boundaries is meaningless,
  # biasing DW toward 2 (anti-conservative for detecting autocorrelation).
  # The AR1 model comparison in Step 3 below is the more reliable diagnostic.
  cat("  Durbin-Watson test on lmer residuals:\n")
  dw_result <- NA_real_
  dw_pvalue <- NA_real_
  tryCatch({
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
    dw_result <<- dw_stat
    cat(sprintf("    DW statistic (ordered by MPA, time) = %.4f\n", dw_stat))
    cat(sprintf("    Interpretation: %s\n",
                ifelse(dw_stat < 1.5, "Substantial positive autocorrelation (DW << 2)",
                       ifelse(dw_stat < 1.8, "Moderate positive autocorrelation (DW < 2)",
                              ifelse(dw_stat > 2.5, "Negative autocorrelation (DW > 2)",
                                     "Weak or no autocorrelation (DW ~ 2)")))))
    # Approximate p-value via permutation test (H0: no autocorrelation)
    # NOTE: DW permutation shuffles residuals globally rather than within-MPA.
    # This is anti-conservative (inflates Type I error) because it breaks
    # within-MPA temporal correlation structure. Interpret permutation p-values
    # as approximate. The parametric DW test results are the primary diagnostic.
    set.seed(42L)  # Reproducibility for DW permutation test
    n_perm <- 999
    dw_perm <- replicate(n_perm, {
      perm_resid <- sample(resid_df$resid)
      sum(diff(perm_resid)^2) / sum(perm_resid^2)
    })
    # One-sided test: positive autocorrelation means DW < 2
    dw_pvalue <<- (sum(dw_perm <= dw_stat) + 1) / (n_perm + 1)
    cat(sprintf("    DW p-value (permutation, H0: no autocorrelation) = %.4f\n", dw_pvalue))
  }, error = function(e) {
    cat("    DW test failed:", conditionMessage(e), "\n")
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

  # NOTE: Attempt 1's random list is interpreted by nlme as nested (source > MPA > resp),
  # not the intended crossed structure. This attempt is expected to fail for most datasets
  # and serves as a "try the richest model first" step before the simpler Attempt 2.
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
  # Check if AR1 model actually converged (returnObject=TRUE may return non-converged fit)
  if (!is.null(temporal_lme_ar1) && inherits(temporal_lme_ar1$apVar, "character")) {
    warning("AR1 model (attempt 1) returned non-positive-definite variance; treating as non-converged")
    temporal_lme_ar1 <- NULL
  }
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
    # Check if AR1 model actually converged (returnObject=TRUE may return non-converged fit)
    if (!is.null(temporal_lme_ar1) && inherits(temporal_lme_ar1$apVar, "character")) {
      warning("AR1 model (attempt 2) returned non-positive-definite variance; treating as non-converged")
      temporal_lme_ar1 <- NULL
    }
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
    # Check if AR1 model actually converged (returnObject=TRUE may return non-converged fit)
    if (!is.null(temporal_lme_ar1) && inherits(temporal_lme_ar1$apVar, "character")) {
      warning("AR1 model (attempt 3) returned non-positive-definite variance; treating as non-converged")
      temporal_lme_ar1 <- NULL
    }
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
    ar1_comparison_out$DW_pvalue <- ifelse(!is.na(dw_pvalue), dw_pvalue, NA_real_)
    ar1_comparison_out$AR1_model <- ar1_model_description
    write.csv(ar1_comparison_out,
              here::here("tables", "table_s_ar1_sensitivity.csv"),
              row.names = FALSE)
    cat("  Saved: tables/table_s_ar1_sensitivity.csv\n")

    # --- Save AR1 fixed effects table (parallel to lmer table) ---
    ar1_fixed_out <- ar1_fixed[, intersect(c("Term", "Estimate", "SE",
                                              "t_value", "DF", "p_value"),
                                            names(ar1_fixed))]
    write.csv(ar1_fixed_out,
              here::here("tables", "table_s_temporal_meta_regression_ar1.csv"),
              row.names = FALSE)
    cat("  Saved: tables/table_s_temporal_meta_regression_ar1.csv\n")

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
      p_value = ifelse(!is.na(dw_pvalue), dw_pvalue, NA_real_),
      note = ifelse(!is.na(dw_result) && dw_result < 1.5,
                    "Positive autocorrelation detected; AR1 model did not converge",
                    "Weak autocorrelation; lmer results likely robust"),
      stringsAsFactors = FALSE
    )
    write.csv(dw_out,
              here::here("tables", "table_s_ar1_sensitivity.csv"),
              row.names = FALSE)
    cat("  Saved: tables/table_s_ar1_sensitivity.csv (DW diagnostic only)\n")
  }

} else {
  cat("  Skipping AR1 analysis — lmer model did not fit.\n")
}

cat("=== AR1 sensitivity analysis complete ===\n")

# RECOMMENDED SENSITIVITY: Extract effect sizes at t = 8, 10, 14 in addition
# to the default t = 11 (EFFECT_SIZE_TIME_YEARS in 00c_analysis_constants.R).
# This would test whether meta-analytic conclusions are robust to the choice
# of extraction time. The current pipeline extracts at t = 11 because it
# corresponds to the age of the youngest MPA in the dataset (MLPA South Coast,
# est. 2012). Extracting at alternative time points would require re-running
# 08_effect_sizes.R with modified EFFECT_SIZE_TIME_YEARS values, then
# re-running 09_meta_analysis.R for each. Not implemented in the current
# pipeline but should be documented as a limitation in the manuscript discussion.


# =============================================================================
# Section D: Analysis 3 -- Trophic Cascade Phase Portrait (Species-Level)
# =============================================================================
# CAVEAT: The set of contributing MPAs changes over time (unbalanced panel).
# Late time points (t > 10) are dominated by the oldest reserves (Channel Islands).
# Phase trajectories at late time points may reflect WHICH MPAs contribute
# rather than within-MPA temporal dynamics. The n >= 2 filter is minimal
# protection. Correlations should be interpreted as exploratory, not confirmatory.
# See figure caption for documentation of this limitation.
#
# Figure S5: Four-panel species-level phase portrait showing how individual
# species effects co-evolve over time. Each panel pairs one species against
# M. pyrifera (kelp) or against an urchin species, revealing whether
# cascade sequencing operates at the species level.
#
# DATA NOTE: Each point is the network-wide average lnRR for a species at a
# given year-since-MPA, averaged across all MPAs with data at that time step
# (two-stage: MPA means first, then cross-MPA average, ≥2 MPAs required).
# Because MPAs were established in different years, the set of contributing
# MPAs changes over time — early time points include most MPAs, late time
# points are dominated by the oldest reserves (Channel Islands, est. ~2003).
# This unbalanced panel means the trajectory may partly reflect compositional
# shifts in which MPAs contribute rather than pure within-MPA dynamics.
# The formal lmer models (Fig 3, S4) handle this via MPA random effects;
# this figure is a complementary descriptive visualization.

if (should_render("fig_s05")) {
cat("\n--- Figure S5: Cascade Phase Portrait (Species-Level) ---\n")

library(patchwork)

# NOTE: Bio and Den response types are pooled here (averaged within MPA-year).
# This means each data point may mix biomass and density lnRR values.
# The formal lmer fits separate Bio/Den models; this descriptive
# aggregation does not account for non-independence between response types.

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

# --- Export phase portrait data ---
write.csv(phase_yearly,
          here::here("tables", "table_s_phase_portrait_yearly_means.csv"),
          row.names = FALSE)
cat("  Saved: tables/table_s_phase_portrait_yearly_means.csv\n")

# Helper: build one phase portrait panel
# q_ arguments label each quadrant: TR = top-right, BL = bottom-left,
# TL = top-left, BR = bottom-right
build_phase_panel <- function(data, x_col, y_col, x_lab, y_lab,
                              q_tr, q_bl, q_tl, q_br) {
  d <- data %>%
    dplyr::filter(!is.na(.data[[x_col]]), !is.na(.data[[y_col]]))
  if (nrow(d) < 3) return(NULL)

  # Compute axis limits that guarantee all four quadrants are visible.
  # Extend at least 25% of the dominant side's extent past zero.
  x_rng <- range(d[[x_col]])
  y_rng <- range(d[[y_col]])
  x_pad <- diff(x_rng) * 0.08
  y_pad <- diff(y_rng) * 0.08
  zero_frac <- 0.35
  x_lo <- min(x_rng[1] - x_pad, -zero_frac * max(abs(x_rng)))
  x_hi <- max(x_rng[2] + x_pad,  zero_frac * max(abs(x_rng)))
  y_lo <- min(y_rng[1] - y_pad, -zero_frac * max(abs(y_rng)))
  y_hi <- max(y_rng[2] + y_pad,  zero_frac * max(abs(y_rng)))

  # Label start and end of trajectory
  d_start <- d %>% dplyr::slice_min(time, n = 1)
  d_end   <- d %>% dplyr::slice_max(time, n = 1)

  ggplot(d, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    # Reference lines
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey75",
               linewidth = 0.35) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey75",
               linewidth = 0.35) +
    # Quadrant labels — all four corners (offset from tags & edges)
    annotate("text", x = Inf, y = Inf, label = q_tr,
             hjust = 1.05, vjust = 1.5, size = 2.2, color = "grey50",
             fontface = "italic") +
    annotate("text", x = -Inf, y = -Inf, label = q_bl,
             hjust = -0.05, vjust = -0.5, size = 2.2, color = "grey50",
             fontface = "italic") +
    annotate("text", x = -Inf, y = Inf, label = q_tl,
             hjust = -0.05, vjust = 3.0, size = 2.2, color = "grey50",
             fontface = "italic") +
    annotate("text", x = Inf, y = -Inf, label = q_br,
             hjust = 1.05, vjust = -0.5, size = 2.2, color = "grey50",
             fontface = "italic") +
    # Temporal path as individual segments with small directional arrows
    geom_segment(
      data = d %>%
        dplyr::mutate(
          xend = dplyr::lead(.data[[x_col]]),
          yend = dplyr::lead(.data[[y_col]])
        ) %>%
        dplyr::filter(!is.na(xend)),
      aes(xend = xend, yend = yend),
      color = "grey50", linewidth = 0.4,
      arrow = arrow(length = unit(1.2, "mm"), type = "closed"),
      lineend = "round"
    ) +
    # Points colored by time
    geom_point(aes(color = time), size = 3, shape = 16) +
    # Year labels at start and end
    ggrepel::geom_text_repel(
      data = d_start,
      aes(label = paste0("t=", time)),
      size = 2.3, color = "grey30", fontface = "bold",
      nudge_x = -0.15, nudge_y = 0.15,
      segment.size = 0.3, segment.color = "grey60",
      min.segment.length = 0, seed = 42, box.padding = 0.3
    ) +
    ggrepel::geom_text_repel(
      data = d_end,
      aes(label = paste0("t=", time)),
      size = 2.3, color = "grey30", fontface = "bold",
      nudge_x = 0.15, nudge_y = -0.15,
      segment.size = 0.3, segment.color = "grey60",
      min.segment.length = 0, seed = 42, box.padding = 0.3
    ) +
    # Sequential gradient using heatmap palette (warm early → cool late)
    scale_color_gradient(
      name = "Years since\nMPA",
      low  = col_heatmap["low"],   # brick red (early)
      high = col_heatmap["high"]   # steel blue (late)
    ) +
    labs(x = x_lab, y = y_lab) +
    coord_cartesian(xlim = c(x_lo, x_hi), ylim = c(y_lo, y_hi),
                    clip = "off") +
    theme_mpa(base_size = 9) +
    theme(
      legend.position = "bottom",
      legend.key.width = unit(10, "mm"),
      legend.key.height = unit(2.5, "mm"),
      panel.grid.major = element_blank(),
      plot.margin = margin(10, 8, 4, 4)
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
  q_tr = "Purples \u2193  Kelp \u2191",
  q_bl = "Purples \u2191  Kelp \u2193",
  q_tl = "Purples \u2191  Kelp \u2191",
  q_br = "Purples \u2193  Kelp \u2193"
)

# Panel (b): M. franciscanus decline vs M. pyrifera increase
pb <- build_phase_panel(
  phase_wide_inv, "M. franciscanus_inv", "M. pyrifera",
  x_lab = expression(italic("M. franciscanus") ~ "lnRR \u00d7 (-1)"),
  y_lab = expression(italic("M. pyrifera") ~ "lnRR"),
  q_tr = "Reds \u2193  Kelp \u2191",
  q_bl = "Reds \u2191  Kelp \u2193",
  q_tl = "Reds \u2191  Kelp \u2191",
  q_br = "Reds \u2193  Kelp \u2193"
)

# Panel (c): P. interruptus increase vs S. purpuratus decline
pc <- build_phase_panel(
  phase_wide_inv, "P. interruptus", "S. purpuratus_inv",
  x_lab = expression(italic("P. interruptus") ~ "lnRR"),
  y_lab = expression(italic("S. purpuratus") ~ "lnRR \u00d7 (-1)"),
  q_tr = "Lobster \u2191  Purples \u2193",
  q_bl = "Lobster \u2193  Purples \u2191",
  q_tl = "Lobster \u2193  Purples \u2193",
  q_br = "Lobster \u2191  Purples \u2191"
)

# Panel (d): S. pulcher increase vs M. franciscanus decline
pd <- build_phase_panel(
  phase_wide_inv, "S. pulcher", "M. franciscanus_inv",
  x_lab = expression(italic("S. pulcher") ~ "lnRR"),
  y_lab = expression(italic("M. franciscanus") ~ "lnRR \u00d7 (-1)"),
  q_tr = "Sheephead \u2191  Reds \u2193",
  q_bl = "Sheephead \u2193  Reds \u2191",
  q_tl = "Sheephead \u2193  Reds \u2193",
  q_br = "Sheephead \u2191  Reds \u2191"
)

# Assemble available panels
panels <- list(pa, pb, pc, pd)
panels <- panels[!vapply(panels, is.null, logical(1))]

if (length(panels) >= 2) {
  fig_s05 <- wrap_plots(panels, ncol = 2) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
    theme(
      plot.tag = element_text(size = 9, face = "bold"),
      legend.position = "bottom",
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.line.x.bottom = element_line(colour = "black", linewidth = 0.5),
      axis.line.y.left = element_line(colour = "black", linewidth = 0.5)
    )

  save_fig(fig_s05, "fig_s05_cascade_phase",
           w = FIG_S05_DIMS["w"], h = FIG_S05_DIMS["h"])
} else {
  cat("  WARNING: Too few species-level pairings for phase portrait.\n")
}

# --- Phase portrait pairwise correlation statistics ---
# Quantify the co-movement shown in the phase portrait figure (S4).
# Each row tests whether two species' yearly mean lnRR values are correlated
# over time (Pearson and Spearman), which provides statistical support for
# the visual trajectories. Urchin columns are inverted (×-1) to match the
# cascade expectation (decline = positive signal).
cat("  Computing phase portrait pairwise correlations...\n")
phase_pairs <- list(
  list(x = "S. purpuratus_inv",  y = "M. pyrifera",         x_lab = "S. purpuratus (inv)", y_lab = "M. pyrifera"),
  list(x = "M. franciscanus_inv", y = "M. pyrifera",        x_lab = "M. franciscanus (inv)", y_lab = "M. pyrifera"),
  list(x = "P. interruptus",     y = "S. purpuratus_inv",   x_lab = "P. interruptus", y_lab = "S. purpuratus (inv)"),
  list(x = "S. pulcher",         y = "M. franciscanus_inv",  x_lab = "S. pulcher", y_lab = "M. franciscanus (inv)")
)

phase_cor_results <- list()
for (pp in phase_pairs) {
  d <- phase_wide_inv %>%
    dplyr::filter(!is.na(.data[[pp$x]]), !is.na(.data[[pp$y]]))
  if (nrow(d) >= 4) {
    ct_pearson  <- cor.test(d[[pp$x]], d[[pp$y]], method = "pearson")
    ct_spearman <- cor.test(d[[pp$x]], d[[pp$y]], method = "spearman", exact = FALSE)
    phase_cor_results[[length(phase_cor_results) + 1]] <- data.frame(
      x_species = pp$x_lab, y_species = pp$y_lab,
      n_years = nrow(d),
      pearson_r = ct_pearson$estimate,
      pearson_p = ct_pearson$p.value,
      pearson_ci_lb = ct_pearson$conf.int[1],
      pearson_ci_ub = ct_pearson$conf.int[2],
      spearman_rho = ct_spearman$estimate,
      spearman_p = ct_spearman$p.value,
      stringsAsFactors = FALSE
    )
  } else {
    phase_cor_results[[length(phase_cor_results) + 1]] <- data.frame(
      x_species = pp$x_lab, y_species = pp$y_lab,
      n_years = nrow(d),
      pearson_r = NA, pearson_p = NA,
      pearson_ci_lb = NA, pearson_ci_ub = NA,
      spearman_rho = NA, spearman_p = NA,
      stringsAsFactors = FALSE
    )
  }
}
phase_cor_df <- do.call(rbind, phase_cor_results)
write.csv(phase_cor_df,
          here::here("tables", "table_s_phase_portrait_correlations.csv"),
          row.names = FALSE)
cat("  Saved: tables/table_s_phase_portrait_correlations.csv\n")

} # end fig_s05


# =============================================================================
# Analysis 4: Species-Level Space-Time Heatmap (Figure S6)
# =============================================================================
# Shows all five species as stacked heatmap panels (replaces old urchin-only version)
# stacked vertically, preserving distinct recovery patterns for each
# predator and urchin species.

if (FALSE) {  # DROPPED: S5 triptych heatmap (redundant with S3+S6)
cat("\n--- Figure S6: Species-Level Heatmap ---\n")

# (FIG_S06_DIMS defined in header constants section)

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
      low = col_heatmap["low"], mid = col_heatmap["mid"], high = col_heatmap["high"],
      midpoint = 0, limits = color_lim, oob = scales::squish,
      na.value = "grey85", name = "lnRR"
    ) +
    ggplot2::scale_x_continuous(breaks = seq(0, 15, by = 3),
                                expand = c(0, 0)) +
    ggplot2::labs(title = sp_abbrev, y = NULL) +
    theme_mpa(base_size = 9) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 8),
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
    tag_levels = "a", tag_prefix = "(", tag_suffix = ")",
    caption = "Grey cells = no data available"
  ) &
  ggplot2::theme(
    plot.tag     = ggplot2::element_text(size = 9, face = "bold"),
    plot.caption = ggplot2::element_text(size = 7, color = "grey40", hjust = 0)
  )

save_fig(fig_s05, "fig_s06_triptych_heatmap", FIG_S06_DIMS["w"], FIG_S06_DIMS["h"])
cat("  Figure S6 saved.\n")

} # end fig_s06 (heatmap)


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
# Slopes computed separately for Bio/Den, then averaged per MPA-Species
# to avoid pooling different response types in a single regression.
slopes_by_resp <- temporal_after %>%
  dplyr::filter(time <= 15) %>%
  dplyr::group_by(CA_MPA_Name_Short, Species, Species_abbrev, resp) %>%
  dplyr::filter(dplyr::n_distinct(time) >= 5) %>%
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

# Average Bio/Den slopes per MPA-Species for cascade consistency analysis
slope_data <- slopes_by_resp %>%
  dplyr::group_by(CA_MPA_Name_Short, Species, Species_abbrev) %>%
  dplyr::summarise(
    slope = mean(slope, na.rm = TRUE),
    slope_se = mean(slope_se, na.rm = TRUE),
    n_years = max(n_years, na.rm = TRUE),
    n_resp = dplyr::n(),
    .groups = "drop"
  )

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

# --- Export per-MPA slopes to CSV ---
write.csv(slope_data,
          here::here("tables", "table_s_per_mpa_slopes.csv"),
          row.names = FALSE)
cat("  Saved: tables/table_s_per_mpa_slopes.csv\n")

# --- Export species-level slope summaries to CSV ---
write.csv(slope_summary,
          here::here("tables", "table_s_slope_summary_by_species.csv"),
          row.names = FALSE)
cat("  Saved: tables/table_s_slope_summary_by_species.csv\n")

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
    width = 0.15, size = 2.0, alpha = 0.65, shape = 16
  ) +
  ggplot2::geom_pointrange(
    data = slope_summary,
    ggplot2::aes(x = Species_abbrev, y = mean_slope,
                 ymin = mean_slope - se_slope,
                 ymax = mean_slope + se_slope),
    size = 0.6, linewidth = 0.8, color = "black",
    position = ggplot2::position_nudge(x = 0.18)
  ) +
  ggplot2::scale_color_manual(values = sp_color_abbrev, guide = "none") +
  ggplot2::coord_flip() +
  ggplot2::labs(
    x = NULL,
    y = expression("Annual rate of change (slope of lnRR" %~% "time)")
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
#
# LIMITATION: Cascade consistency metric
# The cascade consistency score is based on sign agreement of slopes
# (predators up, urchins down, kelp up). This is a descriptive metric
# without a formal null expectation or statistical test. Under the null
# hypothesis of no MPA effect, random sign assignment would yield ~12.5%
# (1/8) chance of all three signs matching the cascade prediction.
# A formal binomial test could assess whether observed consistency exceeds
# this null rate.
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

fill_cascade <- col_cascade  # From centralized palette (00b_color_palette.R)

# Create complete grid so missing combos get a "no data" fill
all_combos <- tidyr::expand_grid(
  MPA_short = levels(cascade_long$MPA_short),
  Species_dir = levels(cascade_long$Species_dir)
)
cascade_complete <- dplyr::left_join(all_combos, cascade_long,
                                      by = c("MPA_short", "Species_dir"))
cascade_complete$fill_val <- dplyr::case_when(
  is.na(cascade_complete$consistent) ~ "NA",
  cascade_complete$consistent         ~ "TRUE",
  TRUE                                ~ "FALSE"
)

p_cascade <- ggplot2::ggplot(
  cascade_complete,
  ggplot2::aes(x = Species_dir, y = MPA_short, fill = fill_val)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.5) +
  ggplot2::scale_fill_manual(
    values = c(fill_cascade, "NA" = "grey90"),
    labels = c("FALSE" = "No", "TRUE" = "Yes", "NA" = "No data"),
    name   = "Expected direction"
  ) +
  ggplot2::labs(x = NULL, y = NULL) +
  theme_mpa(base_size = 9) +
  ggplot2::theme(
    axis.text.y     = ggplot2::element_text(size = 8),
    axis.text.x     = ggplot2::element_text(size = 8.5, face = "italic",
                                            angle = 45, hjust = 1),
    legend.position = "bottom",
    panel.grid      = ggplot2::element_blank()
  )

# ---------------------------------------------------------------------------
# 5d. Assemble 2-panel figure
# ---------------------------------------------------------------------------

fig_s05 <- p_slopes + p_cascade +
  patchwork::plot_layout(widths = c(1, 1.3)) +
  patchwork::plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  ggplot2::theme(
    plot.tag = ggplot2::element_text(size = 9, face = "bold")
  )

save_fig(fig_s05, "fig_s06_slope_comparison", FIG_S06_DIMS["w"], FIG_S06_DIMS["h"])
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
          here::here("tables", "table_s_cascade_consistency.csv"),
          row.names = FALSE)
cat("  Table saved: tables/table_s_cascade_consistency.csv\n")

} # end fig_s06


# =============================================================================
# Completion
# =============================================================================
cat("\n=== 10_temporal_analysis.R complete ===\n")
