# =============================================================================
# 13_additional_analyses.R
# =============================================================================
#
# PURPOSE:
#   Additional supplemental analyses and figures examining how MPA
#   characteristics (protection level, location) relate to effect sizes.
#
# ANALYSES:
#   A1. Combined moderator comparison: SMR vs SMCA + mainland vs Ch. Islands (Fig S9)
#   C1. Formal moderator meta-regression table
#
# INPUTS:
#   - SumStats.Final: Effect sizes (ES at t=11) from 08_effect_sizes.R
#   - Site: MPA metadata from 03_load_harmonized_data.R
#   - Color palette objects from 00b_color_palette.R
#
# OUTPUTS:
#   Figures: fig_s09_moderator_comparisons
#   Tables:  table_s_moderator_meta_regression.csv   (coefficient table)
#            table_s_moderator_diagnostics.csv        (QM, heterogeneity, sample sizes)
#            table_s_moderator_subgroup_effects.csv   (per-subgroup effect sizes & n)
#
# NOTE: Previous versions generated additional supplemental figures (MPA size,
#       regional, rate-of-change, cascade completeness) which were removed to
#       eliminate redundancy.
#
# AUTHORS: Emily Donham & Adrian Stier
# =============================================================================


# =============================================================================
# Section A: Setup
# =============================================================================

dir.create(here::here("plots"), showWarnings = FALSE)

cat("\n=== 13_additional_analyses.R ===\n")

# --- Figure dimension constants ---
FIG_S09_DIMS <- c(w = 17.8, h = 18)   # Combined moderator comparisons (2-row)

# --- Input validation ---
required_objects <- c("SumStats.Final", "Site")
missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), envir = globalenv())
]
if (length(missing_objects) > 0) {
  stop("Missing required data objects: ", paste(missing_objects, collapse = ", "),
       "\nPlease run scripts 00-09 first.")
}

# --- Species name mappings ---
full_to_abbrev <- c(
  "Macrocystis pyrifera"          = "M. pyrifera",
  "Mesocentrotus franciscanus"    = "M. franciscanus",
  "Strongylocentrotus purpuratus" = "S. purpuratus",
  "Panulirus interruptus"         = "P. interruptus",
  "Semicossyphus pulcher"         = "S. pulcher"
)

species_order_full <- c(
  "Panulirus interruptus", "Semicossyphus pulcher",
  "Strongylocentrotus purpuratus", "Mesocentrotus franciscanus",
  "Macrocystis pyrifera"
)
species_order_abbrev <- full_to_abbrev[species_order_full]

# Species colors keyed by abbreviated name
sp_color_abbrev <- setNames(
  unname(col_taxa[species_order_abbrev]),
  species_order_abbrev
)

# --- Merge SumStats.Final with Site metadata ---
# Detect the MPA column name in SumStats.Final
mpa_col <- if ("CA_MPA_Name_Short" %in% names(SumStats.Final)) {
  "CA_MPA_Name_Short"
} else if ("MPA" %in% names(SumStats.Final)) {
  "MPA"
} else {
  stop("Cannot find MPA column in SumStats.Final")
}

# Check which Site columns SumStats.Final already has (from prior joins)
needed_cols <- c("type", "Hectares", "Location", "MPA_Start")
already_has <- intersect(needed_cols, names(SumStats.Final))
still_need  <- setdiff(needed_cols, names(SumStats.Final))

if (length(still_need) > 0) {
  site_meta <- Site %>%
    dplyr::select(CA_MPA_Name_Short, dplyr::any_of(still_need)) %>%
    dplyr::distinct()
  ss_merged <- SumStats.Final %>%
    dplyr::left_join(site_meta, by = setNames("CA_MPA_Name_Short", mpa_col))
} else {
  # All columns already present — just use SumStats.Final directly
  ss_merged <- SumStats.Final
}

# Standardize Taxa to abbreviated names
if (!"Taxa_abbrev" %in% names(ss_merged)) {
  taxa_col_ss <- if ("Taxa" %in% names(ss_merged)) "Taxa"
                 else stop("Required column 'Taxa' not found in ss_merged")
  ss_merged$Taxa_abbrev <- dplyr::case_when(
    ss_merged[[taxa_col_ss]] %in% names(full_to_abbrev) ~
      full_to_abbrev[ss_merged[[taxa_col_ss]]],
    TRUE ~ ss_merged[[taxa_col_ss]]
  )
}

# Detect taxa column
taxa_col_ss <- if ("Taxa" %in% names(ss_merged)) "Taxa" else "Taxa_abbrev"

cat(sprintf("  Merged dataset: %d rows, %d with type info\n",
            nrow(ss_merged),
            sum(!is.na(ss_merged$type))))

cat("=== Setup complete ===\n")


# =============================================================================
# Section A2: Data for Moderator Analysis (consistent with 09_meta_analysis.R)
# =============================================================================
# The primary meta-analysis in script 09 uses the full dataset (no outlier removal).
# Moderator analysis uses the same full data for consistency.
#
# Strategy: if per-taxa full data objects from script 09 exist in the global
# environment, use those directly. Otherwise, use SumStats.Final (which IS the
# full data before any outlier filtering).

if (exists("pertaxa_biomass") && exists("pertaxa_density") &&
    is.list(pertaxa_biomass) && "full_data" %in% names(pertaxa_biomass) &&
    is.list(pertaxa_density) && "full_data" %in% names(pertaxa_density)) {
  # Use the full (no outlier removal) data from script 09 directly
  ss_clean <- dplyr::bind_rows(pertaxa_biomass$full_data,
                                pertaxa_density$full_data)
  cat(sprintf("  Using full dataset from script 09 (no outlier removal): %d rows\n",
              nrow(ss_clean)))
  outlier_source <- "full_data_from_09"
} else {
  # Fallback: use SumStats.Final directly (already the full dataset)
  cat("  Per-taxa full data not found; using SumStats.Final directly.\n")
  ss_clean <- SumStats.Final %>%
    dplyr::filter(!is.na(Mean), !is.na(SE), SE > 0) %>%
    dplyr::mutate(vi = as.numeric(SE)^2, Mean = as.numeric(Mean))
  cat(sprintf("  Full dataset: %d rows\n", nrow(ss_clean)))
  outlier_source <- "sumstats_final_no_filtering"
}

# Re-merge with Site metadata (ss_clean may lack Site columns)
needed_cols_clean <- c("type", "Hectares", "Location", "MPA_Start")
already_has_clean <- intersect(needed_cols_clean, names(ss_clean))
still_need_clean  <- setdiff(needed_cols_clean, names(ss_clean))

if (length(still_need_clean) > 0) {
  site_meta_clean <- Site %>%
    dplyr::select(CA_MPA_Name_Short, dplyr::any_of(still_need_clean)) %>%
    dplyr::distinct()
  # Detect MPA column in clean data
  mpa_col_clean <- if ("CA_MPA_Name_Short" %in% names(ss_clean)) {
    "CA_MPA_Name_Short"
  } else if ("MPA" %in% names(ss_clean)) {
    "MPA"
  } else {
    mpa_col
  }
  ss_clean <- ss_clean %>%
    dplyr::left_join(site_meta_clean, by = setNames("CA_MPA_Name_Short", mpa_col_clean))
}

# Standardize Taxa_abbrev
if (!"Taxa_abbrev" %in% names(ss_clean)) {
  taxa_col_clean <- if ("Taxa" %in% names(ss_clean)) "Taxa" else mpa_col
  ss_clean$Taxa_abbrev <- dplyr::case_when(
    ss_clean[[taxa_col_clean]] %in% names(full_to_abbrev) ~
      full_to_abbrev[ss_clean[[taxa_col_clean]]],
    TRUE ~ ss_clean[[taxa_col_clean]]
  )
}

# Ensure Resp column exists
if (!"Resp" %in% names(ss_clean) && "Response" %in% names(ss_clean)) {
  ss_clean$Resp <- ss_clean$Response
}

# Ensure variance column
if (!"V" %in% names(ss_clean)) {
  ss_clean$V <- as.numeric(ss_clean$SE)^2
}

cat(sprintf("  Outlier-filtered data ready: %d rows, outlier source: %s\n",
            nrow(ss_clean), outlier_source))

# Replace ss_merged with the outlier-filtered version for all downstream analyses
ss_merged <- ss_clean


# =============================================================================
# Section B: Combined Moderator Comparison (Fig S9)
# =============================================================================
# Two-row figure: (top) SMR vs SMCA protection level, (bottom) Ch. Islands vs Mainland

if (should_render("fig_s09")) {
cat("\n--- Figure S9: Combined Moderator Comparisons ---\n")
cat("  Using precision-weighted marginal means from rma/rma.mv models\n")

# RR-scaled y-axis uses shared helper (scale_y_rr) from 01_utils.R

# --- Helper: compute model-predicted marginal means per subgroup ---
# For each Taxa x Resp combination within a moderator, fit a small meta-analytic
# model with mods = ~ moderator_level - 1 (no intercept, cell-means). This gives
# the precision-weighted mean for each level with proper SEs from the model
# (accounting for inverse-variance weighting and between-MPA heterogeneity).
#
# Falls back to simple weighted mean if model fitting fails.
compute_marginal_means <- function(data, moderator_col, mpa_col_name) {
  combos <- unique(data[, c("Taxa_abbrev", "Resp"), drop = FALSE])
  results <- list()
  for (i in seq_len(nrow(combos))) {
    taxon <- combos$Taxa_abbrev[i]
    resp  <- combos$Resp[i]
    sub <- data[data$Taxa_abbrev == taxon & data$Resp == resp, ]
    if (nrow(sub) < 2 || length(unique(sub[[moderator_col]])) < 2) {
      # Not enough data for a model — use simple per-level means as fallback
      for (lev in unique(sub[[moderator_col]])) {
        lev_sub <- sub[sub[[moderator_col]] == lev, ]
        results[[length(results) + 1]] <- data.frame(
          Taxa_abbrev = taxon, Resp = resp,
          level = lev,
          mean_es = mean(as.numeric(lev_sub$Mean), na.rm = TRUE),
          se_es = sd(as.numeric(lev_sub$Mean), na.rm = TRUE) / sqrt(nrow(lev_sub)),
          ci_lo = NA_real_, ci_hi = NA_real_,
          n = nrow(lev_sub), method = "arithmetic_mean",
          stringsAsFactors = FALSE
        )
      }
      next
    }
    # Ensure moderator is a factor for the formula
    sub[[moderator_col]] <- factor(sub[[moderator_col]])
    mod_formula <- as.formula(paste0("~ ", moderator_col, " - 1"))
    re_formula <- as.formula(paste0("~ 1 | ", mpa_col_name))
    # Try rma.mv with MPA random effect; fall back to rma if needed
    mod <- tryCatch({
      if (nrow(sub) >= 5 && length(unique(sub[[mpa_col_name]])) >= 3) {
        metafor::rma.mv(yi = Mean, V = V, mods = mod_formula,
                        random = list(re_formula),
                        data = sub, method = "REML", test = "t")
      } else {
        metafor::rma(yi = Mean, vi = V, mods = mod_formula,
                     data = sub, method = "REML", test = "t")
      }
    }, error = function(e) {
      tryCatch(
        metafor::rma(yi = Mean, vi = V, mods = mod_formula,
                     data = sub, method = "REML"),
        error = function(e2) NULL
      )
    })
    if (!is.null(mod)) {
      # With mods = ~ moderator - 1, each coefficient IS the marginal mean
      # for that moderator level. Extract them.
      ctbl <- coef(summary(mod))
      level_names <- gsub(paste0("^", moderator_col), "", rownames(ctbl))
      for (j in seq_len(nrow(ctbl))) {
        results[[length(results) + 1]] <- data.frame(
          Taxa_abbrev = taxon, Resp = resp,
          level = level_names[j],
          mean_es = ctbl[j, "estimate"],
          se_es = ctbl[j, "se"],
          ci_lo = ctbl[j, "ci.lb"],
          ci_hi = ctbl[j, "ci.ub"],
          n = sum(sub[[moderator_col]] == level_names[j]),
          method = "model_marginal_mean",
          stringsAsFactors = FALSE
        )
      }
    } else {
      # Fallback: simple weighted mean per level
      for (lev in levels(sub[[moderator_col]])) {
        lev_sub <- sub[sub[[moderator_col]] == lev, ]
        w <- 1 / as.numeric(lev_sub$V)
        wm <- sum(w * as.numeric(lev_sub$Mean)) / sum(w)
        wse <- sqrt(1 / sum(w))
        results[[length(results) + 1]] <- data.frame(
          Taxa_abbrev = taxon, Resp = resp,
          level = lev,
          mean_es = wm, se_es = wse,
          ci_lo = wm - 1.96 * wse, ci_hi = wm + 1.96 * wse,
          n = nrow(lev_sub), method = "inverse_variance_weighted",
          stringsAsFactors = FALSE
        )
      }
    }
  }
  dplyr::bind_rows(results)
}

# --- Panel A: SMR vs SMCA ---
ss_type <- ss_merged %>%
  dplyr::filter(!is.na(type), type %in% c("SMR", "SMCA"),
                !is.na(Mean), !is.na(SE), SE > 0) %>%
  dplyr::mutate(V = ifelse(is.na(V), as.numeric(SE)^2, V))

panel_a <- NULL
if (nrow(ss_type) > 0) {
  ss_type$Taxa_abbrev <- factor(ss_type$Taxa_abbrev,
                                 levels = species_order_abbrev)
  ss_type <- ss_type %>% dplyr::filter(!is.na(Taxa_abbrev))

  # Compute model-predicted marginal means per Taxa x Resp x type
  type_summary <- compute_marginal_means(ss_type, "type", mpa_col)
  names(type_summary)[names(type_summary) == "level"] <- "type"
  type_summary$Taxa_abbrev <- factor(type_summary$Taxa_abbrev,
                                      levels = species_order_abbrev)

  panel_a <- ggplot(ss_type, aes(x = type, y = Mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
               linewidth = 0.4) +
    geom_jitter(aes(color = Taxa_abbrev), width = 0.15,
                size = 1.8, alpha = 0.55) +
    geom_pointrange(
      data = type_summary,
      aes(x = type, y = mean_es,
          ymin = ci_lo,
          ymax = ci_hi),
      size = 0.5, linewidth = 0.7, color = "black",
      position = position_nudge(x = 0.25)
    ) +
    facet_grid(Resp ~ Taxa_abbrev, scales = "free_y",
               labeller = labeller(
                 Taxa_abbrev = function(x) x,
                 Resp = c("Bio" = "Biomass", "Den" = "Density")
               )) +
    scale_color_manual(values = sp_color_abbrev, guide = "none") +
    scale_y_rr(-5, 5, name = "MPA / Reference at t = 11") +
    labs(
      x = "MPA protection level",
      y = "MPA / Reference at t = 11"
    ) +
    theme_mpa(base_size = 9) +
    theme(
      strip.text.x = element_text(face = "italic", size = 8),
      strip.text.y = element_text(size = 8),
      panel.grid.major = element_blank(),
      panel.spacing.y = unit(0.8, "lines")
    )

  cat("  Panel A (protection): SMR = ", sum(ss_type$type == "SMR"),
      " obs, SMCA = ", sum(ss_type$type == "SMCA"), " obs\n")
  cat("    Summary method breakdown: ",
      paste(unique(type_summary$method), collapse = ", "), "\n")
}

# --- Panel B: Channel Islands vs Mainland ---
ss_region <- ss_merged %>%
  dplyr::filter(!is.na(Location),
                !is.na(Mean), !is.na(SE), SE > 0) %>%
  dplyr::mutate(V = ifelse(is.na(V), as.numeric(SE)^2, V))

panel_b <- NULL
if (nrow(ss_region) > 0 && length(unique(ss_region$Location)) >= 2) {
  ss_region$Region <- dplyr::case_when(
    ss_region$Location == "C"  ~ "Ch. Islands",
    ss_region$Location == "ML" ~ "Mainland",
    TRUE ~ ss_region$Location
  )

  ss_region$Taxa_abbrev <- factor(ss_region$Taxa_abbrev,
                                   levels = species_order_abbrev)
  ss_region <- ss_region %>% dplyr::filter(!is.na(Taxa_abbrev))

  # Compute model-predicted marginal means per Taxa x Resp x Region
  region_summary <- compute_marginal_means(ss_region, "Region", mpa_col)
  names(region_summary)[names(region_summary) == "level"] <- "Region"
  region_summary$Taxa_abbrev <- factor(region_summary$Taxa_abbrev,
                                        levels = species_order_abbrev)

  panel_b <- ggplot(ss_region, aes(x = Region, y = Mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
               linewidth = 0.4) +
    geom_jitter(aes(color = Taxa_abbrev), width = 0.12,
                size = 1.8, alpha = 0.55) +
    geom_pointrange(
      data = region_summary,
      aes(x = Region, y = mean_es,
          ymin = ci_lo,
          ymax = ci_hi),
      size = 0.5, linewidth = 0.7, color = "black",
      position = position_nudge(x = 0.25)
    ) +
    facet_grid(Resp ~ Taxa_abbrev, scales = "free_y",
               labeller = labeller(
                 Taxa_abbrev = function(x) x,
                 Resp = c("Bio" = "Biomass", "Den" = "Density")
               )) +
    scale_color_manual(values = sp_color_abbrev, guide = "none") +
    scale_y_rr(-5, 5, name = "MPA / Reference at t = 11") +
    labs(
      x = "Region",
      y = "MPA / Reference at t = 11"
    ) +
    theme_mpa(base_size = 9) +
    theme(
      strip.text.x = element_text(face = "italic", size = 8),
      strip.text.y = element_text(size = 8),
      panel.grid.major = element_blank(),
      panel.spacing.y = unit(0.8, "lines")
    )

  cat("  Panel B (region): Ch. Islands = ", sum(ss_region$Region == "Ch. Islands"),
      " obs, Mainland = ", sum(ss_region$Region == "Mainland"), " obs\n")
  cat("    Summary method breakdown: ",
      paste(unique(region_summary$method), collapse = ", "), "\n")
}

# --- Combine panels with patchwork ---
if (!is.null(panel_a) && !is.null(panel_b)) {
  fig_s09 <- panel_a / panel_b +
    patchwork::plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") +
    patchwork::plot_layout(heights = c(1, 1)) &
    theme(plot.tag = element_text(size = 9, face = "bold"))

  save_fig(fig_s09, "fig_s09_moderator_comparisons",
           FIG_S09_DIMS["w"], FIG_S09_DIMS["h"])
} else if (!is.null(panel_a)) {
  save_fig(panel_a, "fig_s09_moderator_comparisons",
           FIG_S09_DIMS["w"], FIG_S09_DIMS["h"] / 2)
  cat("  WARNING: Region panel unavailable. Saved protection panel only.\n")
} else {
  cat("  WARNING: No moderator data available. Skipping Fig S9 (moderator).\n")
}

} # end fig_s09


# =============================================================================
# Section C: Moderator Meta-Regression
# =============================================================================
# Tests whether MPA type or location significantly predict effect size.
#
# LIMITATION: Moderator analysis power
# The moderator comparisons (MPA type, geographic region) have limited power
# due to small and unbalanced subgroup sizes. SMR vs SMCA and Islands vs
# Mainland comparisons should be interpreted as exploratory. The Q_M test
# for moderator significance may be underpowered with k < 20 per subgroup.
# Results are presented for transparency but strong conclusions about
# moderator effects are not warranted.
#
# NOTE: Moderator meta-regression pools Bio and Den effect sizes from the same MPA.
# These are not fully independent (shared MPA environment), but the rma.mv model
# includes (1|MPA) random effects which partially accounts for this clustering.
# A fully multivariate approach (e.g., rma.mv with V matrix specifying within-MPA
# covariance) would be more rigorous but requires assumptions about Bio-Den correlation.

cat("\n--- Moderator Meta-Regression ---\n")

if (requireNamespace("metafor", quietly = TRUE)) {

  # Prepare data with variance column
  ss_mod <- ss_merged %>%
    dplyr::filter(!is.na(Mean), !is.na(SE), SE > 0) %>%
    dplyr::mutate(V = SE^2)

  # Standardize Taxa to abbreviated names for consistent output
  if (!"Taxa_abbrev" %in% names(ss_mod) || all(is.na(ss_mod$Taxa_abbrev))) {
    ss_mod$Taxa_abbrev <- ss_mod[[taxa_col_ss]]
  }

  mod_results <- list()
  mod_diagnostics <- list()  # Model-level test statistics & heterogeneity

  # Helper: extract model-level diagnostics from rma / rma.mv object
  extract_mod_diagnostics <- function(mod, model_label, data_used, mpa_col_name) {
    diag <- data.frame(
      Model         = model_label,
      n_obs         = mod$k,
      n_MPAs        = length(unique(data_used[[mpa_col_name]])),
      QM            = mod$QM,
      QM_df         = ifelse(!is.null(mod$m), mod$m, mod$p - 1),
      QM_pval       = mod$QMp,
      QE            = if (!is.null(mod$QE)) mod$QE else NA_real_,
      QE_df         = if (!is.null(mod$QEp)) mod$k - mod$p else NA_real_,
      QE_pval       = if (!is.null(mod$QEp)) mod$QEp else NA_real_,
      stringsAsFactors = FALSE
    )
    # Heterogeneity: rma.mv stores sigma2 vector; rma stores tau2, I2, H2
    if (inherits(mod, "rma.mv")) {
      # sigma2 components correspond to random effects in order specified
      diag$sigma2_MPA    <- if (length(mod$sigma2) >= 1) mod$sigma2[1] else NA_real_
      diag$sigma2_Source <- if (length(mod$sigma2) >= 2) mod$sigma2[2] else NA_real_
      # Compute pseudo-I2 for rma.mv (Nakagawa & Santos 2012 approximation)
      total_sigma2 <- sum(mod$sigma2)
      typical_v <- tryCatch(1 / mean(1 / data_used$V), error = function(e) NA_real_)
      diag$tau2  <- total_sigma2  # total random effects variance
      diag$I2    <- if (!is.na(typical_v) && (total_sigma2 + typical_v) > 0) {
        round(100 * total_sigma2 / (total_sigma2 + typical_v), 1)
      } else NA_real_
      diag$H2    <- if (!is.na(typical_v) && typical_v > 0) {
        round((total_sigma2 + typical_v) / typical_v, 2)
      } else NA_real_
    } else {
      diag$sigma2_MPA    <- NA_real_
      diag$sigma2_Source <- NA_real_
      diag$tau2  <- if (!is.null(mod$tau2)) mod$tau2 else NA_real_
      diag$I2    <- if (!is.null(mod$I2))   mod$I2   else NA_real_
      diag$H2    <- if (!is.null(mod$H2))   mod$H2   else NA_real_
    }
    diag
  }

  # NOTE: Biomass and density effect sizes are pooled in the moderator analysis,
  # unlike the base meta-analysis (script 09) which fits separate per-response models.
  # This pooling increases statistical power but assumes the moderator effect
  # (e.g., SMR vs SMCA) is consistent across response types. The model includes
  # ~1|MPA and ~1|Source as random effects but does NOT include (1|Resp) — meaning
  # different baselines by response type are not modeled here. Taxa_abbrev is included
  # as a fixed-effect covariate to adjust for species-level differences, but the
  # moderator coefficient represents an average across both biomass and density.
  # This is a deliberate trade-off: splitting by Resp would further reduce an
  # already-small sample size per moderator level.
  #
  # CAUTION: Pooling Bio/Den inflates effective sample size. Standard errors for
  # moderator coefficients may be underestimated. See BUG_BASH_REPORT.md B11.
  # A more rigorous approach would use ~1|MPA:Resp or separate Bio/Den models.

  # --- Model 1: Type moderator (SMR vs SMCA) ---
  ss_type_mod <- ss_mod %>%
    dplyr::filter(!is.na(type), type %in% c("SMR", "SMCA"))

  if (nrow(ss_type_mod) >= 10) {
    re_mpa <- as.formula(paste0("~1 | ", mpa_col))
    mod_type <- tryCatch({
      metafor::rma.mv(
        yi = Mean, V = V,
        mods = ~ Taxa_abbrev + type,
        random = list(re_mpa, ~1 | Source),
        data = ss_type_mod,
        method = "REML",
        test = "t"
      )
    }, error = function(e) {
      tryCatch(
        metafor::rma(yi = Mean, vi = V,
                     mods = ~ Taxa_abbrev + type,
                     data = ss_type_mod,
                     test = "knha"),
        error = function(e2) NULL
      )
    })

    if (!is.null(mod_type)) {
      type_coefs <- data.frame(
        Model = "Type (SMR vs SMCA)",
        Term = rownames(coef(summary(mod_type))),
        coef(summary(mod_type)),
        stringsAsFactors = FALSE
      )
      mod_results[["type"]] <- type_coefs
      mod_diagnostics[["type"]] <- extract_mod_diagnostics(
        mod_type, "Type (SMR vs SMCA)", ss_type_mod, mpa_col
      )
      # QM tests the omnibus null that all moderator coefficients = 0.
      # For a single binary moderator (with taxa covariates), QM tests all
      # fixed effects jointly. The specific type coefficient p-value (below)
      # isolates the SMR vs SMCA contrast, adjusting for taxa.
      cat("  Type moderator (omnibus): QM = ",
          round(mod_type$QM, 2), ", p = ",
          round(mod_type$QMp, 4), "\n")
      # Report the specific moderator coefficient p-value
      type_coef_summary <- coef(summary(mod_type))
      type_row <- grep("^type", rownames(type_coef_summary))
      if (length(type_row) > 0) {
        cat("  Type coefficient (SMR vs SMCA): estimate = ",
            round(type_coef_summary[type_row[1], "estimate"], 4),
            ", SE = ", round(type_coef_summary[type_row[1], "se"], 4),
            ", p = ", round(type_coef_summary[type_row[1], "pval"], 4), "\n")
      }
    }
  } else {
    cat("  Insufficient data for type moderator (n = ",
        nrow(ss_type_mod), ")\n")
  }

  # --- Model 2: Location moderator (mainland vs CI) ---
  ss_loc_mod <- ss_mod %>%
    dplyr::filter(!is.na(Location))

  if (nrow(ss_loc_mod) >= 10 && length(unique(ss_loc_mod$Location)) >= 2) {
    re_mpa <- as.formula(paste0("~1 | ", mpa_col))
    mod_loc <- tryCatch({
      metafor::rma.mv(
        yi = Mean, V = V,
        mods = ~ Taxa_abbrev + Location,
        random = list(re_mpa, ~1 | Source),
        data = ss_loc_mod,
        method = "REML",
        test = "t"
      )
    }, error = function(e) {
      tryCatch(
        metafor::rma(yi = Mean, vi = V,
                     mods = ~ Taxa_abbrev + Location,
                     data = ss_loc_mod,
                     test = "knha"),
        error = function(e2) NULL
      )
    })

    if (!is.null(mod_loc)) {
      loc_coefs <- data.frame(
        Model = "Location (ML vs CI)",
        Term = rownames(coef(summary(mod_loc))),
        coef(summary(mod_loc)),
        stringsAsFactors = FALSE
      )
      mod_results[["location"]] <- loc_coefs
      mod_diagnostics[["location"]] <- extract_mod_diagnostics(
        mod_loc, "Location (ML vs CI)", ss_loc_mod, mpa_col
      )
      # QM tests the omnibus null that all moderator coefficients = 0.
      # For a single binary moderator (with taxa covariates), QM tests all
      # fixed effects jointly. The specific Location coefficient p-value
      # (below) isolates the mainland vs Ch. Islands contrast, adjusting for taxa.
      cat("  Location moderator (omnibus): QM = ",
          round(mod_loc$QM, 2), ", p = ",
          round(mod_loc$QMp, 4), "\n")
      # Report the specific moderator coefficient p-value
      loc_coef_summary <- coef(summary(mod_loc))
      loc_row <- grep("^Location", rownames(loc_coef_summary))
      if (length(loc_row) > 0) {
        cat("  Location coefficient (ML vs CI): estimate = ",
            round(loc_coef_summary[loc_row[1], "estimate"], 4),
            ", SE = ", round(loc_coef_summary[loc_row[1], "se"], 4),
            ", p = ", round(loc_coef_summary[loc_row[1], "pval"], 4), "\n")
      }
    }
  } else {
    cat("  Insufficient data for location moderator (n = ",
        nrow(ss_loc_mod), ")\n")
  }

  # --- Save combined moderator coefficient table ---
  if (length(mod_results) > 0) {
    mod_table <- dplyr::bind_rows(mod_results)
    rownames(mod_table) <- NULL
    write.csv(mod_table,
              here::here("tables", "table_s_moderator_meta_regression.csv"),
              row.names = FALSE)
    cat("  Saved: tables/table_s_moderator_meta_regression.csv\n")
  } else {
    cat("  WARNING: No moderator models converged.\n")
  }

  # --- Save model-level diagnostics (QM, heterogeneity, sample sizes) ---
  if (length(mod_diagnostics) > 0) {
    diag_table <- dplyr::bind_rows(mod_diagnostics)
    rownames(diag_table) <- NULL
    write.csv(diag_table,
              here::here("tables", "table_s_moderator_diagnostics.csv"),
              row.names = FALSE)
    cat("  Saved: tables/table_s_moderator_diagnostics.csv\n")
  }

  # --- Save subgroup sample sizes and effect sizes ---
  # These correspond to the summary statistics plotted in Fig S9
  subgroup_rows <- list()

  # Type subgroups (SMR vs SMCA)
  if (exists("ss_type_mod") && nrow(ss_type_mod) >= 10) {
    type_subgroup <- ss_type_mod %>%
      dplyr::group_by(Taxa_abbrev, type, Resp) %>%
      dplyr::summarise(
        n_obs      = dplyr::n(),
        n_MPAs     = length(unique(.data[[mpa_col]])),
        mean_lnRR  = {
          w <- 1 / (as.numeric(SE)^2)
          weighted.mean(as.numeric(Mean), w, na.rm = TRUE)
        },
        sd_lnRR    = sd(as.numeric(Mean), na.rm = TRUE),
        se_lnRR    = {
          w <- 1 / (as.numeric(SE)^2)
          sqrt(1 / sum(w, na.rm = TRUE))
        },
        ci_lo_lnRR = {
          w <- 1 / (as.numeric(SE)^2)
          wm <- weighted.mean(as.numeric(Mean), w, na.rm = TRUE)
          wse <- sqrt(1 / sum(w, na.rm = TRUE))
          wm - qt(0.975, dplyr::n() - 1) * wse
        },
        ci_hi_lnRR = {
          w <- 1 / (as.numeric(SE)^2)
          wm <- weighted.mean(as.numeric(Mean), w, na.rm = TRUE)
          wse <- sqrt(1 / sum(w, na.rm = TRUE))
          wm + qt(0.975, dplyr::n() - 1) * wse
        },
        mean_RR    = {
          w <- 1 / (as.numeric(SE)^2)
          exp(weighted.mean(as.numeric(Mean), w, na.rm = TRUE))
        },
        .groups    = "drop"
      ) %>%
      dplyr::mutate(Moderator = "Protection level", Subgroup = type) %>%
      dplyr::select(Moderator, Subgroup, Taxa_abbrev, Resp, dplyr::everything(),
                    -type)
    subgroup_rows[["type"]] <- type_subgroup
  }

  # Location subgroups (Channel Islands vs Mainland)
  if (exists("ss_loc_mod") && nrow(ss_loc_mod) >= 10) {
    loc_subgroup <- ss_loc_mod %>%
      dplyr::mutate(Region = dplyr::case_when(
        Location == "C"  ~ "Ch. Islands",
        Location == "ML" ~ "Mainland",
        TRUE ~ Location
      )) %>%
      dplyr::group_by(Taxa_abbrev, Region, Resp) %>%
      dplyr::summarise(
        n_obs      = dplyr::n(),
        n_MPAs     = length(unique(.data[[mpa_col]])),
        mean_lnRR  = {
          w <- 1 / (as.numeric(SE)^2)
          weighted.mean(as.numeric(Mean), w, na.rm = TRUE)
        },
        sd_lnRR    = sd(as.numeric(Mean), na.rm = TRUE),
        se_lnRR    = {
          w <- 1 / (as.numeric(SE)^2)
          sqrt(1 / sum(w, na.rm = TRUE))
        },
        ci_lo_lnRR = {
          w <- 1 / (as.numeric(SE)^2)
          wm <- weighted.mean(as.numeric(Mean), w, na.rm = TRUE)
          wse <- sqrt(1 / sum(w, na.rm = TRUE))
          wm - qt(0.975, dplyr::n() - 1) * wse
        },
        ci_hi_lnRR = {
          w <- 1 / (as.numeric(SE)^2)
          wm <- weighted.mean(as.numeric(Mean), w, na.rm = TRUE)
          wse <- sqrt(1 / sum(w, na.rm = TRUE))
          wm + qt(0.975, dplyr::n() - 1) * wse
        },
        mean_RR    = {
          w <- 1 / (as.numeric(SE)^2)
          exp(weighted.mean(as.numeric(Mean), w, na.rm = TRUE))
        },
        .groups    = "drop"
      ) %>%
      dplyr::mutate(Moderator = "Location", Subgroup = Region) %>%
      dplyr::select(Moderator, Subgroup, Taxa_abbrev, Resp, dplyr::everything(),
                    -Region)
    subgroup_rows[["location"]] <- loc_subgroup
  }

  if (length(subgroup_rows) > 0) {
    subgroup_table <- dplyr::bind_rows(subgroup_rows)
    write.csv(subgroup_table,
              here::here("tables", "table_s_moderator_subgroup_effects.csv"),
              row.names = FALSE)
    cat("  Saved: tables/table_s_moderator_subgroup_effects.csv\n")
  }

} else {
  cat("  WARNING: Package 'metafor' not available. Skipping moderator analysis.\n")
}


# =============================================================================
# Completion
# =============================================================================
cat("\n=== 13_additional_analyses.R complete ===\n")
