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
#   - Site: MPA metadata from 03_data_import.R
#   - Color palette objects from 00b_color_palette.R
#
# OUTPUTS:
#   Figures: fig_s09_moderator_comparisons
#   Tables:  table_s_moderator_meta_regression.csv
#
# NOTE: Previous versions generated S10 (MPA size), S11 (regional, now merged
#       into S9), S12 (rate-of-change, redundant with Fig 5), and S13 (cascade
#       completeness, redundant with S6 panel b). These were removed to
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
FIG_S09_DIMS <- c(w = 18, h = 20)   # Combined moderator comparisons (2-row)

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
needed_cols <- c("type", "Hectares", "Location", "MPA_Start", "ChannelIsland")
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
  taxa_col_ss <- if ("Taxa" %in% names(ss_merged)) "Taxa" else mpa_col
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
# Section B: Combined Moderator Comparison (Fig S9)
# =============================================================================
# Two-row figure: (top) SMR vs SMCA protection level, (bottom) Ch. Islands vs Mainland

if (should_render("fig_s09")) {
cat("\n--- Figure S9: Combined Moderator Comparisons ---\n")

# --- Shared axis breaks for RR scale ---
rr_breaks <- log(c(0.01, 0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8, 20, 100))
rr_labels <- c("0.01", "0.05", "0.1", "0.25", "0.5", "1", "2", "4", "8", "20", "100")

# --- Panel A: SMR vs SMCA ---
ss_type <- ss_merged %>%
  dplyr::filter(!is.na(type), type %in% c("SMR", "SMCA"))

panel_a <- NULL
if (nrow(ss_type) > 0) {
  ss_type$Taxa_abbrev <- factor(ss_type$Taxa_abbrev,
                                 levels = species_order_abbrev)
  ss_type <- ss_type %>% dplyr::filter(!is.na(Taxa_abbrev))

  type_summary <- ss_type %>%
    dplyr::group_by(Taxa_abbrev, type, Resp) %>%
    dplyr::summarise(
      mean_es = mean(Mean, na.rm = TRUE),
      se_es   = sd(Mean, na.rm = TRUE) / sqrt(dplyr::n()),
      n       = dplyr::n(),
      .groups = "drop"
    )

  panel_a <- ggplot(ss_type, aes(x = type, y = Mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
               linewidth = 0.4) +
    geom_jitter(aes(color = Taxa_abbrev), width = 0.15,
                size = 1.5, alpha = 0.5) +
    geom_pointrange(
      data = type_summary,
      aes(x = type, y = mean_es,
          ymin = mean_es - se_es,
          ymax = mean_es + se_es),
      size = 0.5, linewidth = 0.7, color = "black",
      position = position_nudge(x = 0.3)
    ) +
    facet_grid(Resp ~ Taxa_abbrev, scales = "free_y",
               labeller = labeller(
                 Taxa_abbrev = function(x) x,
                 Resp = c("Bio" = "Biomass", "Den" = "Density")
               )) +
    scale_color_manual(values = sp_color_abbrev, guide = "none") +
    scale_y_continuous(breaks = rr_breaks, labels = rr_labels) +
    labs(
      x = "MPA protection level",
      y = "Response ratio at t = 11"
    ) +
    theme_mpa(base_size = 9) +
    theme(
      strip.text.x = element_text(face = "italic", size = 8),
      strip.text.y = element_text(size = 8),
      panel.grid.major = element_blank()
    )

  cat("  Panel A (protection): SMR = ", sum(ss_type$type == "SMR"),
      " obs, SMCA = ", sum(ss_type$type == "SMCA"), " obs\n")
}

# --- Panel B: Channel Islands vs Mainland ---
ss_region <- ss_merged %>%
  dplyr::filter(!is.na(Location))

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

  region_summary <- ss_region %>%
    dplyr::group_by(Taxa_abbrev, Region, Resp) %>%
    dplyr::summarise(
      mean_es = mean(Mean, na.rm = TRUE),
      se_es   = sd(Mean, na.rm = TRUE) / sqrt(dplyr::n()),
      n       = dplyr::n(),
      .groups = "drop"
    )

  panel_b <- ggplot(ss_region, aes(x = Region, y = Mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
               linewidth = 0.4) +
    geom_jitter(aes(color = Taxa_abbrev), width = 0.12,
                size = 1.5, alpha = 0.5) +
    geom_pointrange(
      data = region_summary,
      aes(x = Region, y = mean_es,
          ymin = mean_es - se_es,
          ymax = mean_es + se_es),
      size = 0.5, linewidth = 0.7, color = "black",
      position = position_nudge(x = 0.25)
    ) +
    facet_grid(Resp ~ Taxa_abbrev, scales = "free_y",
               labeller = labeller(
                 Taxa_abbrev = function(x) x,
                 Resp = c("Bio" = "Biomass", "Den" = "Density")
               )) +
    scale_color_manual(values = sp_color_abbrev, guide = "none") +
    scale_y_continuous(breaks = rr_breaks, labels = rr_labels) +
    labs(
      x = "Region",
      y = "Response ratio at t = 11"
    ) +
    theme_mpa(base_size = 9) +
    theme(
      strip.text.x = element_text(face = "italic", size = 8),
      strip.text.y = element_text(size = 8),
      panel.grid.major = element_blank()
    )

  cat("  Panel B (region): Ch. Islands = ", sum(ss_region$Region == "Ch. Islands"),
      " obs, Mainland = ", sum(ss_region$Region == "Mainland"), " obs\n")
}

# --- Combine panels with patchwork ---
if (!is.null(panel_a) && !is.null(panel_b)) {
  fig_s09 <- panel_a / panel_b +
    patchwork::plot_annotation(tag_levels = "a") +
    patchwork::plot_layout(heights = c(1, 1))

  save_fig(fig_s09, "fig_s09_moderator_comparisons",
           FIG_S09_DIMS["w"], FIG_S09_DIMS["h"])
} else if (!is.null(panel_a)) {
  save_fig(panel_a, "fig_s09_moderator_comparisons",
           FIG_S09_DIMS["w"], FIG_S09_DIMS["h"] / 2)
  cat("  WARNING: Region panel unavailable. Saved protection panel only.\n")
} else {
  cat("  WARNING: No moderator data available. Skipping Fig S9.\n")
}

} # end fig_s09


# =============================================================================
# Section C: Moderator Meta-Regression
# =============================================================================
# Tests whether MPA type or location significantly predict effect size.

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
      cat("  Type moderator: QM = ",
          round(mod_type$QM, 2), ", p = ",
          round(mod_type$QMp, 4), "\n")
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
      cat("  Location moderator: QM = ",
          round(mod_loc$QM, 2), ", p = ",
          round(mod_loc$QMp, 4), "\n")
    }
  } else {
    cat("  Insufficient data for location moderator (n = ",
        nrow(ss_loc_mod), ")\n")
  }

  # --- Save combined moderator table ---
  if (length(mod_results) > 0) {
    mod_table <- dplyr::bind_rows(mod_results)
    rownames(mod_table) <- NULL
    write.csv(mod_table,
              here::here("data", "table_s_moderator_meta_regression.csv"),
              row.names = FALSE)
    cat("  Saved: data/table_s_moderator_meta_regression.csv\n")
  } else {
    cat("  WARNING: No moderator models converged.\n")
  }

} else {
  cat("  WARNING: Package 'metafor' not available. Skipping moderator analysis.\n")
}


# =============================================================================
# Completion
# =============================================================================
cat("\n=== 13_additional_analyses.R complete ===\n")
