# =============================================================================
# 10_figures.R
# =============================================================================
#
# PURPOSE:
#   Generate publication-quality figures for the Conservation Letters manuscript
#   on MPA effects on kelp forest trophic cascades.
#
# WHAT THIS SCRIPT DOES:
#   Produces publication figures for the manuscript:
#
#   Figure 1: Map of MPAs with Channel Islands + inset time series
#     - Base map showing Southern California coastline and Channel Islands
#     - MPA site markers (shape indicates data source: NPS-KFM, LTER, PISCO)
#     - Inset time series panels showing kelp biomass at 6 featured MPAs
#     - Scale bar and site labels
#
#   Figure 2: Data processing pipeline example
#     - 4-panel illustration using KFM purple urchin at Scorpion SMR
#     - (a) Raw density time series
#     - (b) Proportion of maximum
#     - (c) Log response ratio
#     - (d) Log response ratio with linear trend fit
#
#   Figure 3: Forest plot of effect sizes
#     - Individual effect sizes by MPA and taxa
#     - Color-coded by response type (density vs biomass)
#     - Shape-coded by data source (PISCO, KFM, LTER, Landsat)
#
#   Figure 4: Mean effect sizes from meta-analysis
#     - Summarizes Table 2 graphically
#     - Shows meta-analytic means with 95% CIs
#     - Individual effect sizes shown as background points
#
#   Figure 5: All taxa time series at example MPAs
#     - Log response ratios for all 5 taxa at 3 example MPAs
#     - Demonstrates the trophic cascade dynamics
#
# DESIGN PRINCIPLES:
#   - Uses colorblind-safe palette from 00b_color_palette.R
#   - Publication-ready sizing for Conservation Letters (80-170mm width)
#   - Consistent theme_mpa() styling across all figures
#   - Exported as both PDF (vector) and PNG (raster at 300 DPI)
#
# INPUTS:
#   - All.RR.sub.trans: Response ratio data
#   - All.Resp.sub: Raw response data
#   - SumStats.Final: Effect size estimates
#   - Table2: Meta-analysis summary
#   - Site: Site metadata
#   - Color palette objects from 00b_color_palette.R
#
# OUTPUTS (saved to plots/ directory):
#   - fig_01_mpa_map_composite.pdf / .png (main Figure 1 with insets)
#   - fig_01_map_only.pdf / .png (map without insets)
#   - fig_02_data_processing.pdf / .png
#   - fig_03_mean_effects.pdf / .png
#   - fig_04_urchin_kelp_scatter.pdf / .png
#   - fig_s01_forest_plot.pdf / .png
#   - fig_s02_all_taxa_timeseries.pdf / .png
#
# DEPENDENCIES:
#   Requires 00-09 scripts to be sourced first
#
# AUTHORS: Emily Donham & Adrian Stier
# PROJECT: CA MPA Kelp Forest pBACIPS Analysis
# =============================================================================

# =============================================================================
# Setup
# =============================================================================

dir.create(here::here("plots"), showWarnings = FALSE)

# Standard factor levels used throughout
taxa_levels  <- c("S. purpuratus", "M. franciscanus", "M. pyrifera",
                   "P. interruptus", "S. pulcher")
source_levels <- c("KFM", "LTER", "PISCO", "Landsat")

cat("=== Starting figure generation ===\n")

# Verify color palette is loaded from 00b_color_palette.R
if (!exists("col_taxa") || !exists("col_response") || !exists("theme_mpa")) {
  stop("Color palette not loaded. Please source 00b_color_palette.R first.")
}
cat("  Color palette verified: col_taxa, col_response, col_site loaded\n")

# =============================================================================
# Input validation: Check required data objects exist and have expected structure
# =============================================================================

# Required data objects from previous scripts
required_objects <- c("All.RR.sub.trans", "All.Resp.sub", "SumStats.Final", "Table2", "Site")
missing_objects <- required_objects[!sapply(required_objects, exists, envir = globalenv())]
if (length(missing_objects) > 0) {
  stop("Missing required data objects: ", paste(missing_objects, collapse = ", "),
       "\nPlease run scripts 00-09 first.")
}

# Validate Table2 structure (meta-analysis summary)
if (!is.data.frame(Table2) || nrow(Table2) == 0) {
  stop("Table2 must be a non-empty dataframe. Check 09_meta_analysis.R output.")
}
required_table2_cols <- c("Taxa", "Response", "Estimate", "CI_lower", "CI_upper")
missing_cols <- required_table2_cols[!required_table2_cols %in% names(Table2)]
if (length(missing_cols) > 0) {
  stop("Table2 missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Validate SumStats.Final structure
if (!is.data.frame(SumStats.Final) || nrow(SumStats.Final) == 0) {
  stop("SumStats.Final must be a non-empty dataframe. Check 08_effect_sizes.R output.")
}
required_sumstats_cols <- c("Taxa", "MPA", "Mean", "SE", "CI", "Source", "Resp")
missing_cols <- required_sumstats_cols[!required_sumstats_cols %in% names(SumStats.Final)]
if (length(missing_cols) > 0) {
  stop("SumStats.Final missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Validate All.RR.sub.trans structure
if (!is.data.frame(All.RR.sub.trans) || nrow(All.RR.sub.trans) == 0) {
  stop("All.RR.sub.trans must be a non-empty dataframe. Check 07_combine_data.R output.")
}

# Validate All.Resp.sub structure
if (!is.data.frame(All.Resp.sub) || nrow(All.Resp.sub) == 0) {
  stop("All.Resp.sub must be a non-empty dataframe. Check 07_combine_data.R output.")
}

cat("  Input data validation passed\n")

# =============================================================================
# Simple save function - PDF and PNG
save_fig <- function(plot, name, w, h) {
  ggsave(here::here("plots", paste0(name, ".pdf")), plot, width = w, height = h, units = "cm")
  ggsave(here::here("plots", paste0(name, ".png")), plot, width = w, height = h, units = "cm", dpi = 300)
  cat("  Saved:", name, "\n")
}

# =============================================================================
# Figure 1: Map of MPAs with Channel Islands + Inset Time Series
# =============================================================================

cat("Building Figure 1: MPA Map with inset time series...\n")

# Load mapping packages
library(maps)
library(mapdata)
library(cowplot)

# --- 1.1 Site data with coordinates ---
sites_map <- read.csv(here::here("data", "Site_List_All.csv"))

# Define which data source covers which MPAs
kfm_sites <- c("Harris Point SMR", "South Point SMR", "Gull Island SMR",
               "Scorpion SMR", "Santa Barbara Island SMR", "Anacapa Island SMR 2003")
lter_sites <- c("Campus Point SMCA", "Naples SMCA")
pisco_sites <- c("Point Vicente SMCA", "Carrington Pt SMR", "Painted Cave SMCA",
                 "Skunk Pt SMR", "Anacapa Island SMCA")

sites_map <- sites_map %>%
  dplyr::mutate(
    data_source = case_when(
      CA_MPA_Name_Short %in% kfm_sites ~ "NPS-KFM",
      CA_MPA_Name_Short %in% lter_sites ~ "LTER",
      CA_MPA_Name_Short %in% pisco_sites ~ "PISCO",
      TRUE ~ "Other"
    )
  ) %>%
  dplyr::filter(data_source != "Other")

# Short labels for sites
label_lookup <- c(
  "Campus Point SMCA" = "CP", "Point Vicente SMCA" = "PV",
  "Harris Point SMR" = "HP", "South Point SMR" = "SP",
  "Gull Island SMR" = "GI", "Santa Barbara Island SMR" = "SB",
  "Scorpion SMR" = "S", "Anacapa Island SMR 2003" = "AI",
  "Naples SMCA" = "NR", "Painted Cave SMCA" = "PC",
  "Carrington Pt SMR" = "AC", "Skunk Pt SMR" = "PD",
  "Anacapa Island SMCA" = "AI"
)
sites_map$label <- label_lookup[sites_map$CA_MPA_Name_Short]

# --- 1.2 Get basemap data ---
ca_state <- map_data("state", region = "california")

# Get Channel Islands
channel_islands <- c(
  "California:Santa Cruz Island", "California:Santa Rosa Island",
  "California:San Miguel Island", "California:Santa Catalina Island",
  "California:San Clemente Island", "California:San Nicolas Island"
)
islands_data <- map_data("worldHires", region = channel_islands)

# Add Santa Barbara Island manually
santa_barbara_island <- data.frame(
  long = c(-119.05, -119.02, -119.00, -119.02, -119.05),
  lat = c(33.46, 33.48, 33.47, 33.45, 33.46),
  group = max(islands_data$group) + 1,
  order = 1:5, region = "Santa Barbara Island", subregion = NA
)
islands_data <- bind_rows(islands_data, santa_barbara_island)

# --- 1.3 Map aesthetics ---
land_color <- "#C8C0B0"
ocean_color <- "#9DBCD4"
coastline_color <- "#555555"
source_colors <- c("NPS-KFM" = "#5B8B9F", "LTER" = "#D4843E", "PISCO" = "#6B9E5A")
source_shapes <- c("NPS-KFM" = 22, "LTER" = 21, "PISCO" = 24)

# --- 1.4 Create base map ---
base_map <- ggplot() +
  geom_rect(aes(xmin = -121.5, xmax = -117.3, ymin = 32.5, ymax = 35.5),
            fill = ocean_color, color = NA) +
  geom_polygon(data = ca_state, aes(x = long, y = lat, group = group),
               fill = land_color, color = coastline_color, linewidth = 0.4) +
  geom_polygon(data = islands_data, aes(x = long, y = lat, group = group),
               fill = land_color, color = coastline_color, linewidth = 0.3) +
  geom_point(data = sites_map,
             aes(x = Lon, y = Lat, shape = data_source, fill = data_source),
             size = 4.5, color = "grey20", stroke = 1) +
  geom_text(data = sites_map, aes(x = Lon, y = Lat, label = label),
            nudge_x = 0.14, nudge_y = 0.08, size = 3.2, fontface = "bold", color = "grey10") +
  scale_shape_manual(name = "Data Source", values = source_shapes,
                     labels = c("NPS-KFM" = "KFM (NPS)", "LTER" = "SBC LTER", "PISCO" = "PISCO")) +
  scale_fill_manual(name = "Data Source", values = source_colors,
                    labels = c("NPS-KFM" = "KFM (NPS)", "LTER" = "SBC LTER", "PISCO" = "PISCO")) +
  coord_fixed(ratio = 1.15, xlim = c(-121.1, -117.5), ylim = c(32.85, 35.15), expand = FALSE) +
  scale_x_continuous(breaks = c(-121, -120, -119, -118),
                     labels = c("121.0\u00B0W", "120.0\u00B0W", "119.0\u00B0W", "118.0\u00B0W")) +
  scale_y_continuous(breaks = c(33, 34, 35),
                     labels = c("33.0\u00B0N", "34.0\u00B0N", "35.0\u00B0N")) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_line(color = alpha("grey50", 0.3), linewidth = 0.15),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = ocean_color, color = coastline_color, linewidth = 0.5),
    legend.position = c(0.12, 0.18),  # Moved slightly for better visibility
    legend.background = element_rect(fill = alpha("white", 0.97), color = "grey40", linewidth = 0.4),
    legend.key.size = unit(0.5, "cm"),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    legend.margin = margin(6, 8, 6, 6),
    axis.text = element_text(size = 8, color = "grey20"),
    plot.margin = margin(3, 3, 3, 3)
  ) +
  guides(
    shape = guide_legend(order = 1, override.aes = list(size = 4)),
    fill = guide_legend(order = 1, override.aes = list(size = 4))
  )

# Add scale bar
km_per_deg <- 92
scale_deg <- 40 / km_per_deg
map_with_scale <- base_map +
  annotate("rect", xmin = -118.1, xmax = -117.55, ymin = 32.88, ymax = 33.06,
           fill = alpha("white", 0.92), color = "grey40", linewidth = 0.25) +
  annotate("rect", xmin = -118.05, xmax = -118.05 + scale_deg/2, ymin = 32.93, ymax = 32.97,
           fill = "black", color = NA) +
  annotate("rect", xmin = -118.05 + scale_deg/2, xmax = -118.05 + scale_deg,
           ymin = 32.93, ymax = 32.97, fill = "white", color = "black", linewidth = 0.25) +
  annotate("text", x = -118.05, y = 33.01, label = "0", size = 2.2, hjust = 0.5) +
  annotate("text", x = -118.05 + scale_deg/2, y = 33.01, label = "20", size = 2.2, hjust = 0.5) +
  annotate("text", x = -118.05 + scale_deg, y = 33.01, label = "40 km", size = 2.2, hjust = 0.5)

# --- 1.5 Create time series insets ---
inset_mpas <- c("Campus Point SMCA", "Point Vicente SMCA", "Harris Point SMR",
                "South Point SMR", "Gull Island SMR", "Santa Barbara Island SMR")
panel_labels <- c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)")
# Status color mapping - map various status labels to Inside/Outside
ts_colors <- c("Inside" = "#D4A84B", "Outside" = "#5A5A78")

# Helper function to standardize status values
standardize_status <- function(status) {
  status <- as.character(status)
  result <- dplyr::case_when(
    is.na(status) ~ NA_character_,
    tolower(status) %in% c("inside", "mpa", "impact", "i") ~ "Inside",
    tolower(status) %in% c("outside", "reference", "control", "ref", "o", "r") ~ "Outside",
    TRUE ~ NA_character_  # Convert unrecognized to NA rather than keeping original

  )
  return(result)
}

create_ts_inset <- function(mpa_name, panel_label) {
  # Get MPA start year from Site table first, then fall back to sites_map
  mpa_start <- Site %>%
    dplyr::filter(CA_MPA_Name_Short == mpa_name) %>%
    dplyr::pull(MPA_Start)
  if (length(mpa_start) == 0 || is.na(mpa_start[1])) {
    mpa_start <- sites_map %>%
      dplyr::filter(CA_MPA_Name_Short == mpa_name) %>%
      dplyr::pull(MPA_Start)
  }
  if (length(mpa_start) == 0 || is.na(mpa_start[1])) {
    mpa_start <- ifelse(grepl("Point Vicente|Campus", mpa_name), 2012, 2003)
  } else {
    mpa_start <- mpa_start[1]
  }

  # Debug: Show available taxon names and resp values for this MPA
  mpa_data <- All.Resp.sub %>% dplyr::filter(CA_MPA_Name_Short == mpa_name)
  cat("    Inset", panel_label, mpa_name, "- available taxon_names:",
      paste(unique(mpa_data$taxon_name), collapse = ", "), "\n")
  cat("      Available resp types:", paste(unique(mpa_data$resp), collapse = ", "), "\n")

  # Try multiple species name variants for kelp biomass (case-insensitive matching)
  kelp_names <- c("Macrocystis pyrifera", "M. pyrifera", "MACPYRAD",
                  "macrocystis pyrifera", "Giant Kelp", "giant kelp")
  bio_types <- c("Bio", "Biomass", "biomass", "bio")
  den_types <- c("Den", "Density", "density", "den")

  d <- All.Resp.sub %>%
    dplyr::filter(
      tolower(taxon_name) %in% tolower(kelp_names),
      tolower(resp) %in% tolower(bio_types),
      CA_MPA_Name_Short == mpa_name
    )

  # If no biomass data, try density data
  if (nrow(d) == 0) {
    d <- All.Resp.sub %>%
      dplyr::filter(
        tolower(taxon_name) %in% tolower(kelp_names),
        tolower(resp) %in% tolower(den_types),
        CA_MPA_Name_Short == mpa_name
      )
  }

  # Print data availability and status values
  cat("      Kelp data rows found:", nrow(d), "\n")
  if (nrow(d) > 0) {
    cat("      Raw status values:", paste(unique(d$status), collapse = ", "), "\n")
  }

  # Standardize status values to Inside/Outside
  if (nrow(d) > 0) {
    d$status <- standardize_status(d$status)
    cat("      Standardized status values:", paste(unique(d$status), collapse = ", "), "\n")
    # Remove rows with NA status (unrecognized values)
    d <- d[!is.na(d$status), ]
    if (nrow(d) > 0) {
      d$status <- factor(d$status, levels = c("Inside", "Outside"))
    }
  }

  if (nrow(d) == 0) {
    # Return informative placeholder with MPA name
    short_name <- gsub(" SMCA| SMR", "", mpa_name)
    return(ggplot() +
           annotate("text", x = 0.5, y = 0.6,
                    label = paste0(panel_label, " ", short_name),
                    size = 3, fontface = "bold") +
           annotate("text", x = 0.5, y = 0.4,
                    label = "No kelp data",
                    size = 2.5, color = "grey50") +
           theme_void() +
           theme(panel.background = element_rect(fill = "grey95", color = "grey60", linewidth = 0.4)))
  }

  short_name <- gsub(" SMCA| SMR", "", mpa_name)
  short_code <- label_lookup[mpa_name]
  if (is.na(short_code)) short_code <- ""

  # Improved inset plot with better typography
  # Use na.value for any remaining NA values
  ggplot(d, aes(x = year, y = value, color = status)) +
    geom_vline(xintercept = mpa_start, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_line(aes(group = interaction(status, source)), linewidth = 0.6, alpha = 0.9) +
    geom_point(size = 1.8, alpha = 0.95) +
    scale_color_manual(values = ts_colors, na.value = "grey50", drop = FALSE, guide = "none") +
    labs(title = paste0(panel_label, " ", short_name),
         x = NULL, y = expression('Biomass (g '~m^{-2}~')')) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    theme_mpa(base_size = 8) +
    theme(
      plot.title = element_text(size = 9, face = "bold", hjust = 0, margin = margin(0,0,3,0)),
      axis.title = element_text(size = 7),
      axis.title.y = element_text(margin = margin(0,4,0,0)),
      axis.text = element_text(size = 6),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.15, color = "grey85"),
      panel.background = element_rect(fill = "white", color = "grey50", linewidth = 0.4),
      plot.margin = margin(4, 6, 4, 4)
    )
}

insets <- lapply(seq_along(inset_mpas), function(i) {
  create_ts_inset(inset_mpas[i], panel_labels[i])
})

# --- 1.6 Create inset legend for Inside/Outside ---
inset_legend_data <- data.frame(
  x = c(1, 2), y = c(1, 1),
  status = factor(c("Inside", "Outside"), levels = c("Inside", "Outside"))
)
inset_legend_plot <- ggplot(inset_legend_data, aes(x = x, y = y, color = status)) +
  geom_point(size = 3) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = ts_colors, name = "Site Status") +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 9),
    legend.key.width = unit(1.0, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0)
  )
inset_legend <- cowplot::get_legend(inset_legend_plot)

# --- 1.7 Compose final figure ---
fig1_composite <- ggdraw() +
  draw_plot(map_with_scale, x = 0, y = 0.05, width = 1, height = 0.95) +
  draw_plot(insets[[1]], x = 0.42, y = 0.75, width = 0.26, height = 0.22) +
  draw_plot(insets[[2]], x = 0.72, y = 0.60, width = 0.26, height = 0.22) +
  draw_plot(insets[[3]], x = 0.02, y = 0.48, width = 0.26, height = 0.22) +
  draw_plot(insets[[4]], x = 0.02, y = 0.18, width = 0.26, height = 0.22) +
  draw_plot(insets[[5]], x = 0.30, y = 0.18, width = 0.26, height = 0.22) +
  draw_plot(insets[[6]], x = 0.58, y = 0.18, width = 0.26, height = 0.22) +
  draw_plot(inset_legend, x = 0.60, y = 0.01, width = 0.40, height = 0.06)

# Conservation Letters: max 170mm double-column width
save_fig(fig1_composite, "fig_01_mpa_map_composite", 17, 14)

# Also save map-only version
save_fig(map_with_scale, "fig_01_map_only", 20, 16)

# =============================================================================
# Figure 2: Data processing example (KFM, S. purpuratus, Scorpion SMR)
# =============================================================================

cat("Building Figure 2: Data processing pipeline...\n")

# MPA implementation year for Scorpion SMR (Channel Islands)
scorpion_start <- Site %>%
  dplyr::filter(CA_MPA_Name_Short == "Scorpion SMR") %>%
  dplyr::pull(MPA_Start)
if (length(scorpion_start) == 0 || is.na(scorpion_start[1])) {
  scorpion_start <- 2003
} else {
  scorpion_start <- scorpion_start[1]
}

# Panel (a): Raw density time series from All.Resp.sub
fig2a_data <- All.Resp.sub %>%
  dplyr::filter(
    source == "KFM",
    CA_MPA_Name_Short == "Scorpion SMR",
    taxon_name %in% c("Strongylocentrotus purpuratus", "S. purpuratus"),
    resp == "Den"
  ) %>%
  dplyr::mutate(year = as.numeric(as.character(year))) %>%  # Ensure year is numeric
  dplyr::group_by(year, status) %>%
  dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(
    status = standardize_status(status),
    status = factor(status, levels = c("Inside", "Outside"))
  )

# Debug: Check status values
cat("  Figure 2 - Status values in fig2a_data:", paste(unique(fig2a_data$status), collapse = ", "), "\n")
cat("  Figure 2 - Year range:", min(fig2a_data$year, na.rm = TRUE), "-", max(fig2a_data$year, na.rm = TRUE), "\n")

p2a <- ggplot(fig2a_data, aes(x = year, y = value, color = status)) +
  geom_vline(xintercept = scorpion_start, linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.5) +
  scale_color_site(name = "Status") +
  labs(title = "(a) Raw density",
       x = NULL, y = expression('Density (ind '~m^{-2}~')')) +
  scale_x_continuous(breaks = seq(2000, 2025, by = 10), limits = c(1995, 2025)) +
  theme_mpa(base_size = 9) +
  theme(legend.position = "none",
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.text.x = element_text(angle = 0, hjust = 0.5))

# Panel (b): Proportion of maximum
fig2b_data <- fig2a_data %>%
  dplyr::group_by(status) %>%
  dplyr::mutate(
    prop = value / max(value, na.rm = TRUE),
    prop = ifelse(is.na(prop) | is.infinite(prop), 0, prop)
  ) %>%
  dplyr::ungroup()

p2b <- ggplot(fig2b_data, aes(x = year, y = prop, color = status)) +
  geom_vline(xintercept = scorpion_start, linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.5) +
  scale_color_site(name = "Status") +
  labs(title = "(b) Standardized",
       x = NULL, y = "Proportion of max") +
  scale_x_continuous(breaks = seq(2000, 2025, by = 10), limits = c(1995, 2025)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5), labels = scales::number_format(accuracy = 0.1)) +
  theme_mpa(base_size = 9) +
  theme(legend.position = "none",
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.text.x = element_text(angle = 0, hjust = 0.5))

# Panel (c) and (d): Log response ratio from All.RR.sub.trans
fig2cd_data <- All.RR.sub.trans %>%
  dplyr::filter(
    source == "KFM",
    CA_MPA_Name_Short == "Scorpion SMR",
    y %in% c("Strongylocentrotus purpuratus", "S. purpuratus"),
    resp == "Den"
  ) %>%
  dplyr::mutate(
    BA = factor(BA, levels = c("Before", "After"))
  )

# If BA column is missing, derive from MPA_Start
if (!"BA" %in% names(fig2cd_data) || all(is.na(fig2cd_data$BA))) {
  fig2cd_data <- fig2cd_data %>%
    dplyr::mutate(BA = factor(
      ifelse(year < scorpion_start, "Before", "After"),
      levels = c("Before", "After")
    ))
}

p2c <- ggplot(fig2cd_data, aes(x = year, y = lnDiff, shape = BA)) +
  geom_vline(xintercept = scorpion_start, linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey70",
             linewidth = 0.4) +
  geom_point(size = 2.5, color = col_taxa["S. purpuratus"]) +
  scale_shape_manual(
    name = "Period",
    values = c("Before" = 1, "After" = 16)
  ) +
  labs(title = "(c) Log response ratio",
       x = NULL, y = "ln(MPA / Reference)") +
  scale_x_continuous(breaks = seq(2000, 2025, by = 10), limits = c(1995, 2025)) +
  theme_mpa(base_size = 9) +
  theme(legend.position = "none",
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.text.x = element_text(angle = 0, hjust = 0.5))

# Panel (d): Same with linear trend in After period
fig2d_after <- dplyr::filter(fig2cd_data, BA == "After")

p2d <- ggplot(fig2cd_data, aes(x = year, y = lnDiff, shape = BA)) +
  geom_vline(xintercept = scorpion_start, linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey70",
             linewidth = 0.4) +
  geom_point(size = 2.5, color = col_taxa["S. purpuratus"]) +
  {
    if (nrow(fig2d_after) >= 3) {
      geom_smooth(
        data = fig2d_after,
        aes(x = year, y = lnDiff),
        method = "lm", se = TRUE,
        color = col_taxa["S. purpuratus"],
        fill = col_taxa["S. purpuratus"],
        alpha = 0.2, linewidth = 0.8,
        inherit.aes = FALSE
      )
    }
  } +
  scale_shape_manual(
    name = "Period",
    values = c("Before" = 1, "After" = 16)
  ) +
  labs(title = "(d) Effect size",
       x = NULL, y = "ln(MPA / Reference)") +
  scale_x_continuous(breaks = seq(2000, 2025, by = 10), limits = c(1995, 2025)) +
  theme_mpa(base_size = 9) +
  theme(legend.position = "none",
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.text.x = element_text(angle = 0, hjust = 0.5))

# Create a shared legend for Inside/Outside status (panels a & b)
# and Before/After period (panels c & d)
legend_status <- ggplot(fig2a_data, aes(x = year, y = value, color = status)) +
  geom_point(size = 3) +
  scale_color_site(name = "Site") +
  theme_mpa() +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 9),
        legend.text = element_text(size = 9))

legend_period <- ggplot(fig2cd_data, aes(x = year, y = lnDiff, shape = BA)) +
  geom_point(size = 3, color = "grey30") +
  scale_shape_manual(name = "Period", values = c("Before" = 1, "After" = 16)) +
  theme_mpa() +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 9),
        legend.text = element_text(size = 9))

# Combine legends
combined_legend <- cowplot::plot_grid(
  cowplot::get_legend(legend_status),
  cowplot::get_legend(legend_period),
  nrow = 1
)

# Arrange panels without duplicate labels (already in titles)
fig2_panels <- ggpubr::ggarrange(
  p2a, p2b, p2c, p2d,
  ncol = 4, nrow = 1,
  align = "hv"
)

# Add shared x-axis label and combine with legend
fig2_with_xlab <- cowplot::ggdraw(fig2_panels) +
  cowplot::draw_label("Year", x = 0.5, y = 0.02, hjust = 0.5, vjust = 0,
                      fontface = "plain", size = 10)

# Add overall title and legend
fig2_titled <- cowplot::plot_grid(
  cowplot::ggdraw() +
    cowplot::draw_label(
      expression(italic("S. purpuratus") * " density at Scorpion SMR"),
      x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5,
      fontface = "bold", size = 11
    ),
  fig2_with_xlab,
  combined_legend,
  ncol = 1,
  rel_heights = c(0.08, 0.82, 0.10)
)

# Conservation Letters: max 170mm double-column width
save_fig(fig2_titled, "fig_02_data_processing", 18, 7)

# =============================================================================
# Figure S1 (Supplemental): Forest plot of effect sizes by MPA and taxa
# Manuscript: Supplemental Figure S1
# =============================================================================

cat("Building Figure S1 (Supplemental): Forest plot...\n")

# Filter SumStats.Final for forest plot
excluded_mpas <- c("Painted Cave SMCA", "San Miguel Island SC",
                   "Arrow Point to Lion Head Point SMCA",
                   "Judith Rk SMR", "Point Conception SMR")

# Debug: Check input data structure
cat("  Figure S1 - Input SumStats.Final MPA column class:", class(SumStats.Final$MPA), "\n")
cat("  Figure S1 - Input MPA levels (if factor):",
    if(is.factor(SumStats.Final$MPA)) paste(head(levels(SumStats.Final$MPA), 5), collapse = ", ") else "N/A", "\n")

fig_s1_data <- SumStats.Final %>%
  dplyr::filter(
    AnalysisType %in% c("pBACIPS", "CI"),
    !(MPA %in% excluded_mpas)
  ) %>%
  dplyr::mutate(
    Taxa = factor(Taxa, levels = taxa_levels),
    Source = factor(Source, levels = source_levels),
    Resp = factor(Resp, levels = c("Den", "Bio")),
    Mean = as.numeric(Mean),
    CI = as.numeric(CI)
  )

# CRITICAL FIX: Convert MPA to character properly, handling factor levels
# If MPA is a factor, as.character() will return the level labels (not codes)
# But if MPA was coerced to factor from numeric, we need the original values
if (is.factor(fig_s1_data$MPA)) {
  fig_s1_data$MPA <- levels(fig_s1_data$MPA)[fig_s1_data$MPA]
}
fig_s1_data$MPA <- as.character(fig_s1_data$MPA)

# Debug: Check MPA values after conversion
cat("  Figure S1 - Unique MPAs after conversion:",
    paste(head(unique(fig_s1_data$MPA), 8), collapse = ", "), "\n")
cat("  Figure S1 - Total rows:", nrow(fig_s1_data), "\n")

# Create shortened MPA names for better readability
fig_s1_data <- fig_s1_data %>%
  dplyr::mutate(
    MPA_short = gsub(" SMCA| SMR| SC| 2003", "", MPA),
    MPA_short = gsub("Anacapa Island", "Anacapa Is.", MPA_short),
    MPA_short = gsub("Santa Barbara Island", "Santa Barbara Is.", MPA_short),
    MPA_short = gsub("San Miguel Island", "San Miguel Is.", MPA_short),
    MPA_short = gsub("Campus Point", "Campus Pt.", MPA_short),
    MPA_short = gsub("Point Vicente", "Pt. Vicente", MPA_short),
    MPA_short = gsub("Harris Point", "Harris Pt.", MPA_short),
    MPA_short = gsub("South Point", "South Pt.", MPA_short),
    MPA_short = gsub("Carrington Pt", "Carrington Pt.", MPA_short)
  )

# Debug: Check MPA_short values
cat("  Figure S1 - MPA_short values:", paste(head(unique(fig_s1_data$MPA_short), 8), collapse = ", "), "\n")

# Order MPAs by mean effect size within each taxa for better visual hierarchy
# Use forcats::fct_reorder to ensure proper factor handling
fig_s1_data <- fig_s1_data %>%
  dplyr::group_by(Taxa) %>%
  dplyr::mutate(
    MPA_order = forcats::fct_reorder(MPA_short, Mean, .fun = mean, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()

fig_s1 <- ggplot(fig_s1_data,
               aes(x = MPA_order, y = Mean, ymin = Mean - CI, ymax = Mean + CI,
                   color = Resp, shape = Source)) +
  # Subtle vertical gridlines at major breaks
  geom_hline(yintercept = seq(-1.5, 1.5, by = 0.5),
             color = "grey92", linewidth = 0.2) +
  # Reference line at zero (no effect) - more prominent
  geom_hline(yintercept = 0, linetype = "solid", color = "grey40",
             linewidth = 0.6) +
  # Confidence intervals as error bars
  geom_errorbarh(aes(xmin = MPA_order, xmax = MPA_order),
                 height = 0, linewidth = 0.6,
                 position = position_dodge(width = 0.7)) +
  geom_pointrange(aes(ymin = Mean - CI, ymax = Mean + CI),
                  size = 0.5, linewidth = 0.6,
                  position = position_dodge(width = 0.7)) +
  geom_point(size = 3.5, position = position_dodge(width = 0.7)) +
  facet_wrap(~ Taxa, ncol = 2, scales = "free_y") +
  scale_color_response(name = "Response",
                       labels = c("Den" = "Density", "Bio" = "Biomass")) +
  scale_shape_source(name = "Source") +
  coord_flip() +
  labs(x = "Marine Protected Area", y = "Effect Size (lnRR)",
       caption = "Error bars = 95% confidence intervals") +
  theme_mpa(base_size = 10) +
  theme(
    strip.text = element_text(face = "italic", size = 10, margin = margin(5, 0, 5, 0)),
    strip.background = element_rect(fill = "grey95", color = "grey60", linewidth = 0.4),
    axis.text.y = element_text(size = 8, color = "grey20"),
    axis.text.x = element_text(size = 9),
    axis.title.y = element_text(size = 10, margin = margin(r = 10)),
    axis.title.x = element_text(size = 10, margin = margin(t = 8)),
    legend.position = "bottom",
    legend.box = "vertical",  # Stack legends vertically to prevent cutoff
    legend.spacing.y = unit(0.3, "cm"),
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 9),
    legend.margin = margin(t = 5, b = 5),
    panel.spacing = unit(1, "lines"),
    panel.grid.major = element_blank(),
    plot.caption = element_text(size = 8, color = "grey50", hjust = 1,
                                 margin = margin(t = 10)),
    plot.margin = margin(10, 15, 10, 10)  # Add right margin for legend
  )

# Conservation Letters: Supplemental figures
save_fig(fig_s1, "fig_s01_forest_plot", 20, 24)

# =============================================================================
# Figure 3: Mean effect sizes by taxa from meta-analysis
# Manuscript: Main Text Figure 3
# =============================================================================

cat("Building Figure 3: Mean effect sizes from meta-analysis...\n")

# Prepare Table2 for plotting
fig3_meta <- Table2 %>%
  dplyr::mutate(
    Taxa = factor(Taxa, levels = taxa_levels),
    Response = factor(Response, levels = c("Density", "Biomass"))
  )

# Prepare individual effect sizes for background points
fig3_individual <- SumStats.Final %>%
  dplyr::filter(
    !(MPA %in% excluded_mpas)
  ) %>%
  dplyr::mutate(
    Taxa = factor(Taxa, levels = taxa_levels),
    Mean = as.numeric(Mean),
    Response = factor(
      ifelse(Resp == "Den", "Density", "Biomass"),
      levels = c("Density", "Biomass")
    )
  )

# Map col_response to long labels for consistency
col_response_fig3 <- c("Density" = unname(col_response["Den"]),
                        "Biomass" = unname(col_response["Bio"]))

# Calculate sample sizes for annotation
fig3_n <- fig3_individual %>%
  dplyr::group_by(Taxa, Response) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop")

fig3 <- ggplot() +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.6) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -0.2, ymax = 0.2),
            fill = "grey93", alpha = 0.5) +
  geom_point(data = fig3_individual, aes(x = Taxa, y = Mean, color = Response),
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7, seed = 42),
             size = 2.5, alpha = 0.5) +
  geom_errorbar(data = fig3_meta, aes(x = Taxa, ymin = CI_lower, ymax = CI_upper, color = Response),
                position = position_dodge(width = 0.7), width = 0.2, linewidth = 1) +
  geom_point(data = fig3_meta, aes(x = Taxa, y = Estimate, color = Response),
             position = position_dodge(width = 0.7), size = 5, shape = 18) +
  scale_color_manual(name = "Response", values = col_response_fig3) +
  coord_cartesian(ylim = c(-6, 3)) +
  labs(x = NULL, y = "Effect Size (lnRR)",
       caption = "Diamonds = meta-analytic means; circles = individual MPA estimates") +
  theme_mpa(base_size = 11) +
  theme(axis.text.x = element_text(face = "italic", size = 10),
        legend.position = "bottom")

save_fig(fig3, "fig_03_mean_effects", 16, 11)

# =============================================================================
# Figure 4: Urchin density vs Kelp biomass scatterplot
# Manuscript: Main Text Figure 4
# Shows relationship between S. purpuratus lnRR and M. pyrifera lnRR
# =============================================================================

cat("Building Figure 4: Urchin vs Kelp scatterplot...\n")

# Prepare data for urchin vs kelp scatterplot
# Need to join urchin (S. purpuratus) and kelp (M. pyrifera) lnRR by MPA and year
fig4_urchin <- All.RR.sub.trans %>%
  dplyr::filter(
    y %in% c("Strongylocentrotus purpuratus", "S. purpuratus", "STRPURAD"),
    resp %in% c("Den", "Density")
  ) %>%
  dplyr::select(CA_MPA_Name_Short, year, source, lnDiff) %>%
  dplyr::rename(lnRR_urchin = lnDiff)

fig4_kelp <- All.RR.sub.trans %>%
  dplyr::filter(
    y %in% c("Macrocystis pyrifera", "M. pyrifera", "MACPYRAD"),
    resp %in% c("Den", "Density")
  ) %>%
  dplyr::select(CA_MPA_Name_Short, year, source, lnDiff) %>%
  dplyr::rename(lnRR_kelp = lnDiff)

fig4_data <- dplyr::inner_join(
  fig4_urchin, fig4_kelp,
  by = c("CA_MPA_Name_Short", "year", "source")
)

# Debug output
cat("  Figure 4 - Urchin rows:", nrow(fig4_urchin), ", Kelp rows:", nrow(fig4_kelp), "\n")
cat("  Figure 4 - Joined data rows:", nrow(fig4_data), "\n")

# Only create figure if we have data
if (nrow(fig4_data) > 0) {
  # Calculate correlation for annotation
  cor_test <- cor.test(fig4_data$lnRR_urchin, fig4_data$lnRR_kelp, method = "pearson")
  cor_label <- sprintf("r = %.2f, p %s",
                       cor_test$estimate,
                       ifelse(cor_test$p.value < 0.001, "< 0.001",
                              sprintf("= %.3f", cor_test$p.value)))

  # Calculate sample size for annotation
  n_points <- nrow(fig4_data)

  fig4 <- ggplot(fig4_data, aes(x = lnRR_urchin, y = lnRR_kelp)) +
    # Reference lines at zero
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    # Quadrant labels for interpretation
    annotate("text", x = -5.5, y = 6, label = "Urchin ↓ Kelp ↑", size = 2.5,
             color = "grey40", fontface = "italic", hjust = 0) +
    annotate("text", x = 4, y = 6, label = "Urchin ↑ Kelp ↑", size = 2.5,
             color = "grey40", fontface = "italic", hjust = 1) +
    annotate("text", x = -5.5, y = -6, label = "Urchin ↓ Kelp ↓", size = 2.5,
             color = "grey40", fontface = "italic", hjust = 0) +
    annotate("text", x = 4, y = -6, label = "Urchin ↑ Kelp ↓", size = 2.5,
             color = "grey40", fontface = "italic", hjust = 1) +
    # Data points - smaller size, more transparent for dense areas
    geom_point(size = 1.8, alpha = 0.35, color = "#2E4E5E",
               position = position_jitter(width = 0.05, height = 0.05, seed = 42)) +
    # Linear regression line with confidence band
    geom_smooth(method = "lm", se = TRUE, color = "#C45B28", fill = "#C45B28",
                linewidth = 1.2, alpha = 0.25) +
    # Correlation and sample size annotation
    annotate("text", x = Inf, y = Inf,
             label = paste0(cor_label, "\nn = ", n_points),
             hjust = 1.1, vjust = 1.3, size = 3.5, fontface = "italic") +
    labs(
      x = expression("Effect Size (lnRR " * italic("S. purpuratus") * " density)"),
      y = expression("Effect Size (lnRR " * italic("M. pyrifera") * " density)")
    ) +
    coord_fixed(ratio = 1, xlim = c(-6, 6), ylim = c(-8, 8)) +
    theme_mpa(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      axis.title = element_text(size = 10),
      plot.margin = margin(10, 10, 10, 10)
    )

  save_fig(fig4, "fig_04_urchin_kelp_scatter", 14, 12)
}

# =============================================================================
# Figure S2 (Supplemental): All taxa log response ratios at example MPAs
# Not in current manuscript - kept as supplemental
# =============================================================================

cat("Building Figure S2 (Supplemental): All taxa time series at example MPAs...\n")

fig_s2_mpas <- c("Naples SMCA", "Scorpion SMR", "Anacapa Island SMR 2003")
fig_s2_starts <- c("Naples SMCA" = 2012, "Scorpion SMR" = 2005,
                  "Anacapa Island SMR 2003" = 2005)

# Look up actual start years from Site table, fall back to hard-coded
for (mpa_name in fig_s2_mpas) {
  site_row <- Site %>%
    dplyr::filter(CA_MPA_Name_Short == mpa_name)
  if (nrow(site_row) > 0 && !is.na(site_row$MPA_Start[1])) {
    fig_s2_starts[mpa_name] <- site_row$MPA_Start[1]
  }
}

# Species to exclude (aggregates)
exclude_species <- c("legal", "sublegal", "All urchins",
                     "Legal", "Sublegal", "all urchins")

fig_s2_data <- All.RR.sub.trans %>%
  dplyr::filter(
    CA_MPA_Name_Short %in% fig_s2_mpas,
    resp == "Den",
    !(tolower(y) %in% tolower(exclude_species))
  ) %>%
  dplyr::mutate(
    y = factor(y, levels = c(
      "Strongylocentrotus purpuratus", "Mesocentrotus franciscanus",
      "Macrocystis pyrifera", "Panulirus interruptus",
      "Semicossyphus pulcher",
      "S. purpuratus", "M. franciscanus", "M. pyrifera",
      "P. interruptus", "S. pulcher"
    )),
    source = factor(source, levels = source_levels),
    CA_MPA_Name_Short = factor(CA_MPA_Name_Short, levels = fig_s2_mpas)
  )

# Build the y (species) column to match col_taxa names if needed
# The col_taxa uses abbreviated names; All.RR.sub.trans may use full or short
# Create a mapping for color matching
species_name_map <- c(
  "Strongylocentrotus purpuratus" = "S. purpuratus",
  "Mesocentrotus franciscanus"    = "M. franciscanus",
  "Macrocystis pyrifera"          = "M. pyrifera",
  "Panulirus interruptus"         = "P. interruptus",
  "Semicossyphus pulcher"         = "S. pulcher",
  "S. purpuratus"                 = "S. purpuratus",
  "M. franciscanus"               = "M. franciscanus",
  "M. pyrifera"                   = "M. pyrifera",
  "P. interruptus"                = "P. interruptus",
  "S. pulcher"                    = "S. pulcher"
)

fig_s2_data <- fig_s2_data %>%
  dplyr::mutate(
    species_short = dplyr::recode(as.character(y), !!!species_name_map),
    species_short = factor(species_short, levels = taxa_levels)
  )

fig_s2_panels <- lapply(fig_s2_mpas, function(mpa_name) {
  d <- dplyr::filter(fig_s2_data, CA_MPA_Name_Short == mpa_name)
  mpa_start_yr <- fig_s2_starts[mpa_name]

  # Clean MPA name for display
  mpa_display <- gsub(" SMCA| SMR| 2003", "", mpa_name)

  ggplot(d, aes(x = year, y = lnDiff,
                color = species_short, shape = source)) +
    # Subtle horizontal gridlines
    geom_hline(yintercept = seq(-1.5, 1.5, by = 0.5),
               color = "grey92", linewidth = 0.2) +
    # Reference line at zero
    geom_hline(yintercept = 0, linetype = "solid", color = "grey50",
               linewidth = 0.5) +
    # MPA implementation vertical line
    geom_vline(xintercept = mpa_start_yr, linetype = "dashed",
               color = "grey40", linewidth = 0.5) +
    # Add MPA label
    annotate("text", x = mpa_start_yr + 0.5, y = max(d$lnDiff, na.rm = TRUE) * 0.9,
             label = "MPA\nStart", hjust = 0, size = 2.5, color = "grey40") +
    # LOESS smoothers for each species (after MPA only)
    geom_smooth(data = dplyr::filter(d, year >= mpa_start_yr),
                aes(group = species_short),
                method = "loess", se = FALSE, span = 0.75,
                linewidth = 0.8, alpha = 0.6) +
    # Data points
    geom_point(size = 2.5, alpha = 0.75) +
    scale_color_taxa(name = "Species") +
    scale_shape_source(name = "Source") +
    scale_x_continuous(breaks = seq(2000, 2020, by = 5)) +
    labs(
      title = mpa_display,
      x = "Year",
      y = "Log response ratio (lnRR)"
    ) +
    theme_mpa(base_size = 9) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0),
      legend.text = element_text(face = "italic", size = 8),
      legend.title = element_text(face = "bold", size = 9),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      panel.grid.major = element_blank()
    )
})

fig_s2 <- ggpubr::ggarrange(
  plotlist = fig_s2_panels,
  ncol = 1, nrow = 3,
  labels = c("(a)", "(b)", "(c)"),
  font.label = list(size = 10, face = "bold"),
  common.legend = TRUE,
  legend = "right"
)

save_fig(fig_s2, "fig_s02_all_taxa_timeseries", 17, 24)

cat("=== All figures saved to:", here::here("plots"), "===\n")
