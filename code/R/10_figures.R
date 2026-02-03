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
ts_colors <- c("Inside" = "#D4A84B", "Outside" = "#CC6633")

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

  # Try multiple species name variants for kelp biomass
  d <- All.Resp.sub %>%
    dplyr::filter(
      taxon_name %in% c("Macrocystis pyrifera", "M. pyrifera", "MACPYRAD"),
      resp %in% c("Bio", "Biomass", "biomass"),
      CA_MPA_Name_Short == mpa_name
    )

  # If no biomass data, try density data
  if (nrow(d) == 0) {
    d <- All.Resp.sub %>%
      dplyr::filter(
        taxon_name %in% c("Macrocystis pyrifera", "M. pyrifera", "MACPYRAD"),
        resp %in% c("Den", "Density", "density"),
        CA_MPA_Name_Short == mpa_name
      )
  }

  # Debug: Print data availability
  cat("    Inset", panel_label, mpa_name, "- rows found:", nrow(d), "\n")

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
  ggplot(d, aes(x = year, y = value, color = status)) +
    geom_vline(xintercept = mpa_start, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_line(aes(group = interaction(status, source)), linewidth = 0.5, alpha = 0.8) +
    geom_point(size = 1.5, alpha = 0.9) +
    scale_color_manual(values = ts_colors, guide = "none") +
    labs(title = paste0(panel_label, " ", short_name),
         x = NULL, y = expression('Biomass (g '~m^{-2}~')')) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
    theme_mpa(base_size = 8) +
    theme(
      plot.title = element_text(size = 8, face = "bold", hjust = 0, margin = margin(0,0,3,0)),
      axis.title = element_text(size = 7),
      axis.title.y = element_text(margin = margin(0,3,0,0)),
      axis.text = element_text(size = 6),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.15, color = "grey85"),
      panel.background = element_rect(fill = "white", color = "grey50", linewidth = 0.4),
      plot.margin = margin(3, 5, 3, 3)
    )
}

insets <- lapply(seq_along(inset_mpas), function(i) {
  create_ts_inset(inset_mpas[i], panel_labels[i])
})

# --- 1.6 Compose final figure ---
fig1_composite <- ggdraw() +
  draw_plot(map_with_scale, x = 0, y = 0, width = 1, height = 1) +
  draw_plot(insets[[1]], x = 0.42, y = 0.72, width = 0.25, height = 0.22) +
  draw_plot(insets[[2]], x = 0.72, y = 0.58, width = 0.25, height = 0.22) +
  draw_plot(insets[[3]], x = 0.02, y = 0.45, width = 0.25, height = 0.22) +
  draw_plot(insets[[4]], x = 0.02, y = 0.15, width = 0.25, height = 0.22) +
  draw_plot(insets[[5]], x = 0.30, y = 0.15, width = 0.25, height = 0.22) +
  draw_plot(insets[[6]], x = 0.58, y = 0.15, width = 0.25, height = 0.22)

# Conservation Letters: max 170mm double-column width
# Reduced from 240mm to 170mm for journal compliance
ggsave(here::here("plots", "fig_01_mpa_map_composite.pdf"), fig1_composite,
       width = 17, height = 14, units = "cm")
ggsave(here::here("plots", "fig_01_mpa_map_composite.png"), fig1_composite,
       width = 17, height = 14, units = "cm", dpi = 300)

# Also save map-only version
ggsave(here::here("plots", "fig_01_map_only.pdf"), map_with_scale,
       width = 20, height = 16, units = "cm")
ggsave(here::here("plots", "fig_01_map_only.png"), map_with_scale,
       width = 20, height = 16, units = "cm", dpi = 300)

cat("  Saved: fig_01_mpa_map_composite.pdf / .png\n")
cat("  Saved: fig_01_map_only.pdf / .png\n")

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
  dplyr::group_by(year, status) %>%
  dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(status = factor(status, levels = c("Inside", "Outside")))

p2a <- ggplot(fig2a_data, aes(x = year, y = value, color = status)) +
  geom_vline(xintercept = scorpion_start, linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2) +
  scale_color_site(name = "Status") +
  labs(x = "Year", y = expression('Density (ind '~m^{-2}~')')) +
  theme_mpa() +
  theme(legend.position = "none")

# Panel (b): Proportion of maximum
fig2b_data <- fig2a_data %>%
  dplyr::group_by(status) %>%
  dplyr::mutate(prop = value / max(value, na.rm = TRUE)) %>%
  dplyr::ungroup()

p2b <- ggplot(fig2b_data, aes(x = year, y = prop, color = status)) +
  geom_vline(xintercept = scorpion_start, linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2) +
  scale_color_site(name = "Status") +
  labs(x = "Year", y = "Proportion of maximum") +
  theme_mpa() +
  theme(legend.position = "none")

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
             linewidth = 0.3) +
  geom_point(size = 2, color = col_taxa["S. purpuratus"]) +
  scale_shape_manual(
    name = "Period",
    values = c("Before" = 1, "After" = 16)
  ) +
  labs(x = "Year", y = "Log response ratio") +
  theme_mpa() +
  theme(legend.position = "none")

# Panel (d): Same with linear trend in After period
fig2d_after <- dplyr::filter(fig2cd_data, BA == "After")

p2d <- ggplot(fig2cd_data, aes(x = year, y = lnDiff, shape = BA)) +
  geom_vline(xintercept = scorpion_start, linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey70",
             linewidth = 0.3) +
  geom_point(size = 2, color = col_taxa["S. purpuratus"]) +
  {
    if (nrow(fig2d_after) >= 3) {
      geom_smooth(
        data = fig2d_after,
        aes(x = year, y = lnDiff),
        method = "lm", se = TRUE,
        color = col_taxa["S. purpuratus"],
        fill = col_taxa["S. purpuratus"],
        alpha = 0.2, linewidth = 0.6,
        inherit.aes = FALSE
      )
    }
  } +
  scale_shape_manual(
    name = "Period",
    values = c("Before" = 1, "After" = 16)
  ) +
  labs(x = "Year", y = "Log response ratio") +
  theme_mpa() +
  theme(legend.position = "none")

fig2 <- ggpubr::ggarrange(
  p2a, p2b, p2c, p2d,
  ncol = 4, nrow = 1,
  labels = c("(a)", "(b)", "(c)", "(d)"),
  font.label = list(size = 10, face = "bold"),
  common.legend = TRUE,
  legend = "bottom"
)

# Conservation Letters: max 170mm double-column width
# Reduced from 280mm to 170mm and adjusted height proportionally
ggsave(here::here("plots", "fig_02_data_processing.pdf"), fig2,
       width = 17, height = 5, units = "cm", device = "pdf")
ggsave(here::here("plots", "fig_02_data_processing.png"), fig2,
       width = 17, height = 5, units = "cm", dpi = 300)
cat("  Saved: fig_02_data_processing.pdf / .png\n")

# =============================================================================
# Figure S1 (Supplemental): Forest plot of effect sizes by MPA and taxa
# Manuscript: Supplemental Figure S1
# =============================================================================

cat("Building Figure S1 (Supplemental): Forest plot...\n")

# Filter SumStats.Final for forest plot
excluded_mpas <- c("Painted Cave SMCA", "San Miguel Island SC",
                   "Arrow Point to Lion Head Point SMCA",
                   "Judith Rk SMR", "Point Conception SMR")

fig_s1_data <- SumStats.Final %>%
  dplyr::filter(
    Type.x %in% c("pBACIPS", "CI"),
    !(MPA %in% excluded_mpas)
  ) %>%
  dplyr::mutate(
    Taxa = factor(Taxa, levels = taxa_levels),
    Source = factor(Source, levels = source_levels),
    Resp = factor(Resp, levels = c("Den", "Bio")),
    Mean = as.numeric(Mean),
    CI = as.numeric(CI),
    # Ensure MPA is character for proper labeling (not numeric factor)
    MPA = as.character(MPA)
  )

# Debug: Check MPA values
cat("  Figure S1 - Unique MPAs:", paste(head(unique(fig_s1_data$MPA), 5), collapse = ", "), "...\n")
cat("  Figure S1 - Total rows:", nrow(fig_s1_data), "\n")

# Create shortened MPA names for better readability
fig_s1_data <- fig_s1_data %>%
  dplyr::mutate(
    MPA_short = gsub(" SMCA| SMR| SC", "", MPA),
    MPA_short = gsub("Anacapa Island", "Anacapa Is.", MPA_short),
    MPA_short = gsub("Santa Barbara Island", "Santa Barbara Is.", MPA_short),
    MPA_short = gsub("San Miguel Island", "San Miguel Is.", MPA_short)
  )

# Order MPAs by mean effect size within each taxa for better visual hierarchy
fig_s1_data <- fig_s1_data %>%
  dplyr::group_by(Taxa) %>%
  dplyr::mutate(MPA_order = reorder(MPA_short, Mean, FUN = mean, na.rm = TRUE)) %>%
  dplyr::ungroup()

fig_s1 <- ggplot(fig_s1_data,
               aes(x = MPA_order, y = Mean, ymin = Mean - CI, ymax = Mean + CI,
                   color = Resp, shape = Source)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
             linewidth = 0.4) +
  geom_pointrange(size = 0.5, linewidth = 0.5,
                  position = position_dodge(width = 0.6)) +
  geom_point(size = 3.5, position = position_dodge(width = 0.6)) +
  facet_wrap(~ Taxa, ncol = 2, scales = "free_y") +
  scale_color_response(name = "Response",
                       labels = c("Den" = "Density", "Bio" = "Biomass")) +
  scale_shape_source(name = "Source") +
  coord_flip() +
  labs(x = NULL, y = "Effect Size (lnRR)") +
  theme_mpa(base_size = 10) +
  theme(
    strip.text = element_text(face = "italic", size = 10, margin = margin(4, 0, 4, 0)),
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 9),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.spacing.x = unit(0.8, "cm"),
    panel.spacing = unit(0.8, "lines")
  )

# Conservation Letters: Supplemental figures can be larger but still reasonable
# Reduced from 280mm to 180mm width for better proportions
ggsave(here::here("plots", "fig_s01_forest_plot.pdf"), fig_s1,
       width = 18, height = 22, units = "cm", device = "pdf")
ggsave(here::here("plots", "fig_s01_forest_plot.png"), fig_s1,
       width = 18, height = 22, units = "cm", dpi = 300)
cat("  Saved: fig_s01_forest_plot.pdf / .png\n")

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

fig3 <- ggplot() +
  # Reference line at zero (no effect)
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
             linewidth = 0.5) +
  # Individual effect sizes as semi-transparent background points
  geom_point(
    data = fig3_individual,
    aes(x = Taxa, y = Mean, color = Response),
    position = position_dodge(width = 0.6),
    size = 2, alpha = 0.25, shape = 16
  ) +
  # Meta-analysis means as large diamonds with error bars
  geom_errorbar(
    data = fig3_meta,
    aes(x = Taxa, ymin = CI_lower, ymax = CI_upper, color = Response),
    position = position_dodge(width = 0.6),
    width = 0.2, linewidth = 0.7
  ) +
  geom_point(
    data = fig3_meta,
    aes(x = Taxa, y = Estimate, color = Response),
    position = position_dodge(width = 0.6),
    size = 5, shape = 18  # Diamond shape for meta-analysis means
  ) +
  scale_color_manual(
    name = "Response",
    values = col_response_fig3,
    labels = c("Density" = "Density", "Biomass" = "Biomass")
  ) +
  scale_y_continuous(
    breaks = seq(-1.5, 1.5, by = 0.5),
    limits = c(-1.5, 1.5)
  ) +
  labs(
    x = NULL,
    y = "Effect Size (lnRR)"
  ) +
  theme_mpa(base_size = 11) +
  theme(
    axis.text.x = element_text(face = "italic", size = 10, margin = margin(t = 5)),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 11, margin = margin(r = 8)),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.5, "cm"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3)
  )

ggsave(here::here("plots", "fig_03_mean_effects.pdf"), fig3,
       width = 15, height = 10, units = "cm", device = "pdf")
ggsave(here::here("plots", "fig_03_mean_effects.png"), fig3,
       width = 15, height = 10, units = "cm", dpi = 300)
cat("  Saved: fig_03_mean_effects.pdf / .png\n")

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

  fig4 <- ggplot(fig4_data, aes(x = lnRR_urchin, y = lnRR_kelp)) +
    # Reference lines at zero
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
    # Data points with slight jitter for overlapping points
    geom_point(size = 2.5, alpha = 0.5, color = "grey30",
               position = position_jitter(width = 0.02, height = 0.02, seed = 42)) +
    # Linear regression line with confidence band
    geom_smooth(method = "lm", se = TRUE, color = "#2E6E7E", fill = "#2E6E7E",
                linewidth = 1, alpha = 0.2) +
    # Correlation annotation
    annotate("text", x = Inf, y = Inf, label = cor_label,
             hjust = 1.1, vjust = 1.5, size = 3.5, fontface = "italic") +
    labs(
      x = expression("Effect Size (lnRR " * italic("S. purpuratus") * " density)"),
      y = expression("Effect Size (lnRR " * italic("M. pyrifera") * " density)")
    ) +
    coord_fixed(ratio = 1) +  # Equal aspect ratio for scatterplot
    theme_mpa(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      axis.title = element_text(size = 10)
    )

  # Save with error handling
  tryCatch({
    ggsave(here::here("plots", "fig_04_urchin_kelp_scatter.pdf"), fig4,
           width = 12, height = 10, units = "cm", device = "pdf")
    ggsave(here::here("plots", "fig_04_urchin_kelp_scatter.png"), fig4,
           width = 12, height = 10, units = "cm", dpi = 300)
    cat("  Saved: fig_04_urchin_kelp_scatter.pdf / .png\n")
  }, error = function(e) {
    cat("  ERROR saving Figure 4:", e$message, "\n")
  })
} else {
  cat("  WARNING: No data available for Figure 4 (urchin-kelp relationship)\n")
  cat("  Check that All.RR.sub.trans contains both S. purpuratus and M. pyrifera density data\n")
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

  ggplot(d, aes(x = year, y = lnDiff,
                color = species_short, shape = source)) +
    geom_vline(xintercept = mpa_start_yr, linetype = "dashed",
               color = "grey50", linewidth = 0.4) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey70",
               linewidth = 0.3) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_taxa(name = "Species") +
    scale_shape_source(name = "Source") +
    labs(
      title = mpa_name,
      x = "Year",
      y = "Log response ratio (lnRR)"
    ) +
    theme_mpa() +
    theme(
      plot.title = element_text(size = 9, face = "bold", hjust = 0),
      legend.text = element_text(face = "italic")
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

# Conservation Letters: Supplemental - reasonable size for multi-panel
ggsave(here::here("plots", "fig_s02_all_taxa_timeseries.pdf"), fig_s2,
       width = 17, height = 24, units = "cm", device = "pdf")
ggsave(here::here("plots", "fig_s02_all_taxa_timeseries.png"), fig_s2,
       width = 17, height = 24, units = "cm", dpi = 300)
cat("  Saved: fig_s02_all_taxa_timeseries.pdf / .png\n")

# =============================================================================
# Done
# =============================================================================

cat("=== All figures saved to:", here::here("plots"), "===\n")
