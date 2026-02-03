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
             size = 4, color = "grey20", stroke = 0.8) +
  geom_text(data = sites_map, aes(x = Lon, y = Lat, label = label),
            nudge_x = 0.12, nudge_y = 0.06, size = 2.8, fontface = "bold", color = "grey10") +
  scale_shape_manual(name = NULL, values = source_shapes) +
  scale_fill_manual(name = NULL, values = source_colors) +
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
    legend.position = c(0.08, 0.12),
    legend.background = element_rect(fill = alpha("white", 0.95), color = "grey50", linewidth = 0.3),
    legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8),
    legend.margin = margin(4, 6, 4, 4), axis.text = element_text(size = 7, color = "grey25"),
    plot.margin = margin(2, 2, 2, 2)
  ) +
  guides(shape = guide_legend(order = 1), fill = "none")

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
  mpa_start <- sites_map %>%
    dplyr::filter(CA_MPA_Name_Short == mpa_name) %>%
    dplyr::pull(MPA_Start)
  if (length(mpa_start) == 0 || is.na(mpa_start[1])) {
    mpa_start <- ifelse(grepl("Point Vicente|Campus", mpa_name), 2012, 2003)
  } else {
    mpa_start <- mpa_start[1]
  }

  d <- All.Resp.sub %>%
    dplyr::filter(taxon_name == "Macrocystis pyrifera", resp == "Bio",
                  CA_MPA_Name_Short == mpa_name)

  if (nrow(d) == 0) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5,
           label = paste(panel_label, gsub(" SM.*", "", mpa_name)), size = 2.5) +
           theme_void() + theme(panel.background = element_rect(fill = "white", color = "grey70")))
  }

  short_name <- gsub(" SMCA| SMR", "", mpa_name)
  short_code <- label_lookup[mpa_name]

  ggplot(d, aes(x = year, y = value, color = status)) +
    geom_vline(xintercept = mpa_start, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_line(aes(group = interaction(status, source)), linewidth = 0.4, alpha = 0.7) +
    geom_point(size = 1.2, alpha = 0.85) +
    scale_color_manual(values = ts_colors, guide = "none") +
    labs(title = paste0(panel_label, " ", short_name, " (", short_code, ")"),
         x = "Years", y = expression('Biomass (g '~m^{-2}~')')) +
    theme_minimal(base_size = 7) +
    theme(
      plot.title = element_text(size = 7, face = "bold", hjust = 0, margin = margin(0,0,2,0)),
      axis.title = element_text(size = 6), axis.title.y = element_text(margin = margin(0,2,0,0)),
      axis.text = element_text(size = 5), panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2),
      panel.background = element_rect(fill = "white", color = "grey60", linewidth = 0.3),
      plot.margin = margin(2, 4, 2, 2)
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

ggsave(here::here("plots", "fig_01_mpa_map_composite.pdf"), fig1_composite,
       width = 24, height = 20, units = "cm")
ggsave(here::here("plots", "fig_01_mpa_map_composite.png"), fig1_composite,
       width = 24, height = 20, units = "cm", dpi = 300)

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

ggsave(here::here("plots", "fig_02_data_processing.pdf"), fig2,
       width = 28, height = 8, units = "cm", device = "pdf")
ggsave(here::here("plots", "fig_02_data_processing.png"), fig2,
       width = 28, height = 8, units = "cm", dpi = 300)
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
    CI = as.numeric(CI)
  )

fig_s1 <- ggplot(fig_s1_data,
               aes(x = MPA, y = Mean, ymin = Mean - CI, ymax = Mean + CI,
                   color = Resp, shape = Source)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60",
             linewidth = 0.3) +
  geom_pointrange(size = 0.4, linewidth = 0.4,
                  position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  facet_wrap(~ Taxa, ncol = 2, scales = "free_y") +
  scale_color_response(name = "Response",
                       labels = c("Den" = "Density", "Bio" = "Biomass")) +
  scale_shape_source(name = "Source") +
  coord_flip() +
  labs(x = "Marine Protected Area", y = "Effect Size (lnRR)") +
  theme_mpa(base_size = 9) +
  theme(
    strip.text = element_text(face = "italic", size = 9),
    axis.text.y = element_text(size = 7),
    legend.position = "bottom",
    legend.box = "horizontal"
  )

ggsave(here::here("plots", "fig_s01_forest_plot.pdf"), fig_s1,
       width = 28, height = 30, units = "cm", device = "pdf")
ggsave(here::here("plots", "fig_s01_forest_plot.png"), fig_s1,
       width = 28, height = 30, units = "cm", dpi = 300)
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
  # Individual effect sizes as semi-transparent background points
  geom_point(
    data = fig3_individual,
    aes(x = Taxa, y = Mean, color = Response),
    position = position_dodge(width = 0.5),
    size = 1.5, alpha = 0.3, shape = 16
  ) +
  # Meta-analysis means as large squares with error bars
  geom_errorbar(
    data = fig3_meta,
    aes(x = Taxa, ymin = CI_lower, ymax = CI_upper, color = Response),
    position = position_dodge(width = 0.5),
    width = 0.15, linewidth = 0.5
  ) +
  geom_point(
    data = fig3_meta,
    aes(x = Taxa, y = Estimate, color = Response),
    position = position_dodge(width = 0.5),
    size = 4, shape = 15
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60",
             linewidth = 0.3) +
  scale_color_manual(name = "Response", values = col_response_fig3) +
  labs(
    x = NULL,
    y = "Effect Size (lnRR)"
  ) +
  theme_mpa(base_size = 10) +
  theme(
    axis.text.x = element_text(face = "italic", size = 9),
    legend.position = "bottom"
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
    y %in% c("Strongylocentrotus purpuratus", "S. purpuratus"),
    resp == "Den"
  ) %>%
  dplyr::select(CA_MPA_Name_Short, year, source, lnDiff) %>%
  dplyr::rename(lnRR_urchin = lnDiff)

fig4_kelp <- All.RR.sub.trans %>%
  dplyr::filter(
    y %in% c("Macrocystis pyrifera", "M. pyrifera"),
    resp == "Den"
  ) %>%
  dplyr::select(CA_MPA_Name_Short, year, source, lnDiff) %>%
  dplyr::rename(lnRR_kelp = lnDiff)

fig4_data <- dplyr::inner_join(
  fig4_urchin, fig4_kelp,
  by = c("CA_MPA_Name_Short", "year", "source")
)

fig4 <- ggplot(fig4_data, aes(x = lnRR_urchin, y = lnRR_kelp)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey70", linewidth = 0.3) +
  geom_point(size = 2, alpha = 0.6, color = "black") +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey80",
              linewidth = 0.8, alpha = 0.3) +
  labs(
    x = expression("Effect Size (lnRR " * italic("S. purpuratus") * ")"),
    y = expression("Effect Size (lnRR " * italic("M. pyrifera") * ")")
  ) +
  theme_mpa(base_size = 10) +
  theme(
    panel.grid.minor = element_blank()
  )

ggsave(here::here("plots", "fig_04_urchin_kelp_scatter.pdf"), fig4,
       width = 12, height = 10, units = "cm", device = "pdf")
ggsave(here::here("plots", "fig_04_urchin_kelp_scatter.png"), fig4,
       width = 12, height = 10, units = "cm", dpi = 300)
cat("  Saved: fig_04_urchin_kelp_scatter.pdf / .png\n")

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

ggsave(here::here("plots", "fig_s02_all_taxa_timeseries.pdf"), fig_s2,
       width = 20, height = 30, units = "cm", device = "pdf")
ggsave(here::here("plots", "fig_s02_all_taxa_timeseries.png"), fig_s2,
       width = 20, height = 30, units = "cm", dpi = 300)
cat("  Saved: fig_s02_all_taxa_timeseries.pdf / .png\n")

# =============================================================================
# Done
# =============================================================================

cat("=== All figures saved to:", here::here("plots"), "===\n")
