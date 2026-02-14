# =============================================================================
# 11_figures.R
# =============================================================================
#
# PURPOSE:
#   Generate publication-quality figures for the Conservation Letters manuscript
#   on MPA effects on kelp forest trophic cascades.
#
# WHAT THIS SCRIPT DOES:
#   Produces publication figures for the manuscript:
#
#   Figure 1: MPA map with bathymetry and time series panels
#     - Ocean bathymetry (depth gradient)
#     - MPA boundaries with monitoring sites
#     - 4 kelp biomass time series panels (Campus Point, Harris Point, South Point, SBI)
#
#   Figure 2: Data processing pipeline example
#     - 4-panel illustration using KFM purple urchin at Scorpion SMR
#     - (a) Raw density time series
#     - (b) Proportion of maximum
#     - (c) Log response ratio
#     - (d) Log response ratio with linear trend fit
#
#   Figure 3: Mean effect sizes from meta-analysis
#     - Summarizes Table 2 graphically
#     - Shows meta-analytic means with 95% CIs (diamonds)
#     - Individual MPA effect sizes shown as background points (circles)
#
#   Figure 4: Urchin density vs Kelp biomass scatterplot
#     - Relationship between MPA effects on urchins and kelp
#     - Includes regression line with confidence band
#     - Quadrant labels indicate trophic cascade dynamics
#
#   Figure 5: Recovery trajectories over time
#     - GAM smooths of annual lnRR for all 5 species
#     - Validates t=11 standardization (vertical reference line)
#     - Shows approximately linear recovery, no saturation
#
#   Figure S1 (Supplemental): Forest plot of effect sizes
#     - Individual effect sizes by MPA and taxa
#     - Color-coded by response type (density vs biomass)
#     - Shape-coded by data source (PISCO, KFM, LTER, Landsat)
#
#   Figure S2 (Supplemental): All taxa time series at example MPAs
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
#   - fig_01_mpa_map.pdf / .png
#   - fig_02_data_processing.pdf / .png
#   - fig_03_mean_effects.pdf / .png
#   - fig_04_trophic_scatter.pdf / .png
#   - fig_05_recovery_curves.pdf / .png
#   - fig_s01_forest_plot.pdf / .png
#   - fig_s02_all_taxa_timeseries.pdf / .png
#   - fig_s07_statistical_transparency.pdf / .png
#   - fig_s08_appendix_*.pdf / .png
#
# DEPENDENCIES:
#   Requires 00-10 scripts to be sourced first
#
# AUTHORS: Emily Donham & Adrian Stier
# PROJECT: CA MPA Kelp Forest pBACIPS Analysis
# =============================================================================

# =============================================================================
# Setup
# =============================================================================

dir.create(here::here("plots"), showWarnings = FALSE)

# --- Selective rendering control ---
# Set RENDER_FIGURES before sourcing to render a subset (e.g., c("fig03", "fig04"))
# Default: render all figures
if (!exists("RENDER_FIGURES", envir = .GlobalEnv)) {
  RENDER_FIGURES <- "all"
}

# Standard factor levels used throughout
taxa_levels  <- c("S. purpuratus", "M. franciscanus", "M. pyrifera",
                   "P. interruptus", "S. pulcher")
source_levels <- c("KFM", "LTER", "PISCO", "Landsat")

# =============================================================================
# Figure dimension constants (Conservation Letters specifications)
# =============================================================================
# Wiley standard widths: single column = 80mm, double column = 180mm
# Figures must be 80-180mm wide; line art at 600 DPI preferred
FIG_WIDTH_SINGLE <- 8    # cm (80mm), for single-column figures
FIG_WIDTH_DOUBLE <- 17   # cm (170mm), for double-column figures
FIG_WIDTH_WIDE   <- 18   # cm (180mm), Wiley max width
FIG_WIDTH_SUPP   <- 17.8 # cm, Conservation Letters max width for supplemental figures

# Figure-specific dimensions (width, height in cm)
# Note: Figure 1 dimensions are defined in analysis/R/fig01_map.R
FIG2_DIMS <- c(w = 17, h = 20)   # Data processing pipeline — 3 vertical panels, double-column
FIG3_DIMS <- c(w = 17, h = 11)   # Mean effects — trophic cascade layout, double-column
FIG4_DIMS <- c(w = 17, h = 17)   # Trophic cascade scatters (2x2 panel)
FIG5_DIMS <- c(w = 17, h = 19)   # Recovery curves (main text, 5 species panels in 3-row trophic layout)
FIG0_DIMS  <- c(w = 8, h = 6)     # Trophic cascade concept diagram (single-column)
FIG_S1_DIMS <- c(w = 17.8, h = 22) # Forest plot (supplemental)
FIG_S2_DIMS <- c(w = 17.8, h = 26) # All taxa time series (supplemental, faceted layout)

# =============================================================================
# Shared variables used across multiple figure sections
# =============================================================================
# MPAs excluded from forest plot and individual effect size displays
excluded_mpas <- c("Painted Cave SMCA", "San Miguel Island SC",
                   "Arrow Point to Lion Head Point SMCA",
                   "Judith Rk SMR", "Point Conception SMR")

# trophic_assignment is defined in 00c_analysis_constants.R (shared with 10_temporal_analysis.R)

# Detect taxa column name in All.RR.sub.trans (used by Fig S3 and S4)
# Deferred until after data validation; set to NULL as placeholder
taxa_col <- NULL

cat("=== Starting figure generation ===\n")

# =============================================================================
# FIX #1: Robust dependency checks with clear error messages
# =============================================================================
# Helper function to check and load required packages
require_pkgs <- function(pkgs, optional = FALSE) {
  missing <- character(0)
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing <- c(missing, pkg)
    }
  }
  if (length(missing) > 0) {
    msg <- paste0("Missing package(s): ", paste(missing, collapse = ", "),
                  "\nInstall with: install.packages(c('",
                  paste(missing, collapse = "', '"), "'))")
    if (optional) {
      warning(msg, call. = FALSE)
      return(FALSE)
    } else {
      stop(msg, call. = FALSE)
    }
  }
  return(TRUE)
}

# Required packages (script will fail without these)
require_pkgs(c("ggplot2", "dplyr", "here", "scales", "forcats", "purrr", "stringr"))

# Optional packages (used if available, with fallbacks)
has_ggpubr  <- require_pkgs("ggpubr", optional = TRUE)
has_cowplot <- require_pkgs("cowplot", optional = TRUE)
has_patchwork <- require_pkgs("patchwork", optional = TRUE)
has_ggrepel <- require_pkgs("ggrepel", optional = TRUE)

# Figure 1 (map) packages - optional but needed for Figure 1
has_fig1_pkgs <- require_pkgs(c("sf", "ggspatial", "marmap", "ggnewscale",
                                 "terra", "tidyterra", "elevatr", "rnaturalearth"),
                               optional = TRUE)
if (!has_fig1_pkgs) {
  cat("  NOTE: Figure 1 packages not available - Figure 1 will be skipped\n")
}

# patchwork is used for multi-panel figure layouts
if (!has_patchwork) {
  stop("Package 'patchwork' is required for figure layouts. Install with: install.packages('patchwork')")
}

# =============================================================================
# FIX #2: Expanded palette/scales validation
# =============================================================================
# Verify color palette is loaded from 00b_color_palette.R
required_palette_objects <- c(
  "col_taxa", "col_response", "col_response_long", "col_site", "theme_mpa"
)
missing_palette <- required_palette_objects[
  !purrr::map_lgl(required_palette_objects, exists, envir = globalenv())
]
if (length(missing_palette) > 0) {
  stop("Color palette objects missing: ", paste(missing_palette, collapse = ", "),
       "\nPlease source 00b_color_palette.R first.")
}

# Verify scale functions exist
required_scale_fns <- c(
  "scale_color_site", "scale_color_response", "scale_shape_source", "scale_color_taxa"
)
missing_scales <- required_scale_fns[
  !purrr::map_lgl(required_scale_fns, exists, envir = globalenv())
]
if (length(missing_scales) > 0) {
  stop("Scale functions missing: ", paste(missing_scales, collapse = ", "),
       "\nPlease ensure 00b_color_palette.R defines all scale_*() helpers.")
}

cat("  Color palette verified: col_taxa, col_response, col_site, theme_mpa loaded\n")
cat("  Scale functions verified: scale_color_site, scale_color_response, scale_shape_source, scale_color_taxa\n")

# =============================================================================
# Input validation: Check required data objects exist and have expected structure
# =============================================================================

# Required data objects from previous scripts
required_objects <- c(
  "All.RR.sub.trans", "All.Resp.sub", "SumStats.Final", "Table2", "Site"
)
missing_objects <- required_objects[
  !purrr::map_lgl(required_objects, exists, envir = globalenv())
]
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

# Set taxa_col now that data is validated
if ("y" %in% names(All.RR.sub.trans)) {
  taxa_col <- "y"
} else if ("Taxa" %in% names(All.RR.sub.trans)) {
  taxa_col <- "Taxa"
} else {
  taxa_col <- names(All.RR.sub.trans)[1]  # fallback
}

# =============================================================================
# Helper function to standardize status values (used in multiple figures)
# =============================================================================
standardize_status <- function(status) {
  status <- as.character(status)
  result <- dplyr::case_when(
    is.na(status) ~ NA_character_,
    tolower(status) %in% c("inside", "mpa", "impact", "i") ~ "Inside",
    tolower(status) %in% c("outside", "reference", "control", "ref", "o", "r") ~ "Outside",
    TRUE ~ NA_character_
  )
  return(result)
}

# =============================================================================
# Consistent MPA implementation annotation style
# =============================================================================
# Standard visual settings for MPA implementation vertical lines
MPA_LINE_COLOR <- "grey40"
MPA_LINE_TYPE <- "dashed"
MPA_LINE_WIDTH <- 0.5
MPA_LABEL_SIZE <- 2.5
MPA_LABEL_COLOR <- "grey40"

# Helper to add MPA implementation line to a ggplot (returns list of geoms)
add_mpa_vline <- function(mpa_year) {
  geom_vline(
    xintercept = mpa_year,
    linetype = MPA_LINE_TYPE,
    color = MPA_LINE_COLOR,
    linewidth = MPA_LINE_WIDTH
  )
}

# =============================================================================
# Reusable theme helpers for consistent legend styling across figures
# =============================================================================

# Theme modifications for bottom-positioned legends (most figures)
theme_legend_bottom <- function(title_size = 10, text_size = 9) {
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(face = "plain", size = title_size),
    legend.text = element_text(size = text_size),
    legend.spacing.x = unit(0.4, "cm")
  )
}

# Theme modifications for right-positioned legends (supplemental time series)
theme_legend_right <- function(title_size = 9, text_size = 8.5, italic = TRUE) {
  theme(
    legend.position = "right",
    legend.title = element_text(face = "plain", size = title_size),
    legend.text = element_text(
      size = text_size,
      face = if (italic) "italic" else "plain"
    )
  )
}

# Small arrow panel used to visually link pipeline panels in Figure 2.
fig2_arrow_panel <- function(label = "\u2192") {
  ggplot() +
    annotate("text", x = 0, y = 0, label = label, size = 10, color = "grey40") +
    xlim(-1, 1) + ylim(-1, 1) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

# =============================================================================
# Figure 0: Trophic Cascade Conceptual Diagram
# =============================================================================
# Schematic showing the three-level trophic cascade hypothesis:
#   MPA Protection -> Predators (+) -> Urchins (-) -> Kelp (+)
# Uses taxa colors from 00b_color_palette.R for species labels.

if (should_render("fig00")) {
cat("Building Figure 0: Trophic Cascade Concept Diagram...\n")

# --- Layout coordinates ---
# Horizontal flow: MPA -> Predators -> Urchins -> Kelp
# All coordinates in [0, 1] normalized space

# Box centers (x positions, evenly spaced)
x_mpa  <- 0.10
x_pred <- 0.35
x_urch <- 0.60
x_kelp <- 0.85

# Box half-widths and half-heights
bw <- 0.10   # half-width
bh <- 0.22   # half-height (tall enough for species text)

# y center for boxes (shifted up to leave room for indirect-effect arc below)
yc <- 0.58

# Arrow y positions (at box vertical center)
arrow_y <- yc

# Species text y positions (relative to yc)
# Row 1: header label at yc + 0.15
# Row 2: first species italic at yc + 0.02
# Row 3: first species common at yc - 0.04
# Row 4: second species italic at yc - 0.10
# Row 5: second species common at yc - 0.16
header_y <- yc + 0.15
sp1_latin_y  <- yc + 0.02
sp1_common_y <- yc - 0.04
sp2_latin_y  <- yc - 0.10
sp2_common_y <- yc - 0.16

# Box colors: use taxa colors with light fill
box_fills   <- c("grey90",
                  col_taxa["S. pulcher"],
                  col_taxa["S. purpuratus"],
                  col_taxa["M. pyrifera"])
box_borders <- c("grey50",
                  col_taxa["S. pulcher"],
                  col_taxa["S. purpuratus"],
                  col_taxa["M. pyrifera"])
box_xc <- c(x_mpa, x_pred, x_urch, x_kelp)

fig0 <- ggplot() +
  # --- Boxes ---
  annotate("rect",
           xmin = box_xc - bw, xmax = box_xc + bw,
           ymin = yc - bh,     ymax = yc + bh,
           fill  = alpha(box_fills, 0.15),
           color = box_borders,
           linewidth = 0.7) +

  # --- Box headers (trophic level labels) ---
  annotate("text", x = x_mpa,  y = header_y, label = "MPA\nProtection",
           fontface = "bold", size = 2.8, color = "grey30", lineheight = 0.85) +
  annotate("text", x = x_pred, y = header_y, label = "Predators",
           fontface = "bold", size = 2.8, color = col_taxa["S. pulcher"]) +
  annotate("text", x = x_urch, y = header_y, label = "Urchins",
           fontface = "bold", size = 2.8, color = col_taxa["S. purpuratus"]) +
  annotate("text", x = x_kelp, y = header_y, label = "Kelp",
           fontface = "bold", size = 2.8, color = col_taxa["M. pyrifera"]) +

  # --- Species names (italic, in taxa colors) ---
  # Predators
  annotate("text", x = x_pred, y = sp1_latin_y,
           label = "S. pulcher",
           fontface = "italic", size = 2.0, color = col_taxa["S. pulcher"]) +
  annotate("text", x = x_pred, y = sp1_common_y,
           label = "(sheephead)", size = 1.7, color = "grey50") +
  annotate("text", x = x_pred, y = sp2_latin_y,
           label = "P. interruptus",
           fontface = "italic", size = 2.0, color = col_taxa["P. interruptus"]) +
  annotate("text", x = x_pred, y = sp2_common_y,
           label = "(lobster)", size = 1.7, color = "grey50") +

  # Urchins
  annotate("text", x = x_urch, y = sp1_latin_y,
           label = "S. purpuratus",
           fontface = "italic", size = 2.0, color = col_taxa["S. purpuratus"]) +
  annotate("text", x = x_urch, y = sp1_common_y,
           label = "(purple urchin)", size = 1.7, color = "grey50") +
  annotate("text", x = x_urch, y = sp2_latin_y,
           label = "M. franciscanus",
           fontface = "italic", size = 2.0, color = col_taxa["M. franciscanus"]) +
  annotate("text", x = x_urch, y = sp2_common_y,
           label = "(red urchin)", size = 1.7, color = "grey50") +

  # Kelp (only one species, centered vertically)
  annotate("text", x = x_kelp, y = (sp1_latin_y + sp2_latin_y) / 2,
           label = "M. pyrifera",
           fontface = "italic", size = 2.0, color = col_taxa["M. pyrifera"]) +
  annotate("text", x = x_kelp, y = (sp1_common_y + sp2_common_y) / 2,
           label = "(giant kelp)", size = 1.7, color = "grey50") +

  # --- Arrows between boxes ---
  # MPA -> Predators (positive)
  annotate("segment",
           x = x_mpa + bw + 0.005, xend = x_pred - bw - 0.005,
           y = arrow_y, yend = arrow_y,
           arrow = arrow(length = unit(1.8, "mm"), type = "closed"),
           linewidth = 0.9, color = "grey30") +
  annotate("text",
           x = (x_mpa + x_pred) / 2, y = arrow_y + 0.045,
           label = "+", fontface = "bold", size = 3.5, color = "#009E73") +

  # Predators -> Urchins (negative)
  annotate("segment",
           x = x_pred + bw + 0.005, xend = x_urch - bw - 0.005,
           y = arrow_y, yend = arrow_y,
           arrow = arrow(length = unit(1.8, "mm"), type = "closed"),
           linewidth = 0.9, color = "grey30") +
  annotate("text",
           x = (x_pred + x_urch) / 2, y = arrow_y + 0.045,
           label = "-", fontface = "bold", size = 4, color = "#D55E00") +

  # Urchins -> Kelp (negative)
  annotate("segment",
           x = x_urch + bw + 0.005, xend = x_kelp - bw - 0.005,
           y = arrow_y, yend = arrow_y,
           arrow = arrow(length = unit(1.8, "mm"), type = "closed"),
           linewidth = 0.9, color = "grey30") +
  annotate("text",
           x = (x_urch + x_kelp) / 2, y = arrow_y + 0.045,
           label = "-", fontface = "bold", size = 4, color = "#D55E00") +

  # --- Indirect effect arc: MPA -> Kelp (+) ---
  # Curved dashed arrow underneath the boxes
  annotate("curve",
           x = x_mpa + bw, xend = x_kelp - bw,
           y = yc - bh - 0.03, yend = yc - bh - 0.03,
           curvature = 0.30,
           arrow = arrow(length = unit(1.8, "mm"), type = "closed"),
           linewidth = 0.6, color = "#009E73", linetype = "dashed") +
  annotate("text",
           x = (x_mpa + x_kelp) / 2, y = yc - bh - 0.12,
           label = "Net indirect effect: +",
           fontface = "bold.italic", size = 2.0, color = "#009E73") +

  # --- Theme and limits ---
  coord_cartesian(xlim = c(-0.03, 0.98), ylim = c(0.18, 0.85), expand = FALSE) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(2, 2, 2, 2, "mm")
  )

save_fig(fig0, "fig_00_trophic_concept", FIG0_DIMS["w"], FIG0_DIMS["h"])
cat("  Figure 0 saved: plots/fig_00_trophic_concept.pdf\n")

} # end fig00


# =============================================================================
# Figure 1: MPA Map with Bathymetry and Time Series Panels
# =============================================================================
# Publication-quality map showing:
#   - Ocean bathymetry (depth gradient)
#   - MPA boundaries
#   - Monitoring sites with paired Inside/Outside design
#   - 4 kelp biomass time series panels

# Load packages needed for figures (outside guards so always available)
library(patchwork)

if (should_render("fig01")) {
cat("Building Figure 1: MPA Map with Bathymetry + Time Series...\n")

if (has_fig1_pkgs) {
  library(sf)
  library(ggspatial)
  library(marmap)
  library(ggnewscale)
  library(terra)
  library(tidyterra)
  library(elevatr)

  # Figure 1 specific constants
  FIG1_PLOT_MARGIN <- ggplot2::margin(2, 2, 2, 2)
  # Reduce overall height to avoid unused vertical whitespace (esp. with square bottom panels).
  FIG1_DIMS <- c(w = 17, h = 14.5)  # cm — taller for balanced map/panel composition

  # --- 1. Define Study Region ---
  BBOX_LONLAT <- c(xmin = -120.75, ymin = 33.28, xmax = -117.65, ymax = 34.50)

  # --- 2. Load Bathymetry Data ---
  bathy_cache <- here::here("data", "cache", "socal_bathy_hires.rds")

  if (file.exists(bathy_cache)) {
    bathy_raw <- readRDS(bathy_cache)
    if (inherits(bathy_raw, "bathy")) {
      bathy_df <- fortify.bathy(bathy_raw)
      names(bathy_df) <- c("lon", "lat", "depth")
    } else {
      bathy_df <- bathy_raw
    }
    cat("  Loaded cached bathymetry\n")
  } else {
    cat("  Downloading bathymetry from NOAA...\n")
    bathy <- getNOAA.bathy(
      lon1 = -121.5, lon2 = -117,
      lat1 = 32.5, lat2 = 35.5,
      resolution = 0.5
    )
    bathy_df <- fortify.bathy(bathy)
    names(bathy_df) <- c("lon", "lat", "depth")
    saveRDS(bathy_df, bathy_cache)
  }

  if (!"lon" %in% names(bathy_df)) names(bathy_df) <- c("lon", "lat", "depth")

	  bathy_ocean <- bathy_df %>%
	    filter(depth < 0,
	           lon >= BBOX_LONLAT["xmin"] - 0.1, lon <= BBOX_LONLAT["xmax"] + 0.1,
	           lat >= BBOX_LONLAT["ymin"] - 0.1, lat <= BBOX_LONLAT["ymax"] + 0.1)

		  # Depth scale should reflect the actual bathymetry range in the plotted layer.
		  # Marmap bathy depths are negative (meters); map legend uses positive meters.
		  depth_max_m <- abs(min(bathy_ocean$depth, na.rm = TRUE))
		  # Use a "nice" step that yields enough labels without overcrowding the compact colorbar.
		  depth_step_m <- if (depth_max_m <= 1000) 500 else 1500
		  depth_breaks_m <- seq(0, floor(depth_max_m / depth_step_m) * depth_step_m, by = depth_step_m)

  # --- 2b. Load/Compute Land Hillshade ---
  hillshade_cache <- here::here("data", "cache", "socal_hillshade.rds")

  if (file.exists(hillshade_cache)) {
    hillshade <- readRDS(hillshade_cache)
    cat("  Loaded cached hillshade\n")
  } else {
    cat("  Downloading elevation and computing hillshade...\n")
    study_bbox <- st_as_sf(
      data.frame(id = 1),
      geometry = st_sfc(st_polygon(list(matrix(c(
        BBOX_LONLAT["xmin"], BBOX_LONLAT["ymin"],
        BBOX_LONLAT["xmax"], BBOX_LONLAT["ymin"],
        BBOX_LONLAT["xmax"], BBOX_LONLAT["ymax"],
        BBOX_LONLAT["xmin"], BBOX_LONLAT["ymax"],
        BBOX_LONLAT["xmin"], BBOX_LONLAT["ymin"]
      ), ncol = 2, byrow = TRUE))), crs = 4326)
    )
    elev_raster <- get_elev_raster(locations = study_bbox, z = 9, clip = "locations")
    dem <- rast(elev_raster)
    dem_land <- classify(dem, cbind(-Inf, 0, NA))
    dem_land <- crop(dem_land, ext(
      BBOX_LONLAT["xmin"] - 0.1, BBOX_LONLAT["xmax"] + 0.1,
      BBOX_LONLAT["ymin"] - 0.1, BBOX_LONLAT["ymax"] + 0.1
    ))
    slope  <- terrain(dem_land, v = "slope",  unit = "radians")
    aspect <- terrain(dem_land, v = "aspect", unit = "radians")
    hill1 <- shade(slope, aspect, angle = 35, direction = 270)
    hill2 <- shade(slope, aspect, angle = 35, direction = 315)
    hill3 <- shade(slope, aspect, angle = 35, direction = 225)
    hill4 <- shade(slope, aspect, angle = 30, direction = 180)
    hillshade <- (hill1 + hill2 + hill3 + hill4) / 4
    names(hillshade) <- "shading"
    saveRDS(hillshade, hillshade_cache)
  }

  terrain_pal <- colorRampPalette(
    c("#5C4A3A", "#8B7D6B", "#C4B9A0", "#E8DCC8", "#F5F0E1", "#FEFCF7")
  )(256)

  # --- 3. Load MPA Boundaries ---
  mpa_shp <- here::here("data", "MPA", "California_Marine_Protected_Areas_[ds582].shp")
  mpa <- st_read(mpa_shp, quiet = TRUE) %>%
    st_transform(4326) %>%
    st_make_valid()

  centroids <- st_coordinates(st_centroid(st_geometry(mpa)))
  in_bbox <- centroids[, 1] >= BBOX_LONLAT["xmin"] & centroids[, 1] <= BBOX_LONLAT["xmax"] &
             centroids[, 2] >= BBOX_LONLAT["ymin"] & centroids[, 2] <= BBOX_LONLAT["ymax"]
  mpa <- mpa[in_bbox, ]

  # Classify MPA types: no-take vs. partial protection (most meaningful distinction)
  mpa <- mpa %>%
    dplyr::mutate(mpa_group = dplyr::case_when(
      Type %in% c("SMR", "FMR", "SMCA (No-Take)")  ~ "No-Take MPA",
      TRUE                                           ~ "Partial Protection"
    ),
    mpa_group = factor(mpa_group, levels = c("No-Take MPA", "Partial Protection")))

  # Distinct fills: darker for no-take, lighter for partial
  mpa_fill_colors <- c(
    "No-Take MPA"         = "#8296A6",
    "Partial Protection"  = "#C5D3DC"
  )

  # --- 4. Load Coastline ---
  coast <- rnaturalearth::ne_states(country = "united states of america", returnclass = "sf") %>%
    filter(name == "California")

  # --- 5. Load Monitoring Sites ---
  site_csv <- here::here("data", "Site_List_All.csv")
  sites_raw <- read.csv(site_csv)

  kfm_sites <- c("Harris Point SMR", "South Point SMR", "Gull Island SMR",
                 "Scorpion SMR", "Santa Barbara Island SMR", "Anacapa Island SMR 2003")
  lter_sites <- c("Campus Point SMCA", "Naples SMCA")
  pisco_sites <- c("Point Vicente SMCA", "Carrington Pt SMR", "Painted Cave SMCA",
                   "Skunk Pt SMR", "Anacapa Island SMCA")

  PANEL_SITES <- c("b" = "Campus Point SMCA", "c" = "Harris Point SMR",
                   "d" = "South Point SMR", "e" = "Santa Barbara Island SMR")
  MPA_YEARS <- c("Campus Point SMCA" = 2012, "Harris Point SMR" = 2003,
                 "South Point SMR" = 2003, "Santa Barbara Island SMR" = 2003)

  sites_base <- sites_raw %>%
    filter(!is.na(Lon) & !is.na(Lat)) %>%
    mutate(program = case_when(
      CA_MPA_Name_Short %in% kfm_sites ~ "KFM",
      CA_MPA_Name_Short %in% lter_sites ~ "LTER",
      CA_MPA_Name_Short %in% pisco_sites ~ "PISCO",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(program)) %>%
    # Jitter Santa Barbara Island marker so the island is visible
    mutate(Lon = ifelse(CA_MPA_Name_Short == "Santa Barbara Island SMR", Lon + 0.08, Lon),
           Lat = ifelse(CA_MPA_Name_Short == "Santa Barbara Island SMR", Lat - 0.06, Lat))

  # Letter labels for ALL monitoring sites on the map.
  # (b)-(e) have time series panels; (f)-(n) are map-only (defined in caption).
  site_abbrev <- c(
    "Campus Point SMCA"       = "b",
    "Harris Point SMR"        = "c",
    "South Point SMR"         = "d",
    "Santa Barbara Island SMR" = "e",
    "Carrington Pt SMR"       = "f",
    "Skunk Pt SMR"            = "g",
    "Painted Cave SMCA"       = "h",
    "Gull Island SMR"         = "i",
    "Scorpion SMR"            = "j",
    "Anacapa Island SMR 2003" = "k",
    "Anacapa Island SMCA"     = "l",
    "Naples SMCA"             = "m",
    "Point Vicente SMCA"      = "n"
  )

  # Manual nudge offsets (lon, lat) for label readability in dense areas
  label_nudge <- list(
    "b" = c( 0.08, -0.05),   # Campus Point — right-below
    "c" = c(-0.14, -0.02),   # Harris Point — left
    "d" = c(-0.14, -0.04),   # South Point — far left-below
    "e" = c(-0.10, -0.07),   # Santa Barbara Is. — left-below
    "f" = c( 0.10,  0.05),   # Carrington Pt — right-above
    "g" = c(-0.12,  0.03),   # Skunk Pt — left-above, scooted
    "h" = c( 0.00,  0.06),   # Painted Cave — above
    "i" = c(-0.12, -0.05),   # Gull Island — further left-below
    "j" = c( 0.10,  0.05),   # Scorpion — right-above
    "k" = c( 0.10,  0.05),   # Anacapa SMR — right-above
    "l" = c(-0.13, -0.04),   # Anacapa SMCA — further left-below
    "m" = c(-0.10,  0.00),   # Naples — left, scooted down
    "n" = c( 0.08, -0.05)    # Point Vicente — right-below
  )

  all_label_df <- tibble::tibble(
    CA_MPA_Name_Short = names(site_abbrev),
    site_abbrev = unname(site_abbrev)
  )
  sites_labels <- sites_base %>%
    inner_join(all_label_df, by = "CA_MPA_Name_Short") %>%
    mutate(
      nudge_lon = sapply(site_abbrev, function(s) label_nudge[[s]][1]),
      nudge_lat = sapply(site_abbrev, function(s) label_nudge[[s]][2]),
      label_x = Lon + nudge_lon,
      label_y = Lat + nudge_lat
    )

  # --- 6. Load Time Series Data ---
  ts_cache <- here::here("data", "cache", "figure_data.rds")
  ts_data <- NULL
  if (file.exists(ts_cache)) {
    fig_data <- readRDS(ts_cache)
    ts_data <- fig_data$All.Resp.sub
  }

  # --- 7. Color Palettes ---
  ocean_colors <- c("#1A4A7A", "#5A8DB8", "#8FBDD6", "#D6E8F0")
  fig1_status_colors <- c(
    "Inside MPA" = if (exists("col_site")) unname(col_site["Inside"]) else "#2A7B8E",
    "Outside MPA" = if (exists("col_site")) unname(col_site["Outside"]) else "#8C7B6A"
  )
  program_shapes <- c("KFM" = 22, "LTER" = 21, "PISCO" = 24)
  fig1_map_colors <- list(
    land = if (exists("col_map")) unname(col_map["land"]) else "#F2EBE1",
    coastline = if (exists("col_map")) unname(col_map["coastline"]) else "#3D3D3D",
    mpa_fill = "#B8C4D0", mpa_border = "#5A6A7A"
  )

	  # --- 8. Build Main Map ---
		  main_map <- ggplot() +
		    # Use positive depth (meters) for legend readability and intuitive left-to-right scale.
		    geom_raster(data = bathy_ocean, aes(x = lon, y = lat, fill = -depth),
		                interpolate = TRUE, alpha = 0.85) +
		    scale_fill_gradientn(
		      colors = rev(ocean_colors),
		      name = "Depth (m)",
		      breaks = depth_breaks_m,
		      labels = depth_breaks_m,
		      limits = c(0, depth_max_m),
		      oob = scales::squish,
		      guide = guide_colorbar(
		        direction = "horizontal",
		        barwidth = unit(2.5, "cm"),
		        barheight = unit(0.3, "cm"),
		        title.position = "left",
		        title.hjust = 0.5,
		        title.theme = element_text(size = 7.5, face = "plain"),
		        frame.colour = "grey60",
		        frame.linewidth = 0.25,
	        ticks = TRUE,
	        ticks.colour = "grey40",
		        ticks.linewidth = 0.35,
		        label.position = "bottom",
		        label.theme = element_text(size = 7, color = "grey25"),
		        order = 3
		      )
		    ) +
	    geom_contour(data = bathy_ocean, aes(x = lon, y = lat, z = depth),
	                 breaks = c(-100, -200, -500, -1000, -2000, -3000), color = "white",
	                 linewidth = 0.20, alpha = 0.25) +
    new_scale_fill() +
    geom_sf(data = mpa, aes(fill = mpa_group), color = fig1_map_colors$mpa_border,
            alpha = 0.45, linewidth = 0.5, inherit.aes = FALSE) +
    scale_fill_manual(name = NULL, values = mpa_fill_colors,
                      guide = guide_legend(order = 2, nrow = 1,
                                           override.aes = list(alpha = 0.6))) +
    new_scale_fill() +
    geom_sf(data = coast, fill = alpha("#F5F0E1", 0.50), color = NA, inherit.aes = FALSE) +
    geom_spatraster(data = hillshade, maxcell = 1e6) +
    scale_fill_gradientn(colors = terrain_pal, na.value = NA, guide = "none") +
    geom_sf(data = coast, fill = NA, color = fig1_map_colors$coastline,
            linewidth = 0.4, inherit.aes = FALSE) +
    # Single point per site (shape by data source)
    geom_point(data = sites_base, aes(x = Lon, y = Lat, shape = program),
               size = 2.8, fill = "grey40", color = "white", stroke = 0.7) +
    scale_shape_manual(name = NULL, values = program_shapes,
                       guide = guide_legend(order = 1, nrow = 1, byrow = TRUE,
                                            override.aes = list(fill = "grey40", size = 2.5))) +
    # Letter labels at nudged positions with leader lines
    geom_segment(data = sites_labels,
                 aes(x = Lon, y = Lat, xend = label_x, yend = label_y),
                 color = "grey50", linewidth = 0.25, alpha = 0.6) +
    geom_label(
      data = sites_labels,
      aes(x = label_x, y = label_y, label = site_abbrev),
      size = 3.5,
      fontface = "bold",
      color = "grey15",
      fill = scales::alpha("white", 0.75),
      label.size = 0,
      label.padding = unit(1.5, "pt"),
      show.legend = FALSE
    ) +
    coord_sf(xlim = c(BBOX_LONLAT["xmin"], BBOX_LONLAT["xmax"]),
             ylim = c(BBOX_LONLAT["ymin"], BBOX_LONLAT["ymax"]), expand = FALSE, crs = 4326) +
    annotation_scale(location = "bl", width_hint = 0.2, pad_x = unit(0.3, "in"),
                     pad_y = unit(0.3, "in"), style = "ticks", text_cex = 0.75, line_width = 0.4) +
    # North arrow — top-right, fancy orienteering style
    annotation_north_arrow(location = "tr", which_north = "true", pad_x = unit(0.15, "in"),
                           pad_y = unit(0.6, "in"), style = north_arrow_fancy_orienteering,
                           height = unit(0.7, "cm"), width = unit(0.7, "cm")) +
    theme_mpa(base_size = 9) +
    labs(tag = "(a)") +
    theme(
      panel.background = element_rect(fill = "#C6DBEF", color = NA),
      panel.border = element_rect(fill = NA, color = "grey35", linewidth = 0.4),
      axis.line.x.bottom = element_blank(),
      axis.line.y.left = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(size = 7.5, color = "grey20"),
      # Legends above the map
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = element_text(size = 7, face = "plain"),
      legend.text = element_text(size = 7),
      legend.key.size = unit(0.3, "cm"),
      legend.spacing.x = unit(4, "mm"),
      legend.spacing.y = unit(0, "mm"),
      legend.key.width = unit(0.5, "cm"),
      legend.key.height = unit(0.3, "cm"),
      legend.background = element_blank(),
      legend.margin = margin(2, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.box.background = element_blank(),
      legend.ticks = element_line(color = "grey40", linewidth = 0.35),
      legend.ticks.length = unit(1.5, "mm"),
      plot.tag = element_text(size = 9, face = "bold"),
      plot.tag.position = "topleft",
      plot.margin = FIG1_PLOT_MARGIN
    )

  # --- 9. California Inset Map ---
  inset_bbox <- st_bbox(BBOX_LONLAT, crs = 4326)
  inset_rect <- st_as_sfc(inset_bbox)

	  inset_map <- ggplot() +
	    geom_sf(data = coast, fill = "grey92", color = "grey60", linewidth = 0.25) +
	    geom_sf(data = inset_rect, fill = NA, color = "#C03A2B", linewidth = 0.7) +
	    coord_sf(xlim = c(-125, -114), ylim = c(32, 42), expand = FALSE) +
	    theme_void() +
	    theme(
		      panel.background = element_rect(fill = scales::alpha("white", 0.70), color = "grey60", linewidth = 0.25),
	      # Keep the inset tight so it can be aligned flush to the main map panel edges.
	      plot.margin = margin(0, 0, 0, 0)
	    )

		  main_map <- main_map +
		    # Align inset to the panel (not the full plot area) so it sits cleanly in the corner.
		    patchwork::inset_element(
		      inset_map,
		      # Requested: top-left corner.
		      # Nudge slightly past the panel bounds to visually "kiss" the top/left panel border.
		      # patchwork clips to the panel, so this avoids tiny anti-aliased gaps.
		      left = -0.01, bottom = 0.825, right = 0.145, top = 1.01,
		      align_to = "panel"
		    )

  # --- 10. Time Series Panels ---
  TS_COLORS <- c("Inside MPA" = "#2A7B8E", "Outside MPA" = "#8C7B6A")

	  build_ts_panel <- function(data, site_name, letter, mpa_year, accent_color = "black") {
	    # Extract MPA protection type before shortening
	    mpa_type <- regmatches(site_name, regexpr("SMCA|SMR", site_name))
	    if (length(mpa_type) == 0) mpa_type <- ""
	    short_name <- gsub(" SMCA| SMR", "", site_name)
	    short_name <- gsub("Point", "Pt.", short_name)
	    short_name <- gsub("Island", "Is.", short_name)
	    short_name <- gsub("Santa Barbara", "S.B.", short_name)
	    # Append MPA type to short name for panel header
	    if (nchar(mpa_type) > 0) short_name <- paste(short_name, mpa_type)
	    site_data <- data %>%
	      filter(CA_MPA_Name_Short == site_name,
	             grepl("pyrifera|Macrocystis", taxon_name, ignore.case = TRUE),
	             resp == "Bio") %>%
      mutate(Status = ifelse(grepl("mpa|inside", status, ignore.case = TRUE),
                             "Inside MPA", "Outside MPA")) %>%
      group_by(year, Status) %>%
      summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")

		    if (nrow(site_data) == 0) {
		      return(
		        ggplot() +
		          theme_void() +
		          labs(title = paste0("(", letter, ") ", short_name))
		      )
		    }

	    # Annotation data for MPA establishment label (all panels)
	    mpa_label_df <- data.frame(x = mpa_year, y = Inf, label = "MPA")

	    p <- ggplot(site_data, aes(x = year, y = mean_value, color = Status,
	                                linetype = Status)) +
	      geom_vline(xintercept = mpa_year, linetype = "dashed", color = "grey45", linewidth = 0.5) +
	      geom_text(
	        data = mpa_label_df, aes(x = x, y = y, label = label),
	        inherit.aes = FALSE, size = 2.8, color = "grey40", hjust = -0.1, vjust = 1.5,
	        fontface = "italic"
	      ) +
	      geom_line(linewidth = 0.6) +
	      geom_point(size = 1.3, alpha = 0.8) +
	      scale_color_manual(values = TS_COLORS, name = NULL) +
	      scale_linetype_manual(values = c("Inside MPA" = "solid", "Outside MPA" = "32"),
	                            guide = "none") +
	      scale_x_continuous(breaks = seq(1990, 2020, by = 10), limits = c(1988, 2024)) +
	      scale_y_continuous(
	        limits = c(0, NA),
	        expand = expansion(mult = c(0, 0.1)),
	        labels = waiver()
	      ) +
			      labs(title = paste0("(", letter, ") ", short_name),
			           x = NULL,
			           y = expression(Kelp~biomass~(g~m^{-2}))) +
	      theme_mpa(base_size = 9) +
		      theme(panel.grid.major.x = element_blank(),
		            panel.grid.major.y = element_blank(),
			            plot.title = element_text(size = 9, face = "bold", hjust = 0,
			                                      color = "black", margin = margin(b = 1)),
			            plot.tag = element_blank(),
		            axis.title = element_text(size = 8, color = "grey20"),
		            axis.title.y = element_text(margin = margin(r = 3)),
		            axis.text = element_text(size = 7.5, color = "grey20"),
		            axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
		            legend.position = "none",
		            plot.margin = margin(2, 3, 1, 3))

	    # Requested: remove y-axis title on panels (b–d) WITHOUT changing widths.
	    # Use transparent text so the grob width is preserved (prevents plots b/c/d looking wider).
	    if (!identical(letter, "b")) {
	      p <- p + theme(axis.title.y = element_text(color = "transparent"))
	    }

	    p
	  }

  # Build panels (accent colors match map labels for perceptual linkage)
  PANEL_ACCENT <- c("b" = "#1B9E77", "c" = "#D95F02", "d" = "#7570B3", "e" = "#E7298A")
  fig1_panels <- list()
  if (!is.null(ts_data)) {
    for (letter in names(PANEL_SITES)) {
      site <- PANEL_SITES[letter]
      fig1_panels[[letter]] <- build_ts_panel(ts_data, site, letter, MPA_YEARS[site],
                                              PANEL_ACCENT[letter])
    }
  }

  # --- 11. Combine Map + Panels (panels in row below map, same width) ---
	  if (length(fig1_panels) == 4) {
	    # Reduce bottom margin of main map
	    main_map <- main_map + theme(plot.margin = margin(1, 1, 0, 1))

    panels_row <- (fig1_panels$b + fig1_panels$c + fig1_panels$d + fig1_panels$e) +
	      plot_layout(nrow = 1, widths = rep(1, 4)) &
	      theme(panel.border = element_blank(),
            axis.line = element_blank(),
            axis.line.x.bottom = element_line(colour = "black", linewidth = 0.3),
            axis.line.y.left = element_line(colour = "black", linewidth = 0.3))

	    fig1 <- main_map / panels_row +
	      plot_layout(heights = c(1.2, 0.8))
	  } else {
	    fig1 <- main_map
	  }

	  fig1 <- fig1 + plot_annotation(theme = theme(plot.margin = margin(1, 8, 2, 2)))

  # --- 12. Save Figure 1 ---
  save_fig(fig1, "fig_01_mpa_map", FIG1_DIMS["w"], FIG1_DIMS["h"], dpi = 600)

} else {
  cat("  Skipping Figure 1 (required packages not available)\n")
}
} # end fig01

# Load additional packages for remaining figures
if (has_cowplot) library(cowplot)
if (has_ggrepel) library(ggrepel)


# =============================================================================
# Figure 2: Conceptual diagram — simulated kelp biomass (3 vertical panels)
# Simulated data modeled on real M. pyrifera patterns at Scorpion SMR
# =============================================================================

if (should_render("fig02")) {
cat("Building Figure 2: Conceptual data processing pipeline (simulated data)...\n")

# --- Simulate kelp biomass data ---
# Modeled on real M. pyrifera patterns at Scorpion SMR (strongest MPA effect)
set.seed(42)  # Reproducibility

mpa_year <- 2003
sim_years <- 1990:2023
n_years <- length(sim_years)

# Noise function: proportional noise with a small additive floor
add_noise <- function(x, sd_frac = 0.15) {
  pmax(0, x + rnorm(length(x), 0, sd_frac * abs(x) + 0.5))
}

# Inside MPA: moderate baseline pre-MPA, then accelerating recovery
inside_before <- 55
inside_trend <- ifelse(sim_years <= mpa_year,
  inside_before,
  55 + pmin((sim_years - mpa_year)^1.6 * 8, 500)  # accelerating recovery
)
inside_bio <- add_noise(inside_trend, sd_frac = 0.20)
inside_bio[sim_years <= mpa_year] <- add_noise(
  rep(inside_before, sum(sim_years <= mpa_year)), sd_frac = 0.20)

# Outside (Reference): matched baseline, gradual decline after MPA
outside_before <- 55
outside_trend <- ifelse(sim_years <= mpa_year,
  outside_before,
  55 - pmin((sim_years - mpa_year) * 1.2, 35)  # decline to ~20
)
outside_bio <- add_noise(outside_trend, sd_frac = 0.20)
outside_bio <- pmax(8, outside_bio)  # floor at 8 to avoid log(0) issues

sim_raw <- data.frame(
  year = rep(sim_years, 2),
  status = factor(rep(c("Inside", "Outside"), each = n_years),
                  levels = c("Inside", "Outside")),
  value = c(inside_bio, outside_bio)
)

# Standardized (proportion of max within each status)
sim_prop <- sim_raw %>%
  dplyr::group_by(status) %>%
  dplyr::mutate(prop = value / max(value, na.rm = TRUE)) %>%
  dplyr::ungroup()

# Log response ratio: ln(Inside / Outside)
sim_lnrr <- data.frame(year = sim_years) %>%
  dplyr::mutate(
    lnRR = log(inside_bio / outside_bio),
    lnRR = ifelse(is.infinite(lnRR) | is.nan(lnRR), NA_real_, lnRR),
    BA = factor(ifelse(year <= mpa_year, "Before", "After"),
                levels = c("Before", "After"))
  ) %>%
  dplyr::filter(!is.na(lnRR))

# X-axis shared settings
sim_x_limits <- c(1988, 2025)
sim_x_breaks <- seq(1990, 2020, by = 10)

# Kelp color
kelp_col <- unname(col_taxa["M. pyrifera"])

# --- Desaturated color variants for supporting panels (a-b) ---
# Lighten/desaturate the teal and taupe for the process-step panels
col_inside_muted  <- adjustcolor(col_site["Inside"], alpha.f = 1,
                                  red.f = 0.90, green.f = 0.95, blue.f = 0.95)
col_outside_muted <- adjustcolor(col_site["Outside"], alpha.f = 1,
                                  red.f = 0.94, green.f = 0.94, blue.f = 0.92)

# --- Theme hierarchy ---
# Supporting panels (a-b): quiet, schematic, process-oriented
fig2_theme_step <- theme_mpa(base_size = 9.5) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 9.5, face = "plain", hjust = 0,
                              margin = margin(b = 1)),
    plot.subtitle = element_text(size = 7.5, color = "grey50",
                                  hjust = 0, margin = margin(b = 4)),
    axis.title.y = element_text(size = 8.5, margin = margin(r = 4)),
    axis.text = element_text(size = 8, color = "grey30"),
    plot.margin = margin(4, 8, 2, 6),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  )

# Causal inference panel (c): dominant, full visual weight
fig2_theme_causal <- theme_mpa(base_size = 10.5) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(face = "plain", size = 8.5),
    legend.text = element_text(size = 8),
    legend.key.size = unit(3.5, "mm"),
    legend.spacing.x = unit(5, "mm"),
    legend.margin = margin(t = 4),
    plot.title = element_text(size = 10.5, face = "plain", hjust = 0,
                              margin = margin(b = 1)),
    plot.subtitle = element_text(size = 7.5, color = "grey50",
                                  hjust = 0, margin = margin(b = 4)),
    axis.title = element_text(size = 9.5, margin = margin(r = 4)),
    axis.text = element_text(size = 8.5, color = "black"),
    plot.margin = margin(6, 8, 4, 6),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  )

# ----------------------------------------------------------
# Panel (a): Raw data — communicate divergence, not noise
# ----------------------------------------------------------
p2a <- ggplot(sim_raw, aes(x = year, y = value, color = status, linetype = status)) +
  add_mpa_vline(mpa_year) +
  annotate("text", x = mpa_year, y = Inf,
           label = "MPA", hjust = -0.15, vjust = 1.5,
           size = 2.5, color = MPA_LABEL_COLOR, fontface = "italic") +
  geom_line(linewidth = 0.5, alpha = 0.6) +
  geom_point(size = 1.1, shape = 21, fill = "white", stroke = 0.4, alpha = 0.35) +
  scale_color_manual(values = c("Inside" = col_inside_muted,
                                 "Outside" = col_outside_muted),
                      name = "Site") +
  scale_linetype_manual(values = c("Inside" = "solid", "Outside" = "32"),
                        guide = "none") +
  labs(title = expression(bold("(a)") ~ "Raw kelp biomass"),
       subtitle = "Example: M. pyrifera biomass at Scorpion SMR (simulated)",
       x = NULL,
       y = expression(Biomass ~ (g ~ m^{-2}))) +
  scale_x_continuous(breaks = sim_x_breaks, limits = sim_x_limits) +
  fig2_theme_step +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# ----------------------------------------------------------
# Panel (b): Standardization — schematic, process-oriented
# ----------------------------------------------------------
p2b <- ggplot(sim_prop, aes(x = year, y = prop, color = status, linetype = status)) +
  add_mpa_vline(mpa_year) +
  annotate("text", x = mpa_year, y = Inf,
           label = "MPA", hjust = -0.15, vjust = 1.5,
           size = 2.5, color = MPA_LABEL_COLOR, fontface = "italic") +
  geom_line(linewidth = 0.5, alpha = 0.6) +
  geom_point(size = 1.1, shape = 21, fill = "white", stroke = 0.4, alpha = 0.35) +
  scale_color_manual(values = c("Inside" = col_inside_muted,
                                 "Outside" = col_outside_muted),
                      name = "Site") +
  scale_linetype_manual(values = c("Inside" = "solid", "Outside" = "32"),
                        guide = "none") +
  labs(title = expression(bold("(b)") ~ "Standardized"),
       subtitle = "Example: M. pyrifera at Scorpion SMR, proportion of site maximum",
       x = NULL,
       y = "Proportion of max") +
  scale_x_continuous(breaks = sim_x_breaks, limits = sim_x_limits) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5)) +
  fig2_theme_step +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# ----------------------------------------------------------
# Panel (c): Causal inference — dominant visual focal point
# ----------------------------------------------------------
sim_lnrr_after <- dplyr::filter(sim_lnrr, BA == "After")

p2c <- ggplot(sim_lnrr, aes(x = year, y = lnRR, shape = BA)) +
  # MPA establishment line (slightly heavier for emphasis)
  geom_vline(xintercept = mpa_year, linetype = "dashed",
             color = "grey30", linewidth = 0.6) +
  annotate("text", x = mpa_year, y = Inf,
           label = "MPA", hjust = -0.15, vjust = 1.5,
           size = 2.8, color = "grey30", fontface = "italic") +
  # Zero-effect reference line (heavier than axes, but restrained)
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey25", linewidth = 0.65) +
  annotate("text", x = sim_x_limits[1] + 1, y = 0,
           label = "No MPA effect", hjust = 0, vjust = -0.6,
           size = 2.5, color = "grey35", fontface = "italic") +
  # Trend line and CI ribbon (drawn before points so points overlay)
  {
    if (nrow(sim_lnrr_after) >= 3) {
      geom_smooth(
        data = sim_lnrr_after,
        aes(x = year, y = lnRR),
        method = "lm", se = TRUE,
        color = kelp_col, fill = kelp_col,
        alpha = 0.18, linewidth = 1.4,
        inherit.aes = FALSE
      )
    }
  } +
  # Data points: After filled and prominent, Before hollow and lighter
  geom_point(aes(alpha = BA), size = 2.3, color = kelp_col) +
  scale_alpha_manual(values = c("Before" = 0.45, "After" = 1.0), guide = "none") +
  scale_shape_manual(
    name = "Period",
    values = c("Before" = 1, "After" = 16),
    guide = guide_legend(override.aes = list(color = kelp_col, size = 2.8,
                                              alpha = c(0.45, 1.0)))
  ) +
  # Causal interpretation annotation (upper right)
  annotate("text",
           x = max(sim_years) - 0.5, y = Inf,
           label = paste0("Positive lnRR after MPA establishment\n",
                          "\u2192 higher kelp biomass inside MPAs"),
           hjust = 1, vjust = 1.5,
           size = 2.8, color = "grey25", lineheight = 1.1,
           fontface = "italic") +
  labs(title = expression(bold("(c)") ~ "Estimated MPA effect"),
       subtitle = "Example: M. pyrifera lnRR at Scorpion SMR with post-MPA trend",
       x = "Year",
       y = "ln(MPA / Reference)") +
  scale_x_continuous(breaks = sim_x_breaks, limits = sim_x_limits) +
  fig2_theme_causal

# ----------------------------------------------------------
# Assemble: 3 vertical panels — (c) is 25% taller
# ----------------------------------------------------------
fig2_final <- (p2a / p2b / p2c) +
  plot_layout(heights = c(1, 1, 1.35), guides = "collect") +
  plot_annotation(
    title = "MPAs increase kelp biomass via post-establishment divergence",
    subtitle = "Illustrative example: M. pyrifera at Scorpion SMR (simulated data modeled on observed patterns)",
    theme = theme(
      plot.title = element_text(size = 10.5, color = "grey15", face = "bold",
                                hjust = 0.5, margin = margin(b = 1)),
      plot.subtitle = element_text(size = 8, color = "grey50", face = "plain",
                                    hjust = 0.5, margin = margin(b = 4)),
      plot.background = element_rect(fill = "white", colour = NA))
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.line.x.bottom = element_line(colour = "black", linewidth = 0.25),
    axis.line.y.left = element_line(colour = "black", linewidth = 0.25),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

save_fig(fig2_final, "fig_02_data_processing", FIG2_DIMS["w"], FIG2_DIMS["h"])
} # end fig02

# =============================================================================
# Figure S1 (Supplemental): Forest plot of effect sizes by MPA and taxa
# Manuscript: Supplemental Figure S1
# =============================================================================

if (should_render("fig_s01")) {
cat("Building Figure S1 (Supplemental): Forest plot...\n")

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

# Convert MPA to character (as.character() on a factor returns level labels, not codes)
fig_s1_data$MPA <- as.character(fig_s1_data$MPA)

cat("  Figure S1 - Unique MPAs after conversion:",
    paste(head(unique(fig_s1_data$MPA), 8), collapse = ", "), "\n")
cat("  Figure S1 - Total rows:", nrow(fig_s1_data), "\n")

# Create shortened MPA names using standardized function
fig_s1_data <- fig_s1_data %>%
  dplyr::mutate(MPA_short = shorten_mpa_name(MPA))

cat("  Figure S1 - MPA_short values:", paste(head(unique(fig_s1_data$MPA_short), 8), collapse = ", "), "\n")

# Order MPAs by mean effect size within each taxa for better visual hierarchy
# Use forcats::fct_reorder to ensure proper factor handling
fig_s1_data <- fig_s1_data %>%
  dplyr::group_by(Taxa) %>%
  dplyr::mutate(
    MPA_order = forcats::fct_reorder(MPA_short, Mean, .fun = mean, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()

# FIX: Calculate consistent x-axis limits across all taxa
x_range_all <- range(c(fig_s1_data$Mean - fig_s1_data$CI,
                       fig_s1_data$Mean + fig_s1_data$CI), na.rm = TRUE)
x_limit <- max(abs(x_range_all)) * 1.1
# RR-labelled breaks: data stays on lnRR, labels show back-transformed RR
# Use sparse log-decade breaks for wide-range forest plots (every 2 decades)
rr_candidates_s1 <- c(0.01, 1, 100)
rr_in_range_s1 <- rr_candidates_s1[abs(log(rr_candidates_s1)) <= x_limit]
x_breaks <- log(rr_in_range_s1)
x_labels <- as.character(rr_in_range_s1)

pd_s1 <- position_dodge(width = 0.55)
fig_s1 <- ggplot(
  fig_s1_data,
  aes(
    y = MPA_order,
    x = Mean,
    xmin = Mean - CI,
    xmax = Mean + CI,
    color = Resp,
    shape = Source
  )
) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey40", linewidth = 0.4) +
  geom_errorbarh(height = 0, linewidth = 0.6, position = pd_s1) +
  geom_point(size = 2.4, position = pd_s1) +
  # Key layout fix: allow each taxa to show only the MPAs that have data
  # (reduces repeated y-label lists and eliminates excess whitespace).
  facet_wrap(~ Taxa, ncol = 2, scales = "free_y", axes = "all_x") +
  scale_color_response(
    name = "Response",
    labels = c("Den" = "Density", "Bio" = "Biomass")
  ) +
  scale_shape_source(name = "Source") +
  scale_x_continuous(limits = c(-x_limit, x_limit),
                     breaks = x_breaks, labels = x_labels) +
  labs(x = "Response ratio (MPA / Reference)", y = NULL) +
  theme_mpa(base_size = 9) +
  theme(
    strip.text = element_text(face = "italic", size = 10, margin = margin(3, 0, 3, 0)),
    strip.background = element_blank(),
    axis.text.y = element_text(size = 8.5, color = "grey15"),
    axis.text.x = element_text(size = 8),
    axis.title.x = element_text(size = 9, margin = margin(t = 6)),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.direction = "horizontal",
    legend.box.spacing = unit(0.25, "cm"),
    legend.spacing.x = unit(0.3, "cm"),
    legend.title = element_text(face = "plain", size = 8.5),
    legend.text = element_text(size = 8),
    legend.margin = margin(t = 4, b = 2),
    legend.key.width = unit(0.55, "cm"),
    legend.box.margin = margin(0, 0, 0, 0),
    panel.spacing.x = unit(0.9, "lines"),
    panel.spacing.y = unit(0.5, "lines"),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    plot.margin = margin(6, 6, 6, 6)
  )

# Conservation Letters: Supplemental figures - slightly wider to fit MPA names
save_fig(fig_s1, "fig_s01_forest_plot", FIG_S1_DIMS["w"], FIG_S1_DIMS["h"])
} # end fig_s01

# =============================================================================
# Figure 3: Mean effect sizes by taxa from meta-analysis
# Manuscript: Main Text Figure 3
# =============================================================================

if (should_render("fig03")) {
cat("Building Figure 3: Mean effect sizes from meta-analysis...\n")

# --- Trophic-level ordering (Predators → Urchins → Kelp) ---
# Use numeric x-axis for precise annotation placement (separators, headers)
fig3_trophic_levels <- c("S. pulcher", "P. interruptus",
                          "S. purpuratus", "M. franciscanus",
                          "M. pyrifera")
fig3_taxa_x <- setNames(seq_along(fig3_trophic_levels), fig3_trophic_levels)

# Common names for x-axis labels (Latin + common)
fig3_common_names <- c(
  "S. pulcher"      = "Sheephead",
  "P. interruptus"  = "Spiny lobster",
  "S. purpuratus"   = "Purple urchin",
  "M. franciscanus" = "Red urchin",
  "M. pyrifera"     = "Giant kelp"
)

# Prepare Table2 for plotting with numeric x positions
fig3_meta <- Table2 %>%
  dplyr::mutate(
    Response = factor(Response, levels = c("Density", "Biomass")),
    x_pos = fig3_taxa_x[as.character(Taxa)]
  ) %>%
  dplyr::filter(!is.na(x_pos))

# Prepare individual effect sizes for background points
fig3_individual <- SumStats.Final %>%
  dplyr::filter(!(MPA %in% excluded_mpas)) %>%
  dplyr::mutate(
    Mean = as.numeric(Mean),
    Response = factor(
      ifelse(Resp == "Den", "Density", "Biomass"),
      levels = c("Density", "Biomass")
    ),
    x_pos = fig3_taxa_x[as.character(Taxa)]
  ) %>%
  dplyr::filter(!is.na(x_pos))

# Calculate sample sizes for annotation
fig3_n <- fig3_individual %>%
  dplyr::group_by(Taxa, Response) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop")

# Add sample sizes to meta data (ensure Taxa types match for join)
fig3_meta <- fig3_meta %>%
  dplyr::mutate(Taxa = as.character(Taxa))
fig3_n <- fig3_n %>%
  dplyr::mutate(Taxa = as.character(Taxa))
fig3_meta <- fig3_meta %>%
  dplyr::left_join(fig3_n, by = c("Taxa", "Response"))

# FDR significance stars: * p<0.05, ** p<0.01, *** p<0.001
fig3_meta <- fig3_meta %>%
  dplyr::mutate(sig_star = dplyr::case_when(
    pval_fdr < 0.001 ~ "***",
    pval_fdr < 0.01  ~ "**",
    pval_fdr < 0.05  ~ "*",
    TRUE ~ ""
  ))

# Dynamic y-limits
fig3_y_min <- min(c(fig3_individual$Mean, fig3_meta$CI_lower), na.rm = TRUE)
fig3_y_max <- max(c(fig3_individual$Mean, fig3_meta$CI_upper), na.rm = TRUE)
fig3_label_y <- fig3_y_min - 0.4

# Dodge offset for Density vs Biomass (manual, since x is numeric)
# Center taxa that have only one response type (e.g., M. pyrifera = Biomass only)
fig3_dodge <- 0.3
fig3_resp_count <- fig3_meta %>%
  dplyr::count(Taxa) %>%
  dplyr::rename(n_resp = n)
fig3_meta <- fig3_meta %>%
  dplyr::left_join(fig3_resp_count, by = "Taxa") %>%
  dplyr::mutate(x_dodge = ifelse(n_resp == 1, x_pos,
                  x_pos + ifelse(Response == "Density", -fig3_dodge/2, fig3_dodge/2)))
fig3_ind_resp_count <- fig3_individual %>%
  dplyr::distinct(Taxa, Response) %>%
  dplyr::count(Taxa) %>%
  dplyr::rename(n_resp = n)
set.seed(42)
fig3_individual <- fig3_individual %>%
  dplyr::left_join(fig3_ind_resp_count, by = "Taxa") %>%
  dplyr::mutate(
    x_dodge = ifelse(n_resp == 1, x_pos,
                x_pos + ifelse(Response == "Density", -fig3_dodge/2, fig3_dodge/2)),
    x_jitter = x_dodge + runif(dplyr::n(), -0.08, 0.08)
  )

# Shape mapping: Density = diamond (18), Biomass = circle (16) for grayscale safety
fig3_shape_map <- c("Density" = 18, "Biomass" = 16)

# Trophic group separator and header positions
# Predators (1,2) | Urchins (3,4) | Kelp (5)
fig3_sep_x <- c(2.5, 4.5)
fig3_header_y <- fig3_y_max + 1.2

fig3 <- ggplot() +
  # Faint vertical separators between trophic groups
  geom_vline(xintercept = fig3_sep_x, color = "grey85", linewidth = 0.3) +
  # Zero-effect reference line
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.5, linetype = "dashed") +
  # Layer 1: Individual MPA effect sizes (subdued background, shape by Response)
  geom_point(data = fig3_individual,
             aes(x = x_jitter, y = Mean, color = Response, shape = Response),
             size = 1.6, alpha = 0.25) +
  # Layer 2: White halo behind focal points
  geom_point(data = fig3_meta,
             aes(x = x_dodge, y = Estimate, shape = Response),
             size = 8.5, color = "white", show.legend = FALSE) +
  # Layer 3: 95% CIs
  geom_errorbar(data = fig3_meta,
                aes(x = x_dodge, ymin = CI_lower, ymax = CI_upper, color = Response),
                width = 0.10, linewidth = 1.15) +
  # Layer 4: Meta-analytic means (shape by Response for grayscale safety)
  geom_point(data = fig3_meta,
             aes(x = x_dodge, y = Estimate, color = Response, shape = Response),
             size = 7.0) +
  # Sample size labels
  geom_text(data = fig3_meta,
            aes(x = x_dodge, y = fig3_label_y, label = paste0("k=", k),
                color = Response),
            size = 2.9, show.legend = FALSE) +
  # FDR significance stars above CI
  geom_text(data = fig3_meta %>% dplyr::filter(sig_star != ""),
            aes(x = x_dodge, y = CI_upper + 0.25, label = sig_star),
            size = 4.0, color = "grey20", show.legend = FALSE) +
  # Trophic group headers (above plot area)
  annotate("text", x = c(1.5, 3.5, 5), y = fig3_header_y,
           label = c("PREDATORS", "URCHINS", "KELP"),
           size = 3.3, color = "grey35", fontface = "bold") +
  scale_color_manual(name = NULL, values = col_response_long,
                     guide = guide_legend(override.aes = list(size = 3.5, alpha = 1))) +
  scale_shape_manual(name = NULL, values = fig3_shape_map) +
  scale_x_continuous(
    breaks = seq_along(fig3_trophic_levels),
    labels = paste0("*", fig3_trophic_levels, "*<br>(", fig3_common_names[fig3_trophic_levels], ")"),
    expand = expansion(mult = 0.08)
  ) +
  # RR-labelled y-axis: data stays on lnRR scale, labels show back-transformed RR
  # Dynamic breaks spanning the full data + CI range
  {
    rr_pool <- c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8, 20, 100, 500)
    vis_lo <- fig3_y_min - 0.9; vis_hi <- fig3_y_max + 0.7
    rr_vis <- rr_pool[log(rr_pool) >= vis_lo & log(rr_pool) <= vis_hi]
    scale_y_continuous(
      breaks = log(rr_vis), labels = as.character(rr_vis),
      name   = "Response ratio (MPA / Reference)"
    )
  } +
  coord_cartesian(ylim = c(fig3_y_min - 0.9, fig3_y_max + 0.7),
                  clip = "off") +
  labs(x = NULL) +
  theme_mpa(base_size = 10) +
  theme(
    axis.text.x = ggtext::element_markdown(size = 8.5, lineheight = 1.2),
    axis.title.y = element_text(size = 10),
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    legend.key.width = unit(0.8, "cm"),
    legend.margin = margin(t = 2),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(18, 10, 6, 6),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

save_fig(fig3, "fig_03_mean_effects", FIG3_DIMS["w"], FIG3_DIMS["h"])
} # end fig03

# =============================================================================
# Figure 4: Trophic cascade scatterplots (4-panel: biomass top, density bottom)
# Uses metafor::rma() meta-regression matching Table 3 methodology from
# 09_meta_analysis.R.  Individual species (not combined predator/urchin averages).
# Top row (biomass):  (a) S. pulcher → S. purpuratus  (b) S. purpuratus → M. pyrifera
# Bottom row (density→biomass): (c) S. pulcher → S. purpuratus  (d) S. purpuratus → M. pyrifera
# =============================================================================

if (should_render("fig04")) {
cat("Building Figure 4: Trophic cascade scatterplots (4-panel)...\n")

# ---------------------------------------------------------------------------
# Build wide-format MPA-level data (same methodology as 09_meta_analysis.R)
# One row per MPA; columns = {Taxa}_{Resp} for effect sizes and SEs
# ---------------------------------------------------------------------------
# NOTE: If an MPA has multiple effect sizes for the same Taxa-Resp (e.g., from
# different data sources), pivot_wider averages them. This is documented and
# intentional — the meta-regression in create_cascade_panel() weights by SE.
fig4_wide_mean <- SumStats.Final %>%
  dplyr::select(Taxa, MPA, Mean, Resp) %>%
  dplyr::mutate(Mean = as.numeric(Mean)) %>%
  tidyr::unite("Taxa_Resp", Taxa, Resp, sep = "_") %>%
  tidyr::pivot_wider(names_from = Taxa_Resp, values_from = Mean, values_fn = mean)

fig4_wide_se <- SumStats.Final %>%
  dplyr::select(Taxa, MPA, SE, Resp) %>%
  dplyr::mutate(SE = as.numeric(SE)) %>%
  tidyr::unite("Taxa_Resp", Taxa, Resp, sep = "_SE_") %>%
  tidyr::pivot_wider(names_from = Taxa_Resp, values_from = SE, values_fn = mean)

fig4_wide <- dplyr::left_join(fig4_wide_mean, fig4_wide_se, by = "MPA")
names(fig4_wide) <- gsub(" ", "_", names(fig4_wide))
names(fig4_wide) <- gsub("\\.", "_", names(fig4_wide))

cat("  Wide-format columns:", paste(names(fig4_wide), collapse = ", "), "\n")

# (Raw annual scatter removed — figure shows only the MPA-level effect sizes
#  that enter the rma() meta-regression, matching the actual analysis.)

# ---------------------------------------------------------------------------
# Helper: create one trophic cascade panel using rma() meta-regression
# (matches 09_meta_analysis.R Table 3 methodology)
# ---------------------------------------------------------------------------
create_cascade_panel <- function(wide_data,
                                  x_col, y_col, se_x_col, se_y_col,
                                  x_lab, y_lab, point_color,
                                  x_short, y_short) {
  # Subset to MPAs with complete data for required columns
  req_cols <- intersect(c(x_col, y_col, se_y_col), names(wide_data))
  if (length(req_cols) < 3) {
    cat("    WARNING: Missing columns for panel. Available:", paste(names(wide_data), collapse = ", "), "\n")
    return(ggplot() + theme_void() +
      annotate("text", x = 0.5, y = 0.5, label = "Missing data columns", size = 4))
  }
  check_cols <- req_cols
  if (se_x_col %in% names(wide_data)) check_cols <- c(check_cols, se_x_col)
  mpa_df <- wide_data[complete.cases(wide_data[, check_cols]), ]
  mpa_df[[x_col]] <- as.numeric(mpa_df[[x_col]])
  mpa_df[[y_col]] <- as.numeric(mpa_df[[y_col]])
  mpa_df[[se_y_col]] <- as.numeric(mpa_df[[se_y_col]])

  cat("    x=", x_col, " y=", y_col, " k=", nrow(mpa_df), "\n")

  if (nrow(mpa_df) < 3) {
    return(ggplot() + theme_void() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste0("Insufficient data (k=", nrow(mpa_df), ")"), size = 4))
  }

  # ---- Meta-regression: rma(yi = y, vi = SE_y^2, mods = ~ x) ----
  rma_df <- data.frame(
    yi    = mpa_df[[y_col]],
    vi    = mpa_df[[se_y_col]]^2,
    x_mod = mpa_df[[x_col]]
  )

  meta_mod <- tryCatch(
    metafor::rma(yi = yi, vi = vi, mods = ~ x_mod, data = rma_df, method = "REML"),
    error = function(e) {
      cat("    rma() error:", conditionMessage(e), "\n")
      NULL
    }
  )

  if (is.null(meta_mod)) {
    return(ggplot() + theme_void() +
      annotate("text", x = 0.5, y = 0.5, label = "Model failed", size = 4))
  }

  # Extract slope statistics
  coef_tbl <- coef(summary(meta_mod))
  slope_est <- coef_tbl[2, "estimate"]
  slope_p   <- coef_tbl[2, "pval"]
  r_sq      <- ifelse(is.null(meta_mod$R2), 0, meta_mod$R2)

  stat_label <- sprintf(
    "\u03b2 = %.2f, p %s\nk = %d MPAs, R\u00b2 = %.0f%%",
    slope_est,
    ifelse(slope_p < 0.001, "< 0.001", sprintf("= %.3f", slope_p)),
    meta_mod$k, r_sq
  )

  # ---- Predicted regression line with 95% CI from rma model ----
  x_seq <- seq(min(rma_df$x_mod), max(rma_df$x_mod), length.out = 100)
  preds <- predict(meta_mod, newmods = x_seq)
  reg_df <- data.frame(x = x_seq, y = preds$pred,
                        ci_lb = preds$ci.lb, ci_ub = preds$ci.ub)

  # ---- Symmetric axis limits centred on zero (tight but balanced) ----
  x_vals <- rma_df$x_mod
  y_vals <- rma_df$yi
  se_x_vals <- if (se_x_col %in% names(mpa_df)) as.numeric(mpa_df[[se_x_col]]) else rep(0, nrow(mpa_df))
  se_y_vals <- as.numeric(mpa_df[[se_y_col]])
  # Symmetric per-axis: extend to the furthest data point + SE, then pad 25%
  max_x <- max(abs(c(x_vals - se_x_vals, x_vals + se_x_vals)), na.rm = TRUE) * 1.25
  max_y <- max(abs(c(y_vals - se_y_vals, y_vals + se_y_vals,
                      reg_df$ci_lb, reg_df$ci_ub)), na.rm = TRUE) * 1.25
  lim_x <- c(-max_x, max_x)
  lim_y <- c(-max_y, max_y)

  # ---- Quadrant labels (ecological interpretation) ----
  # Position labels at 55% of axis limit (centre of each quadrant)
  qx <- max_x * 0.55
  qy <- max_y * 0.55
  quad_labels <- data.frame(
    x = c(-qx, qx, -qx, qx),
    y = c( qy, qy, -qy, -qy),
    label = c(
      paste0(x_short, " \u2193\n", y_short, " \u2191"),   # top-left
      paste0(x_short, " \u2191\n", y_short, " \u2191"),   # top-right
      paste0(x_short, " \u2193\n", y_short, " \u2193"),   # bottom-left
      paste0(x_short, " \u2191\n", y_short, " \u2193")    # bottom-right
    ),
    stringsAsFactors = FALSE
  )

  # ---- Build plot ----
  p <- ggplot() +
    # Quadrant labels (behind everything)
    geom_text(data = quad_labels, aes(x = x, y = y, label = label),
              size = 2.4, color = "grey70", lineheight = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.3)

  # Vertical error bars (y SE)
  p <- p + geom_errorbar(data = mpa_df,
                aes(x = .data[[x_col]],
                    ymin = .data[[y_col]] - .data[[se_y_col]],
                    ymax = .data[[y_col]] + .data[[se_y_col]]),
                width = 0, linewidth = 0.4, color = "grey30", alpha = 0.6)

  # Horizontal error bars (x SE) — only if column exists
  if (se_x_col %in% names(mpa_df)) {
    mpa_df[[se_x_col]] <- as.numeric(mpa_df[[se_x_col]])
    p <- p + geom_errorbarh(data = mpa_df,
                   aes(y = .data[[y_col]],
                       xmin = .data[[x_col]] - .data[[se_x_col]],
                       xmax = .data[[x_col]] + .data[[se_x_col]]),
                   height = 0, linewidth = 0.4, color = "grey30", alpha = 0.6)
  }

  # Solid line for significant (p < 0.05), dashed for non-significant
  reg_linetype <- ifelse(slope_p < 0.05, "solid", "dashed")

  p <- p +
    # Meta-regression CI ribbon + line
    geom_ribbon(data = reg_df, aes(x = x, ymin = ci_lb, ymax = ci_ub),
                fill = point_color, alpha = ifelse(slope_p < 0.05, 0.12, 0.05)) +
    geom_line(data = reg_df, aes(x = x, y = y),
              color = point_color, linewidth = 0.8,
              alpha = ifelse(slope_p < 0.05, 0.7, 0.4),
              linetype = reg_linetype) +
    # MPA-level points (prominent — the actual data)
    geom_point(data = mpa_df, aes(x = .data[[x_col]], y = .data[[y_col]]),
               size = 3.5, alpha = 0.90, color = point_color, shape = 16) +
    # Statistics annotation (bottom-left to avoid quadrant labels)
    annotate("label",
      x = lim_x[1] + 0.02 * diff(lim_x),
      y = lim_y[1] + 0.02 * diff(lim_y),
      label = stat_label,
      hjust = 0, vjust = 0,
      size = 2.8, color = "grey25",
      label.size = 0, fill = scales::alpha("white", 0.75)
    ) +
    labs(x = x_lab, y = y_lab) +
    coord_cartesian(xlim = lim_x, ylim = lim_y, expand = FALSE) +
    theme_mpa(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      plot.margin = margin(4, 6, 2, 6)
    )

  p
}

# ---------------------------------------------------------------------------
# Build all 4 panels (each matches a Table 3 meta-regression model)
# Column naming: after gsub(" ","_") and gsub("\\.","_") on {Taxa}_{Resp}
# ---------------------------------------------------------------------------

# Detect column naming convention (handles both S_pulcher and S__pulcher)
fig4_find_col <- function(pattern) {
  matches <- grep(pattern, names(fig4_wide), value = TRUE)
  if (length(matches) == 1) return(matches)
  if (length(matches) > 1) return(matches[1])
  return(NA_character_)
}

# Panel column mappings
fig4_cols <- list(
  sp_bio    = fig4_find_col("pulcher_Bio$"),
  sp_den    = fig4_find_col("pulcher_Den$"),
  spur_bio  = fig4_find_col("purpuratus_Bio$"),
  spur_den  = fig4_find_col("purpuratus_Den$"),
  mp_bio    = fig4_find_col("pyrifera_Bio$"),
  sp_se_bio = fig4_find_col("pulcher_SE_Bio$"),
  sp_se_den = fig4_find_col("pulcher_SE_Den$"),
  spur_se_bio = fig4_find_col("purpuratus_SE_Bio$"),
  spur_se_den = fig4_find_col("purpuratus_SE_Den$"),
  mp_se_bio   = fig4_find_col("pyrifera_SE_Bio$")
)

cat("  Column mapping:\n")
for (nm in names(fig4_cols)) cat("    ", nm, "=", fig4_cols[[nm]], "\n")

cat("  Building Panel A: S. pulcher biomass -> S. purpuratus biomass...\n")
panel_A <- create_cascade_panel(
  wide_data = fig4_wide,
  x_col = fig4_cols$sp_bio, y_col = fig4_cols$spur_bio,
  se_x_col = fig4_cols$sp_se_bio, se_y_col = fig4_cols$spur_se_bio,
  x_lab = expression(italic("S. pulcher")~"biomass effect (lnRR)"),
  y_lab = expression(italic("S. purpuratus")~"biomass effect (lnRR)"),
  point_color = col_taxa["S. pulcher"],
  x_short = "Sheephead", y_short = "Urchin"
)

cat("  Building Panel B: S. purpuratus biomass -> M. pyrifera biomass...\n")
panel_B <- create_cascade_panel(
  wide_data = fig4_wide,
  x_col = fig4_cols$spur_bio, y_col = fig4_cols$mp_bio,
  se_x_col = fig4_cols$spur_se_bio, se_y_col = fig4_cols$mp_se_bio,
  x_lab = expression(italic("S. purpuratus")~"biomass effect (lnRR)"),
  y_lab = expression(italic("M. pyrifera")~"biomass effect (lnRR)"),
  point_color = col_taxa["S. purpuratus"],
  x_short = "Urchin", y_short = "Kelp"
)

cat("  Building Panel C: S. pulcher density -> S. purpuratus density...\n")
panel_C <- create_cascade_panel(
  wide_data = fig4_wide,
  x_col = fig4_cols$sp_den, y_col = fig4_cols$spur_den,
  se_x_col = fig4_cols$sp_se_den, se_y_col = fig4_cols$spur_se_den,
  x_lab = expression(italic("S. pulcher")~"density effect (lnRR)"),
  y_lab = expression(italic("S. purpuratus")~"density effect (lnRR)"),
  point_color = col_taxa["S. pulcher"],
  x_short = "Sheephead", y_short = "Urchin"
)

cat("  Building Panel D: S. purpuratus density -> M. pyrifera biomass...\n")
panel_D <- create_cascade_panel(
  wide_data = fig4_wide,
  x_col = fig4_cols$spur_den, y_col = fig4_cols$mp_bio,
  se_x_col = fig4_cols$spur_se_den, se_y_col = fig4_cols$mp_se_bio,
  x_lab = expression(italic("S. purpuratus")~"density effect (lnRR)"),
  y_lab = expression(italic("M. pyrifera")~"biomass effect (lnRR)"),
  point_color = col_taxa["S. purpuratus"],
  x_short = "Urchin", y_short = "Kelp"
)

# =============================================================================
# Assemble Figure 4: 4-panel trophic cascade (biomass top, density bottom)
# =============================================================================

cat("Assembling Figure 4 (4-panel trophic cascade)...\n")

fig4 <- (panel_A | panel_B) / (panel_C | panel_D) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")",
                  theme = theme(plot.tag = element_text(face = "plain", size = 10))) +
  plot_layout(heights = c(1, 1)) &  # Equal height for both rows
  theme(
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.line.x.bottom = element_line(colour = "black", linewidth = 0.3),
    axis.line.y.left = element_line(colour = "black", linewidth = 0.3)
  )

save_fig(fig4, "fig_04_trophic_scatter", FIG4_DIMS["w"], FIG4_DIMS["h"])
} # end fig04


# =============================================================================
# Figure 5 (Main Text): Recovery Trajectories Over Time
# =============================================================================
# Validates the t=11 standardization by showing approximately linear GAM
# trajectories for all five species. The linearity confirms that evaluating
# MPAs at a common time horizon captures the same rate-based process.
#
# This is a simplified main-text version of Figure S3 (which shows full
# MPA-level spaghetti). Here we emphasize the population-level smooth and
# mark the t=11 standardization point.

if (should_render("fig05")) {
cat("Building Figure 5 (Main Text): Recovery trajectories over time...\n")

# trophic_assignment defined in 00c_analysis_constants.R
# taxa_col set during input validation above

# --- Full-name to abbreviated-name mapping ---
fig5_full_to_abbrev <- c(
  "Macrocystis pyrifera"          = "M. pyrifera",
  "Mesocentrotus franciscanus"    = "M. franciscanus",
  "Strongylocentrotus purpuratus" = "S. purpuratus",
  "Panulirus interruptus"         = "P. interruptus",
  "Semicossyphus pulcher"         = "S. pulcher"
)

# Species order: predators -> urchins -> kelp (trophic cascade top-down)
fig5_species_order <- c(
  "Panulirus interruptus", "Semicossyphus pulcher",
  "Strongylocentrotus purpuratus", "Mesocentrotus franciscanus",
  "Macrocystis pyrifera"
)

# Species colors keyed by full name
fig5_species_colors <- setNames(
  col_taxa[fig5_full_to_abbrev],
  names(fig5_full_to_abbrev)
)

# Build temporal dataset (After period, time >= 0)
# NOTE: Includes both biomass and density observations. Linear trends are
# shown for consistency with the formal temporal meta-regression (lmer with
# Species x time interaction, random slopes per MPA, random intercepts for
# source and response type). See 10_temporal_analysis.R and Figure S3 for
# non-parametric GAM smooths.
fig5_data <- All.RR.sub.trans %>%
  dplyr::filter(BA == "After", time >= 0, time <= 15) %>%
  dplyr::mutate(
    Species = .data[[taxa_col]],
    time = as.numeric(time)
  ) %>%
  dplyr::filter(Species %in% fig5_species_order)

# Factor species in cascade order for facet ordering
fig5_data$Species <- factor(fig5_data$Species, levels = fig5_species_order)

# Count MPAs per species for panel annotation
fig5_n_mpas <- fig5_data %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(n_mpa = dplyr::n_distinct(CA_MPA_Name_Short), .groups = "drop")

# --- Extract formal lmer slopes from temporal meta-regression CSV ---
# These slopes come from the multilevel model in 10_temporal_analysis.R which
# accounts for MPA-level clustering, source effects, and Bio/Den pooling.
fig5_lmer_csv <- here::here("data", "table_s_temporal_meta_regression.csv")
if (file.exists(fig5_lmer_csv)) {
  fig5_coefs <- read.csv(fig5_lmer_csv, stringsAsFactors = FALSE)
  # Extract reference intercept and slope (P. interruptus)
  ref_intercept <- fig5_coefs$Estimate[fig5_coefs$Term == "(Intercept)"]
  ref_slope <- fig5_coefs$Estimate[fig5_coefs$Term == "time"]
  # Build per-species slopes and intercepts
  fig5_slopes <- setNames(numeric(5), fig5_species_order)
  fig5_intercepts <- setNames(numeric(5), fig5_species_order)
  fig5_slopes["Panulirus interruptus"] <- ref_slope
  fig5_intercepts["Panulirus interruptus"] <- ref_intercept
  for (sp in fig5_species_order[-1]) {
    # Match species interaction term (e.g., "SpeciesSemicossyphus pulcher:time")
    sp_int <- paste0("Species", sp, ":time")
    sp_main <- paste0("Species", sp)
    int_row <- fig5_coefs[fig5_coefs$Term == sp_int, ]
    main_row <- fig5_coefs[fig5_coefs$Term == sp_main, ]
    fig5_slopes[sp] <- ref_slope + ifelse(nrow(int_row) > 0, int_row$Estimate, 0)
    fig5_intercepts[sp] <- ref_intercept + ifelse(nrow(main_row) > 0, main_row$Estimate, 0)
  }
  cat("  Loaded formal lmer slopes from table_s_temporal_meta_regression.csv\n")
} else {
  # Fallback: will use OLS slopes from geom_smooth only
  fig5_slopes <- NULL
  fig5_intercepts <- NULL
  cat("  NOTE: Temporal meta-regression CSV not found; using OLS slopes\n")
}

# --- Compute shared y-axis limits within trophic pairs ---
# Predators share one range, herbivores share another, kelp gets its own
fig5_ylims_raw <- fig5_data %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(
    q_lo = quantile(lnDiff, 0.025, na.rm = TRUE),
    q_hi = quantile(lnDiff, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# Trophic group assignments
fig5_trophic <- c(
  "Panulirus interruptus"         = "predator",
  "Semicossyphus pulcher"         = "predator",
  "Strongylocentrotus purpuratus" = "herbivore",
  "Mesocentrotus franciscanus"    = "herbivore",
  "Macrocystis pyrifera"          = "producer"
)

fig5_ylims_raw$group <- fig5_trophic[as.character(fig5_ylims_raw$Species)]

# Compute shared limits per trophic group
fig5_group_lims <- fig5_ylims_raw %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(
    g_lo = min(q_lo),
    g_hi = max(q_hi),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    pad = pmax((g_hi - g_lo) * 0.15, 0.5),
    ymin = g_lo - pad,
    ymax = g_hi + pad
  )

# Cap producer (kelp) panel y-range to reduce compression from outlier
# spaghetti lines. This clips extreme individual MPA trajectories while
# keeping the trend line and bulk of data fully visible.
fig5_group_lims <- fig5_group_lims %>%
  dplyr::mutate(
    ymin = ifelse(group == "producer", pmax(ymin, -6), ymin),
    ymax = ifelse(group == "producer", pmin(ymax, 8), ymax)
  )

# Merge back to species level
fig5_ylims <- fig5_ylims_raw %>%
  dplyr::left_join(fig5_group_lims, by = "group")

# --- Helper: build one species panel ---
make_fig5_panel <- function(sp, tag_label, show_xlab = FALSE) {
  sp_data <- fig5_data[fig5_data$Species == sp, ]
  sp_color <- fig5_species_colors[[sp]]
  sp_abbrev <- fig5_full_to_abbrev[[sp]]
  n_mpa <- fig5_n_mpas$n_mpa[fig5_n_mpas$Species == sp]
  ylims <- fig5_ylims[fig5_ylims$Species == sp, ]
  sp_slope <- if (!is.null(fig5_slopes)) fig5_slopes[[sp]] else NA
  sp_int <- if (!is.null(fig5_intercepts)) fig5_intercepts[[sp]] else NA

  # Formal lmer prediction at t=11 for point marker
  lnRR_at_11 <- if (!is.na(sp_int)) sp_int + sp_slope * 11 else NA

  p <- ggplot(sp_data, aes(x = time, y = lnDiff)) +
    # Individual MPA trajectories — species-tinted, subtle behind trend
    geom_line(aes(group = CA_MPA_Name_Short),
              color = sp_color, alpha = 0.12, linewidth = 0.25) +
    # Reference line at zero effect (solid, subtle)
    geom_hline(yintercept = 0, linetype = "solid", color = "grey45",
               linewidth = 0.3) +
    # Vertical reference at t=11 (standardized comparison point)
    geom_vline(xintercept = 11, linetype = "dotted", color = "grey50",
               linewidth = 0.4) +
    annotate("text", x = 11.3, y = Inf, label = "t = 11",
             hjust = 0, vjust = 1.3, size = 2.3, color = "grey40") +
    # Population-level linear trend with 95% CI (OLS for visualization)
    geom_smooth(color = sp_color, fill = sp_color,
                method = "lm", formula = y ~ x,
                linewidth = 1.5, alpha = 0.20,
                se = TRUE, level = 0.95) +
    # MPA count annotation (top-left below title)
    annotate("text", x = 0.3, y = Inf,
             label = paste0("n = ", n_mpa, " MPAs"),
             hjust = 0, vjust = 1.5, size = 2.2, color = "grey40")

  # Add formal lmer slope annotation if available
  if (!is.na(sp_slope)) {
    slope_label <- paste0(
      ifelse(sp_slope >= 0, "+", "\u2212"),
      sprintf("%.2f", abs(sp_slope)),
      " lnRR yr\u207b\u00b9"
    )
    p <- p + annotate("text", x = 0.3, y = Inf,
                       label = slope_label,
                       hjust = 0, vjust = 3.0, size = 2.2, color = "grey30")
  }

  # Add t=11 point marker from lmer prediction
  if (!is.na(lnRR_at_11)) {
    p <- p + annotate("point", x = 11, y = lnRR_at_11,
                       color = sp_color, fill = "white",
                       shape = 21, size = 2.5, stroke = 0.8)
  }

  p <- p +
    # Shared y-axis within trophic group (clips outlier spaghetti)
    coord_cartesian(ylim = c(ylims$ymin, ylims$ymax),
                    clip = "off") +
    scale_x_continuous(breaks = seq(0, 15, by = 5),
                       limits = c(0, 16), expand = c(0.02, 0)) +
    scale_y_continuous(expand = expansion(mult = 0.02)) +
    labs(
      title = sp_abbrev,
      x = if (show_xlab) "Years since MPA implementation" else NULL,
      y = NULL
    ) +
    theme_mpa(base_size = 9) +
    theme(
      plot.title = element_text(face = "italic", size = 9, hjust = 0,
                                margin = margin(b = 1)),
      plot.tag = element_text(size = 10, face = "bold"),
      plot.tag.position = c(0, 1),
      plot.margin = margin(3, 6, 3, 4, "pt")
    ) +
    labs(tag = tag_label)

  p
}

# --- Build individual panels in trophic order ---
# Rows: Predators | Herbivores | Kelp (only bottom row gets x-axis label)
p_a <- make_fig5_panel("Panulirus interruptus", "a")
p_b <- make_fig5_panel("Semicossyphus pulcher", "b")
p_c <- make_fig5_panel("Strongylocentrotus purpuratus", "c")
p_d <- make_fig5_panel("Mesocentrotus franciscanus", "d")
p_e <- make_fig5_panel("Macrocystis pyrifera", "e", show_xlab = TRUE)

# --- Assemble with patchwork: trophic cascade layout ---
# Row 1: Predators (lobster, sheephead)
# Row 2: Herbivores (purple urchin, red urchin)
# Row 3: Primary producer (kelp, centered — spacers match half-row width)
fig5 <- (p_a | p_b) /
        (p_c | p_d) /
        (plot_spacer() + p_e + plot_spacer() +
           plot_layout(widths = c(0.5, 1, 0.5))) +
  plot_layout(heights = c(1, 1, 1.08))

# Add shared y-axis label via cowplot draw_label
fig5_final <- cowplot::ggdraw() +
  cowplot::draw_label("Log response ratio (lnRR)", x = 0.015, y = 0.52,
                      angle = 90, size = 9, fontfamily = "Helvetica") +
  cowplot::draw_plot(fig5, x = 0.04, y = 0, width = 0.96, height = 1)

save_fig(fig5_final, "fig_05_recovery_curves", FIG5_DIMS["w"], FIG5_DIMS["h"])
cat("  Figure 5 saved: plots/fig_05_recovery_curves.pdf\n")

} # end fig05


# =============================================================================
# Figure S2 (Supplemental): All taxa log response ratios at example MPAs
# Not in current manuscript - kept as supplemental
# =============================================================================

if (should_render("fig_s02")) {
cat("Building Figure S2 (Supplemental): All taxa time series at example MPAs...\n")

fig_s2_mpas <- c("Naples SMCA", "Scorpion SMR", "Anacapa Island SMR 2003")

# Look up start years from Site table with fallback defaults
# Using tidyverse pattern: create defaults, then update from Site data
fig_s2_starts <- tibble::tibble(
  CA_MPA_Name_Short = fig_s2_mpas,
  MPA_Start = c(2012, 2005, 2005)  # fallback defaults
) %>%
  dplyr::rows_update(
    Site %>%
      dplyr::filter(CA_MPA_Name_Short %in% fig_s2_mpas) %>%
      dplyr::select(CA_MPA_Name_Short, MPA_Start) %>%
      dplyr::filter(!is.na(MPA_Start)),
    by = "CA_MPA_Name_Short"
  ) %>%
  tibble::deframe()  # Convert to named vector

# Species to exclude (aggregates) - lowercase for case-insensitive matching
exclude_species <- c("legal", "sublegal", "all urchins")

fig_s2_data <- All.RR.sub.trans %>%
  dplyr::filter(
    CA_MPA_Name_Short %in% fig_s2_mpas,
    resp == "Den",
    !stringr::str_to_lower(y) %in% exclude_species
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

# Build species_short column to match col_taxa abbreviated names
# Uses forcats::fct_recode for tidyverse-consistent factor recoding
fig_s2_data <- fig_s2_data %>%
  dplyr::mutate(
    species_short = forcats::fct_recode(
      factor(y),
      "S. purpuratus"   = "Strongylocentrotus purpuratus",
      "M. franciscanus" = "Mesocentrotus franciscanus",
      "M. pyrifera"     = "Macrocystis pyrifera",
      "P. interruptus"  = "Panulirus interruptus",
      "S. pulcher"      = "Semicossyphus pulcher"
    ) %>%
      forcats::fct_relevel(taxa_levels)
  )

# Create short MPA display names for facet column labels
fig_s2_data <- fig_s2_data %>%
  dplyr::mutate(
    site_display = shorten_mpa_name(as.character(CA_MPA_Name_Short)),
    site_display = factor(site_display,
                          levels = shorten_mpa_name(fig_s2_mpas))
  )

# Build a per-site MPA start year lookup for the vertical line geom.
# Each facet column shares the same MPA start year for its site.
fig_s2_vline_data <- fig_s2_data %>%
  dplyr::distinct(site_display, CA_MPA_Name_Short) %>%
  dplyr::mutate(
    mpa_start_yr = fig_s2_starts[as.character(CA_MPA_Name_Short)]
  )

# Join MPA start years into main data for filtering LOESS to post-MPA period
fig_s2_data_with_start <- fig_s2_data %>%
  dplyr::left_join(
    tibble::enframe(fig_s2_starts,
                    name = "CA_MPA_Name_Short",
                    value = "mpa_start_yr"),
    by = "CA_MPA_Name_Short"
  )

# Faceted figure: species_short (rows) x site_display (columns)
# Each facet shows one species at one site, eliminating overplotting.
# Data source shape encoding removed (already shown in Figure S1 forest plot).
fig_s2 <- ggplot(fig_s2_data, aes(x = year, y = lnDiff,
                                   color = species_short)) +
  # Reference line at zero — prominent for orientation
  geom_hline(yintercept = 0, linetype = "solid", color = "grey30",
             linewidth = 0.7) +
  # MPA implementation vertical line (per-site, drawn via vline data)
  geom_vline(data = fig_s2_vline_data,
             aes(xintercept = mpa_start_yr),
             linetype = MPA_LINE_TYPE,
             color = MPA_LINE_COLOR,
             linewidth = MPA_LINE_WIDTH,
             inherit.aes = FALSE) +
  # LOESS smoothers (post-MPA only)
  geom_smooth(data = fig_s2_data_with_start %>%
                dplyr::filter(year >= mpa_start_yr),
              method = "loess", se = FALSE, span = 0.75,
              linewidth = 0.7, alpha = 0.55) +
  # Data points — consistent filled circle (shape 16), ~30% smaller than old 1.8
  geom_point(size = 1.25, alpha = 0.6, shape = 16, na.rm = TRUE) +
  scale_color_taxa(name = "Species") +
  facet_grid(species_short ~ site_display, scales = "free_x") +
  scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 1)) +
  labs(
    x = "Year",
    y = "Log response ratio (lnRR)"
  ) +
  theme_mpa(base_size = 9) +
  theme(
    # Legend not needed: species identity is encoded in row strip labels + color
    legend.position = "none",
    strip.text.x = element_text(size = 9, face = "bold",
                                margin = margin(b = 3, t = 3)),
    strip.text.y = element_text(size = 8, face = "italic",
                                margin = margin(l = 3, r = 3)),
    axis.title = element_text(size = 8.5),
    axis.text = element_text(size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(fill = NA, colour = "grey80", linewidth = 0.4),
    panel.spacing = unit(3, "mm"),
    plot.margin = margin(6, 6, 6, 6)
  )

save_fig(fig_s2, "fig_s02_all_taxa_timeseries", FIG_S2_DIMS["w"], FIG_S2_DIMS["h"])
} # end fig_s02

# (Old Figures S3 and S4 removed — superseded by species-level versions
#  in 10_temporal_analysis.R: fig_s03 through fig_s06)


# =============================================================================
# Figure S7: Model Selection & Heterogeneity
# =============================================================================
# Statistical transparency: variance components and model selection

if (should_render("fig_s07")) {
cat("\n--- Figure S7: Statistical Transparency ---\n")

FIG_S7_DIMS <- c(w = 17, h = 9)  # Conservation Letters max width (reduced height to tighten layout)

# Panel A: Model selection distribution by taxa
model_dist <- SumStats.Final %>%
  group_by(Taxa, Model) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Taxa) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

model_dist$Taxa <- factor(model_dist$Taxa, levels = taxa_levels)

# Model colors
col_model <- c(
  "Step" = "#4A90A4",
  "Linear" = "#7B68A6",
  "Asymptotic" = "#D4933B",
  "Sigmoid" = "#B85A4C",
  "Mean" = "#8C7B6A"
)

# Pre-compute cumulative positions for Sigmoid labels
sigmoid_labels <- model_dist %>%
  arrange(Taxa, desc(Model)) %>%
  group_by(Taxa) %>%
  mutate(ymax = cumsum(prop), ymin = ymax - prop) %>%
  ungroup() %>%
  filter(Model == "Sigmoid", prop > 0) %>%
  mutate(ymid = (ymin + ymax) / 2,
         label = paste0(round(prop * 100), "%"))

panel_S7A <- ggplot(model_dist, aes(x = Taxa, y = prop, fill = Model)) +
  geom_col(position = "stack", width = 0.7, color = "white", linewidth = 0.3) +
  # Label small Sigmoid slivers so readers know they are intentional

  {if (nrow(sigmoid_labels) > 0) {
    geom_text(
      data = sigmoid_labels,
      aes(x = Taxa, y = ymid, label = label),
      inherit.aes = FALSE,
      size = 2.5, hjust = -0.3, color = "grey30"
    )
  }} +
  scale_fill_manual(
    values = col_model,
    name = "Effect size method",
    guide = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  labs(
    x = NULL,
    y = "Proportion of MPAs"
  ) +
  theme_mpa(base_size = 10) +
  theme(
    axis.text.x = element_text(face = "italic", angle = 25, hjust = 1),
    legend.position = "bottom",
    legend.direction = "horizontal",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Panel B: Variance components (if meta-analysis objects exist)
if (exists("meta_biomass") && exists("meta_density")) {

  var_comp <- data.frame(
    Response = rep(c("Biomass", "Density"), each = 2),
    Component = rep(c("MPA", "Source"), 2),
    tau2 = c(
      meta_biomass$sigma2[1], meta_biomass$sigma2[2],
      meta_density$sigma2[1], meta_density$sigma2[2]
    )
  )
  var_comp$Response <- factor(var_comp$Response, levels = c("Density", "Biomass"))

  panel_S7B <- ggplot(var_comp, aes(x = tau2, y = Response, fill = Component)) +
    geom_col(position = position_dodge(width = 0.6), width = 0.6, color = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2f", tau2)), position = position_dodge(width = 0.6), hjust = -0.15, size = 2.8, color = "grey20") +
    scale_fill_manual(
      values = c("MPA" = "#2A7B8E", "Source" = "#D4933B"),
      name = "Variance component",
      guide = guide_legend(nrow = 1, byrow = TRUE)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(
      x = expression(tau^2),
      y = NULL
    ) +
    theme_mpa(base_size = 10) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  fig_s7 <- panel_S7A / panel_S7B +
    # Put the collected guides into an explicit guide row to avoid the large
    # whitespace patchwork can leave when legends are collected implicitly.
    patchwork::guide_area() +
    plot_layout(heights = c(1, 0.6, 0.10), guides = "collect") +
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      # Two legends: keep them side-by-side to reduce overall figure height.
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(1, "mm"),
      legend.spacing.y = unit(1, "mm"),
      legend.spacing.x = unit(4, "mm"),
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.line.x.bottom = element_line(colour = "black", linewidth = 0.3),
      axis.line.y.left = element_line(colour = "black", linewidth = 0.3)
    )

} else {
  fig_s7 <- panel_S7A +
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
}

save_fig(fig_s7, "fig_s07_statistical_transparency", FIG_S7_DIMS["w"], FIG_S7_DIMS["h"])
} # end fig_s07


# =============================================================================
# Figure S8: Comprehensive Site-Level Appendix
# =============================================================================
# Individual lnRR time series for ALL taxa at ALL sites, including sites
# excluded from the final analysis. Provides transparency and allows readers
# to inspect site-level variation underlying the meta-analytic summaries.

if (should_render("fig_s08")) {
cat("\n--- Figure S8: Site-Level Appendix ---\n")

# Determine which MPAs are in the final analysis
final_mpas <- if (exists("SumStats.Final")) unique(SumStats.Final$MPA) else character(0)

# Map full scientific names to abbreviations used in color palette
taxa_name_map <- c(
  "Strongylocentrotus purpuratus" = "S. purpuratus",
  "Mesocentrotus franciscanus" = "M. franciscanus",
  "Macrocystis pyrifera" = "M. pyrifera",
  "Panulirus interruptus" = "P. interruptus",
  "Semicossyphus pulcher" = "S. pulcher"
)

# Trophic ordering for consistent plotting
taxa_plot_order <- c("P. interruptus", "S. pulcher",
                     "S. purpuratus", "M. franciscanus",
                     "M. pyrifera")

# Use All.RR.sub.trans — the combined response ratio dataset with all sites
appendix_data <- All.RR.sub.trans %>%
  mutate(
    Taxa_Short = ifelse(y %in% names(taxa_name_map),
                        taxa_name_map[y], y),
    year_num = as.numeric(year)
  ) %>%
  filter(Taxa_Short %in% taxa_plot_order)

# Get MPA implementation years from Site metadata
if (exists("Site")) {
  mpa_impl <- Site %>%
    dplyr::select(CA_MPA_Name_Short, MPA_Start) %>%
    distinct()
} else {
  mpa_impl <- data.frame(CA_MPA_Name_Short = character(0),
                          MPA_Start = numeric(0))
}

# Build one figure per taxa
for (taxa_i in taxa_plot_order) {

  cat("  Building appendix for:", taxa_i, "\n")

  taxa_dat <- appendix_data %>%
    filter(Taxa_Short == taxa_i) %>%
    left_join(mpa_impl, by = "CA_MPA_Name_Short") %>%
    mutate(
      In_Final = CA_MPA_Name_Short %in% final_mpas,
      Status = ifelse(In_Final, "Included in analysis", "Excluded from analysis")
    )

  if (nrow(taxa_dat) < 5) {
    cat("    Skipping", taxa_i, "- insufficient data\n")
    next
  }

  # Order MPAs: included first, then excluded, alphabetically within each
  mpa_levels <- taxa_dat %>%
    distinct(CA_MPA_Name_Short, In_Final) %>%
    arrange(desc(In_Final), CA_MPA_Name_Short) %>%
    pull(CA_MPA_Name_Short)

  taxa_dat$CA_MPA_Name_Short <- factor(taxa_dat$CA_MPA_Name_Short,
                                         levels = mpa_levels)

  n_mpas <- length(mpa_levels)
  n_cols <- 4
  n_rows <- ceiling(n_mpas / n_cols)

  # Compute annual means per MPA (aggregate across response types within this taxa)
  annual_means <- taxa_dat %>%
    group_by(CA_MPA_Name_Short, year_num, resp, Status, MPA_Start) %>%
    summarise(
      mean_lnRR = mean(lnDiff, na.rm = TRUE),
      se_lnRR = sd(lnDiff, na.rm = TRUE) / sqrt(n()),
      n = n(),
      .groups = "drop"
    ) %>%
    filter(is.finite(mean_lnRR))

  # Use density as primary response where available; fall back to biomass
  primary_resp <- if (taxa_i == "M. pyrifera") "Bio" else "Den"
  plot_dat <- annual_means %>% filter(resp == primary_resp)
  # If no data for primary response, try the other

  if (nrow(plot_dat) < 5) {
    alt_resp <- ifelse(primary_resp == "Den", "Bio", "Den")
    plot_dat <- annual_means %>% filter(resp == alt_resp)
    resp_label <- ifelse(alt_resp == "Bio", "Biomass lnRR", "Density lnRR")
  } else {
    resp_label <- ifelse(primary_resp == "Bio", "Biomass lnRR", "Density lnRR")
  }

  if (nrow(plot_dat) < 5) {
    cat("    Skipping", taxa_i, "- insufficient plot data\n")
    next
  }

  # Get the taxa color
  taxa_color <- col_taxa[taxa_i]

  # Use consistent axes across all facets for this taxa (improves readability
  # and makes small-multiple comparisons possible).
  year_rng <- range(plot_dat$year_num, na.rm = TRUE)
  if (!all(is.finite(year_rng))) year_rng <- c(2000, 2022)
  year_span <- diff(year_rng)
  by_val <- if (year_span > 30) 20 else 10
  year_min <- floor(year_rng[1] / by_val) * by_val
  year_max <- min(ceiling(year_rng[2] / by_val) * by_val, year_rng[2] + 2)
  year_breaks <- pretty(year_min:year_max, n = 4)

  y_lim <- max(
    abs(plot_dat$mean_lnRR - plot_dat$se_lnRR),
    abs(plot_dat$mean_lnRR + plot_dat$se_lnRR),
    na.rm = TRUE
  ) * 1.15
  if (!is.finite(y_lim)) y_lim <- 2
  y_lim <- max(2, ceiling(y_lim))
  y_by <- if (y_lim > 4) 2 else 1

  fig_appendix <- ggplot(plot_dat,
                          aes(x = year_num, y = mean_lnRR)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
               linewidth = 0.3) +
    # MPA implementation year (vertical line)
    geom_vline(aes(xintercept = MPA_Start),
               linetype = "dotted", color = "grey40", linewidth = 0.4) +
    # Ribbon for SE
    geom_ribbon(aes(ymin = mean_lnRR - se_lnRR,
                    ymax = mean_lnRR + se_lnRR),
                fill = taxa_color, alpha = 0.15) +
    geom_line(color = taxa_color, linewidth = 0.6, alpha = 0.8) +
    geom_point(color = taxa_color, size = 1, alpha = 0.7) +
    scale_x_continuous(limits = c(year_min, year_max), breaks = year_breaks) +
    coord_cartesian(ylim = c(-y_lim, y_lim)) +
    scale_y_continuous(breaks = seq(-y_lim, y_lim, by = y_by)) +
    facet_wrap(~ CA_MPA_Name_Short, ncol = n_cols, scales = "fixed") +
    labs(
      title = taxa_i,
      x = "Year",
      y = resp_label
    ) +
    theme_mpa(base_size = 9) +
    theme(
      plot.title = element_text(size = 10, face = "italic",
                                 margin = margin(b = 6)),
      strip.text = element_text(size = 8, face = "plain"),
      strip.background = element_blank(),
      axis.text = element_text(size = 7.5),
      panel.grid.major = element_blank(),
      panel.spacing = unit(0.4, "lines"),
      # Reinforce L-shaped axes
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.line.x.bottom = element_line(colour = "black", linewidth = 0.3),
      axis.line.y.left = element_line(colour = "black", linewidth = 0.3)
    )

  # Mark excluded sites with a subtle indicator in the strip
  # (ggplot2 doesn't natively support per-facet strip styling, so we add
  # a small annotation instead)
  excluded_sites <- mpa_levels[!mpa_levels %in% final_mpas]
  if (length(excluded_sites) > 0) {
    # Add a dagger symbol to excluded MPA names in the factor labels
    new_labels <- setNames(
      ifelse(mpa_levels %in% final_mpas,
             as.character(mpa_levels),
             paste0(as.character(mpa_levels), " \u2020")),
      mpa_levels
    )
    fig_appendix <- fig_appendix +
      facet_wrap(~ CA_MPA_Name_Short, ncol = n_cols, scales = "fixed",
                 labeller = labeller(CA_MPA_Name_Short = new_labels))
  }

  # Dimensions scale with number of MPAs (capped at Conservation Letters max)
  fig_w <- min(17, n_cols * 4)  # Conservation Letters max width: 17cm
  fig_h <- max(10, n_rows * 4)

  # Create clean filename from taxa name
  taxa_filename <- tolower(gsub("\\. ", "", taxa_i))
  save_fig(fig_appendix,
           paste0("fig_s08_appendix_", taxa_filename),
           fig_w, fig_h)
}

cat("  Site-level appendix complete. Dagger (\u2020) marks sites excluded from analysis.\n")
} # end fig_s08

cat("\n")
cat("========================================================================\n")
if (identical(RENDER_FIGURES, "all")) {
  cat("  ALL MANUSCRIPT FIGURES GENERATED SUCCESSFULLY\n")
  cat("  Main text: Fig 1-4\n")
  cat("  Supplemental: Fig S1-S2, S7-S8\n")
} else {
  cat("  SELECTED FIGURES GENERATED:", paste(RENDER_FIGURES, collapse = ", "), "\n")
}
cat("========================================================================\n")
cat("\n=== Figures saved to:", here::here("plots"), "===\n")
