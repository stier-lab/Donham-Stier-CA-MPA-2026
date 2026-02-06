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

# =============================================================================
# Figure dimension constants (Conservation Letters specifications)
# =============================================================================
# Journal max widths: single column = 80mm, double column = 170mm
FIG_WIDTH_SINGLE <- 12   # cm, for single-column figures
FIG_WIDTH_DOUBLE <- 17   # cm, for double-column figures
FIG_WIDTH_WIDE   <- 18   # cm, for figures needing extra width
FIG_WIDTH_SUPP   <- 22   # cm, for supplemental figures with legends

# Figure-specific dimensions (width, height in cm)
# Note: Figure 1 dimensions are defined in analysis/R/fig01_map.R
FIG2_DIMS <- c(w = 18, h = 7)    # Data processing pipeline
FIG3_DIMS <- c(w = 16, h = 11)   # Mean effects
FIG4_DIMS <- c(w = 14, h = 12)   # Urchin-kelp scatter
FIG_S1_DIMS <- c(w = 22, h = 26) # Forest plot (supplemental)
FIG_S2_DIMS <- c(w = 20, h = 26) # All taxa time series (supplemental) - wider for legend, taller for title

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

# =============================================================================
# Simple save function - PDF and PNG
save_fig <- function(plot, name, w, h) {
  pdf_path <- here::here("plots", paste0(name, ".pdf"))
  png_path <- here::here("plots", paste0(name, ".png"))

  # Prefer cairo for cleaner text rendering when available.
  if (capabilities("cairo")) {
    ggsave(pdf_path, plot, width = w, height = h, units = "cm",
           device = cairo_pdf, bg = "white", limitsize = FALSE)
  } else {
    ggsave(pdf_path, plot, width = w, height = h, units = "cm",
           device = "pdf", bg = "white", limitsize = FALSE)
  }
  ggsave(png_path, plot, width = w, height = h, units = "cm",
         dpi = 300, bg = "white", limitsize = FALSE)
  cat("  Saved:", name, "\n")
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
# Helper function to standardize MPA names for display (used in multiple figures)
# =============================================================================
shorten_mpa_name <- function(mpa_name) {

  # Replacement patterns: names are regex patterns, values are replacements
  # Using stringr::str_replace_all with named vector for maintainability
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
    legend.title = element_text(face = "bold", size = title_size),
    legend.text = element_text(size = text_size),
    legend.spacing.x = unit(0.4, "cm")
  )
}

# Theme modifications for right-positioned legends (supplemental time series)
theme_legend_right <- function(title_size = 9, text_size = 8.5, italic = TRUE) {
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = title_size),
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
# Figure 1: MPA Map with Bathymetry and Time Series Panels
# =============================================================================
# Publication-quality map showing:
#   - Ocean bathymetry (depth gradient)
#   - MPA boundaries
#   - Monitoring sites with paired Inside/Outside design
#   - 4 kelp biomass time series panels

cat("Building Figure 1: MPA Map with Bathymetry + Time Series...\n")

# Load packages needed for figures
library(patchwork)

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
  FIG1_DIMS <- c(w = 40, h = 32)  # cm (taller to accommodate panels below map)

  # --- 1. Define Study Region ---
  BBOX_LONLAT <- c(xmin = -120.75, ymin = 33.42, xmax = -117.65, ymax = 34.50)

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

  PANEL_SITES <- c("a" = "Campus Point SMCA", "b" = "Harris Point SMR",
                   "c" = "South Point SMR", "d" = "Santa Barbara Island SMR")
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
    filter(!is.na(program))

  offset <- 0.012
  sites_inside <- sites_base %>%
    mutate(status = "Inside MPA", Lon = Lon - offset / 2, Lat = Lat + offset / 2)
  sites_outside <- sites_base %>%
    mutate(status = "Outside MPA", Lon = Lon + offset / 2, Lat = Lat - offset / 2)

  fig1_sites <- bind_rows(sites_inside, sites_outside)
  sites_pairs <- sites_base %>%
    mutate(x_in = Lon - offset / 2, y_in = Lat + offset / 2,
           x_out = Lon + offset / 2, y_out = Lat - offset / 2)

  # Create short site names for labels
  site_short_names <- c(
    "Campus Point SMCA" = "Campus Pt.",
    "Harris Point SMR" = "Harris Pt.",
    "South Point SMR" = "South Pt.",
    "Santa Barbara Island SMR" = "SB Island"
  )

  panel_label_df <- tibble::tibble(
    CA_MPA_Name_Short = unname(PANEL_SITES),
    panel_letter = names(PANEL_SITES),
    site_label = paste0("(", names(PANEL_SITES), ") ", site_short_names[unname(PANEL_SITES)])
  )
  sites_labels <- sites_base %>% inner_join(panel_label_df, by = "CA_MPA_Name_Short")

  # --- 6. Load Time Series Data ---
  ts_cache <- here::here("data", "cache", "figure_data.rds")
  ts_data <- NULL
  if (file.exists(ts_cache)) {
    fig_data <- readRDS(ts_cache)
    ts_data <- fig_data$All.Resp.sub
  }

  # --- 7. Color Palettes ---
  ocean_colors <- c("#08306B", "#2171B5", "#6BAED6", "#C6DBEF")
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
    geom_raster(data = bathy_ocean, aes(x = lon, y = lat, fill = depth),
                interpolate = TRUE, alpha = 0.85) +
    scale_fill_gradientn(
      colors = ocean_colors,
      values = scales::rescale(c(-2000, -500, -100, 0)),
      name = "Depth (m)",
      breaks = c(-1500, -1000, -500, -100),
      labels = c("1500", "1000", "500", "100"),
      limits = c(min(bathy_ocean$depth, na.rm = TRUE), 0),
      guide = guide_colorbar(barwidth = 1.0, barheight = 6, title.position = "top",
                             title.hjust = 0.5, frame.colour = "grey50", ticks.colour = "grey50")
    ) +
    geom_contour(data = bathy_ocean, aes(x = lon, y = lat, z = depth),
                 breaks = c(-100, -200, -500, -1000), color = "white",
                 linewidth = 0.15, alpha = 0.35) +
    geom_sf(data = mpa, fill = fig1_map_colors$mpa_fill, color = fig1_map_colors$mpa_border,
            alpha = 0.30, linewidth = 0.5, inherit.aes = FALSE) +
    new_scale_fill() +
    geom_spatraster(data = hillshade, maxcell = 5e5) +
    scale_fill_gradientn(colors = terrain_pal, na.value = NA, guide = "none") +
    geom_sf(data = coast, fill = alpha("#F5F0E1", 0.15), color = NA, inherit.aes = FALSE) +
    geom_sf(data = coast, fill = NA, color = fig1_map_colors$coastline,
            linewidth = 0.5, inherit.aes = FALSE) +
    geom_segment(data = sites_pairs, aes(x = x_in, y = y_in, xend = x_out, yend = y_out),
                 color = "grey60", linewidth = 0.8, alpha = 0.7) +
    new_scale_fill() +
    geom_point(data = fig1_sites, aes(x = Lon, y = Lat, fill = status, shape = program),
               size = 3.5, color = "white", stroke = 0.8) +
    scale_fill_manual(name = "MPA Status", values = fig1_status_colors,
                      guide = guide_legend(order = 1, override.aes = list(shape = 21, size = 4))) +
    scale_shape_manual(name = "Program", values = program_shapes,
                       guide = guide_legend(order = 2, override.aes = list(fill = "grey50", size = 4))) +
    geom_text(data = sites_labels, aes(x = Lon, y = Lat - 0.08, label = site_label),
              size = 3.0, fontface = "bold", color = "#1A1A1A") +
    coord_sf(xlim = c(BBOX_LONLAT["xmin"], BBOX_LONLAT["xmax"]),
             ylim = c(BBOX_LONLAT["ymin"], BBOX_LONLAT["ymax"]), expand = FALSE, crs = 4326) +
    annotation_scale(location = "bl", width_hint = 0.2, pad_x = unit(0.3, "in"),
                     pad_y = unit(0.3, "in"), style = "ticks", text_cex = 0.75, line_width = 0.4) +
    annotation_north_arrow(location = "tr", which_north = "true", pad_x = unit(0.25, "in"),
                           pad_y = unit(0.25, "in"), style = north_arrow_minimal,
                           height = unit(0.7, "cm"), width = unit(0.7, "cm")) +
    theme_minimal(base_size = 10) +
    theme(
      panel.background = element_rect(fill = "#C6DBEF", color = NA),
      panel.grid.major = element_line(color = "white", linewidth = 0.15),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey30", linewidth = 0.5),
      axis.title = element_blank(),
      axis.text = element_text(size = 8, color = "grey30"),
      legend.position = c(0.95, 0.65),
      legend.justification = c(1, 0.5),
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.45, "cm"),
      legend.background = element_rect(fill = alpha("white", 0.90), color = "grey50", linewidth = 0.3),
      legend.margin = margin(5, 5, 5, 5),
      legend.box.background = element_rect(fill = alpha("white", 0.90), color = "grey50", linewidth = 0.3),
      plot.margin = FIG1_PLOT_MARGIN
    )

  # --- 9. California Inset Map ---
  inset_bbox <- st_bbox(BBOX_LONLAT, crs = 4326)
  inset_rect <- st_as_sfc(inset_bbox)

  inset_map <- ggplot() +
    geom_sf(data = coast, fill = "grey90", color = "grey60", linewidth = 0.25) +
    geom_sf(data = inset_rect, fill = NA, color = "#C03A2B", linewidth = 0.7) +
    coord_sf(xlim = c(-125, -114), ylim = c(32, 42), expand = FALSE) +
    theme_void() +
    theme(panel.background = element_rect(fill = "aliceblue", color = "grey50", linewidth = 0.4),
          plot.margin = margin(2, 2, 2, 2))

  main_map <- main_map +
    patchwork::inset_element(inset_map, left = 0.76, bottom = 0.05, right = 0.99, top = 0.28, align_to = "full")

  # --- 10. Time Series Panels ---
  TS_COLORS <- c("Inside MPA" = "#2A7B8E", "Outside MPA" = "#8C7B6A")

  build_ts_panel <- function(data, site_name, letter, mpa_year) {
    short_name <- gsub(" SMCA| SMR", "", site_name)
    site_data <- data %>%
      filter(CA_MPA_Name_Short == site_name,
             grepl("pyrifera|Macrocystis", taxon_name, ignore.case = TRUE),
             resp == "Bio") %>%
      mutate(Status = ifelse(grepl("mpa|inside", status, ignore.case = TRUE),
                             "Inside MPA", "Outside MPA")) %>%
      group_by(year, Status) %>%
      summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")

    if (nrow(site_data) == 0) {
      return(ggplot() + theme_void() + labs(title = paste0("(", letter, ") ", short_name)))
    }

    ggplot(site_data, aes(x = year, y = mean_value, color = Status)) +
      geom_vline(xintercept = mpa_year, linetype = "dashed", color = "grey50", linewidth = 0.5) +
      geom_line(linewidth = 0.9) +
      geom_point(size = 2, shape = 21, fill = "white", stroke = 0.8) +
      scale_color_manual(values = TS_COLORS) +
      scale_x_continuous(breaks = seq(1990, 2020, by = 10), limits = c(1988, 2024)) +
      scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
      labs(title = paste0("(", letter, ") ", short_name), x = NULL,
           y = expression(Kelp~biomass~(g~m^{-2}))) +
      theme_mpa(base_size = 10) +
      theme(panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(color = "grey92", linewidth = 0.25),
            plot.title = element_text(size = 11, face = "bold", hjust = 0, margin = margin(t = 6, b = 4)),
            axis.title = element_text(size = 9, color = "grey20"),
            axis.title.y = element_text(margin = margin(r = 6)),
            axis.text = element_text(size = 8, color = "grey30"),
            legend.position = "none",
            plot.margin = FIG1_PLOT_MARGIN)
  }

  # Build panels
  fig1_panels <- list()
  if (!is.null(ts_data)) {
    for (letter in names(PANEL_SITES)) {
      site <- PANEL_SITES[letter]
      fig1_panels[[letter]] <- build_ts_panel(ts_data, site, letter, MPA_YEARS[site])
    }
  }

  # --- 11. Combine Map + Panels (panels in row below map, same width) ---
  if (length(fig1_panels) == 4) {
    # Style time series panels for consistency, with reduced top margin
    for (letter in names(fig1_panels)) {
      fig1_panels[[letter]] <- fig1_panels[[letter]] +
        theme(legend.position = "none",
              panel.background = element_rect(fill = "white", color = "grey50", linewidth = 0.4),
              plot.margin = margin(2, 4, 4, 4))  # Reduced top margin
    }

    # Reduce bottom margin of main map
    main_map <- main_map + theme(plot.margin = margin(2, 2, 2, 2))

    # Create row of panels using patchwork
    panels_row <- fig1_panels$a + fig1_panels$b + fig1_panels$c + fig1_panels$d +
      plot_layout(nrow = 1)

    # Stack map on top, panels below (map takes ~70% height, panels ~30%)
    fig1 <- main_map / panels_row +
      plot_layout(heights = c(2.5, 1))
  } else {
    fig1 <- main_map
  }

  fig1 <- fig1 + plot_annotation(theme = theme(plot.margin = margin(4, 8, 8, 8)))

  # --- 12. Save Figure 1 ---
  save_fig(fig1, "fig_01_mpa_map", FIG1_DIMS["w"], FIG1_DIMS["h"])

} else {
  cat("  Skipping Figure 1 (required packages not available)\n")
}

# Load additional packages for remaining figures
if (has_cowplot) library(cowplot)
if (has_ggrepel) library(ggrepel)


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

# =============================================================================
# FIX #4: Standardize status BEFORE aggregation
# =============================================================================
# Previously the code did group_by(status) then standardize, which could
# create separate groups for "control"/"reference"/"outside" before collapsing.
# Now we standardize first, then group.

# Panel (a) data: Raw density - STANDARDIZE FIRST
fig2_raw <- All.Resp.sub %>%
  dplyr::filter(
    source == "KFM",
    CA_MPA_Name_Short == "Scorpion SMR",
    taxon_name %in% c("Strongylocentrotus purpuratus", "S. purpuratus"),
    resp == "Den"
  ) %>%
  dplyr::mutate(
    year = as.numeric(as.character(year)),
    # CRITICAL: Standardize status BEFORE grouping/aggregation
    status = standardize_status(status)
  )

# Warn about dropped rows with unmapped status values
n_unmapped_fig2 <- sum(is.na(fig2_raw$status))
if (n_unmapped_fig2 > 0) {
  warning(sprintf("Figure 2: Dropped %d rows with unmapped status values", n_unmapped_fig2))
}

fig2_raw <- fig2_raw %>%
  dplyr::filter(!is.na(status)) %>%  # Remove any that didn't map
  dplyr::group_by(year, status) %>%
  dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(
    status = factor(status, levels = c("Inside", "Outside"))
  )

# FIX: Validate sufficient data remains after filtering
if (nrow(fig2_raw) == 0) {
  stop("Figure 2: No valid data after status standardization. ",
       "Check All.Resp.sub$status values for KFM/Scorpion SMR data.")
}
if (length(unique(fig2_raw$status)) < 2) {
  stop("Figure 2: Missing Inside or Outside data after standardization. ",
       "Found only: ", paste(unique(fig2_raw$status), collapse = ", "))
}

# Debug: Check status values
cat("  Figure 2 - Status values:", paste(unique(fig2_raw$status), collapse = ", "), "\n")
cat("  Figure 2 - Year range:", min(fig2_raw$year, na.rm = TRUE), "-", max(fig2_raw$year, na.rm = TRUE), "\n")

# Panel (b) data: Proportion of maximum
fig2_prop <- fig2_raw %>%
  dplyr::group_by(status) %>%
  dplyr::mutate(
    value = value / max(value, na.rm = TRUE),
    value = ifelse(is.na(value) | is.infinite(value), 0, value)
  ) %>%
  dplyr::ungroup()

# Panel (c) & (d) data: Log response ratio
fig2_lnrr <- All.RR.sub.trans %>%
  dplyr::filter(
    source == "KFM",
    CA_MPA_Name_Short == "Scorpion SMR",
    y %in% c("Strongylocentrotus purpuratus", "S. purpuratus"),
    resp == "Den"
  )

# Ensure BA column exists
if (!"BA" %in% names(fig2_lnrr) || all(is.na(fig2_lnrr$BA))) {
  fig2_lnrr <- fig2_lnrr %>%
    dplyr::mutate(BA = factor(
      ifelse(year < scorpion_start, "Before", "After"),
      levels = c("Before", "After")
    ))
} else {
  fig2_lnrr <- fig2_lnrr %>%
    dplyr::mutate(BA = factor(BA, levels = c("Before", "After")))
}

# FIX: Calculate dynamic x-axis limits from actual data range
fig2_year_range <- range(c(fig2_raw$year, fig2_lnrr$year), na.rm = TRUE)
fig2_x_limits <- c(floor(fig2_year_range[1] / 5) * 5 - 2,  # Round down to 5, add padding
                   ceiling(fig2_year_range[2] / 5) * 5 + 2)  # Round up to 5, add padding
fig2_x_breaks <- seq(fig2_x_limits[1], fig2_x_limits[2], by = 10)
fig2_x_breaks <- fig2_x_breaks[fig2_x_breaks >= fig2_x_limits[1] & fig2_x_breaks <= fig2_x_limits[2]]

# Build each panel individually for better control, then combine with patchwork
# Common theme settings for consistent panel appearance
fig2_theme <- theme_mpa(base_size = 10) +
  theme(legend.position = "none",
        plot.title = element_text(size = 9.5, face = "bold", hjust = 0,
                                  margin = margin(b = 8)),
        axis.title.y = element_text(size = 9, margin = margin(r = 4)),
        axis.text = element_text(size = 8, color = "black"),
        plot.margin = margin(4, 4, 4, 4))

p2a <- ggplot(fig2_raw, aes(x = year, y = value, color = status)) +
  add_mpa_vline(scorpion_start) +
  # TUFTE: Annotate the MPA implementation line
  annotate("text", x = scorpion_start, y = Inf,
           label = "MPA\nestablished", hjust = 1.1, vjust = 1.3,
           size = 2.3, color = MPA_LABEL_COLOR, lineheight = 0.9) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 0.8) +
  scale_color_site(name = "Site") +
  # TUFTE: Declarative title - concise to fit panel width
  labs(title = "(a) MPA > Reference",
       x = NULL, y = expression(Density~(ind~m^{-2}))) +
  scale_x_continuous(breaks = fig2_x_breaks, limits = fig2_x_limits) +
  fig2_theme

p2b <- ggplot(fig2_prop, aes(x = year, y = value, color = status)) +
  add_mpa_vline(scorpion_start) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 0.8) +
  scale_color_site(name = "Site") +
  # TUFTE: Concise title
  labs(title = "(b) Standardized", x = NULL, y = "Proportion of max") +
  scale_x_continuous(breaks = fig2_x_breaks, limits = fig2_x_limits) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5)) +
  fig2_theme

p2c <- ggplot(fig2_lnrr, aes(x = year, y = lnDiff, shape = BA)) +
  add_mpa_vline(scorpion_start) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey60", linewidth = 0.5) +
  geom_point(size = 2.5, color = col_taxa["S. purpuratus"]) +
  scale_shape_manual(name = "Period", values = c("Before" = 1, "After" = 16)) +
  # TUFTE: Concise title
  labs(title = "(c) MPA effect", x = NULL, y = "ln(MPA / Reference)") +
  scale_x_continuous(breaks = fig2_x_breaks, limits = fig2_x_limits) +
  fig2_theme

# Panel (d): Same with linear trend in After period
fig2d_after <- dplyr::filter(fig2_lnrr, BA == "After")

p2d <- ggplot(fig2_lnrr, aes(x = year, y = lnDiff, shape = BA)) +
  add_mpa_vline(scorpion_start) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey60", linewidth = 0.5) +
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
  scale_shape_manual(name = "Period", values = c("Before" = 1, "After" = 16)) +
  # TUFTE: Concise title with finding
  labs(title = "(d) Declining trend", x = NULL, y = NULL) +
  scale_x_continuous(breaks = fig2_x_breaks, limits = fig2_x_limits) +
  fig2_theme

# FIX: Use patchwork to collect legends into single unified row
# Enable legends on first panel of each type (Site from p2a, Period from p2c)
p2a <- p2a + theme(legend.position = "bottom")
p2c <- p2c + theme(legend.position = "bottom")

# Combine panels using patchwork with collected guides
arrow <- fig2_arrow_panel()
fig2_panels <- (p2a + arrow + p2b + arrow + p2c + arrow + p2d) +
  plot_layout(
    ncol = 7,
    widths = c(1, 0.08, 1, 0.08, 1, 0.08, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold", size = 10, color = "black"),
    legend.text = element_text(size = 9, color = "black"),
    legend.spacing.x = unit(8, "mm")
  )

# TUFTE: Add declarative overall title using patchwork annotation
# States the key finding rather than just describing what's shown
fig2_titled <- fig2_panels +
  plot_annotation(
    title = expression("Purple urchin density declined after MPA protection at Scorpion SMR"),
    theme = theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5,
                                margin = margin(b = 8))
    )
  )

# Conservation Letters: max 170mm double-column width
save_fig(fig2_titled, "fig_02_data_processing", FIG2_DIMS["w"], FIG2_DIMS["h"])

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

# Convert MPA to character (as.character() on a factor returns level labels, not codes)
fig_s1_data$MPA <- as.character(fig_s1_data$MPA)

# Debug: Check MPA values after conversion
cat("  Figure S1 - Unique MPAs after conversion:",
    paste(head(unique(fig_s1_data$MPA), 8), collapse = ", "), "\n")
cat("  Figure S1 - Total rows:", nrow(fig_s1_data), "\n")

# Create shortened MPA names using standardized function
fig_s1_data <- fig_s1_data %>%
  dplyr::mutate(MPA_short = shorten_mpa_name(MPA))

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

# FIX: Calculate consistent x-axis limits across all taxa
x_range_all <- range(c(fig_s1_data$Mean - fig_s1_data$CI,
                       fig_s1_data$Mean + fig_s1_data$CI), na.rm = TRUE)
x_limit <- max(abs(x_range_all)) * 1.1
x_breaks <- seq(-10, 10, by = 5)
x_breaks <- x_breaks[abs(x_breaks) <= x_limit]

fig_s1 <- ggplot(fig_s1_data,
               aes(x = MPA_order, y = Mean, ymin = Mean - CI, ymax = Mean + CI,
                   color = Resp, shape = Source)) +
  # Reference line at zero (no effect)
  geom_hline(yintercept = 0, linetype = "solid", color = "grey40",
             linewidth = 0.6) +
  # Confidence intervals (inherits ymin/ymax from main aes)
  geom_linerange(linewidth = 0.7, position = position_dodge(width = 0.6)) +
  # Points with shape/color encoding
  geom_point(size = 2.8, position = position_dodge(width = 0.6)) +
  # Free the discrete MPA axis within each taxa panel, keep effect-size axis comparable.
  facet_wrap(~ Taxa, ncol = 2, scales = "free_x") +
  scale_color_response(name = "Response",
                       labels = c("Den" = "Density", "Bio" = "Biomass")) +
  scale_shape_source(name = "Source") +
  scale_y_continuous(breaks = x_breaks) +
  coord_flip(ylim = c(-x_limit, x_limit), clip = "off") +
  # TUFTE: Declarative title and informative subtitle
  labs(title = "MPA effects vary across sites but show consistent patterns by taxa",
       subtitle = "Urchins generally decline while predators increase across California MPAs",
       x = NULL, y = "Effect Size (lnRR)",
       caption = "Error bars: 95% CI. Horizontal line at zero indicates no MPA effect.") +
  theme_mpa(base_size = 10) +
  theme(
    # TUFTE: Title and subtitle styling
    plot.title = element_text(face = "bold", size = 12, hjust = 0,
                              margin = margin(b = 4)),
    plot.subtitle = element_text(size = 10, hjust = 0, color = "grey40",
                                 margin = margin(b = 12)),
    # FIX: Consistent strip background across all panels
    strip.text = element_text(face = "italic", size = 11, margin = margin(4, 0, 4, 0)),
    strip.background = element_rect(fill = "grey95", color = "grey70", linewidth = 0.4),
    axis.text.y = element_text(size = 7.5, color = "grey20", hjust = 1),
    axis.text.x = element_text(size = 9),
    axis.title.x = element_text(size = 10, margin = margin(t = 8)),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.box.spacing = unit(0.8, "cm"),
    legend.spacing.x = unit(0.4, "cm"),
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8.5),
    legend.margin = margin(t = 8, b = 5),
    legend.key.width = unit(0.8, "cm"),
    panel.spacing.x = unit(1.2, "lines"),
    panel.spacing.y = unit(0.8, "lines"),
    # FIX: Consistent panel backgrounds
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.25),
    plot.caption = element_text(size = 7.5, color = "grey50", hjust = 1,
                                 margin = margin(t = 10)),
    plot.margin = margin(8, 12, 8, 8)
  )

# Conservation Letters: Supplemental figures - slightly wider to fit MPA names
save_fig(fig_s1, "fig_s01_forest_plot", FIG_S1_DIMS["w"], FIG_S1_DIMS["h"])

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

# Use existing col_response_long palette from 00b_color_palette.R
# (Has "Density" and "Biomass" keys mapped to the same colors as col_response)

# Calculate sample sizes for annotation
# Ensure Taxa is regular factor (not ordered) to avoid join incompatibility
fig3_n <- fig3_individual %>%
  dplyr::mutate(Taxa = factor(as.character(Taxa), levels = taxa_levels)) %>%
  dplyr::group_by(Taxa, Response) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop")

# Define unified dodge width for perfect alignment between raw data and summary stats
pd <- position_dodge(width = 0.6)

# Add sample sizes to meta data for annotation
# Ensure Taxa types match before join
fig3_meta <- fig3_meta %>%
  dplyr::mutate(Taxa = factor(as.character(Taxa), levels = taxa_levels)) %>%
  dplyr::left_join(fig3_n, by = c("Taxa", "Response"))

# Dynamic y-limits and sample-size label position (robust to data updates).
fig3_y_min <- min(c(fig3_individual$Mean, fig3_meta$CI_lower), na.rm = TRUE)
fig3_y_max <- max(c(fig3_individual$Mean, fig3_meta$CI_upper), na.rm = TRUE)
fig3_pad <- 0.6
fig3_label_y <- fig3_y_min - fig3_pad

fig3 <- ggplot() +
  # Reference line at zero (no effect) - lighter dotted line
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.4, linetype = "dashed") +
  # Layer 1: Raw data points (circles) - FIX: more jitter to reduce overlap
  geom_point(data = fig3_individual,
             aes(x = Taxa, y = Mean, color = Response),
             position = position_jitterdodge(jitter.width = 0.22,
                                              dodge.width = 0.6,
                                              seed = 42),
             size = 1.7, alpha = 0.28, shape = 16) +
  # Layer 2: White outline behind diamonds for better visibility
  geom_point(data = fig3_meta,
             aes(x = Taxa, y = Estimate),
             position = pd, size = 8, shape = 18, color = "white") +
  # Layer 3: Confidence intervals (error bars) - FIX: thicker lines
  geom_errorbar(data = fig3_meta,
                aes(x = Taxa, ymin = CI_lower, ymax = CI_upper, color = Response),
                position = pd, width = 0.15, linewidth = 1.0) +
  # Layer 4: Meta-analytic means (diamonds) - FIX: LARGER diamonds
  geom_point(data = fig3_meta,
             aes(x = Taxa, y = Estimate, color = Response),
             position = pd, size = 6.5, shape = 18) +
  # FIX: Add sample size labels below x-axis
  geom_text(data = fig3_meta,
            aes(x = Taxa, y = fig3_label_y, label = paste0("n=", n), color = Response),
            position = pd, size = 2.8, show.legend = FALSE) +
  scale_color_manual(name = "Response", values = col_response_long) +
  coord_cartesian(ylim = c(fig3_y_min - 1.2, fig3_y_max + 0.6)) +
  # TUFTE: Declarative title - states the key finding
  labs(title = "Urchins decreased while kelp and predators increased in MPAs",
       x = NULL, y = "Effect Size (lnRR)",
       caption = "Diamonds = meta-analytic means; circles = individual MPA estimates") +
  theme_mpa(base_size = 11) +
  theme(
    # TUFTE: Bold, centered title
    plot.title = element_text(face = "bold", size = 11, hjust = 0.5,
                              margin = margin(b = 10)),
    axis.text.x = element_text(face = "italic", size = 10),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 10),
    legend.key.width = unit(1.2, "cm"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.25)
  )

save_fig(fig3, "fig_03_mean_effects", FIG3_DIMS["w"], FIG3_DIMS["h"])

# =============================================================================
# Figure 4: Urchin density vs Kelp BIOMASS scatterplot
# Manuscript: Main Text Figure 4
# =============================================================================
# FIX #5: The original code filtered kelp as "density" but the figure title and
# comments said "biomass". This fix implements option A: use kelp BIOMASS.
# If biomass data is unavailable, fall back to density with updated labels.
# =============================================================================

cat("Building Figure 4: Urchin density vs Kelp biomass scatterplot...\n")

# Prepare data for urchin vs kelp scatterplot
# Urchin: always use density
fig4_urchin <- All.RR.sub.trans %>%
  dplyr::filter(
    y %in% c("Strongylocentrotus purpuratus", "S. purpuratus", "STRPURAD"),
    resp %in% c("Den", "Density")
  ) %>%
  dplyr::select(CA_MPA_Name_Short, year, source, lnDiff) %>%
  dplyr::rename(lnRR_urchin = lnDiff)

# FIX: Try BIOMASS first for kelp (to match figure description)
fig4_kelp_bio <- All.RR.sub.trans %>%
  dplyr::filter(
    y %in% c("Macrocystis pyrifera", "M. pyrifera", "MACPYRAD"),
    resp %in% c("Bio", "Biomass")
  ) %>%
  dplyr::select(CA_MPA_Name_Short, year, source, lnDiff) %>%
  dplyr::rename(lnRR_kelp = lnDiff)

# Require biomass data for kelp (this figure is meant to show biomass relationship)
# If biomass is unavailable, fail explicitly rather than silently falling back to density
if (nrow(fig4_kelp_bio) == 0) {
  stop(paste0(
    "Figure 4 ERROR: Kelp biomass data is required but unavailable.\n",
    "  This figure is designed to show the trophic cascade (urchin density vs kelp biomass).\n",
    "  Check that kelp biomass response ('Bio' or 'Biomass') is present in All.RR.sub.trans.\n",
    "  If density is acceptable, modify the figure code explicitly."
  ))
}
fig4_kelp <- fig4_kelp_bio
kelp_metric <- "biomass"
kelp_y_label <- expression("MPA Effect on " * italic("M. pyrifera") * " biomass (lnRR)")
cat("  Using kelp BIOMASS data (", nrow(fig4_kelp), " rows)\n")

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

  # Build scatterplot showing relationship between urchin and kelp lnRR

  # Symmetric limits centered on 0 so quadrants are visually balanced.
  max_abs <- max(abs(c(fig4_data$lnRR_urchin, fig4_data$lnRR_kelp)), na.rm = TRUE)
  max_abs <- max(6, ceiling(max_abs))
  x_lim <- c(-max_abs, max_abs)
  y_lim <- c(-max_abs, max_abs)

  # Quadrant label positions (slightly inside plot area)
  q_x_left <- x_lim[1] + 0.3
  q_x_right <- x_lim[2] - 0.3
  q_y_top <- y_lim[2] - 0.5
  q_y_bot <- y_lim[1] + 0.5

  fig4 <- ggplot(fig4_data, aes(x = lnRR_urchin, y = lnRR_kelp)) +
    # Reference lines at zero (drawn first, underneath everything)
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
    # Quadrant labels - subtle, positioned dynamically in corners
    annotate("text", x = q_x_left, y = q_y_top,
             label = "Urchin \u2193  Kelp \u2191", size = 2.5,
             color = "grey45", fontface = "italic", hjust = 0, vjust = 1) +
    annotate("text", x = q_x_right, y = q_y_top,
             label = "Urchin \u2191  Kelp \u2191", size = 2.5,
             color = "grey45", fontface = "italic", hjust = 1, vjust = 1) +
    annotate("text", x = q_x_left, y = q_y_bot,
             label = "Urchin \u2193  Kelp \u2193", size = 2.5,
             color = "grey45", fontface = "italic", hjust = 0, vjust = 0) +
    annotate("text", x = q_x_right, y = q_y_bot,
             label = "Urchin \u2191  Kelp \u2193", size = 2.5,
             color = "grey45", fontface = "italic", hjust = 1, vjust = 0) +
    # Layer 1: Contour lines for point-density structure
    geom_density_2d(
      color = "grey35",
      linewidth = 0.35,
      alpha = 0.45,
      bins = 6
    ) +
    # Layer 2: Data points - FIX: more alpha transparency to reduce overplotting
    geom_point(
      position = position_jitter(width = 0.12, height = 0.12, seed = 42),
      size = 1.9,
      alpha = 0.28,
      color = col_taxa["S. purpuratus"],
      stroke = 0
    ) +
    # Layer 3: Linear regression line with confidence band
    geom_smooth(
      method = "lm",
      se = TRUE,
      color = col_response["Bio"],
      fill = col_response["Bio"],
      linewidth = 1.0,
      alpha = 0.20,
      linetype = "solid"
    ) +
    # Correlation and sample size annotation - positioned in top-right
    annotate(
      "label",
      x = x_lim[2],
      y = y_lim[2],
      label = paste0(cor_label, "\nn = ", n_points),
      hjust = 1.05,
      vjust = 1.1,
      size = 3.2,
      fontface = "plain",
      fill = "white",
      alpha = 0.85,
      label.size = 0.3,
      label.padding = unit(0.2, "lines")
    ) +
    # TUFTE: Declarative title - states the trophic cascade finding
    labs(
      title = "MPA protection reduces urchins and increases kelp biomass",
      x = expression("MPA Effect on " * italic("S. purpuratus") * " density (lnRR)"),
      y = kelp_y_label  # Uses biomass or density depending on data availability
    ) +
    coord_fixed(ratio = 1, xlim = x_lim, ylim = y_lim, expand = FALSE) +
    theme_mpa(base_size = 11) +
    theme(
      # TUFTE: Bold, centered title
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5,
                                margin = margin(b = 10)),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9),
      plot.margin = margin(12, 12, 8, 8)
    )

  save_fig(fig4, "fig_04_urchin_kelp_scatter", FIG4_DIMS["w"], FIG4_DIMS["h"])
} else {
  cat("  WARNING: No overlapping urchin/kelp data found for Figure 4\n")
}

# =============================================================================
# Figure S2 (Supplemental): All taxa log response ratios at example MPAs
# Not in current manuscript - kept as supplemental
# =============================================================================

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

# Helper function to build individual time series panels for Figure S2
build_s2_panel <- function(mpa_name, data, mpa_starts) {
  d <- dplyr::filter(data, CA_MPA_Name_Short == mpa_name)
  mpa_start_yr <- mpa_starts[mpa_name]

  # Clean MPA name for display using standardized function
  mpa_display <- shorten_mpa_name(mpa_name)

  year_rng <- range(d$year, na.rm = TRUE)
  if (!all(is.finite(year_rng))) year_rng <- c(2000, 2022)
  year_breaks <- seq(floor(year_rng[1] / 5) * 5, ceiling(year_rng[2] / 5) * 5, by = 5)

  y_annot <- suppressWarnings(max(d$lnDiff, na.rm = TRUE))
  if (!is.finite(y_annot)) y_annot <- 0
  y_annot <- min(1.8, y_annot)

  ggplot(d, aes(x = year, y = lnDiff,
                color = species_short, shape = source)) +
    # Reference line at zero
    geom_hline(yintercept = 0, linetype = "solid", color = "grey50",
               linewidth = 0.5) +
    # MPA implementation vertical line (using standardized style)
    add_mpa_vline(mpa_start_yr) +
    # Add MPA label (using standardized style)
    annotate("text", x = mpa_start_yr + 0.5, y = y_annot,
             label = "MPA\nStart", hjust = 0, size = MPA_LABEL_SIZE, color = MPA_LABEL_COLOR) +
    # LOESS smoothers for each species (after MPA only)
    geom_smooth(data = dplyr::filter(d, year >= mpa_start_yr),
                aes(group = species_short),
                method = "loess", se = FALSE, span = 0.75,
                linewidth = 0.8, alpha = 0.6) +
    # Data points
    geom_point(size = 2.1, alpha = 0.65) +
    scale_color_taxa(name = "Species") +
    scale_shape_source(name = "Source") +
    scale_x_continuous(breaks = year_breaks) +
    scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 1)) +
    labs(
      title = mpa_display,
      x = "Year",
      y = "Log response ratio (lnRR)"
    ) +
    theme_mpa(base_size = 9) +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0,
                                 margin = margin(b = 6)),
      legend.text = element_text(face = "italic", size = 8.5),
      legend.title = element_text(face = "bold", size = 9),
      legend.key.height = unit(0.5, "cm"),  # More space between legend items
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.25),
      plot.margin = margin(6, 8, 6, 6)
    )
}

# Build panels using purrr::map for tidyverse consistency
fig_s2_panels <- purrr::map(fig_s2_mpas, build_s2_panel,
                             data = fig_s2_data,
                             mpa_starts = fig_s2_starts)

# TUFTE: Use patchwork for better title control (preferred over ggpubr)
if (has_patchwork) {
  fig_s2 <- (fig_s2_panels[[1]] / fig_s2_panels[[2]] / fig_s2_panels[[3]]) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Trophic cascade: urchins decline, kelp and predators increase post-MPA",
      theme = theme(
        plot.title = element_text(face = "bold", size = 11, hjust = 0.5,
                                  margin = margin(b = 10))
      )
    ) &
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 8.5, face = "italic"),
      legend.key.width = unit(0.9, "cm"),
      legend.spacing.x = unit(0.35, "cm")
    )
} else if (has_ggpubr) {
  fig_s2 <- ggpubr::ggarrange(
    plotlist = fig_s2_panels,
    ncol = 1, nrow = 3,
    labels = c("a", "b", "c"),
    font.label = list(size = 11, face = "bold"),
    label.x = 0.02, label.y = 0.98,
    common.legend = TRUE,
    legend = "right"
  )
} else {
  # Basic fallback
  fig_s2 <- fig_s2_panels[[1]] / fig_s2_panels[[2]] / fig_s2_panels[[3]]
}

# Slightly wider to accommodate right legend
save_fig(fig_s2, "fig_s02_all_taxa_timeseries", FIG_S2_DIMS["w"], FIG_S2_DIMS["h"])

# =============================================================================
# All required manuscript figures complete
# =============================================================================
cat("\n")
cat("========================================================================\n")
cat("  ALL MANUSCRIPT FIGURES GENERATED SUCCESSFULLY\n")
cat("  Main text: Fig 1-4\n")
cat("  Supplemental: Fig S1-S2\n")
cat("========================================================================\n")
cat("\n")

# Exit successfully - remaining figures are exploratory only
if (FALSE) {
  # NOTE: The code below (Figures 5-6) is retained for reference but not
  # included in the manuscript. Set if(FALSE) to if(TRUE) to generate.

# =============================================================================
# Figure 5: Trophic Cascade Pathway (Path Analysis)
# =============================================================================
# Shows the complete causal chain: Predators -> Urchins -> Kelp
# Three panels: (a) Pred vs Urch, (b) Urch vs Kelp, (c) Path diagram

cat("\n--- Figure 5: Trophic Cascade Pathway ---\n")

FIG5_DIMS <- c(w = 18, h = 16)

# Prepare cascade data: pivot to wide format with one row per MPA
cascade_data <- SumStats.Final %>%
  filter(AnalysisType %in% c("pBACIPS", "CI")) %>%
  mutate(
    Taxa_Resp = paste(Taxa, Resp, sep = "_"),
    Mean = as.numeric(Mean),
    SE = as.numeric(SE)
  ) %>%
  dplyr::select(MPA, Taxa_Resp, Mean) %>%
  tidyr::pivot_wider(names_from = Taxa_Resp, values_from = Mean, values_fn = mean)

# Clean column names for formula use
names(cascade_data) <- gsub(" ", "_", names(cascade_data))
names(cascade_data) <- gsub("\\.", "_", names(cascade_data))

# Calculate predator and urchin indices (mean of species within trophic level)
cascade_data <- cascade_data %>%
  rowwise() %>%
  mutate(
    Predator_Den = mean(c_across(any_of(c("P_interruptus_Den", "S_pulcher_Den"))), na.rm = TRUE),
    Urchin_Den = mean(c_across(any_of(c("S_purpuratus_Den", "M_franciscanus_Den"))), na.rm = TRUE)
  ) %>%
  ungroup()

# Filter to MPAs with complete cascade data
cascade_complete <- cascade_data %>%
  filter(
    is.finite(Predator_Den),
    is.finite(Urchin_Den),
    is.finite(M_pyrifera_Bio)
  )

cat("  MPAs with complete cascade data:", nrow(cascade_complete), "\n")

# Only build figure if we have enough data
if (nrow(cascade_complete) >= 5) {

  # Panel A: Predators predict Urchins
  cor_pred_urch <- cor.test(cascade_complete$Predator_Den,
                             cascade_complete$Urchin_Den,
                             method = "pearson")
  cor_label_A <- sprintf("r = %.2f, p = %.3f",
                          cor_pred_urch$estimate,
                          cor_pred_urch$p.value)

  panel_5A <- ggplot(cascade_complete,
                      aes(x = Predator_Den, y = Urchin_Den)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
    annotate("text", x = Inf, y = -Inf, label = "Predators up\nUrchins down",
             hjust = 1.1, vjust = -0.1, size = 2.5, color = "grey50", fontface = "italic") +
    geom_point(size = 3.5, alpha = 0.7, color = col_taxa["P. interruptus"]) +
    geom_smooth(method = "lm", se = TRUE, color = col_taxa["P. interruptus"],
                fill = col_taxa["P. interruptus"], alpha = 0.2, linewidth = 1) +
    annotate("label", x = Inf, y = Inf, label = cor_label_A,
             hjust = 1.1, vjust = 1.3, size = 3, fill = "white", alpha = 0.9,
             label.size = 0.3) +
    labs(
      title = "(a) Predator recovery predicts urchin decline",
      x = "Predator effect (lnRR)",
      y = "Urchin effect (lnRR)"
    ) +
    coord_cartesian(clip = "off") +
    theme_mpa(base_size = 10) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0),
      panel.grid.major = element_line(color = "grey95", linewidth = 0.25)
    )

  # Panel B: Urchins predict Kelp
  cor_urch_kelp <- cor.test(cascade_complete$Urchin_Den,
                             cascade_complete$M_pyrifera_Bio,
                             method = "pearson")
  cor_label_B <- sprintf("r = %.2f, p = %.3f",
                          cor_urch_kelp$estimate,
                          cor_urch_kelp$p.value)

  panel_5B <- ggplot(cascade_complete,
                      aes(x = Urchin_Den, y = M_pyrifera_Bio)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
    annotate("text", x = -Inf, y = Inf, label = "Urchins down\nKelp up",
             hjust = -0.1, vjust = 1.2, size = 2.5, color = "grey50", fontface = "italic") +
    geom_point(size = 3.5, alpha = 0.7, color = col_taxa["S. purpuratus"]) +
    geom_smooth(method = "lm", se = TRUE, color = col_taxa["M. pyrifera"],
                fill = col_taxa["M. pyrifera"], alpha = 0.2, linewidth = 1) +
    annotate("label", x = Inf, y = Inf, label = cor_label_B,
             hjust = 1.1, vjust = 1.3, size = 3, fill = "white", alpha = 0.9,
             label.size = 0.3) +
    labs(
      title = "(b) Urchin decline predicts kelp recovery",
      x = "Urchin effect (lnRR)",
      y = "Kelp biomass effect (lnRR)"
    ) +
    coord_cartesian(clip = "off") +
    theme_mpa(base_size = 10) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0),
      panel.grid.major = element_line(color = "grey95", linewidth = 0.25)
    )

  # Panel C: Path diagram summary
  lm_pred_urch <- lm(Urchin_Den ~ Predator_Den, data = cascade_complete)
  lm_urch_kelp <- lm(M_pyrifera_Bio ~ Urchin_Den, data = cascade_complete)

  path_a <- coef(lm_pred_urch)[2]
  path_b <- coef(lm_urch_kelp)[2]
  indirect_effect <- path_a * path_b

  # Node positions for path diagram
  nodes_df <- tibble(
    label = c("PREDATORS\n(lobster,\nsheephead)", "URCHINS\n(purple, red)", "KELP\n(M. pyrifera)"),
    x = c(0, 1, 2),
    y = c(0, 0, 0),
    color = c(col_taxa["P. interruptus"], col_taxa["S. purpuratus"], col_taxa["M. pyrifera"])
  )

  panel_5C <- ggplot() +
    # Arrows between nodes
    annotate("segment", x = 0.22, xend = 0.78, y = 0, yend = 0,
             arrow = arrow(length = unit(4, "mm"), type = "closed"),
             linewidth = 2, color = "grey30") +
    annotate("segment", x = 1.22, xend = 1.78, y = 0, yend = 0,
             arrow = arrow(length = unit(4, "mm"), type = "closed"),
             linewidth = 2, color = "grey30") +
    # Path coefficients on arrows
    annotate("label", x = 0.5, y = 0.18,
             label = sprintf("β = %.2f", path_a),
             size = 4, fill = "white", label.size = 0.3, fontface = "bold") +
    annotate("label", x = 1.5, y = 0.18,
             label = sprintf("β = %.2f", path_b),
             size = 4, fill = "white", label.size = 0.3, fontface = "bold") +
    # Nodes as large colored circles
    geom_point(data = nodes_df, aes(x = x, y = y),
               size = 22, color = nodes_df$color, fill = nodes_df$color, shape = 21) +
    # Node labels
    geom_text(data = nodes_df, aes(x = x, y = y, label = label),
              size = 2.8, color = "white", fontface = "bold", lineheight = 0.9) +
    # Indirect effect annotation
    annotate("text", x = 1, y = -0.38,
             label = sprintf("Indirect effect: %.2f × %.2f = %.2f", path_a, path_b, indirect_effect),
             size = 4, fontface = "bold") +
    annotate("text", x = 1, y = -0.55,
             label = "MPA protection restores top-down control of kelp forests",
             size = 3.5, color = "grey40", fontface = "italic") +
    coord_cartesian(xlim = c(-0.4, 2.4), ylim = c(-0.7, 0.4)) +
    labs(title = "(c) Trophic cascade pathway with effect sizes") +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5,
                                margin = margin(b = 10))
    )

  # Combine panels
  fig5 <- (panel_5A | panel_5B) / panel_5C +
    plot_layout(heights = c(1, 0.65)) +
    plot_annotation(
      title = "MPA protection restores trophic cascades through predator recovery",
      theme = theme(
        plot.title = element_text(face = "bold", size = 12, hjust = 0.5)
      )
    )

  save_fig(fig5, "fig_05_trophic_cascade_pathway", FIG5_DIMS["w"], FIG5_DIMS["h"])

} else {
  cat("  WARNING: Not enough MPAs with complete cascade data for Figure 5\n")
}


# =============================================================================
# Figure S3: Temporal Dynamics of Trophic Cascade
# =============================================================================
# Shows how effect sizes unfold over time since MPA implementation

cat("\n--- Figure S3: Temporal Dynamics ---\n")

FIG_S3_DIMS <- c(w = 18, h = 18)

# Assign trophic levels to taxa
trophic_assignment <- c(
  "Panulirus interruptus" = "Predators",
  "Semicossyphus pulcher" = "Predators",
  "P. interruptus" = "Predators",
  "S. pulcher" = "Predators",
  "Strongylocentrotus purpuratus" = "Urchins",
  "Mesocentrotus franciscanus" = "Urchins",
  "S. purpuratus" = "Urchins",
  "M. franciscanus" = "Urchins",
  "Macrocystis pyrifera" = "Kelp",
  "M. pyrifera" = "Kelp"
)

# Trophic level colors (blend of species colors within each level)
trophic_colors <- c(
  "Predators" = "#5C8A70",  # Blend of teal + amber

  "Urchins"   = "#956079",  # Blend of purple + rust
  "Kelp"      = col_taxa["M. pyrifera"]
)

# Calculate mean effect by year and trophic level (After period only)
if ("y" %in% names(All.RR.sub.trans)) {
  taxa_col <- "y"
} else if ("Taxa" %in% names(All.RR.sub.trans)) {
  taxa_col <- "Taxa"
} else {
  taxa_col <- names(All.RR.sub.trans)[1]  # fallback
}

temporal_data <- All.RR.sub.trans %>%
  filter(BA == "After", time >= 0) %>%
  mutate(
    Trophic_Level = trophic_assignment[.data[[taxa_col]]],
    time = as.numeric(time)
  ) %>%
  filter(!is.na(Trophic_Level)) %>%
  group_by(Trophic_Level, time) %>%
  summarise(
    mean_lnRR = mean(lnDiff, na.rm = TRUE),
    se_lnRR = sd(lnDiff, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(n >= 3, is.finite(mean_lnRR))

temporal_data$Trophic_Level <- factor(
  temporal_data$Trophic_Level,
  levels = c("Predators", "Urchins", "Kelp")
)

if (nrow(temporal_data) >= 10) {

  # Panel A: Time series of effect sizes by trophic level
  panel_S3A <- ggplot(temporal_data,
                       aes(x = time, y = mean_lnRR, color = Trophic_Level, fill = Trophic_Level)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    annotate("text", x = 0.5, y = Inf, label = "Increase in MPA",
             hjust = 0, vjust = 1.5, size = 2.5, color = "grey50", fontface = "italic") +
    annotate("text", x = 0.5, y = -Inf, label = "Decrease in MPA",
             hjust = 0, vjust = -0.5, size = 2.5, color = "grey50", fontface = "italic") +
    geom_ribbon(aes(ymin = mean_lnRR - 1.96 * se_lnRR,
                    ymax = mean_lnRR + 1.96 * se_lnRR),
                alpha = 0.2, color = NA) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5, shape = 21, stroke = 0.8, color = "white") +
    scale_color_manual(values = trophic_colors, name = "Trophic Level") +
    scale_fill_manual(values = trophic_colors, name = "Trophic Level") +
    scale_x_continuous(breaks = seq(0, 20, by = 5)) +
    labs(
      title = "(a) Effect sizes over time since MPA implementation",
      x = "Years since MPA implementation",
      y = "Mean log response ratio (lnRR)"
    ) +
    theme_mpa(base_size = 10) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      panel.grid.major = element_line(color = "grey95", linewidth = 0.25)
    )

  # Panel B: Cumulative mean trajectories
  cumulative_data <- All.RR.sub.trans %>%
    filter(BA == "After", time >= 0) %>%
    mutate(
      Trophic_Level = trophic_assignment[.data[[taxa_col]]],
      time = as.numeric(time)
    ) %>%
    filter(!is.na(Trophic_Level)) %>%
    arrange(Trophic_Level, time) %>%
    group_by(Trophic_Level) %>%
    mutate(cumulative_mean = cummean(lnDiff)) %>%
    ungroup() %>%
    group_by(Trophic_Level, time) %>%
    summarise(cumulative_mean = last(cumulative_mean), .groups = "drop")

  cumulative_data$Trophic_Level <- factor(
    cumulative_data$Trophic_Level,
    levels = c("Predators", "Urchins", "Kelp")
  )

  panel_S3B <- ggplot(cumulative_data,
                       aes(x = time, y = cumulative_mean, color = Trophic_Level)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_line(linewidth = 1, alpha = 0.8) +
    geom_point(size = 1.5, alpha = 0.6) +
    scale_color_manual(values = trophic_colors, name = "Trophic Level") +
    scale_x_continuous(breaks = seq(0, 20, by = 5)) +
    labs(
      title = "(b) Cumulative mean effect stabilizes over time",
      x = "Years since MPA implementation",
      y = "Cumulative mean lnRR"
    ) +
    theme_mpa(base_size = 10) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0),
      legend.position = "none",
      panel.grid.major = element_line(color = "grey95", linewidth = 0.25)
    )

  # Panel C: Rate of change by trophic level
  slope_data <- All.RR.sub.trans %>%
    filter(BA == "After", time >= 0) %>%
    mutate(
      Trophic_Level = trophic_assignment[.data[[taxa_col]]],
      time = as.numeric(time)
    ) %>%
    filter(!is.na(Trophic_Level)) %>%
    group_by(Trophic_Level, CA_MPA_Name_Short) %>%
    filter(n() >= 3) %>%
    summarise(
      slope = tryCatch(coef(lm(lnDiff ~ time))[2], error = function(e) NA_real_),
      n_years = n(),
      .groups = "drop"
    ) %>%
    filter(is.finite(slope))

  slope_data$Trophic_Level <- factor(
    slope_data$Trophic_Level,
    levels = c("Predators", "Urchins", "Kelp")
  )

  slope_summary <- slope_data %>%
    group_by(Trophic_Level) %>%
    summarise(
      mean_slope = mean(slope, na.rm = TRUE),
      se_slope = sd(slope, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )

  panel_S3C <- ggplot(slope_data, aes(x = Trophic_Level, y = slope, color = Trophic_Level)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 2.5) +
    geom_pointrange(
      data = slope_summary,
      aes(x = Trophic_Level, y = mean_slope,
          ymin = mean_slope - 1.96 * se_slope,
          ymax = mean_slope + 1.96 * se_slope),
      size = 1, linewidth = 1, color = "black",
      position = position_nudge(x = 0.25)
    ) +
    scale_color_manual(values = trophic_colors) +
    labs(
      title = "(c) Rate of effect change by trophic level",
      x = "Trophic Level",
      y = expression("Effect trajectory slope (lnRR" ~ yr^{-1} * ")")
    ) +
    theme_mpa(base_size = 10) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0),
      legend.position = "none",
      axis.text.x = element_text(face = "bold")
    )

  # Combine panels
  fig_s3 <- panel_S3A / (panel_S3B | panel_S3C) +
    plot_layout(heights = c(1.2, 1)) +
    plot_annotation(
      title = "Temporal dynamics of trophic cascade recovery",
      subtitle = "Predator effects emerge first, followed by urchin decline, then kelp recovery",
      theme = theme(
        plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 10, color = "grey40", hjust = 0.5)
      )
    )

  save_fig(fig_s3, "fig_s03_temporal_dynamics", FIG_S3_DIMS["w"], FIG_S3_DIMS["h"])

} else {
  cat("  WARNING: Not enough temporal data for Figure S3\n")
}


# =============================================================================
# Figure S4: Space-Time Effect Heatmap
# =============================================================================
# Matrix showing effect accumulation across MPAs and years

cat("\n--- Figure S4: Space-Time Heatmap ---\n")

FIG_S4_DIMS <- c(w = 20, h = 16)

# Get MPA implementation years
mpa_years <- Site %>%
  dplyr::select(CA_MPA_Name_Short, MPA_Start_Year) %>%
  distinct() %>%
  filter(!is.na(MPA_Start_Year))

# Prepare heatmap data for urchins
heatmap_data <- All.RR.sub.trans %>%
  filter(BA == "After") %>%
  mutate(
    Trophic_Level = trophic_assignment[.data[[taxa_col]]],
    time = as.numeric(time)
  ) %>%
  filter(Trophic_Level == "Urchins", time >= 0, time <= 15) %>%
  group_by(CA_MPA_Name_Short, time) %>%
  summarise(mean_lnRR = mean(lnDiff, na.rm = TRUE), .groups = "drop") %>%
  left_join(mpa_years, by = "CA_MPA_Name_Short") %>%
  filter(!is.na(MPA_Start_Year))

# Order MPAs by implementation year (newer at bottom)
mpa_order <- heatmap_data %>%
  distinct(CA_MPA_Name_Short, MPA_Start_Year) %>%
  arrange(desc(MPA_Start_Year)) %>%
  pull(CA_MPA_Name_Short)

heatmap_data$CA_MPA_Name_Short <- factor(heatmap_data$CA_MPA_Name_Short, levels = mpa_order)

if (nrow(heatmap_data) >= 20) {

  fig_s4 <- ggplot(heatmap_data, aes(x = time, y = CA_MPA_Name_Short, fill = mean_lnRR)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient2(
      low = col_taxa["M. pyrifera"],  # Green for urchin decrease (kelp good)
      mid = "white",
      high = col_taxa["S. purpuratus"],  # Purple for urchin increase
      midpoint = 0,
      name = "Urchin\nEffect\n(lnRR)",
      limits = c(-3, 3),
      oob = scales::squish
    ) +
    scale_x_continuous(breaks = seq(0, 15, by = 3), expand = c(0, 0)) +
    labs(
      title = "Urchin response to MPA protection over time",
      subtitle = "Negative values (green) = urchin decline inside MPA relative to reference",
      x = "Years since MPA implementation",
      y = "MPA (ordered by implementation year)"
    ) +
    theme_mpa(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, color = "grey40"),
      axis.text.y = element_text(size = 8),
      legend.position = "right",
      panel.grid = element_blank()
    )

  save_fig(fig_s4, "fig_s04_spacetime_heatmap", FIG_S4_DIMS["w"], FIG_S4_DIMS["h"])

} else {
  cat("  WARNING: Not enough data for Figure S4 heatmap\n")
}


# =============================================================================
# Figure S5: Model Selection & Heterogeneity
# =============================================================================
# Statistical transparency: variance components and model selection

cat("\n--- Figure S5: Statistical Transparency ---\n")

FIG_S5_DIMS <- c(w = 18, h = 12)

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

panel_S5A <- ggplot(model_dist, aes(x = Taxa, y = prop, fill = Model)) +
  geom_col(position = "stack", width = 0.7, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = col_model, name = "Best Model") +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  labs(
    title = "(a) Model selection distribution by taxa",
    subtitle = "Proportion of MPAs selecting each pBACIPS model form",
    x = NULL,
    y = "Proportion of MPAs"
  ) +
  theme_mpa(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    plot.subtitle = element_text(size = 8, color = "grey40"),
    axis.text.x = element_text(face = "italic", angle = 25, hjust = 1),
    legend.position = "right"
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

  panel_S5B <- ggplot(var_comp, aes(x = tau2, y = Response, fill = Component)) +
    geom_col(position = "dodge", width = 0.6, color = "white", linewidth = 0.3) +
    scale_fill_manual(
      values = c("MPA" = col_taxa["S. pulcher"], "Source" = col_taxa["P. interruptus"]),
      name = "Variance\nComponent"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(
      title = "(b) Meta-analysis variance components",
      subtitle = expression("Between-group heterogeneity ("*tau^2*")"),
      x = expression(tau^2),
      y = NULL
    ) +
    theme_mpa(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 8, color = "grey40"),
      legend.position = "right"
    )

  fig_s5 <- panel_S5A / panel_S5B +
    plot_layout(heights = c(1, 0.8)) +
    plot_annotation(
      title = "Statistical transparency: Model selection and heterogeneity",
      theme = theme(
        plot.title = element_text(face = "bold", size = 12, hjust = 0.5)
      )
    )

} else {
  fig_s5 <- panel_S5A +
    plot_annotation(
      title = "Model selection distribution across taxa",
      theme = theme(
        plot.title = element_text(face = "bold", size = 12, hjust = 0.5)
      )
    )
}

save_fig(fig_s5, "fig_s05_statistical_transparency", FIG_S5_DIMS["w"], FIG_S5_DIMS["h"])


# =============================================================================
# Figure 6: Trophic Cascade Flow Diagram (Conceptual)
# =============================================================================
# Sankey-style visualization showing effect flow through food web

cat("\n--- Figure 6: Trophic Cascade Flow ---\n")

FIG6_DIMS <- c(w = 18, h = 10)

# Get meta-analytic means for each taxa from Table2
if (exists("Table2")) {

  # Extract density effects for the flow diagram
  meta_means <- Table2 %>%
    filter(Response == "Density") %>%
    dplyr::select(Taxa, Estimate, CI_lower, CI_upper) %>%
    mutate(
      Trophic = case_when(
        Taxa %in% c("P. interruptus", "S. pulcher") ~ "Predators",
        Taxa %in% c("S. purpuratus", "M. franciscanus") ~ "Urchins",
        Taxa == "M. pyrifera" ~ "Kelp",
        TRUE ~ NA_character_
      ),
      Direction = ifelse(Estimate > 0, "Increase", "Decrease"),
      Significant = (CI_lower > 0) | (CI_upper < 0)
    )

  # Create conceptual flow diagram
  # Layer 1: Trophic levels as colored bars
  # Layer 2: Arrows showing direction of effects

  trophic_summary <- meta_means %>%
    group_by(Trophic) %>%
    summarise(
      mean_effect = mean(Estimate, na.rm = TRUE),
      n_taxa = n(),
      .groups = "drop"
    ) %>%
    mutate(
      x = case_when(Trophic == "Predators" ~ 1, Trophic == "Urchins" ~ 2, Trophic == "Kelp" ~ 3),
      direction_label = ifelse(mean_effect > 0, "↑ INCREASE", "↓ DECREASE"),
      y = 0
    )

  # Arrow data
  arrows_df <- data.frame(
    x_start = c(1.3, 2.3),
    x_end = c(1.7, 2.7),
    y = c(0, 0),
    label = c("Predation\npressure", "Grazing\nrelease")
  )

  fig6 <- ggplot() +
    # Background rectangles for trophic levels
    annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -0.5, ymax = 0.5,
             fill = col_taxa["P. interruptus"], alpha = 0.3) +
    annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -0.5, ymax = 0.5,
             fill = col_taxa["S. purpuratus"], alpha = 0.3) +
    annotate("rect", xmin = 2.5, xmax = 3.5, ymin = -0.5, ymax = 0.5,
             fill = col_taxa["M. pyrifera"], alpha = 0.3) +
    # Arrows between levels
    geom_segment(data = arrows_df,
                 aes(x = x_start, xend = x_end, y = y, yend = y),
                 arrow = arrow(length = unit(5, "mm"), type = "closed"),
                 linewidth = 3, color = "grey30") +
    # Arrow labels
    geom_text(data = arrows_df, aes(x = (x_start + x_end)/2, y = 0.3, label = label),
              size = 3, fontface = "italic", color = "grey30", lineheight = 0.9) +
    # Trophic level labels
    annotate("text", x = 1, y = 0.15, label = "PREDATORS", fontface = "bold", size = 4) +
    annotate("text", x = 1, y = -0.05, label = "(lobster, sheephead)", size = 3, color = "grey40") +
    annotate("text", x = 1, y = -0.25, label = "↑ INCREASE", size = 4, fontface = "bold",
             color = col_taxa["P. interruptus"]) +
    annotate("text", x = 2, y = 0.15, label = "URCHINS", fontface = "bold", size = 4) +
    annotate("text", x = 2, y = -0.05, label = "(purple, red)", size = 3, color = "grey40") +
    annotate("text", x = 2, y = -0.25, label = "↓ DECREASE", size = 4, fontface = "bold",
             color = col_taxa["S. purpuratus"]) +
    annotate("text", x = 3, y = 0.15, label = "KELP", fontface = "bold", size = 4) +
    annotate("text", x = 3, y = -0.05, label = "(M. pyrifera)", size = 3, color = "grey40", fontface = "italic") +
    annotate("text", x = 3, y = -0.25, label = "↑ INCREASE", size = 4, fontface = "bold",
             color = col_taxa["M. pyrifera"]) +
    # MPA protection arrow at left
    annotate("segment", x = 0.3, xend = 0.3, y = -0.4, yend = 0.4,
             arrow = arrow(length = unit(4, "mm"), type = "closed", ends = "first"),
             linewidth = 2, color = col_site["Inside"]) +
    annotate("text", x = 0.15, y = 0, label = "MPA\nPROTECTION",
             size = 3, fontface = "bold", color = col_site["Inside"], angle = 90) +
    coord_cartesian(xlim = c(0, 3.7), ylim = c(-0.6, 0.5)) +
    labs(
      title = "Trophic cascade mechanism in California MPA kelp forests",
      subtitle = "MPA protection triggers a three-level cascade: predator recovery → urchin decline → kelp recovery"
    ) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5, margin = margin(b = 5)),
      plot.subtitle = element_text(size = 10, color = "grey40", hjust = 0.5, margin = margin(b = 10))
    )

  save_fig(fig6, "fig_06_trophic_flow", FIG6_DIMS["w"], FIG6_DIMS["h"])

} else {
  cat("  WARNING: Table2 not found - skipping Figure 6\n")
}

} # End if(FALSE) - exploratory figures not included in manuscript

cat("\n=== All required manuscript figures saved to:", here::here("plots"), "===\n")
