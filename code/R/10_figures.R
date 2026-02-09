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
# Journal max widths: single column = 85mm, double column = 178mm
FIG_WIDTH_SINGLE <- 12   # cm, for single-column figures
FIG_WIDTH_DOUBLE <- 17   # cm, for double-column figures
FIG_WIDTH_WIDE   <- 18   # cm, for figures needing extra width
FIG_WIDTH_SUPP   <- 17.8 # cm, Conservation Letters max width for supplemental figures

# Figure-specific dimensions (width, height in cm)
# Note: Figure 1 dimensions are defined in analysis/R/fig01_map.R
FIG2_DIMS <- c(w = 17, h = 7.5)  # Data processing pipeline — double-column width
FIG3_DIMS <- c(w = 14, h = 10)   # Mean effects — single-column+
FIG4_DIMS <- c(w = 17, h = 9)    # Trophic cascade scatters (2-panel: pred→urchin | urchin→kelp)
FIG_S1_DIMS <- c(w = 17.8, h = 22) # Forest plot (supplemental)
FIG_S2_DIMS <- c(w = 17.8, h = 22) # All taxa time series (supplemental, max CL width)

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
save_fig <- function(plot, name, w, h, dpi = 300) {
  pdf_path <- here::here("plots", paste0(name, ".pdf"))
  png_path <- here::here("plots", paste0(name, ".png"))

  # Detect if this is a patchwork plot (complex multi-panel)
  is_patchwork <- inherits(plot, "patchwork") || inherits(plot, "gg") && !is.null(plot$patches)

  # Save PDF with error handling and fallback strategies
  pdf_success <- FALSE

  # Strategy 1: For patchwork plots, use pdf() device (most reliable)
  if (is_patchwork && !pdf_success) {
    pdf_success <- tryCatch({
      pdf(pdf_path, width = w / 2.54, height = h / 2.54, bg = "white")
      print(plot)
      dev.off()
      TRUE
    }, error = function(e) {
      if (dev.cur() > 1) dev.off()  # Close any open devices
      FALSE
    })
  }

  # Strategy 2: Try cairo_pdf via ggsave
  if (!pdf_success && capabilities("cairo")) {
    pdf_success <- tryCatch({
      ggsave(pdf_path, plot, width = w, height = h, units = "cm",
             device = cairo_pdf, bg = "white", limitsize = FALSE)
      TRUE
    }, error = function(e) { FALSE })
  }

  # Strategy 3: Try standard pdf device via ggsave
  if (!pdf_success) {
    pdf_success <- tryCatch({
      ggsave(pdf_path, plot, width = w, height = h, units = "cm",
             device = "pdf", bg = "white", limitsize = FALSE)
      TRUE
    }, error = function(e) { FALSE })
  }

  # Strategy 4: Last resort - use pdf() device directly
  if (!pdf_success) {
    pdf_success <- tryCatch({
      pdf(pdf_path, width = w / 2.54, height = h / 2.54, bg = "white")
      print(plot)
      dev.off()
      TRUE
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
    cat(sprintf("  ✓ Saved: %s (PDF: %s KB, PNG: %s KB @ %d DPI)\n",
                name, pdf_size, png_size, dpi))
  } else if (!png_exists) {
    # PNG is critical - always fail if PNG doesn't exist
    stop("CRITICAL: Failed to create PNG for ", name,
         " at ", png_path, "\n",
         "  Check disk space and write permissions.")
  } else if (!pdf_exists) {
    # PDF failed but PNG succeeded - log warning and continue
    png_size <- format(file.size(png_path) / 1024, digits = 1, nsmall = 1)
    cat(sprintf("  ⚠ Saved: %s (PNG: %s KB @ %d DPI) - PDF generation failed\n",
                name, png_size, dpi))
    warning("PDF generation failed for ", name, " but PNG was created successfully. ",
            "This is a known issue with some complex ggplot2 figures. ",
            "PNG can be converted to PDF externally if needed.")
  }

  invisible(list(pdf = pdf_path, png = png_path))
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
  FIG1_DIMS <- c(w = 17, h = 16)  # cm — Conservation Letters double-column max (17cm)

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
    geom_raster(data = bathy_ocean, aes(x = lon, y = lat, fill = depth),
                interpolate = TRUE, alpha = 0.85) +
    scale_fill_gradientn(
      colors = ocean_colors,
      values = scales::rescale(c(-2000, -500, -100, 0)),
      name = "Depth (m)",
      breaks = c(-1500, -1000, -500, -100),
      labels = c("1500", "1000", "500", "100"),
      limits = c(min(bathy_ocean$depth, na.rm = TRUE), 0),
      guide = guide_colorbar(barwidth = 0.5, barheight = 3, title.position = "top",
                             title.hjust = 0.5, frame.colour = "grey50", ticks.colour = "grey50",
                             order = 3)
    ) +
    geom_contour(data = bathy_ocean, aes(x = lon, y = lat, z = depth),
                 breaks = c(-100, -200, -500, -1000), color = "white",
                 linewidth = 0.12, alpha = 0.25) +
    geom_sf(data = mpa, fill = fig1_map_colors$mpa_fill, color = fig1_map_colors$mpa_border,
            alpha = 0.22, linewidth = 0.4, inherit.aes = FALSE) +
    new_scale_fill() +
    geom_spatraster(data = hillshade, maxcell = 5e5) +
    scale_fill_gradientn(colors = terrain_pal, na.value = NA, guide = "none") +
    geom_sf(data = coast, fill = alpha("#F5F0E1", 0.15), color = NA, inherit.aes = FALSE) +
    geom_sf(data = coast, fill = NA, color = fig1_map_colors$coastline,
            linewidth = 0.4, inherit.aes = FALSE) +
    geom_segment(data = sites_pairs, aes(x = x_in, y = y_in, xend = x_out, yend = y_out),
                 color = "grey65", linewidth = 0.45, alpha = 0.7) +
    new_scale_fill() +
    geom_point(data = fig1_sites, aes(x = Lon, y = Lat, fill = status, shape = program),
               size = 3.6, color = "white", stroke = 0.9) +
    scale_fill_manual(name = "Inside vs Outside MPA", values = fig1_status_colors,
                      guide = guide_legend(order = 1, override.aes = list(shape = 21, size = 3))) +
    scale_shape_manual(name = "Monitoring program", values = program_shapes,
                       guide = guide_legend(order = 2, override.aes = list(fill = "grey60", size = 3))) +
    # Label only the panel sites with (a)–(d) at the true site location (not offset points)
    geom_text(data = sites_labels, aes(x = Lon, y = Lat - 0.04, label = paste0("(", panel_letter, ")")),
              size = 3.3, fontface = "plain", color = "#1A1A1A") +
    # Compact key so readers can map (a)-(d) to site names without cluttering the coast.
    annotate(
      "label",
      # Put the (a)-(d) key in the top-left so it never competes with the inset
      # (bottom-right) and keeps the bottom-left free for scale + compass.
      x = BBOX_LONLAT["xmin"] + 0.10, y = BBOX_LONLAT["ymax"] - 0.06,
      label = paste(panel_label_df$site_label, collapse = "\n"),
      hjust = 0, vjust = 1,
      size = 2.8, label.size = 0.25,
      label.padding = unit(0.18, "lines"),
      fill = scales::alpha("white", 0.85), color = "grey40"
    ) +
    coord_sf(xlim = c(BBOX_LONLAT["xmin"], BBOX_LONLAT["xmax"]),
             ylim = c(BBOX_LONLAT["ymin"], BBOX_LONLAT["ymax"]), expand = FALSE, crs = 4326) +
    annotation_scale(location = "bl", width_hint = 0.2, pad_x = unit(0.3, "in"),
                     pad_y = unit(0.3, "in"), style = "ticks", text_cex = 0.75, line_width = 0.4) +
    # Compass moved away from legends and inset; lift it above the scalebar.
    annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.3, "in"),
                           pad_y = unit(1.1, "in"), style = north_arrow_minimal,
                           height = unit(0.7, "cm"), width = unit(0.7, "cm")) +
    theme_mpa(base_size = 10) +
    theme(
      panel.background = element_rect(fill = "#C6DBEF", color = NA),
      panel.border = element_rect(fill = NA, color = "grey35", linewidth = 0.4),
      axis.line.x.bottom = element_blank(),
      axis.line.y.left = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(size = 8, color = "grey30"),
      # Put legends top-right and keep them aligned as a single vertical box.
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.title = element_text(size = 8, face = "plain"),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.35, "cm"),
      legend.spacing.x = unit(2, "mm"),
      legend.spacing.y = unit(1, "mm"),
      legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
      legend.margin = margin(2, 4, 2, 4),
      legend.box.background = element_blank(),
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
      panel.background = element_rect(fill = scales::alpha("white", 0.85), color = "grey60", linewidth = 0.35),
      plot.margin = margin(2, 2, 2, 2)
    )

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
      return(ggplot() + theme_void() +
        labs(tag = paste0("(", letter, ")")))
    }

    ggplot(site_data, aes(x = year, y = mean_value, color = Status)) +
      geom_vline(xintercept = mpa_year, linetype = "dashed", color = "grey40", linewidth = 0.7) +
      geom_line(linewidth = 0.6) +
      geom_point(size = 1.2, alpha = 0.7, shape = 21, fill = "white", stroke = 0.5) +
      scale_color_manual(values = TS_COLORS, name = NULL) +
      scale_x_continuous(breaks = seq(1990, 2020, by = 10), limits = c(1988, 2024)) +
      scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
      labs(tag = paste0("(", letter, ")"),
           title = short_name, x = NULL,
           y = expression(Kelp~biomass~(g~m^{-2}))) +
      theme_mpa(base_size = 8) +
      theme(panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_blank(),
            plot.title = element_text(size = 8, face = "plain", hjust = 0.5,
                                      margin = margin(b = 2)),
            plot.tag = element_text(size = 9, face = "plain"),
            plot.tag.position = c(0.02, 0.98),
            axis.title = element_text(size = 7.5, color = "grey20"),
            axis.title.y = element_text(margin = margin(r = 3)),
            axis.text = element_text(size = 6.5, color = "grey30"),
            legend.position = "bottom",
            legend.text = element_text(size = 7),
            legend.key.width = unit(0.5, "cm"),
            legend.margin = margin(0, 0, 0, 0),
            plot.margin = margin(2, 3, 2, 3))
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
    # Reduce bottom margin of main map
    main_map <- main_map + theme(plot.margin = margin(2, 2, 0, 2))

    # Create row of panels using patchwork with collected legend
    panels_row <- (fig1_panels$a + fig1_panels$b + fig1_panels$c + fig1_panels$d) +
      plot_layout(nrow = 1, guides = "collect") &
      theme(legend.position = "bottom",
            legend.justification = "center",
            panel.border = element_blank(),
            axis.line = element_blank(),
            axis.line.x.bottom = element_line(colour = "black", linewidth = 0.3),
            axis.line.y.left = element_line(colour = "black", linewidth = 0.3))

    # Stack map on top, panels below — give panels more vertical space
    fig1 <- main_map / panels_row +
      plot_layout(heights = c(2, 1.1))
  } else {
    fig1 <- main_map
  }

  fig1 <- fig1 + plot_annotation(theme = theme(plot.margin = margin(2, 4, 4, 4)))

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
fig2_x_breaks <- seq(fig2_x_limits[1], fig2_x_limits[2], by = 20)
fig2_x_breaks <- fig2_x_breaks[fig2_x_breaks >= fig2_x_limits[1] & fig2_x_breaks <= fig2_x_limits[2]]

# Build each panel individually for better control, then combine with patchwork
# Common theme settings for consistent panel appearance
fig2_theme <- theme_mpa(base_size = 10) +
  theme(legend.position = "none",
        plot.title = element_text(size = 9.5, face = "plain", hjust = 0,
                                  margin = margin(b = 8)),
        axis.title.y = element_text(size = 9, margin = margin(r = 4)),
        axis.text = element_text(size = 8, color = "black"),
        plot.margin = margin(4, 4, 4, 4))

p2a <- ggplot(fig2_raw, aes(x = year, y = value, color = status)) +
  add_mpa_vline(scorpion_start) +
  # TUFTE: Annotate the MPA implementation line
  annotate("text", x = scorpion_start, y = Inf,
           label = "MPA\nestablished", hjust = -0.1, vjust = 1.8,
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
        alpha = 0.25, linewidth = 1.0,
        inherit.aes = FALSE
      )
    }
  } +
  scale_shape_manual(name = "Period", values = c("Before" = 1, "After" = 16)) +
  # TUFTE: Concise title with finding
  labs(title = "(d) Declining trend", x = NULL, y = NULL) +
  scale_x_continuous(breaks = fig2_x_breaks, limits = fig2_x_limits) +
  fig2_theme

# Enable legends on first panel of each type for collection
p2a <- p2a + theme(legend.position = "bottom")
p2c <- p2c + theme(legend.position = "bottom")

# Combine 4 data panels directly — no arrow spacers
fig2_final <- (p2a + p2b + p2c + p2d) +
  plot_layout(ncol = 4, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(face = "plain", size = 9, color = "black"),
    legend.text = element_text(size = 8, color = "black"),
    legend.spacing.x = unit(6, "mm"),
    legend.margin = margin(t = 2),
    # Reinforce L-shaped axes through patchwork & operator
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.line.x.bottom = element_line(colour = "black", linewidth = 0.3),
    axis.line.y.left = element_line(colour = "black", linewidth = 0.3)
  )

# Conservation Letters: max 170mm double-column width
save_fig(fig2_final, "fig_02_data_processing", FIG2_DIMS["w"], FIG2_DIMS["h"])

# =============================================================================
# Figure S1 (Supplemental): Forest plot of effect sizes by MPA and taxa
# Manuscript: Supplemental Figure S1
# =============================================================================

cat("Building Figure S1 (Supplemental): Forest plot...\n")

# Filter SumStats.Final for forest plot
excluded_mpas <- c("Painted Cave SMCA", "San Miguel Island SC",
                   "Arrow Point to Lion Head Point SMCA",
                   "Judith Rk SMR", "Point Conception SMR")

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
x_breaks <- seq(-10, 10, by = 5)
x_breaks <- x_breaks[abs(x_breaks) <= x_limit]

fig_s1 <- ggplot(fig_s1_data,
               aes(x = MPA_order, y = Mean, ymin = Mean - CI, ymax = Mean + CI,
                   color = Resp, shape = Source)) +
  # Reference line at zero
  geom_hline(yintercept = 0, linetype = "solid", color = "grey40",
             linewidth = 0.4) +
  # 95% CIs
  geom_linerange(linewidth = 0.6, position = position_dodge(width = 0.55)) +
  # Effect size points
  geom_point(size = 2.5, position = position_dodge(width = 0.55)) +
  facet_wrap(~ Taxa, ncol = 2, scales = "free_x") +
  scale_color_response(name = "Response",
                       labels = c("Den" = "Density", "Bio" = "Biomass")) +
  scale_shape_source(name = "Source") +
  scale_y_continuous(breaks = x_breaks) +
  coord_flip(ylim = c(-x_limit, x_limit), clip = "off") +
  labs(x = NULL, y = "Effect size (lnRR)") +
  theme_mpa(base_size = 9) +
  theme(
    strip.text = element_text(face = "italic", size = 10,
                              margin = margin(3, 0, 3, 0)),
    strip.background = element_blank(),
    # Larger MPA name font for readability
    axis.text.y = element_text(size = 8, color = "grey15", hjust = 1),
    axis.text.x = element_text(size = 8),
    axis.title.x = element_text(size = 9, margin = margin(t = 6)),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.box.spacing = unit(0.5, "cm"),
    legend.spacing.x = unit(0.3, "cm"),
    legend.title = element_text(face = "plain", size = 8.5),
    legend.text = element_text(size = 8),
    legend.margin = margin(t = 4, b = 2),
    legend.key.width = unit(0.5, "cm"),
    legend.box.margin = margin(0, 0, 0, 0),
    panel.spacing.x = unit(1, "lines"),
    panel.spacing.y = unit(0.4, "lines"),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(6, 10, 6, 6)
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
pd <- position_dodge(width = 0.75)

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
  # Reference line at zero (no effect)
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3, linetype = "dashed") +
  # Layer 1: Individual MPA effect sizes (background)
  geom_point(data = fig3_individual,
             aes(x = Taxa, y = Mean, color = Response),
             position = position_jitterdodge(jitter.width = 0.2,
                                              dodge.width = 0.6,
                                              seed = 42),
             size = 1.8, alpha = 0.30, shape = 16) +
  # Layer 2: White halo behind diamonds
  geom_point(data = fig3_meta,
             aes(x = Taxa, y = Estimate),
             position = pd, size = 7, shape = 18, color = "white") +
  # Layer 3: 95% CIs
  geom_errorbar(data = fig3_meta,
                aes(x = Taxa, ymin = CI_lower, ymax = CI_upper, color = Response),
                position = pd, width = 0.12, linewidth = 1.0) +
  # Layer 4: Meta-analytic means (diamonds)
  geom_point(data = fig3_meta,
             aes(x = Taxa, y = Estimate, color = Response),
             position = pd, size = 6.0, shape = 18) +
  # Sample size labels below x-axis (wider dodge to prevent overlap)
  geom_text(data = fig3_meta,
            aes(x = Taxa, y = fig3_label_y, label = paste0("n=", n), color = Response),
            position = position_dodge(width = 0.75), size = 2.8, show.legend = FALSE) +
  scale_color_manual(name = NULL, values = col_response_long,
                     guide = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  coord_cartesian(ylim = c(fig3_y_min - 1.2, fig3_y_max + 0.6)) +
  labs(x = NULL, y = "Effect size (lnRR)") +
  theme_mpa(base_size = 10) +
  theme(
    axis.text.x = element_text(face = "italic", size = 9),
    axis.title.y = element_text(size = 10),
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    legend.key.width = unit(0.8, "cm"),
    legend.margin = margin(t = 2),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(8, 8, 6, 6)
  )

save_fig(fig3, "fig_03_mean_effects", FIG3_DIMS["w"], FIG3_DIMS["h"])

# =============================================================================
# Figure 4: Trophic cascade scatterplots (4-panel: biomass top, density bottom)
# Raw annual data (background) + MPA-level meta-analytic means (foreground)
# Top row (biomass): (a) Predator → Urchin  (b) Urchin → Kelp
# Bottom row (density): (c) Predator → Urchin  (d) Urchin → Kelp
# =============================================================================

cat("Building Figure 4: Trophic cascade scatterplots (4-panel)...\n")

# ---------------------------------------------------------------------------
# Prepare MPA-level summary data from SumStats.Final (for foreground overlay)
# ---------------------------------------------------------------------------

# BIOMASS panels (top row)
# Panel A: Predator biomass → Urchin biomass
fig4_mpa_pred_bio <- SumStats.Final %>%
  dplyr::filter(
    Taxa %in% c("Panulirus interruptus", "P. interruptus",
                "Semicossyphus pulcher", "S. pulcher"),
    Resp %in% c("Bio", "Biomass"),
    AnalysisType %in% c("pBACIPS", "CI")
  ) %>%
  dplyr::mutate(Mean = as.numeric(Mean), SE = as.numeric(SE)) %>%
  dplyr::group_by(MPA) %>%
  dplyr::summarise(
    lnRR_predator = mean(Mean, na.rm = TRUE),
    SE_predator   = sqrt(sum(SE^2, na.rm = TRUE)) / dplyr::n(),
    .groups = "drop"
  )

fig4_mpa_urch_bio <- SumStats.Final %>%
  dplyr::filter(
    Taxa %in% c("Strongylocentrotus purpuratus", "S. purpuratus",
                "Mesocentrotus franciscanus", "M. franciscanus"),
    Resp %in% c("Bio", "Biomass"),
    AnalysisType %in% c("pBACIPS", "CI")
  ) %>%
  dplyr::mutate(Mean = as.numeric(Mean), SE = as.numeric(SE)) %>%
  dplyr::group_by(MPA) %>%
  dplyr::summarise(
    lnRR_urchin = mean(Mean, na.rm = TRUE),
    SE_urchin   = mean(SE, na.rm = TRUE),
    .groups = "drop"
  )

fig4_mpa_A <- dplyr::inner_join(fig4_mpa_pred_bio, fig4_mpa_urch_bio, by = "MPA")

# Panel B: Urchin biomass → Kelp biomass
fig4_mpa_kelp_bio <- SumStats.Final %>%
  dplyr::filter(
    Taxa %in% c("Macrocystis pyrifera", "M. pyrifera"),
    Resp %in% c("Bio", "Biomass"),
    AnalysisType %in% c("pBACIPS", "CI")
  ) %>%
  dplyr::mutate(Mean = as.numeric(Mean), SE = as.numeric(SE)) %>%
  dplyr::group_by(MPA) %>%
  dplyr::summarise(
    lnRR_kelp = mean(Mean, na.rm = TRUE),
    SE_kelp   = mean(SE, na.rm = TRUE),
    .groups = "drop"
  )

fig4_mpa_B <- dplyr::inner_join(fig4_mpa_urch_bio, fig4_mpa_kelp_bio, by = "MPA")

# DENSITY panels (bottom row)
# Panel C: predator and urchin density means per MPA
fig4_mpa_pred <- SumStats.Final %>%
  dplyr::filter(
    Taxa %in% c("Panulirus interruptus", "P. interruptus",
                "Semicossyphus pulcher", "S. pulcher"),
    Resp %in% c("Den", "Density"),
    AnalysisType %in% c("pBACIPS", "CI")
  ) %>%
  dplyr::mutate(Mean = as.numeric(Mean), SE = as.numeric(SE)) %>%
  dplyr::group_by(MPA) %>%
  dplyr::summarise(
    lnRR_predator = mean(Mean, na.rm = TRUE),
    SE_predator   = sqrt(sum(SE^2, na.rm = TRUE)) / dplyr::n(),
    .groups = "drop"
  )

fig4_mpa_urch <- SumStats.Final %>%
  dplyr::filter(
    Taxa %in% c("Strongylocentrotus purpuratus", "S. purpuratus"),
    Resp %in% c("Den", "Density"),
    AnalysisType %in% c("pBACIPS", "CI")
  ) %>%
  dplyr::mutate(Mean = as.numeric(Mean), SE = as.numeric(SE)) %>%
  dplyr::group_by(MPA) %>%
  dplyr::summarise(
    lnRR_urchin = mean(Mean, na.rm = TRUE),
    SE_urchin   = mean(SE, na.rm = TRUE),
    .groups = "drop"
  )

fig4_mpa_C <- dplyr::inner_join(fig4_mpa_pred, fig4_mpa_urch, by = "MPA")

# Panel D: urchin density and kelp biomass means per MPA
# Note: Kelp has no density metric, so this uses biomass (same as panel B but with density urchins)
fig4_mpa_D <- dplyr::inner_join(fig4_mpa_urch, fig4_mpa_kelp_bio, by = "MPA")

cat("  MPA-level overlay: panel A =", nrow(fig4_mpa_A),
    "MPAs, panel B =", nrow(fig4_mpa_B),
    "MPAs, panel C =", nrow(fig4_mpa_C),
    "MPAs, panel D =", nrow(fig4_mpa_D), "MPAs\n")

# ---------------------------------------------------------------------------
# Prepare raw annual data (for background points)
# ---------------------------------------------------------------------------

# BIOMASS raw data (top row)
# Panel A raw: predator biomass
figA_raw_pred <- All.RR.sub.trans %>%
  dplyr::filter(
    y %in% c("Panulirus interruptus", "P. interruptus",
             "Semicossyphus pulcher", "S. pulcher"),
    resp %in% c("Bio", "Biomass")
  ) %>%
  dplyr::group_by(CA_MPA_Name_Short, year, source) %>%
  dplyr::summarise(lnRR_predator = mean(lnDiff, na.rm = TRUE), .groups = "drop")

figA_raw_urch <- All.RR.sub.trans %>%
  dplyr::filter(
    y %in% c("Strongylocentrotus purpuratus", "S. purpuratus",
             "Mesocentrotus franciscanus", "M. franciscanus"),
    resp %in% c("Bio", "Biomass")
  ) %>%
  dplyr::select(CA_MPA_Name_Short, year, source, lnDiff) %>%
  dplyr::rename(lnRR_urchin = lnDiff)

figA_raw <- dplyr::inner_join(
  figA_raw_pred, figA_raw_urch,
  by = c("CA_MPA_Name_Short", "year", "source")
)

# Panel B raw: urchin biomass vs kelp biomass
figB_raw_urch <- All.RR.sub.trans %>%
  dplyr::filter(
    y %in% c("Strongylocentrotus purpuratus", "S. purpuratus",
             "Mesocentrotus franciscanus", "M. franciscanus", "STRPURAD"),
    resp %in% c("Bio", "Biomass")
  ) %>%
  dplyr::select(CA_MPA_Name_Short, year, source, lnDiff) %>%
  dplyr::rename(lnRR_urchin = lnDiff)

figB_raw_kelp <- All.RR.sub.trans %>%
  dplyr::filter(
    y %in% c("Macrocystis pyrifera", "M. pyrifera", "MACPYRAD"),
    resp %in% c("Bio", "Biomass")
  ) %>%
  dplyr::select(CA_MPA_Name_Short, year, source, lnDiff) %>%
  dplyr::rename(lnRR_kelp = lnDiff)

figB_raw <- dplyr::inner_join(
  figB_raw_urch, figB_raw_kelp,
  by = c("CA_MPA_Name_Short", "year", "source")
)

# DENSITY raw data (bottom row)
# Panel C raw: predator density (lobster + sheephead combined per MPA-year)
figC_raw_pred <- All.RR.sub.trans %>%
  dplyr::filter(
    y %in% c("Panulirus interruptus", "P. interruptus",
             "Semicossyphus pulcher", "S. pulcher"),
    resp %in% c("Den", "Density")
  ) %>%
  dplyr::group_by(CA_MPA_Name_Short, year, source) %>%
  dplyr::summarise(lnRR_predator = mean(lnDiff, na.rm = TRUE), .groups = "drop")

figC_raw_urch <- All.RR.sub.trans %>%
  dplyr::filter(
    y %in% c("Strongylocentrotus purpuratus", "S. purpuratus"),
    resp %in% c("Den", "Density")
  ) %>%
  dplyr::select(CA_MPA_Name_Short, year, source, lnDiff) %>%
  dplyr::rename(lnRR_urchin = lnDiff)

figC_raw <- dplyr::inner_join(
  figC_raw_pred, figC_raw_urch,
  by = c("CA_MPA_Name_Short", "year", "source")
)

# Panel 4b raw: urchin density vs kelp biomass
figD_raw_urch <- All.RR.sub.trans %>%
  dplyr::filter(
    y %in% c("Strongylocentrotus purpuratus", "S. purpuratus", "STRPURAD"),
    resp %in% c("Den", "Density")
  ) %>%
  dplyr::select(CA_MPA_Name_Short, year, source, lnDiff) %>%
  dplyr::rename(lnRR_urchin = lnDiff)

figD_raw_kelp <- All.RR.sub.trans %>%
  dplyr::filter(
    y %in% c("Macrocystis pyrifera", "M. pyrifera", "MACPYRAD"),
    resp %in% c("Bio", "Biomass")
  ) %>%
  dplyr::select(CA_MPA_Name_Short, year, source, lnDiff) %>%
  dplyr::rename(lnRR_kelp = lnDiff)

if (nrow(figD_raw_kelp) == 0) {
  stop("Figure 4b: Kelp biomass data required but unavailable in All.RR.sub.trans.")
}

figD_raw <- dplyr::inner_join(
  figD_raw_urch, figD_raw_kelp,
  by = c("CA_MPA_Name_Short", "year", "source")
)

cat("  Raw annual data: panel A (bio pred→urch) =", nrow(figA_raw),
    "obs, panel B (bio urch→kelp) =", nrow(figB_raw),
    "obs, panel C (den pred→urch) =", nrow(figC_raw),
    "obs, panel D (den urch→kelp) =", nrow(figD_raw), "obs\n")

# Diagnostic: Show data ranges for all 4 panels
if (nrow(figA_raw) > 0) {
  cat("    Panel A (biomass) ranges: Predator [",
      paste(round(range(figA_raw$lnRR_predator, na.rm=TRUE), 1), collapse=" to "),
      "], Urchin [",
      paste(round(range(figA_raw$lnRR_urchin, na.rm=TRUE), 1), collapse=" to "),
      "]\n", sep="")
}
if (nrow(figB_raw) > 0) {
  cat("    Panel B (biomass) ranges: Urchin [",
      paste(round(range(figB_raw$lnRR_urchin, na.rm=TRUE), 1), collapse=" to "),
      "], Kelp [",
      paste(round(range(figB_raw$lnRR_kelp, na.rm=TRUE), 1), collapse=" to "),
      "]\n", sep="")
}
if (nrow(figC_raw) > 0) {
  cat("    Panel C (density) ranges: Predator [",
      paste(round(range(figC_raw$lnRR_predator, na.rm=TRUE), 1), collapse=" to "),
      "], Urchin [",
      paste(round(range(figC_raw$lnRR_urchin, na.rm=TRUE), 1), collapse=" to "),
      "]\n", sep="")
}
if (nrow(figD_raw) > 0) {
  cat("    Panel D (density) ranges: Urchin [",
      paste(round(range(figD_raw$lnRR_urchin, na.rm=TRUE), 1), collapse=" to "),
      "], Kelp [",
      paste(round(range(figD_raw$lnRR_kelp, na.rm=TRUE), 1), collapse=" to "),
      "]\n", sep="")
}

# ---------------------------------------------------------------------------
# Helper function to create a trophic cascade panel
# ---------------------------------------------------------------------------
create_cascade_panel <- function(raw_data, mpa_data,
                                  x_var, y_var, x_lab, y_lab,
                                  point_color, panel_label) {
  # Check if sufficient data
  if (nrow(raw_data) < 5 || nrow(mpa_data) < 3) {
    return(ggplot() + theme_void() +
      annotate("text", x = 0.5, y = 0.5, label = "Insufficient data", size = 4))
  }

  # Correlation on MPA-level means
  cor_test <- cor.test(mpa_data[[x_var]], mpa_data[[y_var]], method = "pearson")
  cor_label <- sprintf("r = %.2f, p %s\nn = %d MPAs",
                       cor_test$estimate,
                       ifelse(cor_test$p.value < 0.001, "< 0.001",
                              sprintf("= %.3f", cor_test$p.value)),
                       nrow(mpa_data))

  # Symmetric axis limits based on raw data range with 15% buffer
  max_abs_x <- max(abs(raw_data[[x_var]]), na.rm = TRUE) * 1.15
  max_abs_y <- max(abs(raw_data[[y_var]]), na.rm = TRUE) * 1.15
  lim_x <- c(-max_abs_x, max_abs_x)
  lim_y <- c(-max_abs_y, max_abs_y)

  # Create plot
  ggplot() +
    # Reference lines
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.3) +
    # Background: raw annual observations
    geom_point(data = raw_data, aes(x = .data[[x_var]], y = .data[[y_var]]),
               position = position_jitter(width = 0.08, height = 0.08, seed = 42),
               size = 1.3, alpha = 0.20, color = point_color, stroke = 0) +
    # Foreground: MPA-level means with error bars
    geom_errorbar(data = mpa_data,
                  aes(x = .data[[x_var]],
                      ymin = .data[[y_var]] - .data[[paste0("SE_", sub("lnRR_", "", y_var))]],
                      ymax = .data[[y_var]] + .data[[paste0("SE_", sub("lnRR_", "", y_var))]]),
                  width = 0, linewidth = 0.4, color = "grey30", alpha = 0.6) +
    geom_errorbarh(data = mpa_data,
                   aes(y = .data[[y_var]],
                       xmin = .data[[x_var]] - .data[[paste0("SE_", sub("lnRR_", "", x_var))]],
                       xmax = .data[[x_var]] + .data[[paste0("SE_", sub("lnRR_", "", x_var))]]),
                   height = 0, linewidth = 0.4, color = "grey30", alpha = 0.6) +
    geom_point(data = mpa_data, aes(x = .data[[x_var]], y = .data[[y_var]]),
               size = 3.5, alpha = 0.90, color = point_color, shape = 16) +
    # Regression
    geom_smooth(data = mpa_data, aes(x = .data[[x_var]], y = .data[[y_var]]),
                method = "lm", se = TRUE, color = point_color, fill = point_color,
                linewidth = 1.2, alpha = 0.20) +
    # Correlation label
    annotate("text", x = lim_x[2], y = lim_y[2], label = cor_label,
             hjust = 1.05, vjust = 1.1, size = 3.5, color = "grey25") +
    labs(x = x_lab, y = y_lab) +
    coord_fixed(ratio = 1, xlim = lim_x, ylim = lim_y, expand = FALSE) +
    theme_mpa(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      plot.margin = margin(4, 6, 2, 6)
    )
}

# ---------------------------------------------------------------------------
# Build all 4 panels
# ---------------------------------------------------------------------------

cat("  Building Panel A (biomass): Predator vs Urchin...\n")
panel_A <- create_cascade_panel(
  figA_raw, fig4_mpa_A,
  "lnRR_predator", "lnRR_urchin",
  "Predator biomass effect (lnRR)",
  "Urchin biomass effect (lnRR)",
  col_taxa["P. interruptus"],
  "A"
)

cat("  Building Panel B (biomass): Urchin vs Kelp...\n")
panel_B <- create_cascade_panel(
  figB_raw, fig4_mpa_B,
  "lnRR_urchin", "lnRR_kelp",
  "Urchin biomass effect (lnRR)",
  "Kelp biomass effect (lnRR)",
  col_taxa["S. purpuratus"],
  "B"
)

cat("  Building Panel C (density): Predator vs Urchin...\n")
panel_C <- create_cascade_panel(
  figC_raw, fig4_mpa_C,
  "lnRR_predator", "lnRR_urchin",
  "Predator density effect (lnRR)",
  expression(italic("S. purpuratus") * " density effect (lnRR)"),
  col_taxa["P. interruptus"],
  "C"
)

cat("  Building Panel D (density): Urchin vs Kelp...\n")
panel_D <- create_cascade_panel(
  figD_raw, fig4_mpa_D,
  "lnRR_urchin", "lnRR_kelp",
  expression(italic("S. purpuratus") * " density effect (lnRR)"),
  "Kelp biomass effect (lnRR)",
  col_taxa["S. purpuratus"],
  "D"
)

# ---------------------------------------------------------------------------
# Old Panel 4a code - DELETE THIS SECTION
# ---------------------------------------------------------------------------


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

# Update dimensions for 4-panel figure
FIG4_NEW_DIMS <- c(w = 17, h = 17)  # Square-ish for 2x2 layout
save_fig(fig4, "fig_04_trophic_scatter", FIG4_NEW_DIMS["w"], FIG4_NEW_DIMS["h"])

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

  # --- Fix duplicate "Source" legend in patchwork ---
  # Patchwork's guides="collect" only merges legends whose key glyphs are
  # identical.  When panels contain different subsets of source levels, the
  # rendered glyphs differ and two "Source" rows appear.  We inject one
  # invisible dummy row per missing source level so every panel draws the
  # same four shapes in its legend, making them merge correctly.
  present_sources <- unique(as.character(d$source))
  missing_sources <- setdiff(source_levels, present_sources)
  if (length(missing_sources) > 0) {
    # Use the first species_short level present as a placeholder (its color
    # won't matter because alpha = 0).
    first_species <- levels(d$species_short)[1]
    dummy <- tibble::tibble(
      year           = NA_real_,
      lnDiff         = NA_real_,
      source         = factor(missing_sources, levels = source_levels),
      species_short  = factor(first_species, levels = levels(d$species_short)),
      CA_MPA_Name_Short = mpa_name
    )
    d <- dplyr::bind_rows(d, dummy)
  }
  mpa_start_yr <- mpa_starts[mpa_name]

  # Clean MPA name for display using standardized function
  mpa_display <- shorten_mpa_name(mpa_name)

  year_rng <- range(d$year, na.rm = TRUE)
  if (!all(is.finite(year_rng))) year_rng <- c(2000, 2022)
  # Use pretty breaks to avoid x-axis label crowding on wide time ranges
  # (e.g., Scorpion panel spans ~45 years; fixed 5-year intervals crowd)
  year_span <- diff(year_rng)
  by_val <- if (year_span > 30) 10 else 5
  year_breaks <- seq(floor(year_rng[1] / by_val) * by_val,
                     ceiling(year_rng[2] / by_val) * by_val, by = by_val)

  y_annot <- suppressWarnings(max(d$lnDiff, na.rm = TRUE))
  if (!is.finite(y_annot)) y_annot <- 0
  y_annot <- min(1.8, y_annot)

  ggplot(d, aes(x = year, y = lnDiff,
                color = species_short, shape = source)) +
    # Reference line at zero
    geom_hline(yintercept = 0, linetype = "solid", color = "grey50",
               linewidth = 0.4) +
    # MPA implementation vertical line
    add_mpa_vline(mpa_start_yr) +
    annotate("text", x = mpa_start_yr + 0.5, y = y_annot,
             label = "MPA\nStart", hjust = 0, size = 2.3, color = MPA_LABEL_COLOR,
             lineheight = 0.9) +
    # LOESS smoothers (post-MPA only)
    geom_smooth(data = dplyr::filter(d, year >= mpa_start_yr),
                aes(group = species_short),
                method = "loess", se = FALSE, span = 0.75,
                linewidth = 0.7, alpha = 0.55) +
    # Data points (na.rm suppresses warnings from invisible dummy rows)
    geom_point(size = 1.8, alpha = 0.6, na.rm = TRUE) +
    scale_color_taxa(name = "Species") +
    scale_shape_source(name = "Source", drop = FALSE) +
    scale_x_continuous(breaks = year_breaks) +
    scale_y_continuous(limits = c(-2, 2), breaks = seq(-2, 2, by = 1)) +
    labs(
      title = mpa_display,
      x = "Year",
      y = "Log response ratio (lnRR)"
    ) +
    theme_mpa(base_size = 9) +
    theme(
      plot.title = element_text(size = 10, face = "plain", hjust = 0,
                                 margin = margin(b = 4)),
      legend.text = element_text(face = "italic", size = 8),
      legend.title = element_text(face = "plain", size = 8.5),
      legend.key.height = unit(0.45, "cm"),
      axis.title = element_text(size = 8.5),
      axis.text = element_text(size = 7.5),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.border = element_blank(),
      plot.margin = margin(4, 6, 4, 4)
    )
}

# Build panels using purrr::map for tidyverse consistency
fig_s2_panels <- purrr::map(fig_s2_mpas, build_s2_panel,
                             data = fig_s2_data,
                             mpa_starts = fig_s2_starts)

# Apply legend styling to each panel individually (avoid patchwork & theme()
# which corrupts axis.line hierarchy in ggplot2 4.x)
s2_legend_theme <- theme(
  legend.position = "bottom",
  legend.direction = "horizontal",
  legend.box = "vertical",
  legend.justification = "center",
  legend.title = element_text(face = "plain", size = 8.5),
  legend.text = element_text(size = 7.5, face = "italic"),
  legend.key.width = unit(0.4, "cm"),
  legend.key.height = unit(0.35, "cm"),
  legend.spacing.x = unit(0.15, "cm"),
  legend.box.margin = margin(0, 0, 0, 0),
  legend.margin = margin(t = 2, b = 2)
)
fig_s2_panels <- purrr::map(fig_s2_panels, function(p) p + s2_legend_theme)

# Combine panels with collected legend — NO & theme() to avoid axis corruption
cat("  [DEBUG] Fig S2: has_patchwork =", has_patchwork, ", has_ggpubr =", has_ggpubr, "\n")
if (has_patchwork) {
  fig_s2 <- (fig_s2_panels[[1]] / fig_s2_panels[[2]] / fig_s2_panels[[3]]) +
    plot_layout(guides = "collect")
} else if (has_ggpubr) {
  fig_s2 <- ggpubr::ggarrange(
    plotlist = fig_s2_panels,
    ncol = 1, nrow = 3,
    labels = c("a", "b", "c"),
    font.label = list(size = 10, face = "plain"),
    label.x = 0.02, label.y = 0.98,
    common.legend = TRUE,
    legend = "bottom"
  )
} else {
  fig_s2 <- fig_s2_panels[[1]] / fig_s2_panels[[2]] / fig_s2_panels[[3]]
}

# Slightly wider to accommodate right legend
save_fig(fig_s2, "fig_s02_all_taxa_timeseries", FIG_S2_DIMS["w"], FIG_S2_DIMS["h"])

# (Figure 5 removed — trophic cascade now shown in Figure 4 with dual layers)
# The following block is commented out but preserved for reference.


# =============================================================================
# Figure S3: Temporal Dynamics of Trophic Cascade
# =============================================================================
# Shows how effect sizes unfold over time since MPA implementation

cat("\n--- Figure S3: Temporal Dynamics ---\n")

FIG_S3_DIMS <- c(w = 17, h = 17)  # Conservation Letters max width

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
  "Kelp"      = unname(col_taxa["M. pyrifera"])
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
             hjust = 0, vjust = 1.5, size = 2.5, color = "grey50") +
    annotate("text", x = 0.5, y = -Inf, label = "Decrease in MPA",
             hjust = 0, vjust = -0.5, size = 2.5, color = "grey50") +
    geom_ribbon(aes(ymin = mean_lnRR - 1.96 * se_lnRR,
                    ymax = mean_lnRR + 1.96 * se_lnRR),
                alpha = 0.2, color = NA) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5, shape = 21, stroke = 0.8, color = "white",
               show.legend = FALSE) +
    scale_color_manual(values = trophic_colors, name = "Trophic Level", drop = FALSE) +
    scale_fill_manual(values = trophic_colors, name = "Trophic Level", drop = FALSE,
                      guide = "none") +
    # Ensure all trophic levels (including Kelp) appear in legend;
    # show colored lines as legend keys (suppress point glyphs and fill legend)
    guides(
      color = guide_legend(override.aes = list(linewidth = 1.2))
    ) +
    scale_x_continuous(breaks = seq(0, 20, by = 5)) +
    labs(
      x = "Years since MPA implementation",
      y = "Mean log response ratio (lnRR)"
    ) +
    theme_mpa(base_size = 10) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "plain"),
      panel.grid.major = element_blank()
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
      x = "Years since MPA implementation",
      y = "Cumulative mean lnRR"
    ) +
    theme_mpa(base_size = 10) +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank()
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
      x = "Trophic Level",
      y = expression("Effect trajectory slope (lnRR" ~ yr^{-1} * ")")
    ) +
    theme_mpa(base_size = 10) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(face = "plain")
    )

  # Combine panels
  fig_s3 <- panel_S3A / (panel_S3B | panel_S3C) +
    plot_layout(heights = c(1.2, 1)) +
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
    theme(
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.line.x.bottom = element_line(colour = "black", linewidth = 0.3),
      axis.line.y.left = element_line(colour = "black", linewidth = 0.3)
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

FIG_S4_DIMS <- c(w = 17, h = 14)  # Conservation Letters max width (adjusted height to maintain aspect)

# Get MPA implementation years
mpa_years <- Site %>%
  dplyr::select(CA_MPA_Name_Short, MPA_Start) %>%
  distinct() %>%
  filter(!is.na(MPA_Start))

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
  filter(!is.na(MPA_Start))

# Shorten MPA names (remove designation suffixes) for cleaner y-axis labels
heatmap_data$CA_MPA_Name_Short <- shorten_mpa_name(heatmap_data$CA_MPA_Name_Short)

# Order MPAs by implementation year (newer at bottom)
mpa_order <- heatmap_data %>%
  distinct(CA_MPA_Name_Short, MPA_Start) %>%
  arrange(desc(MPA_Start)) %>%
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
      x = "Years since MPA implementation",
      y = "MPA (ordered by implementation year)"
    ) +
    theme_mpa(base_size = 10) +
    theme(
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

FIG_S5_DIMS <- c(w = 17, h = 11)  # Conservation Letters max width (adjusted height to maintain aspect)

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
  scale_fill_manual(values = col_model, name = "Effect size\nmethod") +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  labs(
    x = NULL,
    y = "Proportion of MPAs"
  ) +
  theme_mpa(base_size = 10) +
  theme(
    axis.text.x = element_text(face = "italic", angle = 25, hjust = 1),
    legend.position = "right",
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

  panel_S5B <- ggplot(var_comp, aes(x = tau2, y = Response, fill = Component)) +
    geom_col(position = "dodge", width = 0.6, color = "white", linewidth = 0.3) +
    scale_fill_manual(
      values = c("MPA" = "#2A7B8E", "Source" = "#D4933B"),
      name = "Variance\nComponent"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(
      x = expression(tau^2),
      y = NULL
    ) +
    theme_mpa(base_size = 10) +
    theme(
      legend.position = "right",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  fig_s5 <- panel_S5A / panel_S5B +
    plot_layout(heights = c(1, 0.8)) +
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
    theme(
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.line.x.bottom = element_line(colour = "black", linewidth = 0.3),
      axis.line.y.left = element_line(colour = "black", linewidth = 0.3)
    )

} else {
  fig_s5 <- panel_S5A +
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
}

save_fig(fig_s5, "fig_s05_statistical_transparency", FIG_S5_DIMS["w"], FIG_S5_DIMS["h"])


# =============================================================================
# Figure S6: Comprehensive Site-Level Appendix
# =============================================================================
# Individual lnRR time series for ALL taxa at ALL sites, including sites
# excluded from the final analysis. Provides transparency and allows readers
# to inspect site-level variation underlying the meta-analytic summaries.

cat("\n--- Figure S6: Site-Level Appendix ---\n")

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
    facet_wrap(~ CA_MPA_Name_Short, ncol = n_cols, scales = "free_x") +
    labs(
      title = taxa_i,
      x = "Year",
      y = resp_label
    ) +
    theme_mpa(base_size = 8) +
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
      facet_wrap(~ CA_MPA_Name_Short, ncol = n_cols, scales = "free_x",
                 labeller = labeller(CA_MPA_Name_Short = new_labels))
  }

  # Dimensions scale with number of MPAs (capped at Conservation Letters max)
  fig_w <- min(17, n_cols * 4)  # Conservation Letters max width: 17cm
  fig_h <- max(10, n_rows * 4)

  # Create clean filename from taxa name
  taxa_filename <- tolower(gsub("\\. ", "", taxa_i))
  save_fig(fig_appendix,
           paste0("fig_s06_appendix_", taxa_filename),
           fig_w, fig_h)
}

cat("  Site-level appendix complete. Dagger (\u2020) marks sites excluded from analysis.\n")


cat("\n")
cat("========================================================================\n")
cat("  ALL MANUSCRIPT FIGURES GENERATED SUCCESSFULLY\n")
cat("  Main text: Fig 1-4\n")
cat("  Supplemental: Fig S1-S6\n")
cat("========================================================================\n")
cat("\n=== All figures saved to:", here::here("plots"), "===\n")
