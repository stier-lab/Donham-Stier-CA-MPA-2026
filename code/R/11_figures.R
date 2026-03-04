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
#   Figure 1: MPA map with bathymetry
#     - Ocean bathymetry (depth gradient)
#     - MPA boundaries with monitoring sites
#
#   Figure 2: Trophic cascade case studies (before/after)
#     - 3×3 grid: trophic rows × 3 exemplar MPAs (Scorpion, Gull Is., SB Is.)
#     - Before/after MPA implementation, with linear trends on after period
#     - Demonstrates predators↑, urchins↓, kelp↑ unfolding over time
#
#   Figure 3: Mean effect sizes from meta-analysis
#     - Summarizes Table 2 graphically
#     - Shows meta-analytic means with 95% CIs (diamonds)
#     - Individual MPA effect sizes shown as background points (circles)
#
#   Figure 4: Recovery trajectories over time (biomass)
#     - 3×2 trophic grid: predators/herbivores/producer rows
#     - lmer prediction lines with 95% CI
#
#   Figure S1 (Supplemental): Data processing pipeline example
#     - 4-panel illustration using KFM purple urchin at Scorpion SMR
#
#   Figure S2 (Supplemental): Forest plot of effect sizes
#     - Individual effect sizes by MPA and taxa
#     - Color-coded by response type (density vs biomass)
#     - Shape-coded by data source (PISCO, KFM, LTER, Landsat)
#
#   Dropped from SI (code retained, files still generated on disk):
#     - Figure S3 (old): All taxa time series at example MPAs
#     - Figure S7 (old): Statistical transparency
#
#   Figure S7a-e (Supplemental): Comprehensive site-level appendix
#     - Individual lnRR time series for ALL taxa at ALL sites
#     - Includes sites excluded from final analysis (marked with dagger)
#     - One file per taxa: fig_s08_appendix_*.pdf
#
# DESIGN PRINCIPLES:
#   - Uses colorblind-safe palette from 00b_color_palette.R
#   - Publication-ready sizing for Conservation Letters (80-170mm width)
#   - Consistent theme_mpa() styling across all figures
#   - Exported as both PDF (vector) and PNG (raster at 600 DPI)
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
#   Main text:
#   - fig_01_mpa_map.pdf / .png
#   - fig_02_cascade_case_studies.pdf / .png
#   - fig_03_mean_effects.pdf / .png
#   - fig_04_recovery_curves.pdf / .png
#   Supplemental:
#   - fig_s01_data_processing.pdf / .png
#   - fig_s02_forest_plot.pdf / .png
#   - fig_s08_appendix_*.pdf / .png  (one per taxa, = SI Fig S7a-e)
#   Dropped from SI (still generated on disk):
#   - fig_s03_all_taxa_timeseries.pdf / .png
#   - fig_s07_statistical_transparency.pdf / .png
#   - fig_s10_recovery_bio_den.pdf / .png
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
# Set RENDER_FIGURES before sourcing to render a subset (e.g., c("fig02", "fig03"))
# Default: render all figures
if (!exists("RENDER_FIGURES", envir = .GlobalEnv)) {
  RENDER_FIGURES <- "all"
}

# Standard factor levels: canonical trophic order (Predators → Urchins → Kelp)
taxa_levels  <- c("P. interruptus", "S. pulcher",
                   "S. purpuratus", "M. franciscanus",
                   "M. pyrifera")
source_levels <- c("KFM", "LTER", "PISCO", "Landsat")

# =============================================================================
# Figure dimension constants (Conservation Letters specifications)
# =============================================================================
# Wiley standard widths: single column = 80mm, double column = 180mm
# Figures must be 80-180mm wide; line art at 600 DPI preferred
#
# FONT SIZE NOTE (Conservation Letters requirement):
#   Minimum 8pt at final printed size. All figures in this script are
#   double-column (170-180mm) and use base_size = 9 or 10 in theme_mpa(),
#   with the smallest axis text at 7pt (forest plot MPA names). At 170mm
#   width this exceeds the 8pt minimum when accounting for print scaling.
#   If any figure were reduced to single-column (80mm), axis text sizes
#   would need to be verified against the 8pt floor.
FIG_WIDTH_SINGLE <- 8    # cm (80mm), for single-column figures
FIG_WIDTH_DOUBLE <- 17   # cm (170mm), for double-column figures
FIG_WIDTH_WIDE   <- 18   # cm (180mm), Wiley max width
FIG_WIDTH_SUPP   <- 17.8 # cm, Conservation Letters max width for supplemental figures

# Figure-specific dimensions (width, height in cm)
# Note: Figure 1 dimensions are defined inside the should_render("fig01") block below
FIG_S01_DIMS <- c(w = 17, h = 20)   # Data processing pipeline — 3 vertical panels, double-column
FIG3_DIMS <- c(w = 17, h = 11)   # Mean effects — trophic cascade layout, double-column
FIG4_DIMS <- c(w = 17, h = 12)  # Recovery curves — 3×2 trophic grid, 170mm double-column
FIG_S10_DIMS <- c(w = 17, h = 22)  # Recovery curves — 5×2 grid (biomass + density) + legend
# FIG_S11_DIMS <- c(w = 17, h = 9) # Recovery curves — density [DROPPED: merged concept into Fig 4]
FIG_S2_DIMS <- c(w = 17.8, h = 23.5) # Forest plot (supplemental)
FIG_S3_DIMS <- c(w = 17.8, h = 26) # All taxa time series (supplemental, faceted layout)
FIG_S16_DIMS <- c(w = 17, h = 18)  # Outlier removal rationale — 2×2 panel layout
FIG_S17_DIMS <- c(w = 17.8, h = 25) # Temporal trajectories with outlier status — 5×2 faceted
FIG_S18_DIMS <- c(w = 17.8, h = 25) # Raw data trajectories with outlier status — 5×2 faceted

# =============================================================================
# Shared variables used across multiple figure sections
# =============================================================================
# MPAs excluded from forest plot and individual effect size displays
excluded_mpas <- EXCLUDED_MPAS  # Use canonical list from 01_utils.R

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
# ggtext is needed for element_markdown() in Figure 3 italic x-axis labels
require_pkgs(c("ggplot2", "dplyr", "here", "scales", "forcats", "purrr", "stringr", "ggtext"))

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
  stop("All.RR.sub.trans must be a non-empty dataframe. Check 03_load_harmonized_data.R.")
}

# Validate All.Resp.sub structure
if (!is.data.frame(All.Resp.sub) || nrow(All.Resp.sub) == 0) {
  stop("All.Resp.sub must be a non-empty dataframe. Check 03_load_harmonized_data.R.")
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


# =============================================================================
# Figure 1: MPA Map with Bathymetry
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
  FIG1_DIMS <- c(w = 17, h = 10)  # cm — map only (time series moved to Fig 2)

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
    "No-Take MPA"         = "#3A5A6C",
    "Partial Protection"  = "#D4DFE8"
  )

  # --- 4. Load Coastline ---
  coast <- rnaturalearth::ne_states(country = "united states of america", returnclass = "sf") %>%
    filter(name == "California")

  # --- 5. Load Monitoring Sites ---
  site_csv <- here::here("data", "Site_List_All.csv")
  if (!file.exists(site_csv)) {
    warning("data/Site_List_All.csv not found — skipping Fig 1 site overlay")
    sites_base <- data.frame(Lon = numeric(0), Lat = numeric(0), program = character(0),
                             CA_MPA_Name_Short = character(0))
    sites_labels <- data.frame(Lon = numeric(0), Lat = numeric(0),
                               label_x = numeric(0), label_y = numeric(0),
                               site_abbrev = character(0))
    PANEL_SITES <- c("b" = "Campus Point SMCA", "c" = "Harris Point SMR",
                     "d" = "South Point SMR", "e" = "Santa Barbara Island SMR")
    MPA_YEARS <- c("Campus Point SMCA" = 2012, "Harris Point SMR" = 2003,
                   "South Point SMR" = 2003, "Santa Barbara Island SMR" = 2003)
  } else {
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
    "f" = c( 0.12,  0.07),   # Carrington Pt — right-above (separated from h)
    "g" = c(-0.12,  0.03),   # Skunk Pt — left-above, scooted
    "h" = c(-0.05,  0.08),   # Painted Cave — left-above (separated from f)
    "i" = c(-0.12, -0.05),   # Gull Island — further left-below
    "j" = c( 0.10,  0.05),   # Scorpion — right-above
    "k" = c( 0.12,  0.08),   # Anacapa SMR — right-above (separated from j)
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
  } # end Site_List_All.csv else block

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
		        title.theme = element_text(size = 8, face = "plain"),
		        frame.colour = "grey60",
		        frame.linewidth = 0.25,
	        ticks = TRUE,
	        ticks.colour = "grey40",
		        ticks.linewidth = 0.35,
		        label.position = "bottom",
		        label.theme = element_text(size = 8, color = "grey25"),
		        order = 3
		      )
		    ) +
	    geom_contour(data = bathy_ocean, aes(x = lon, y = lat, z = depth),
	                 breaks = c(-100, -200, -500, -1000, -2000, -3000), color = "white",
	                 linewidth = 0.20, alpha = 0.25) +
    new_scale_fill() +
    geom_sf(data = mpa, aes(fill = mpa_group), color = fig1_map_colors$mpa_border,
            alpha = 0.65, linewidth = 0.5, inherit.aes = FALSE) +
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
                 color = "grey30", linewidth = 0.5, alpha = 1.0) +
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
                     pad_y = unit(0.3, "in"), style = "ticks", text_cex = 0.90, line_width = 0.5) +
    # North arrow — top-right, fancy orienteering style
    annotation_north_arrow(location = "tr", which_north = "true", pad_x = unit(0.15, "in"),
                           pad_y = unit(0.4, "in"), style = north_arrow_fancy_orienteering,
                           height = unit(0.7, "cm"), width = unit(0.7, "cm")) +
    theme_mpa(base_size = 9) +
    labs(tag = "(a)") +
    theme(
      panel.background = element_rect(fill = "#C6DBEF", color = NA),
      panel.border = element_rect(fill = NA, color = "grey35", linewidth = 0.4),
      axis.line.x.bottom = element_blank(),
      axis.line.y.left = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(size = 8, color = "grey20"),
      # Legends above the map
      legend.position = "top",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = element_text(size = 8, face = "plain"),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.3, "cm"),
      legend.spacing.x = unit(6, "mm"),
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

  # --- 9. California Inset Map (removed per author request) ---

  # --- 10. Time Series Panels (removed — time series now in Figure 2) ---

  # --- 11. Finalize map ---
  fig1 <- main_map + plot_annotation(theme = theme(plot.margin = margin(1, 8, 2, 2)))

  # --- 12. Save Figure 1 (full) and map-only (panel a) ---
  save_fig(fig1, "fig_01_mpa_map", FIG1_DIMS["w"], FIG1_DIMS["h"], dpi = 600)
  # fig_01a_map_only removed from figure suite (orphan)

} else {
  cat("  Skipping Figure 1 (required packages not available)\n")
}
} # end fig01

# Load additional packages for remaining figures
if (has_cowplot) library(cowplot)
if (has_ggrepel) library(ggrepel)


# =============================================================================
# Figure S1: Conceptual diagram — simulated kelp biomass (3 vertical panels)
# Simulated data modeled on real M. pyrifera patterns at Scorpion SMR
# =============================================================================

if (should_render("fig_s01")) {
cat("Building Figure S1 (Supplemental): Data processing pipeline (simulated data)...\n")

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
    plot.subtitle = element_text(size = 8, color = "grey50",
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
    plot.subtitle = element_text(size = 8, color = "grey50",
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
           size = 2.8, color = MPA_LABEL_COLOR, fontface = "italic") +
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
           size = 2.8, color = MPA_LABEL_COLOR, fontface = "italic") +
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
           size = 2.8, color = "grey35", fontface = "italic") +
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

save_fig(fig2_final, "fig_s01_data_processing", FIG_S01_DIMS["w"], FIG_S01_DIMS["h"])
} # end fig_s01

# =============================================================================
# Figure S2 (Supplemental): Forest plot of effect sizes by MPA and taxa
# Manuscript: Supplemental Figure S2
# =============================================================================

if (should_render("fig_s02")) {
cat("Building Figure S2 (Supplemental): Forest plot...\n")

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

fig_s1_data$MPA <- as.character(fig_s1_data$MPA)

# Shortened MPA names for display
fig_s1_data <- fig_s1_data %>%
  dplyr::mutate(MPA_short = shorten_mpa_name(MPA))

# Order MPAs by mean effect size within each taxa
fig_s1_data <- fig_s1_data %>%
  dplyr::group_by(Taxa) %>%
  dplyr::mutate(
    MPA_order = forcats::fct_reorder(MPA_short, Mean, .fun = mean, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()

# --- Meta-analytic summary from Table2 (dashed lines + CI bands per panel) ---
meta_summary_s1 <- Table2 %>%
  dplyr::transmute(
    Taxa = factor(Taxa, levels = taxa_levels),
    Resp = factor(
      dplyr::recode(Response, "Density" = "Den", "Biomass" = "Bio"),
      levels = c("Den", "Bio")
    ),
    Mean = Estimate,
    CI_lo = CI_lower,
    CI_hi = CI_upper
  )

# Dynamic RR-labelled breaks: each free-scaled panel picks its own ticks.
# Data is on the lnRR scale internally; labels show back-transformed RR values
# where RR = 1 (lnRR = 0) means no MPA effect.
rr_breaks_fn <- function(lim) {
  rr_pool <- c(0.01, 0.05, 0.1, 0.25, 0.5, 1, 2, 4, 10, 20, 100)
  candidates <- log(rr_pool)
  in_range <- sort(candidates[candidates >= lim[1] & candidates <= lim[2]])
  if (length(in_range) == 0) return(pretty(lim, n = 3))
  # Thin to ~4 ticks max
  if (length(in_range) > 4) {
    idx <- round(seq(1, length(in_range), length.out = 4))
    in_range <- unique(in_range[idx])
  }
  # Always include RR = 1 (lnRR = 0) — the key no-effect reference value.
  # If it crowds a neighbor, drop that neighbor instead.
  if (0 >= lim[1] && 0 <= lim[2] && !(0 %in% in_range)) {
    range_width <- diff(lim)
    too_close <- which(abs(in_range - 0) / range_width < 0.10)
    if (length(too_close) > 0) in_range <- in_range[-too_close]
    in_range <- sort(c(in_range, 0))
  }
  in_range
}
rr_labels_fn <- function(x) {
  rr <- exp(x)
  ifelse(rr >= 1, as.character(as.integer(round(rr))), as.character(round(rr, 2)))
}

# Add panel tags to facet strip labels
fig_s1_data$Taxa <- factor(fig_s1_data$Taxa,
  levels = taxa_levels,
  labels = paste0("(", letters[seq_along(taxa_levels)], ")  ", taxa_levels))
meta_summary_s1$Taxa <- factor(meta_summary_s1$Taxa,
  levels = taxa_levels,
  labels = paste0("(", letters[seq_along(taxa_levels)], ")  ", taxa_levels))

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
  # Meta-analytic CI as a subtle vertical band behind everything
  geom_rect(
    data = meta_summary_s1,
    aes(xmin = CI_lo, xmax = CI_hi, ymin = -Inf, ymax = Inf, fill = Resp),
    inherit.aes = FALSE, alpha = 0.08, show.legend = FALSE
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.4) +
  # Meta-analytic mean as dashed vertical line per response type
  geom_vline(
    data = meta_summary_s1,
    aes(xintercept = Mean, color = Resp),
    linetype = "dashed", linewidth = 0.5, alpha = 0.7,
    show.legend = FALSE
  ) +
  geom_errorbar(aes(y = MPA_order, xmin = Mean - CI, xmax = Mean + CI),
                width = 0, linewidth = 0.5, position = pd_s1, alpha = 0.7) +
  geom_point(size = 2.5, position = pd_s1) +
  facet_wrap(~ Taxa, ncol = 2, scales = "free") +
  scale_color_response(
    name = "Response",
    labels = c("Den" = "Density", "Bio" = "Biomass")
  ) +
  scale_fill_response(
    name = "Response",
    labels = c("Den" = "Density", "Bio" = "Biomass"),
    guide = "none"
  ) +
  scale_shape_source(name = "Source") +
  # RR-labelled x-axis: data on lnRR scale, labels show back-transformed RR
  # where RR = 1 (lnRR = 0) means no MPA effect (free scale per panel)
  scale_x_continuous(breaks = rr_breaks_fn, labels = rr_labels_fn) +
  labs(x = "MPA / Reference", y = NULL) +
  theme_mpa(base_size = 9) +
  theme(
    strip.text = element_text(face = "bold.italic", size = 10, margin = margin(3, 0, 3, 0)),
    strip.background = element_blank(),
    # 8pt minimum per Conservation Letters guidelines; height increased to compensate
    axis.text.y = element_text(size = 8, color = "grey15"),
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

save_fig(fig_s1, "fig_s02_forest_plot", FIG_S2_DIMS["w"], FIG_S2_DIMS["h"])
} # end fig_s02

# =============================================================================
# Figure 3: Mean effect sizes by taxa from meta-analysis
# Manuscript: Main Text Figure 3
# =============================================================================

if (should_render("fig03")) {
cat("Building Figure 3: Mean effect sizes from meta-analysis...\n")

# --- Trophic-level ordering (Predators → Urchins → Kelp) ---
# Use numeric x-axis for precise annotation placement (separators, headers)
fig2_trophic_levels <- c("P. interruptus", "S. pulcher",
                          "S. purpuratus", "M. franciscanus",
                          "M. pyrifera")
fig2_taxa_x <- setNames(seq_along(fig2_trophic_levels), fig2_trophic_levels)

# Common names for x-axis labels (Latin + common)
fig2_common_names <- c(
  "S. pulcher"      = "Sheephead",
  "P. interruptus"  = "Spiny lobster",
  "S. purpuratus"   = "Purple urchin",
  "M. franciscanus" = "Red urchin",
  "M. pyrifera"     = "Giant kelp"
)

# Prepare Table2 for plotting with numeric x positions
fig2_meta <- Table2 %>%
  dplyr::mutate(
    Response = factor(Response, levels = c("Density", "Biomass")),
    x_pos = fig2_taxa_x[as.character(Taxa)]
  ) %>%
  dplyr::filter(!is.na(x_pos))

# Prepare individual effect sizes for background points
fig2_individual <- SumStats.Final %>%
  dplyr::filter(!(MPA %in% excluded_mpas)) %>%
  dplyr::mutate(
    Mean = as.numeric(Mean),
    Response = factor(
      ifelse(Resp == "Den", "Density", "Biomass"),
      levels = c("Density", "Biomass")
    ),
    x_pos = fig2_taxa_x[as.character(Taxa)]
  ) %>%
  dplyr::filter(!is.na(x_pos))

# Calculate sample sizes for annotation
fig2_n <- fig2_individual %>%
  dplyr::group_by(Taxa, Response) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop")

# Add sample sizes to meta data (ensure Taxa types match for join)
fig2_meta <- fig2_meta %>%
  dplyr::mutate(Taxa = as.character(Taxa))
fig2_n <- fig2_n %>%
  dplyr::mutate(Taxa = as.character(Taxa))
fig2_meta <- fig2_meta %>%
  dplyr::left_join(fig2_n, by = c("Taxa", "Response"))

# FDR significance stars: * p<0.05, ** p<0.01, *** p<0.001
fig2_meta <- fig2_meta %>%
  dplyr::mutate(sig_star = dplyr::case_when(
    pval_fdr < 0.001 ~ "***",
    pval_fdr < 0.01  ~ "**",
    pval_fdr < 0.05  ~ "*",
    TRUE ~ ""
  ))

# Dynamic y-limits
fig2_y_min <- min(c(fig2_individual$Mean, fig2_meta$CI_lower), na.rm = TRUE)
fig2_y_max <- max(c(fig2_individual$Mean, fig2_meta$CI_upper), na.rm = TRUE)
fig2_label_y <- fig2_y_min - 0.4

# Dodge offset for Density vs Biomass (manual, since x is numeric)
# Center taxa that have only one response type (e.g., M. pyrifera = Biomass only)
fig2_dodge <- 0.3
fig2_resp_count <- fig2_meta %>%
  dplyr::count(Taxa) %>%
  dplyr::rename(n_resp = n)
fig2_meta <- fig2_meta %>%
  dplyr::left_join(fig2_resp_count, by = "Taxa") %>%
  dplyr::mutate(x_dodge = ifelse(n_resp == 1, x_pos,
                  x_pos + ifelse(Response == "Density", -fig2_dodge/2, fig2_dodge/2)))
fig2_ind_resp_count <- fig2_individual %>%
  dplyr::distinct(Taxa, Response) %>%
  dplyr::count(Taxa) %>%
  dplyr::rename(n_resp = n)
set.seed(42)
fig2_individual <- fig2_individual %>%
  dplyr::left_join(fig2_ind_resp_count, by = "Taxa") %>%
  dplyr::mutate(
    x_dodge = ifelse(n_resp == 1, x_pos,
                x_pos + ifelse(Response == "Density", -fig2_dodge/2, fig2_dodge/2)),
    x_jitter = x_dodge + runif(dplyr::n(), -0.08, 0.08)
  )

# Shape mapping: Density = diamond (18), Biomass = circle (16) for grayscale safety
fig2_shape_map <- c("Density" = 18, "Biomass" = 16)

# Trophic group separator and header positions
# Predators (1,2) | Urchins (3,4) | Kelp (5)
fig2_sep_x <- c(2.5, 4.5)
fig2_header_y <- fig2_y_max + 1.2

fig2 <- ggplot() +
  # Faint vertical separators between trophic groups
  geom_vline(xintercept = fig2_sep_x, color = "grey85", linewidth = 0.3) +
  # Zero-effect reference line
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.5, linetype = "dashed") +
  # Layer 1: Individual MPA effect sizes (subdued background, shape by Response)
  geom_point(data = fig2_individual,
             aes(x = x_jitter, y = Mean, color = Response, shape = Response),
             size = 1.6, alpha = 0.30) +
  # Layer 2: White halo behind focal points
  geom_point(data = fig2_meta,
             aes(x = x_dodge, y = Estimate, shape = Response),
             size = 5.5, color = "white", show.legend = FALSE) +
  # Layer 3: 95% CIs
  geom_errorbar(data = fig2_meta,
                aes(x = x_dodge, ymin = CI_lower, ymax = CI_upper, color = Response),
                width = 0.12, linewidth = 0.9) +
  # Layer 4: Meta-analytic means (shape by Response for grayscale safety)
  geom_point(data = fig2_meta,
             aes(x = x_dodge, y = Estimate, color = Response, shape = Response),
             size = 4.5) +
  # Sample size labels
  geom_text(data = fig2_meta,
            aes(x = x_dodge, y = fig2_label_y, label = paste0("k=", k),
                color = Response),
            size = 3.2, show.legend = FALSE) +
  # FDR significance stars above CI
  geom_text(data = fig2_meta %>% dplyr::filter(sig_star != ""),
            aes(x = x_dodge, y = CI_upper + 0.25, label = sig_star),
            size = 4.0, color = "grey20", show.legend = FALSE) +
  # Trophic group headers (above plot area)
  annotate("text", x = c(1.5, 3.5, 5), y = fig2_header_y,
           label = c("PREDATORS", "URCHINS", "KELP"),
           size = 3.3, color = "grey35", fontface = "bold") +
  scale_color_manual(name = NULL, values = col_response_long,
                     guide = guide_legend(override.aes = list(size = 3.5, alpha = 1))) +
  scale_shape_manual(name = NULL, values = fig2_shape_map) +
  scale_x_continuous(
    breaks = seq_along(fig2_trophic_levels),
    labels = paste0("*", fig2_trophic_levels, "*<br>(", fig2_common_names[fig2_trophic_levels], ")"),
    expand = expansion(mult = 0.08)
  ) +
  # RR-labelled y-axis: data stays on lnRR (log response ratio) scale internally;
  # tick labels show back-transformed RR values. RR = 1 (lnRR = 0) means no MPA
  # effect; RR > 1 = higher inside MPA, RR < 1 = lower inside MPA.
  scale_y_rr(fig2_y_min - 0.9, fig2_y_max + 0.7) +
  coord_cartesian(ylim = c(fig2_y_min - 0.9, fig2_y_max + 0.7),
                  clip = "off") +
  labs(x = NULL) +
  theme_mpa(base_size = 10) +
  theme(
    axis.text.x = ggtext::element_markdown(size = 8.5, lineheight = 1.3),
    axis.title.y = element_text(size = 10),
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    legend.key.width = unit(0.8, "cm"),
    legend.margin = margin(t = 2),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.line.x.bottom = element_line(colour = "black", linewidth = 0.5),
    axis.line.y.left = element_line(colour = "black", linewidth = 0.5),
    plot.margin = margin(18, 10, 6, 6),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

save_fig(fig2, "fig_03_mean_effects", FIG3_DIMS["w"], FIG3_DIMS["h"])
} # end fig03

# =============================================================================
# Figure S11 (Supplemental): Trophic cascade meta-regression scatterplots
# Uses metafor::rma() meta-regression with Knapp-Hartung adjustment.
# 12-panel layout: 3 columns (Sheephead→Urchin, Lobster→Urchin, Urchin→Kelp)
#                  4 rows (Purple Urchin Bio, Purple Urchin Den, Red Urchin Bio, Red Urchin Den)
# Also exports comprehensive cascade results table (20 pathways) to outputs/.
# =============================================================================

# ---------------------------------------------------------------------------
# Cascade data prep (always runs — needed for table_cascade_analysis.csv)
# ---------------------------------------------------------------------------
cat("Building cascade wide-format data for table export...\n")

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
  tidyr::pivot_wider(names_from = Taxa_Resp, values_from = SE,
                     values_fn = function(se_vec) {
                       # Inverse-variance pooling (correct for SE aggregation)
                       w <- 1 / se_vec^2
                       sqrt(1 / sum(w))
                     })

fig4_wide <- dplyr::left_join(fig4_wide_mean, fig4_wide_se, by = "MPA")
names(fig4_wide) <- gsub(" ", "_", names(fig4_wide))
names(fig4_wide) <- gsub("\\.", "_", names(fig4_wide))

cat("  Wide-format columns:", paste(names(fig4_wide), collapse = ", "), "\n")

# Detect column naming convention (handles both S_pulcher and S__pulcher)
fig4_find_col <- function(pattern) {
  matches <- grep(pattern, names(fig4_wide), value = TRUE)
  if (length(matches) == 1) return(matches)
  if (length(matches) > 1) return(matches[1])
  return(NA_character_)
}

# Panel column mappings — all 3 predator/prey pairs
fig4_cols <- list(
  # Sheephead
  sp_bio      = fig4_find_col("pulcher_Bio$"),
  sp_den      = fig4_find_col("pulcher_Den$"),
  sp_se_bio   = fig4_find_col("pulcher_SE_Bio$"),
  sp_se_den   = fig4_find_col("pulcher_SE_Den$"),
  # Lobster
  pi_bio      = fig4_find_col("interruptus_Bio$"),
  pi_den      = fig4_find_col("interruptus_Den$"),
  pi_se_bio   = fig4_find_col("interruptus_SE_Bio$"),
  pi_se_den   = fig4_find_col("interruptus_SE_Den$"),
  # Purple urchin
  spur_bio    = fig4_find_col("purpuratus_Bio$"),
  spur_den    = fig4_find_col("purpuratus_Den$"),
  spur_se_bio = fig4_find_col("purpuratus_SE_Bio$"),
  spur_se_den = fig4_find_col("purpuratus_SE_Den$"),
  # Red urchin (M. franciscanus)
  mf_bio      = fig4_find_col("franciscanus_Bio$"),
  mf_den      = fig4_find_col("franciscanus_Den$"),
  mf_se_bio   = fig4_find_col("franciscanus_SE_Bio$"),
  mf_se_den   = fig4_find_col("franciscanus_SE_Den$"),
  # Kelp
  mp_bio      = fig4_find_col("pyrifera_Bio$"),
  mp_se_bio   = fig4_find_col("pyrifera_SE_Bio$")
)

cat("  Column mapping:\n")
for (nm in names(fig4_cols)) cat("    ", nm, "=", fig4_cols[[nm]], "\n")

# (Raw annual scatter removed — figure shows only the MPA-level effect sizes
#  that enter the rma() meta-regression, matching the actual analysis.)

if (FALSE) {  # DROPPED: cascade scatter figure (weak cross-sectional evidence; "fig_s11" namespace now used by DHARMa diagnostics)
cat("Building cascade scatter (Supplemental): Trophic cascade meta-regression scatterplots...\n")
# ---------------------------------------------------------------------------
# Helper: create one trophic cascade panel using rma() meta-regression
# (matches 09_meta_analysis.R Table 3 methodology)
# ---------------------------------------------------------------------------
create_cascade_panel <- function(wide_data,
                                  x_col, y_col, se_x_col, se_y_col,
                                  x_lab, y_lab, point_color,
                                  x_short, y_short,
                                  base_size = 9,
                                  point_size = 3.5,
                                  stat_size = 2.8,
                                  quad_size = 2.4,
                                  panel_bg = "white",
                                  show_regression = TRUE) {
  # Subset to MPAs with complete data for required columns
  req_cols <- intersect(c(x_col, y_col, se_y_col), names(wide_data))
  if (length(req_cols) < 3) {
    cat("    WARNING: Missing columns for panel. Available:", paste(names(wide_data), collapse = ", "), "\n")
    return(list(
      plot = ggplot() + theme_void() +
        annotate("text", x = 0.5, y = 0.5, label = "Missing data columns", size = 4),
      results = data.frame(x_var = x_col, y_var = y_col, x_resp = NA, y_resp = NA,
                           beta = NA, se = NA, tval = NA,
                           ci_lb = NA, ci_ub = NA,
                           p_raw = NA, k = 0, R2 = NA, I2 = NA, tau2 = NA,
                           QE = NA, QEp = NA,
                           stringsAsFactors = FALSE)
    ))
  }
  check_cols <- req_cols
  if (se_x_col %in% names(wide_data)) check_cols <- c(check_cols, se_x_col)
  mpa_df <- wide_data[complete.cases(wide_data[, check_cols]), ]
  mpa_df[[x_col]] <- as.numeric(mpa_df[[x_col]])
  mpa_df[[y_col]] <- as.numeric(mpa_df[[y_col]])
  mpa_df[[se_y_col]] <- as.numeric(mpa_df[[se_y_col]])

  cat("    x=", x_col, " y=", y_col, " k=", nrow(mpa_df), "\n")

  if (nrow(mpa_df) < 3) {
    return(list(
      plot = ggplot() + theme_void() +
        annotate("text", x = 0.5, y = 0.5,
                 label = paste0("Insufficient data (k=", nrow(mpa_df), ")"), size = 4),
      results = data.frame(x_var = x_col, y_var = y_col, x_resp = NA, y_resp = NA,
                           beta = NA, se = NA, tval = NA,
                           ci_lb = NA, ci_ub = NA,
                           p_raw = NA, k = nrow(mpa_df), R2 = NA, I2 = NA, tau2 = NA,
                           QE = NA, QEp = NA,
                           stringsAsFactors = FALSE)
    ))
  }

  # ---- Meta-regression: rma(yi = y, vi = SE_y^2, mods = ~ x) ----
  rma_df <- data.frame(
    yi    = mpa_df[[y_col]],
    vi    = mpa_df[[se_y_col]]^2,
    x_mod = mpa_df[[x_col]]
  )

  meta_mod <- tryCatch(
    metafor::rma(yi = yi, vi = vi, mods = ~ x_mod, data = rma_df,
                 method = "REML", test = "knha"),
    error = function(e) {
      cat("    rma() error:", conditionMessage(e), "\n")
      NULL
    }
  )

  if (is.null(meta_mod)) {
    return(list(
      plot = ggplot() + theme_void() +
        annotate("text", x = 0.5, y = 0.5, label = "Model failed", size = 4),
      results = data.frame(x_var = x_col, y_var = y_col, x_resp = NA, y_resp = NA,
                           beta = NA, se = NA, tval = NA,
                           ci_lb = NA, ci_ub = NA,
                           p_raw = NA, k = nrow(mpa_df), R2 = NA, I2 = NA, tau2 = NA,
                           QE = NA, QEp = NA,
                           stringsAsFactors = FALSE)
    ))
  }

  # Extract slope statistics
  coef_tbl <- coef(summary(meta_mod))
  slope_est <- coef_tbl[2, "estimate"]
  slope_se  <- coef_tbl[2, "se"]
  slope_p   <- coef_tbl[2, "pval"]
  # tval column present with Knapp-Hartung; zval otherwise
  slope_tval <- if ("tval" %in% colnames(coef_tbl)) coef_tbl[2, "tval"] else
                if ("zval" %in% colnames(coef_tbl)) coef_tbl[2, "zval"] else NA_real_
  r_sq      <- ifelse(is.null(meta_mod$R2), 0, meta_mod$R2)
  tau_sq    <- meta_mod$tau2
  i_sq      <- meta_mod$I2
  QE_stat   <- meta_mod$QE
  QE_p      <- meta_mod$QEp

  # Display raw Knapp-Hartung p-value (FDR applied only in supplemental table)
  display_p <- slope_p

  stat_label <- sprintf(
    "\u03b2 = %.2f, p %s\nk = %d MPAs, R\u00b2 = %.0f%%",
    slope_est,
    ifelse(display_p < 0.001, "< 0.001", sprintf("= %.3f", display_p)),
    meta_mod$k, r_sq
  )

  # ---- Predicted regression line with 95% CI from rma model ----
  x_seq <- seq(min(rma_df$x_mod), max(rma_df$x_mod), length.out = 100)
  preds <- predict(meta_mod, newmods = x_seq)
  reg_df <- data.frame(x = x_seq, y = preds$pred,
                        ci_lb = preds$ci.lb, ci_ub = preds$ci.ub)

  # ---- Axis limits: data-fitted, always include zero ----
  x_vals <- rma_df$x_mod
  y_vals <- rma_df$yi
  se_x_vals <- if (se_x_col %in% names(mpa_df)) as.numeric(mpa_df[[se_x_col]]) else rep(0, nrow(mpa_df))
  se_y_vals <- as.numeric(mpa_df[[se_y_col]])
  # Data extent including SE and regression CI
  x_lo_raw <- min(c(x_vals - se_x_vals, 0), na.rm = TRUE)
  x_hi_raw <- max(c(x_vals + se_x_vals, 0), na.rm = TRUE)
  y_lo_raw <- min(c(y_vals - se_y_vals, reg_df$ci_lb, 0), na.rm = TRUE)
  y_hi_raw <- max(c(y_vals + se_y_vals, reg_df$ci_ub, 0), na.rm = TRUE)
  # 15% padding on each side
  pad_frac <- 0.15
  x_rng <- x_hi_raw - x_lo_raw
  y_rng <- y_hi_raw - y_lo_raw
  lim_x <- c(x_lo_raw - x_rng * pad_frac, x_hi_raw + x_rng * pad_frac)
  lim_y <- c(y_lo_raw - y_rng * pad_frac, y_hi_raw + y_rng * pad_frac)
  # Ensure at least 2.0 lnRR on each side of zero so all four quadrants
  # are visible with room for two-line labels in each corner
  lim_x <- c(min(lim_x[1], -2.0), max(lim_x[2], 2.0))
  lim_y <- c(min(lim_y[1], -2.0), max(lim_y[2], 2.0))

  # ---- Quadrant labels (ecological interpretation) ----
  # Single-line labels anchored at the outer edges of each quadrant.
  # Only show labels where the quadrant is wide enough to avoid overlap
  # with the adjacent label across the reference line.
  pad_x <- (lim_x[2] - lim_x[1]) * 0.03
  pad_y <- (lim_y[2] - lim_y[1]) * 0.03
  quad_labels <- data.frame(
    x = c(lim_x[1] + pad_x, lim_x[2] - pad_x,
           lim_x[1] + pad_x, lim_x[2] - pad_x),
    y = c(lim_y[2] - pad_y, lim_y[2] - pad_y,
           lim_y[1] + pad_y, lim_y[1] + pad_y),
    hj = c(0, 1, 0, 1),   # left-align left cols, right-align right cols
    vj = c(1, 1, 0, 0),   # top-align top row, bottom-align bottom row
    label = c(
      paste0(x_short, " \u2193\n", y_short, " \u2191"),   # top-left
      paste0(x_short, " \u2191\n", y_short, " \u2191"),   # top-right
      paste0(x_short, " \u2193\n", y_short, " \u2193"),   # bottom-left
      paste0(x_short, " \u2191\n", y_short, " \u2193")    # bottom-right
    ),
    stringsAsFactors = FALSE
  )
  # All four quadrants are now guaranteed ≥1.5 lnRR on each side,
  # so all four labels will fit without overlap.

  # ---- Build plot ----
  p <- ggplot() +
    # Quadrant labels (behind everything)
    geom_text(data = quad_labels, aes(x = x, y = y, label = label,
              hjust = hj, vjust = vj),
              size = quad_size, color = "grey70", lineheight = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.5)

  # Vertical error bars (y SE)
  p <- p + geom_errorbar(data = mpa_df,
                aes(x = .data[[x_col]],
                    ymin = .data[[y_col]] - .data[[se_y_col]],
                    ymax = .data[[y_col]] + .data[[se_y_col]]),
                width = 0, linewidth = 0.5, color = "grey30", alpha = 0.6)

  # Horizontal error bars (x SE) — only if column exists
  if (se_x_col %in% names(mpa_df)) {
    mpa_df[[se_x_col]] <- as.numeric(mpa_df[[se_x_col]])
    p <- p + geom_errorbarh(data = mpa_df,
                   aes(y = .data[[y_col]],
                       xmin = .data[[x_col]] - .data[[se_x_col]],
                       xmax = .data[[x_col]] + .data[[se_x_col]]),
                   height = 0, linewidth = 0.5, color = "grey30", alpha = 0.6)
  }

  # Regression line + stats (optional — suppressed for descriptive main figure)
  if (show_regression) {
    sig <- display_p < 0.05
    reg_linetype <- ifelse(sig, "solid", "dashed")
    p <- p +
      geom_ribbon(data = reg_df, aes(x = x, ymin = ci_lb, ymax = ci_ub),
                  fill = point_color, alpha = ifelse(sig, 0.12, 0.05)) +
      geom_line(data = reg_df, aes(x = x, y = y),
                color = point_color, linewidth = 0.8,
                alpha = ifelse(sig, 0.7, 0.4),
                linetype = reg_linetype)
  }

  # MPA-level points (always shown)
  p <- p +
    geom_point(data = mpa_df, aes(x = .data[[x_col]], y = .data[[y_col]]),
               size = point_size, alpha = 0.90, color = point_color, shape = 16)

  # Annotation: full stats when regression shown, just k when not
  if (show_regression) {
    p <- p +
      annotate("label", x = mean(lim_x), y = lim_y[1] + pad_y * 2,
               label = stat_label, hjust = 0.5, vjust = 0,
               size = stat_size, color = "grey20", lineheight = 0.85,
               fill = alpha("white", 0.85), label.size = 0,
               label.padding = unit(1.5, "pt"))
  } else {
    p <- p +
      annotate("label", x = lim_x[2] - pad_x, y = lim_y[1] + pad_y,
               label = sprintf("k = %d", nrow(mpa_df)),
               hjust = 1, vjust = 0,
               size = stat_size, color = "grey40",
               fill = alpha("white", 0.85), label.size = 0,
               label.padding = unit(1.5, "pt"))
  }

  p <- p +
    labs(x = x_lab, y = y_lab) +
    # RR-scaled axes: data on lnRR, ticks labelled with back-transformed RR
    # Use sparse ticks for small panels (base_size < 9)
    scale_x_rr(lim_x[1], lim_x[2], name = x_lab, sparse = (base_size <= 9)) +
    scale_y_rr(lim_y[1], lim_y[2], name = y_lab, sparse = (base_size <= 9)) +
    coord_cartesian(xlim = lim_x, ylim = lim_y, expand = FALSE, clip = "off") +
    theme_mpa(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_rect(fill = panel_bg, colour = NA),
      axis.title = element_text(size = base_size - 1),
      axis.text = element_text(size = base_size - 1),
      plot.margin = margin(6, 6, 2, 4)
    )

  # Return both the plot and model results for FDR correction and table export
  list(
    plot = p,
    results = data.frame(
      x_var    = x_col,
      y_var    = y_col,
      x_resp   = sub(".*_SE_", "", se_x_col),
      y_resp   = sub(".*_SE_", "", se_y_col),
      beta     = slope_est,
      se       = slope_se,
      tval     = slope_tval,
      ci_lb    = coef_tbl[2, "ci.lb"],
      ci_ub    = coef_tbl[2, "ci.ub"],
      p_raw    = slope_p,
      k        = meta_mod$k,
      R2       = r_sq,
      I2       = i_sq,
      tau2     = tau_sq,
      QE       = QE_stat,
      QEp      = QE_p,
      stringsAsFactors = FALSE
    )
  )
}

# Shared sizing for 6-panel layout (smaller than 4-panel)
f4_bs   <- 9     # base_size (≥8pt axis text at base_size - 1)
f4_pt   <- 2.8   # point size
f4_stat <- 2.8   # stat label size (≥8pt)
f4_quad <- 2.85  # quadrant label size (≥8pt: 2.85/0.3528 = 8.08pt)
# No special background for any column
f4_shared_bg <- "white"

# ---------------------------------------------------------------------------
# Define ALL 12 cascade pathway specs for figures (same-response pairs)
# Layout: 3 columns (Sheephead→Urchin, Lobster→Urchin, Urchin→Kelp)
#         4 rows (Purple Bio, Purple Den, Red Bio, Red Den)
# ---------------------------------------------------------------------------

fig4_all_specs <- list(
  # Row 1: Purple urchin, Biomass
  pr_bio_a = list(x = "sp_bio",   y = "spur_bio", sex = "sp_se_bio",   sey = "spur_se_bio",
                  xlab = expression(italic("S. pulcher")~"biomass (RR)"),
                  ylab = expression(italic("S. purpuratus")~"biomass (RR)"),
                  color = "S. pulcher", xs = "Sheep.", ys = "P. Urchin",
                  col_header = "Sheephead \u2192 Urchin", row_label = "Purple Urchin\nBiomass"),
  pr_bio_b = list(x = "pi_bio",   y = "spur_bio", sex = "pi_se_bio",   sey = "spur_se_bio",
                  xlab = expression(italic("P. interruptus")~"biomass (RR)"),
                  ylab = expression(italic("S. purpuratus")~"biomass (RR)"),
                  color = "P. interruptus", xs = "Lobster", ys = "P. Urchin",
                  col_header = "Lobster \u2192 Urchin", row_label = NA),
  pr_bio_c = list(x = "spur_bio", y = "mp_bio",   sex = "spur_se_bio", sey = "mp_se_bio",
                  xlab = expression(italic("S. purpuratus")~"biomass (RR)"),
                  ylab = expression(italic("M. pyrifera")~"biomass (RR)"),
                  color = "S. purpuratus", xs = "P. Urchin", ys = "Kelp",
                  col_header = "Urchin \u2192 Kelp", row_label = NA),
  # Row 2: Purple urchin, Density
  pr_den_a = list(x = "sp_den",   y = "spur_den", sex = "sp_se_den",   sey = "spur_se_den",
                  xlab = expression(italic("S. pulcher")~"density (RR)"),
                  ylab = expression(italic("S. purpuratus")~"density (RR)"),
                  color = "S. pulcher", xs = "Sheep.", ys = "P. Urchin",
                  col_header = NA, row_label = "Purple Urchin\nDensity"),
  pr_den_b = list(x = "pi_den",   y = "spur_den", sex = "pi_se_den",   sey = "spur_se_den",
                  xlab = expression(italic("P. interruptus")~"density (RR)"),
                  ylab = expression(italic("S. purpuratus")~"density (RR)"),
                  color = "P. interruptus", xs = "Lobster", ys = "P. Urchin",
                  col_header = NA, row_label = NA),
  pr_den_c = list(x = "spur_den", y = "mp_bio",   sex = "spur_se_den", sey = "mp_se_bio",
                  xlab = expression(italic("S. purpuratus")~"density (RR)"),
                  ylab = expression(italic("M. pyrifera")~"biomass (RR)"),
                  color = "S. purpuratus", xs = "P. Urchin", ys = "Kelp",
                  col_header = NA, row_label = NA),
  # Row 3: Red urchin, Biomass
  rd_bio_a = list(x = "sp_bio",   y = "mf_bio",   sex = "sp_se_bio",   sey = "mf_se_bio",
                  xlab = expression(italic("S. pulcher")~"biomass (RR)"),
                  ylab = expression(italic("M. franciscanus")~"biomass (RR)"),
                  color = "S. pulcher", xs = "Sheep.", ys = "R. Urchin",
                  col_header = NA, row_label = "Red Urchin\nBiomass"),
  rd_bio_b = list(x = "pi_bio",   y = "mf_bio",   sex = "pi_se_bio",   sey = "mf_se_bio",
                  xlab = expression(italic("P. interruptus")~"biomass (RR)"),
                  ylab = expression(italic("M. franciscanus")~"biomass (RR)"),
                  color = "P. interruptus", xs = "Lobster", ys = "R. Urchin",
                  col_header = NA, row_label = NA),
  rd_bio_c = list(x = "mf_bio",   y = "mp_bio",   sex = "mf_se_bio",   sey = "mp_se_bio",
                  xlab = expression(italic("M. franciscanus")~"biomass (RR)"),
                  ylab = expression(italic("M. pyrifera")~"biomass (RR)"),
                  color = "M. franciscanus", xs = "R. Urchin", ys = "Kelp",
                  col_header = NA, row_label = NA),
  # Row 4: Red urchin, Density
  rd_den_a = list(x = "sp_den",   y = "mf_den",   sex = "sp_se_den",   sey = "mf_se_den",
                  xlab = expression(italic("S. pulcher")~"density (RR)"),
                  ylab = expression(italic("M. franciscanus")~"density (RR)"),
                  color = "S. pulcher", xs = "Sheep.", ys = "R. Urchin",
                  col_header = NA, row_label = "Red Urchin\nDensity"),
  rd_den_b = list(x = "pi_den",   y = "mf_den",   sex = "pi_se_den",   sey = "mf_se_den",
                  xlab = expression(italic("P. interruptus")~"density (RR)"),
                  ylab = expression(italic("M. franciscanus")~"density (RR)"),
                  color = "P. interruptus", xs = "Lobster", ys = "R. Urchin",
                  col_header = NA, row_label = NA),
  rd_den_c = list(x = "mf_den",   y = "mp_bio",   sex = "mf_se_den",   sey = "mp_se_bio",
                  xlab = expression(italic("M. franciscanus")~"density (RR)"),
                  ylab = expression(italic("M. pyrifera")~"biomass (RR)"),
                  color = "M. franciscanus", xs = "R. Urchin", ys = "Kelp",
                  col_header = NA, row_label = NA)
)

# ---------------------------------------------------------------------------
# Build all 12 panels and collect model results
# ---------------------------------------------------------------------------
cat("  Building all 12 cascade panels...\n")
fig4_all_panels <- list()
fig4_all_results <- list()
for (pid in names(fig4_all_specs)) {
  s <- fig4_all_specs[[pid]]
  result <- create_cascade_panel(
    wide_data = fig4_wide,
    x_col = fig4_cols[[s$x]], y_col = fig4_cols[[s$y]],
    se_x_col = fig4_cols[[s$sex]], se_y_col = fig4_cols[[s$sey]],
    x_lab = s$xlab, y_lab = s$ylab,
    point_color = col_taxa[s$color],
    x_short = s$xs, y_short = s$ys,
    base_size = f4_bs, point_size = f4_pt, stat_size = f4_stat,
    quad_size = f4_quad, panel_bg = "white"
  )
  fig4_all_panels[[pid]] <- result$plot
  fig4_all_results[[pid]] <- result$results
  fig4_all_results[[pid]]$panel_id <- pid
}

# Collect results for reporting
fig4_results_df <- do.call(rbind, fig4_all_results)
cat("  Panel results (Knapp-Hartung p-values):\n")
for (i in seq_len(nrow(fig4_results_df))) {
  r <- fig4_results_df[i, ]
  sig_flag <- if (!is.na(r$p_raw) && r$p_raw < 0.05) " *" else ""
  cat("    ", r$panel_id, ": beta=",
      ifelse(is.na(r$beta), "NA", sprintf("%.3f", r$beta)),
      " p=", ifelse(is.na(r$p_raw), "NA", sprintf("%.4f", r$p_raw)),
      " k=", r$k, sig_flag, "\n")
}

# =============================================================================
# Supplemental Figure: All 12 cascade panels (3 cols x 4 rows)
# =============================================================================
cat("Assembling supplemental 12-panel cascade figure...\n")

# Add column headers to top-row panels only
fig4_all_panels[["pr_bio_a"]] <- fig4_all_panels[["pr_bio_a"]] +
  ggtitle("Sheephead \u2192 Urchin") +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5, margin = margin(b = 2)))
fig4_all_panels[["pr_bio_b"]] <- fig4_all_panels[["pr_bio_b"]] +
  ggtitle("Lobster \u2192 Urchin") +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5, margin = margin(b = 2)))
fig4_all_panels[["pr_bio_c"]] <- fig4_all_panels[["pr_bio_c"]] +
  ggtitle("Urchin \u2192 Kelp") +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5, margin = margin(b = 2)))

# Add panel tags (a)-(l)
panel_tags <- letters[1:12]
panel_order <- c("pr_bio_a", "pr_bio_b", "pr_bio_c",
                 "pr_den_a", "pr_den_b", "pr_den_c",
                 "rd_bio_a", "rd_bio_b", "rd_bio_c",
                 "rd_den_a", "rd_den_b", "rd_den_c")
for (i in seq_along(panel_order)) {
  fig4_all_panels[[panel_order[i]]] <- fig4_all_panels[[panel_order[i]]] +
    labs(tag = paste0("(", panel_tags[i], ")"))
}

# Assemble 3x4 grid
fig4_supp_grid <- (fig4_all_panels[["pr_bio_a"]] | fig4_all_panels[["pr_bio_b"]] | fig4_all_panels[["pr_bio_c"]]) /
                  (fig4_all_panels[["pr_den_a"]] | fig4_all_panels[["pr_den_b"]] | fig4_all_panels[["pr_den_c"]]) /
                  (fig4_all_panels[["rd_bio_a"]] | fig4_all_panels[["rd_bio_b"]] | fig4_all_panels[["rd_bio_c"]]) /
                  (fig4_all_panels[["rd_den_a"]] | fig4_all_panels[["rd_den_b"]] | fig4_all_panels[["rd_den_c"]]) +
  plot_annotation(theme = theme(plot.tag = element_text(face = "bold", size = 9))) +
  plot_layout(widths = c(1, 1, 1), heights = c(1, 1, 1, 1)) &
  theme(
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.line.x.bottom = element_line(colour = "black", linewidth = 0.5),
    axis.line.y.left = element_line(colour = "black", linewidth = 0.5)
  )

# Overlay row labels via cowplot
FIG_S11_DIMS <- c(w = 18, h = 24)  # 3x4 panel layout
fig4_supp <- cowplot::ggdraw(fig4_supp_grid) +
  cowplot::draw_label("Purple Urchin\nBiomass",  x = 0.015, y = 0.875, angle = 90,
                      fontface = "bold", size = 8, hjust = 0.5, lineheight = 0.9) +
  cowplot::draw_label("Purple Urchin\nDensity",  x = 0.015, y = 0.625, angle = 90,
                      fontface = "bold", size = 8, hjust = 0.5, lineheight = 0.9) +
  cowplot::draw_label("Red Urchin\nBiomass",     x = 0.015, y = 0.375, angle = 90,
                      fontface = "bold", size = 8, hjust = 0.5, lineheight = 0.9) +
  cowplot::draw_label("Red Urchin\nDensity",     x = 0.015, y = 0.125, angle = 90,
                      fontface = "bold", size = 8, hjust = 0.5, lineheight = 0.9)

save_fig(fig4_supp, "fig_s11_cascade_scatter", FIG_S11_DIMS["w"], FIG_S11_DIMS["h"])

} # end S10 figure rendering (DROPPED)

# ---------------------------------------------------------------------------
# Comprehensive cascade table: all 20 pathways (including cross-response)
# Includes same-response (12 above) plus cross-response combinations.
# FDR correction applied across all testable pathways in the table.
# Always runs (supports Table 3 in manuscript).
# ---------------------------------------------------------------------------
cat("Running comprehensive cascade analysis (all pathways for table)...\n")

cascade_pathways <- list(
  # --- Purple urchin same-response ---
  list(x = "sp_bio",  y = "spur_bio", sex = "sp_se_bio",  sey = "spur_se_bio",
       predator = "S. pulcher", prey = "S. purpuratus", x_resp = "Bio", y_resp = "Bio",
       pathway = "Predator \u2192 Urchin", link = "Sheephead Bio \u2192 Purple Urchin Bio"),
  list(x = "sp_den",  y = "spur_den", sex = "sp_se_den",  sey = "spur_se_den",
       predator = "S. pulcher", prey = "S. purpuratus", x_resp = "Den", y_resp = "Den",
       pathway = "Predator \u2192 Urchin", link = "Sheephead Den \u2192 Purple Urchin Den"),
  list(x = "pi_bio",  y = "spur_bio", sex = "pi_se_bio",  sey = "spur_se_bio",
       predator = "P. interruptus", prey = "S. purpuratus", x_resp = "Bio", y_resp = "Bio",
       pathway = "Predator \u2192 Urchin", link = "Lobster Bio \u2192 Purple Urchin Bio"),
  list(x = "pi_den",  y = "spur_den", sex = "pi_se_den",  sey = "spur_se_den",
       predator = "P. interruptus", prey = "S. purpuratus", x_resp = "Den", y_resp = "Den",
       pathway = "Predator \u2192 Urchin", link = "Lobster Den \u2192 Purple Urchin Den"),
  list(x = "spur_bio", y = "mp_bio", sex = "spur_se_bio", sey = "mp_se_bio",
       predator = "S. purpuratus", prey = "M. pyrifera", x_resp = "Bio", y_resp = "Bio",
       pathway = "Urchin \u2192 Kelp", link = "Purple Urchin Bio \u2192 Kelp Bio"),
  list(x = "spur_den", y = "mp_bio", sex = "spur_se_den", sey = "mp_se_bio",
       predator = "S. purpuratus", prey = "M. pyrifera", x_resp = "Den", y_resp = "Bio",
       pathway = "Urchin \u2192 Kelp", link = "Purple Urchin Den \u2192 Kelp Bio"),
  # --- Purple urchin cross-response ---
  list(x = "sp_bio",  y = "spur_den", sex = "sp_se_bio",  sey = "spur_se_den",
       predator = "S. pulcher", prey = "S. purpuratus", x_resp = "Bio", y_resp = "Den",
       pathway = "Predator \u2192 Urchin", link = "Sheephead Bio \u2192 Purple Urchin Den"),
  list(x = "sp_den",  y = "spur_bio", sex = "sp_se_den",  sey = "spur_se_bio",
       predator = "S. pulcher", prey = "S. purpuratus", x_resp = "Den", y_resp = "Bio",
       pathway = "Predator \u2192 Urchin", link = "Sheephead Den \u2192 Purple Urchin Bio"),
  list(x = "pi_bio",  y = "spur_den", sex = "pi_se_bio",  sey = "spur_se_den",
       predator = "P. interruptus", prey = "S. purpuratus", x_resp = "Bio", y_resp = "Den",
       pathway = "Predator \u2192 Urchin", link = "Lobster Bio \u2192 Purple Urchin Den"),
  list(x = "pi_den",  y = "spur_bio", sex = "pi_se_den",  sey = "spur_se_bio",
       predator = "P. interruptus", prey = "S. purpuratus", x_resp = "Den", y_resp = "Bio",
       pathway = "Predator \u2192 Urchin", link = "Lobster Den \u2192 Purple Urchin Bio"),
  # --- Red urchin same-response ---
  list(x = "sp_bio",  y = "mf_bio",  sex = "sp_se_bio",  sey = "mf_se_bio",
       predator = "S. pulcher", prey = "M. franciscanus", x_resp = "Bio", y_resp = "Bio",
       pathway = "Predator \u2192 Urchin", link = "Sheephead Bio \u2192 Red Urchin Bio"),
  list(x = "sp_den",  y = "mf_den",  sex = "sp_se_den",  sey = "mf_se_den",
       predator = "S. pulcher", prey = "M. franciscanus", x_resp = "Den", y_resp = "Den",
       pathway = "Predator \u2192 Urchin", link = "Sheephead Den \u2192 Red Urchin Den"),
  list(x = "pi_bio",  y = "mf_bio",  sex = "pi_se_bio",  sey = "mf_se_bio",
       predator = "P. interruptus", prey = "M. franciscanus", x_resp = "Bio", y_resp = "Bio",
       pathway = "Predator \u2192 Urchin", link = "Lobster Bio \u2192 Red Urchin Bio"),
  list(x = "pi_den",  y = "mf_den",  sex = "pi_se_den",  sey = "mf_se_den",
       predator = "P. interruptus", prey = "M. franciscanus", x_resp = "Den", y_resp = "Den",
       pathway = "Predator \u2192 Urchin", link = "Lobster Den \u2192 Red Urchin Den"),
  list(x = "mf_bio",  y = "mp_bio", sex = "mf_se_bio",  sey = "mp_se_bio",
       predator = "M. franciscanus", prey = "M. pyrifera", x_resp = "Bio", y_resp = "Bio",
       pathway = "Urchin \u2192 Kelp", link = "Red Urchin Bio \u2192 Kelp Bio"),
  list(x = "mf_den",  y = "mp_bio", sex = "mf_se_den",  sey = "mp_se_bio",
       predator = "M. franciscanus", prey = "M. pyrifera", x_resp = "Den", y_resp = "Bio",
       pathway = "Urchin \u2192 Kelp", link = "Red Urchin Den \u2192 Kelp Bio"),
  # --- Red urchin cross-response ---
  list(x = "sp_bio",  y = "mf_den",  sex = "sp_se_bio",  sey = "mf_se_den",
       predator = "S. pulcher", prey = "M. franciscanus", x_resp = "Bio", y_resp = "Den",
       pathway = "Predator \u2192 Urchin", link = "Sheephead Bio \u2192 Red Urchin Den"),
  list(x = "sp_den",  y = "mf_bio",  sex = "sp_se_den",  sey = "mf_se_bio",
       predator = "S. pulcher", prey = "M. franciscanus", x_resp = "Den", y_resp = "Bio",
       pathway = "Predator \u2192 Urchin", link = "Sheephead Den \u2192 Red Urchin Bio"),
  list(x = "pi_bio",  y = "mf_den",  sex = "pi_se_bio",  sey = "mf_se_den",
       predator = "P. interruptus", prey = "M. franciscanus", x_resp = "Bio", y_resp = "Den",
       pathway = "Predator \u2192 Urchin", link = "Lobster Bio \u2192 Red Urchin Den"),
  list(x = "pi_den",  y = "mf_bio",  sex = "pi_se_den",  sey = "mf_se_bio",
       predator = "P. interruptus", prey = "M. franciscanus", x_resp = "Den", y_resp = "Bio",
       pathway = "Predator \u2192 Urchin", link = "Lobster Den \u2192 Red Urchin Bio")
)

# Run all pathways
all_cascade_results <- list()
for (i in seq_along(cascade_pathways)) {
  pw <- cascade_pathways[[i]]
  x_col_name <- fig4_cols[[pw$x]]
  y_col_name <- fig4_cols[[pw$y]]
  se_y_name  <- fig4_cols[[pw$sey]]

  if (is.na(x_col_name) || is.na(y_col_name) || is.na(se_y_name)) {
    all_cascade_results[[i]] <- data.frame(
      pathway = pw$pathway, link = pw$link,
      predator = pw$predator, prey = pw$prey,
      x_response = pw$x_resp, y_response = pw$y_resp,
      beta = NA, se = NA, tval = NA,
      ci_lb = NA, ci_ub = NA,
      p_raw = NA, k = 0, R2 = NA, I2 = NA, tau2 = NA,
      QE = NA, QEp = NA,
      note = "Missing columns", stringsAsFactors = FALSE)
    next
  }

  req <- c(x_col_name, y_col_name, se_y_name)
  se_x_name <- fig4_cols[[pw$sex]]
  if (!is.na(se_x_name) && se_x_name %in% names(fig4_wide)) req <- c(req, se_x_name)
  sub_df <- fig4_wide[complete.cases(fig4_wide[, intersect(req, names(fig4_wide))]), ]
  k_val <- nrow(sub_df)

  if (k_val < 3) {
    all_cascade_results[[i]] <- data.frame(
      pathway = pw$pathway, link = pw$link,
      predator = pw$predator, prey = pw$prey,
      x_response = pw$x_resp, y_response = pw$y_resp,
      beta = NA, se = NA, tval = NA,
      ci_lb = NA, ci_ub = NA,
      p_raw = NA, k = k_val, R2 = NA, I2 = NA, tau2 = NA,
      QE = NA, QEp = NA,
      note = paste0("Insufficient data (k=", k_val, ")"),
      stringsAsFactors = FALSE)
    next
  }

  rma_df <- data.frame(
    yi = as.numeric(sub_df[[y_col_name]]),
    vi = as.numeric(sub_df[[se_y_name]])^2,
    x_mod = as.numeric(sub_df[[x_col_name]]))

  mod <- tryCatch(
    metafor::rma(yi = yi, vi = vi, mods = ~ x_mod, data = rma_df,
                 method = "REML", test = "knha"),
    error = function(e) NULL)

  if (is.null(mod)) {
    all_cascade_results[[i]] <- data.frame(
      pathway = pw$pathway, link = pw$link,
      predator = pw$predator, prey = pw$prey,
      x_response = pw$x_resp, y_response = pw$y_resp,
      beta = NA, se = NA, tval = NA,
      ci_lb = NA, ci_ub = NA,
      p_raw = NA, k = k_val, R2 = NA, I2 = NA, tau2 = NA,
      QE = NA, QEp = NA,
      note = "Model failed", stringsAsFactors = FALSE)
    next
  }

  ct <- coef(summary(mod))
  # Extract tval (Knapp-Hartung) or zval (standard)
  ct_tval <- if ("tval" %in% colnames(ct)) ct[2, "tval"] else
             if ("zval" %in% colnames(ct)) ct[2, "zval"] else NA_real_
  all_cascade_results[[i]] <- data.frame(
    pathway = pw$pathway, link = pw$link,
    predator = pw$predator, prey = pw$prey,
    x_response = pw$x_resp, y_response = pw$y_resp,
    beta = ct[2, "estimate"], se = ct[2, "se"],
    tval = ct_tval,
    ci_lb = ct[2, "ci.lb"], ci_ub = ct[2, "ci.ub"],
    p_raw = ct[2, "pval"], k = mod$k,
    R2 = ifelse(is.null(mod$R2), 0, mod$R2),
    I2 = mod$I2, tau2 = mod$tau2,
    QE = mod$QE, QEp = mod$QEp,
    note = "",
    stringsAsFactors = FALSE)
}

# Combine and apply FDR across all testable pathways (table only)
cascade_table <- do.call(rbind, all_cascade_results)
cascade_table$p_fdr <- NA_real_
testable <- !is.na(cascade_table$p_raw)
cascade_table$p_fdr[testable] <- p.adjust(cascade_table$p_raw[testable], method = "fdr")
# Expected signs based on trophic cascade ecology:
# Predator->Urchin: negative (predators suppress urchins)
# Urchin->Kelp: negative (urchins suppress kelp)
# Predator->Kelp: positive (indirect positive effect through urchin suppression)
cascade_table$expected_sign <- ifelse(
  grepl("Predator.*Kelp", cascade_table$pathway),
  "positive", "negative"
)
cascade_table$sign_match <- ifelse(is.na(cascade_table$beta), NA_character_,
  ifelse(
    (cascade_table$expected_sign == "negative" & cascade_table$beta < 0) |
    (cascade_table$expected_sign == "positive" & cascade_table$beta > 0),
    "yes", "no"
  ))
cascade_table$sig_raw <- ifelse(!is.na(cascade_table$p_raw) & cascade_table$p_raw < 0.05, "*", "")
cascade_table$sig_fdr <- ifelse(!is.na(cascade_table$p_fdr) & cascade_table$p_fdr < 0.05, "*", "")
cascade_table$method <- "rma(REML, Knapp-Hartung)"

# Round for publication
cascade_table$beta  <- round(cascade_table$beta, 3)
cascade_table$se    <- round(cascade_table$se, 3)
cascade_table$tval  <- round(cascade_table$tval, 3)
cascade_table$ci_lb <- round(cascade_table$ci_lb, 3)
cascade_table$ci_ub <- round(cascade_table$ci_ub, 3)
cascade_table$p_raw <- round(cascade_table$p_raw, 4)
cascade_table$p_fdr <- round(cascade_table$p_fdr, 4)
cascade_table$R2    <- round(cascade_table$R2, 1)
cascade_table$I2    <- round(cascade_table$I2, 1)
cascade_table$tau2  <- round(cascade_table$tau2, 4)
cascade_table$QE    <- round(cascade_table$QE, 3)
cascade_table$QEp   <- round(cascade_table$QEp, 4)

cascade_table <- cascade_table[, c(
  "pathway", "link", "predator", "prey",
  "x_response", "y_response", "expected_sign",
  "beta", "se", "tval", "ci_lb", "ci_ub", "sign_match",
  "p_raw", "p_fdr", "sig_raw", "sig_fdr",
  "k", "R2", "I2", "tau2", "QE", "QEp", "method", "note")]

cascade_csv_path <- here::here("outputs", "table_cascade_analysis.csv")
write.csv(cascade_table, cascade_csv_path, row.names = FALSE)
cat("Cascade table saved to:", cascade_csv_path, "\n")
cat("  Tested:", sum(testable), " | Sig (raw):", sum(cascade_table$sig_raw == "*", na.rm = TRUE),
    " | Sig (FDR):", sum(cascade_table$sig_fdr == "*", na.rm = TRUE),
    " | Sign match:", sum(cascade_table$sign_match == "yes", na.rm = TRUE),
    "of", sum(!is.na(cascade_table$sign_match)), "\n")


# =============================================================================
# Figure 4 (Main Text): Recovery Trajectories Over Time
# =============================================================================
# Validates the t=11 standardization by showing approximately linear
# trajectories for all five species. The linearity confirms that evaluating
# MPAs at a common time horizon captures the same rate-based process.
#
# This is a simplified main-text version of Figure S3 (which shows full
# MPA-level spaghetti). Here we emphasize the population-level smooth and
# mark the t=11 standardization point.

# --- Shared data prep for Fig 4 and S10 (runs if either is rendered) ---
if (should_render("fig04") || should_render("fig_s10")) {
cat("Preparing recovery trajectory data (Fig 4 / S10)...\n")

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

# Build temporal datasets split by response type (biomass vs density)
# NOTE: Separate columns allow response-specific lmer slopes from
# 10_temporal_analysis.R (bio and den models). M. pyrifera has biomass only.
fig5_data_all <- All.RR.sub.trans %>%
  dplyr::filter(BA == "After", time >= 0, time <= 15) %>%
  dplyr::mutate(
    Species = .data[[taxa_col]],
    time = as.numeric(time)
  ) %>%
  dplyr::filter(Species %in% fig5_species_order)

fig5_data_all$Species <- factor(fig5_data_all$Species, levels = fig5_species_order)

fig5_data_bio <- fig5_data_all %>% dplyr::filter(resp == "Bio")
fig5_data_den <- fig5_data_all %>% dplyr::filter(resp == "Den")

# Count MPAs per species x response type for panel annotation
fig5_n_mpas_bio <- fig5_data_bio %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(n_mpa = dplyr::n_distinct(CA_MPA_Name_Short), .groups = "drop")
fig5_n_mpas_den <- fig5_data_den %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(n_mpa = dplyr::n_distinct(CA_MPA_Name_Short), .groups = "drop")

# --- Helper: extract lmer prediction parameters from vcov-derived CSV ---
# Reads table_s_lmer_prediction_params.csv (exported by 10_temporal_analysis.R)
# which contains per-species intercept, slope, their SEs, covariance, and df.
# The CSV has a `resp` column ("Bio" or "Den") — pass resp_filter to select.
# Falls back to the old coefficient CSV (approximate, independence-assumed SEs)
# with a warning if the new CSV is not yet available.
# Returns a list with named vectors: slopes, intercepts, slope_ses,
#   intercept_ses, cov_int_slopes, df_slopes (keyed by full species name).
extract_lmer_slopes <- function(pred_params_path, resp_filter, fallback_csv_path) {
  null_result <- list(slopes = NULL, intercepts = NULL, slope_ses = NULL,
                      intercept_ses = NULL, cov_int_slopes = NULL,
                      df_slopes = NULL)

  # --- Try new prediction params CSV first ---
  if (file.exists(pred_params_path)) {
    params <- read.csv(pred_params_path, stringsAsFactors = FALSE)
    required_cols <- c("resp", "Species", "intercept", "intercept_se", "slope",
                       "slope_se", "cov_int_slope", "df_slope")
    if (!all(required_cols %in% names(params))) {
      warning("Prediction params CSV missing columns: ",
              paste(setdiff(required_cols, names(params)), collapse = ", "),
              ". Falling back to approximate method.", call. = FALSE)
    } else {
      # Filter to requested response type
      params <- params[params$resp == resp_filter, ]
      if (nrow(params) == 0) {
        warning("No rows for resp='", resp_filter, "' in prediction params CSV.",
                " Falling back to approximate method.", call. = FALSE)
      } else {
        # Build named vectors keyed by full species name
        slopes         <- setNames(params$slope,          params$Species)
        intercepts     <- setNames(params$intercept,      params$Species)
        slope_ses      <- setNames(params$slope_se,       params$Species)
        intercept_ses  <- setNames(params$intercept_se,   params$Species)
        cov_int_slopes <- setNames(params$cov_int_slope,  params$Species)
        df_slopes      <- setNames(params$df_slope,       params$Species)
        cat("  Loaded vcov-derived prediction params (", resp_filter,
            ") from: ", basename(pred_params_path), "\n")
        return(list(slopes = slopes, intercepts = intercepts,
                    slope_ses = slope_ses, intercept_ses = intercept_ses,
                    cov_int_slopes = cov_int_slopes, df_slopes = df_slopes))
      }
    }
  }

  # --- Fallback: read old coefficient CSV (approximate SEs) ---
  if (!file.exists(fallback_csv_path)) {
    return(null_result)
  }
  warning("Prediction params CSV not found at:\n  ", pred_params_path,
          "\n  Falling back to approximate CI (assumes independent slope/intercept SEs).",
          call. = FALSE)
  coefs <- read.csv(fallback_csv_path, stringsAsFactors = FALSE)
  ref_intercept <- coefs$Estimate[coefs$Term == "(Intercept)"]
  ref_slope <- coefs$Estimate[coefs$Term == "time"]
  if (length(ref_intercept) == 0 || length(ref_slope) == 0) {
    return(null_result)
  }
  slopes <- setNames(numeric(5), fig5_species_order)
  intercepts <- setNames(numeric(5), fig5_species_order)
  slopes["Panulirus interruptus"] <- ref_slope
  intercepts["Panulirus interruptus"] <- ref_intercept
  ref_slope_se <- coefs$SE[coefs$Term == "time"]
  slope_ses <- setNames(numeric(5), fig5_species_order)
  slope_ses["Panulirus interruptus"] <- ref_slope_se
  for (sp in fig5_species_order[-1]) {
    sp_int_term <- paste0("Species", sp, ":time")
    sp_main_term <- paste0("Species", sp)
    int_row <- coefs[coefs$Term == sp_int_term, ]
    main_row <- coefs[coefs$Term == sp_main_term, ]
    slopes[sp] <- ref_slope + ifelse(nrow(int_row) > 0, int_row$Estimate, 0)
    intercepts[sp] <- ref_intercept + ifelse(nrow(main_row) > 0, main_row$Estimate, 0)
    int_se <- ifelse(nrow(int_row) > 0, int_row$SE, 0)
    slope_ses[sp] <- sqrt(ref_slope_se^2 + int_se^2)
  }
  # Fallback: no intercept SEs or covariance available
  list(slopes = slopes, intercepts = intercepts, slope_ses = slope_ses,
       intercept_ses = NULL, cov_int_slopes = NULL, df_slopes = NULL)
}

# --- Load response-specific lmer prediction parameters ---
# Single prediction params CSV with resp column; fall back to separate coefficient CSVs.
fig5_pred_params_path <- here::here("tables", "table_s_lmer_prediction_params.csv")
fig5_bio_lmer <- extract_lmer_slopes(
  pred_params_path   = fig5_pred_params_path,
  resp_filter        = "Bio",
  fallback_csv_path  = here::here("tables", "table_s_temporal_meta_regression_bio.csv")
)
fig5_den_lmer <- extract_lmer_slopes(
  pred_params_path   = fig5_pred_params_path,
  resp_filter        = "Den",
  fallback_csv_path  = here::here("tables", "table_s_temporal_meta_regression_den.csv")
)
if (!is.null(fig5_bio_lmer$slopes)) {
  cat("  Biomass lmer parameters loaded.\n")
} else {
  cat("  NOTE: Biomass temporal meta-regression CSVs not found; using OLS slopes\n")
}
if (!is.null(fig5_den_lmer$slopes)) {
  cat("  Density lmer parameters loaded.\n")
} else {
  cat("  NOTE: Density temporal meta-regression CSVs not found; using OLS slopes\n")
}

# --- Compute per-species y-axis limits (independent per panel) ---
# Biomass-only limits for Figure 4 (main text, biomass panels only)
fig5_ylims_bio <- fig5_data_bio %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(
    q_lo = min(lnDiff, na.rm = TRUE),
    q_hi = max(lnDiff, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    pad = pmax((q_hi - q_lo) * 0.05, 0.2),
    ymin = q_lo - pad,
    ymax = q_hi + pad
  )

# Combined bio+den limits for Figure S10 (supplemental, both response types)
fig5_ylims_all <- fig5_data_all %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(
    q_lo = min(lnDiff, na.rm = TRUE),
    q_hi = max(lnDiff, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    pad = pmax((q_hi - q_lo) * 0.05, 0.2),
    ymin = q_lo - pad,
    ymax = q_hi + pad
  )

# --- Helper: build one species x response panel ---
# sp_data: data filtered to species AND response type
# sp_lmer: list from extract_lmer_slopes() with slopes, intercepts,
#   slope_ses, intercept_ses, cov_int_slopes, df_slopes (named vectors)
# n_mpas_df: data.frame with Species, n_mpa (response-specific)
# ylims_df: data.frame with Species, ymin, ymax for axis limits
make_fig5_panel <- function(sp, tag_label, show_xlab = FALSE,
                            sp_data, sp_lmer,
                            n_mpas_df, ylims_df) {
  panel_data <- sp_data[sp_data$Species == sp, ]
  sp_color <- fig5_species_colors[[sp]]
  sp_abbrev <- fig5_full_to_abbrev[[sp]]
  n_mpa_row <- n_mpas_df[n_mpas_df$Species == sp, ]
  n_mpa <- if (nrow(n_mpa_row) > 0) n_mpa_row$n_mpa else 0
  ylims <- ylims_df[ylims_df$Species == sp, ]
  sp_slope <- if (!is.null(sp_lmer$slopes)) sp_lmer$slopes[[sp]] else NA
  sp_int   <- if (!is.null(sp_lmer$intercepts)) sp_lmer$intercepts[[sp]] else NA
  sp_slope_se   <- if (!is.null(sp_lmer$slope_ses)) sp_lmer$slope_ses[[sp]] else NA
  sp_int_se     <- if (!is.null(sp_lmer$intercept_ses)) sp_lmer$intercept_ses[[sp]] else NA
  sp_cov        <- if (!is.null(sp_lmer$cov_int_slopes)) sp_lmer$cov_int_slopes[[sp]] else NA
  sp_df         <- if (!is.null(sp_lmer$df_slopes)) sp_lmer$df_slopes[[sp]] else NA

  # Formal lmer prediction at t=11 for point marker
  lnRR_at_11 <- if (!is.na(sp_int)) sp_int + sp_slope * 11 else NA

  p <- ggplot(panel_data, aes(x = time, y = lnDiff)) +
    # Individual MPA trajectories — species-tinted behind trend
    geom_line(aes(group = CA_MPA_Name_Short),
              color = sp_color, alpha = 0.18, linewidth = 0.4) +
    # Reference line at zero effect (solid, subtle)
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey30",
               linewidth = 0.5) +
    # Vertical reference at t=11 (standardized comparison point)
    geom_vline(xintercept = 11, linetype = "dotted", color = "grey50",
               linewidth = 0.5) +
    annotate("text", x = 11.3, y = Inf, label = "t = 11",
             hjust = 0, vjust = 1.3, size = 2.85, color = "grey40")

  # lmer prediction line with 95% CI ribbon (matches formal temporal
  # meta-regression). Falls back to OLS geom_smooth if lmer unavailable.
  if (!is.na(sp_int) && !is.na(sp_slope) && !is.na(sp_slope_se)) {
    # Build prediction data.frame from lmer fixed effects
    pred_df <- data.frame(
      time_seq = seq(0, 15, length.out = 50)
    )
    pred_df$pred <- sp_int + sp_slope * pred_df$time_seq

    # Correct prediction variance: Var(pred) = Var(int) + t^2*Var(slope) + 2*t*Cov(int,slope)
    # Use vcov-derived params if available; fall back to slope-only SE otherwise.
    if (!is.na(sp_int_se) && !is.na(sp_cov) && !is.na(sp_df)) {
      # Full variance formula with intercept-slope covariance
      pred_se <- sqrt(sp_int_se^2 +
                      pred_df$time_seq^2 * sp_slope_se^2 +
                      2 * pred_df$time_seq * sp_cov)
      t_crit <- qt(0.975, sp_df)
    } else {
      # Approximate: slope SE only (collapses to zero at t=0)
      pred_se <- sp_slope_se * pred_df$time_seq
      t_crit <- 1.96
    }
    pred_df$ci_lb <- pred_df$pred - t_crit * pred_se
    pred_df$ci_ub <- pred_df$pred + t_crit * pred_se

    p <- p +
      geom_ribbon(data = pred_df,
                  aes(x = time_seq, ymin = ci_lb, ymax = ci_ub),
                  fill = sp_color, alpha = 0.20, inherit.aes = FALSE) +
      geom_line(data = pred_df, aes(x = time_seq, y = pred),
                color = sp_color, linewidth = 1.3, alpha = 0.85, inherit.aes = FALSE)
  } else {
    # Fallback: OLS trend when lmer coefficients unavailable
    if (nrow(panel_data) > 2) {
      p <- p +
        geom_smooth(color = sp_color, fill = sp_color,
                    method = "lm", formula = y ~ x,
                    linewidth = 1.5, alpha = 0.20,
                    se = TRUE, level = 0.95)
    }
  }

  # Add t=11 point marker from lmer prediction
  if (!is.na(lnRR_at_11)) {
    p <- p + annotate("point", x = 11, y = lnRR_at_11,
                       color = sp_color, fill = "white",
                       shape = 21, size = 2.5, stroke = 0.8)
  }

  p <- p +
    # Shared y-axis within trophic group — clip ON to prevent spaghetti overflow
    coord_cartesian(ylim = c(ylims$ymin, ylims$ymax),
                    clip = "on") +
    scale_x_continuous(breaks = seq(0, 15, by = 5),
                       limits = c(0, 16), expand = c(0.02, 0)) +
    # RR-scaled y-axis: data is stored on lnRR (log response ratio) scale
    # internally, but tick labels show back-transformed RR values via
    # scale_y_rr(). RR = 1 (lnRR = 0) means no MPA effect. Shared y-axis
    # label "MPA / Reference" is added via cowplot::draw_label() below.
    # Sparse ticks to avoid crowding in compact panels.
    scale_y_rr(ylims$ymin, ylims$ymax, name = NULL, sparse = TRUE) +
    labs(
      title = sp_abbrev,
      x = if (show_xlab) "Years since MPA implementation" else NULL,
      y = NULL
    ) +
    theme_mpa(base_size = 9) +
    theme(
      plot.title = element_text(face = "italic", size = 9, hjust = 0,
                                margin = margin(b = 1)),
      plot.tag = element_text(size = 9, face = "bold"),
      plot.tag.position = c(0, 1),
      plot.margin = margin(3, 6, 3, 4, "pt")
    ) +
    labs(tag = paste0("(", tag_label, ")"))

  p
}

if (should_render("fig04")) {
cat("Building Figure 4 (Main Text): Recovery trajectories — 3×2 trophic grid...\n")
# --- Build 3×2 trophic grid: rows = trophic groups, columns = species pairs ---
# Row 1 (Predators):  P. interruptus (a) | S. pulcher (b)
# Row 2 (Herbivores): S. purpuratus (c)  | M. franciscanus (d)
# Row 3 (Producer):   M. pyrifera (e)    | [empty]

p_fig4_a <- make_fig5_panel("Panulirus interruptus", "a",
                             sp_data = fig5_data_bio,
                             sp_lmer = fig5_bio_lmer,
                             n_mpas_df = fig5_n_mpas_bio,
                             ylims_df = fig5_ylims_bio)
p_fig4_b <- make_fig5_panel("Semicossyphus pulcher", "b",
                             sp_data = fig5_data_bio,
                             sp_lmer = fig5_bio_lmer,
                             n_mpas_df = fig5_n_mpas_bio,
                             ylims_df = fig5_ylims_bio)
p_fig4_c <- make_fig5_panel("Strongylocentrotus purpuratus", "c",
                             sp_data = fig5_data_bio,
                             sp_lmer = fig5_bio_lmer,
                             n_mpas_df = fig5_n_mpas_bio,
                             ylims_df = fig5_ylims_bio)
p_fig4_d <- make_fig5_panel("Mesocentrotus franciscanus", "d",
                             sp_data = fig5_data_bio,
                             sp_lmer = fig5_bio_lmer,
                             n_mpas_df = fig5_n_mpas_bio,
                             ylims_df = fig5_ylims_bio)
p_fig4_e <- make_fig5_panel("Macrocystis pyrifera", "e", show_xlab = TRUE,
                             sp_data = fig5_data_bio,
                             sp_lmer = fig5_bio_lmer,
                             n_mpas_df = fig5_n_mpas_bio,
                             ylims_df = fig5_ylims_bio)

# Bottom-right: legend in the empty spacer cell
p_fig4_spacer <- cowplot::ggdraw() +
  # Population mean
  cowplot::draw_line(x = c(0.08, 0.22), y = c(0.78, 0.78),
                     color = "grey30", linewidth = 1.2) +
  cowplot::draw_label("Population mean\n(lmer \u00b1 95% CI)", x = 0.25, y = 0.78,
                      size = 8, fontfamily = "Helvetica", hjust = 0,
                      lineheight = 0.9) +
  # Individual MPA
  cowplot::draw_line(x = c(0.08, 0.22), y = c(0.45, 0.45),
                     color = "grey55", linewidth = 0.4) +
  cowplot::draw_label("Individual MPA", x = 0.25, y = 0.45,
                      size = 8, fontfamily = "Helvetica", hjust = 0) +
  # t = 11
  cowplot::draw_line(x = c(0.15, 0.15), y = c(0.14, 0.26),
                     color = "grey50", linewidth = 0.5, linetype = "dotted") +
  cowplot::draw_label("t = 11 years", x = 0.25, y = 0.20,
                      size = 8, fontfamily = "Helvetica", hjust = 0)

# --- Assemble 3×2 grid ---
fig4_grid <- (p_fig4_a | p_fig4_b) /
             (p_fig4_c | p_fig4_d) /
             (p_fig4_e | p_fig4_spacer) +
  plot_layout(heights = c(1, 1, 1.08))

# Add shared y-axis label and trophic row labels via cowplot
fig4_main <- cowplot::ggdraw() +
  cowplot::draw_plot(fig4_grid, x = 0.07, y = 0, width = 0.93, height = 0.96) +
  # Shared y-axis label
  cowplot::draw_label("MPA / Reference", x = 0.025, y = 0.50,
                      angle = 90, size = 9, fontfamily = "Helvetica") +
  # Trophic row labels
  cowplot::draw_label("Predators", x = 0.045, y = 0.82,
                      angle = 90, size = 8, fontface = "bold",
                      fontfamily = "Helvetica", hjust = 0.5) +
  cowplot::draw_label("Herbivores", x = 0.045, y = 0.50,
                      angle = 90, size = 8, fontface = "bold",
                      fontfamily = "Helvetica", hjust = 0.5) +
  cowplot::draw_label("Producers", x = 0.045, y = 0.17,
                      angle = 90, size = 8, fontface = "bold",
                      fontfamily = "Helvetica", hjust = 0.5)

# Legend is embedded in the bottom-right spacer cell; no separate legend row needed
fig4_final <- fig4_main

save_fig(fig4_final, "fig_04_recovery_curves", FIG4_DIMS["w"], FIG4_DIMS["h"])
cat("  Figure 4 saved: plots/fig_04_recovery_curves.pdf\n")

} # end fig04 inner block

# =============================================================================
# Figure S10 (Supplemental): Biomass + density recovery trajectories (5×2 grid)
# Companion to main text Figure 4 (biomass only). Shows both response types
# side-by-side. M. pyrifera has biomass only (no density metric).
# =============================================================================

if (should_render("fig_s10")) {
cat("Building Figure S10 (Supplemental): Biomass + density recovery grid...\n")

# Reuse data, slopes, and panel function from Fig 4 section above.
# fig5_data_bio, fig5_data_den, fig5_bio_lmer, fig5_den_lmer already computed.

# Helper to strip species title and y-axis from right (density) column panels
strip_for_right <- function(p) {
  p + labs(title = NULL) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
}

# Biomass column (left)
p_s9_bio_1 <- make_fig5_panel("Panulirus interruptus", "a",
                               sp_data = fig5_data_bio,
                               sp_lmer = fig5_bio_lmer,
                               n_mpas_df = fig5_n_mpas_bio,
                               ylims_df = fig5_ylims_all)
p_s9_bio_2 <- make_fig5_panel("Semicossyphus pulcher", "b",
                               sp_data = fig5_data_bio,
                               sp_lmer = fig5_bio_lmer,
                               n_mpas_df = fig5_n_mpas_bio,
                               ylims_df = fig5_ylims_all)
p_s9_bio_3 <- make_fig5_panel("Strongylocentrotus purpuratus", "c",
                               sp_data = fig5_data_bio,
                               sp_lmer = fig5_bio_lmer,
                               n_mpas_df = fig5_n_mpas_bio,
                               ylims_df = fig5_ylims_all)
p_s9_bio_4 <- make_fig5_panel("Mesocentrotus franciscanus", "d",
                               sp_data = fig5_data_bio,
                               sp_lmer = fig5_bio_lmer,
                               n_mpas_df = fig5_n_mpas_bio,
                               ylims_df = fig5_ylims_all)
p_s9_bio_5 <- make_fig5_panel("Macrocystis pyrifera", "e", show_xlab = TRUE,
                               sp_data = fig5_data_bio,
                               sp_lmer = fig5_bio_lmer,
                               n_mpas_df = fig5_n_mpas_bio,
                               ylims_df = fig5_ylims_all)

# Density column (right)
p_s9_den_1 <- strip_for_right(make_fig5_panel("Panulirus interruptus", "f",
                               sp_data = fig5_data_den,
                               sp_lmer = fig5_den_lmer,
                               n_mpas_df = fig5_n_mpas_den,
                               ylims_df = fig5_ylims_all))
p_s9_den_2 <- strip_for_right(make_fig5_panel("Semicossyphus pulcher", "g",
                               sp_data = fig5_data_den,
                               sp_lmer = fig5_den_lmer,
                               n_mpas_df = fig5_n_mpas_den,
                               ylims_df = fig5_ylims_all))
p_s9_den_3 <- strip_for_right(make_fig5_panel("Strongylocentrotus purpuratus", "h",
                               sp_data = fig5_data_den,
                               sp_lmer = fig5_den_lmer,
                               n_mpas_df = fig5_n_mpas_den,
                               ylims_df = fig5_ylims_all))
p_s9_den_4 <- strip_for_right(make_fig5_panel("Mesocentrotus franciscanus", "i",
                               show_xlab = TRUE,
                               sp_data = fig5_data_den,
                               sp_lmer = fig5_den_lmer,
                               n_mpas_df = fig5_n_mpas_den,
                               ylims_df = fig5_ylims_all))

# Empty spacer for bottom-right (M. pyrifera has no density)
p_s9_spacer <- ggplot() +
  annotate("text", x = 0.5, y = 0.5,
           label = "No density metric\nfor M. pyrifera",
           size = 3, color = "grey70", fontface = "italic",
           family = "Helvetica", hjust = 0.5, vjust = 0.5) +
  xlim(0, 1) + ylim(0, 1) +
  theme_void() +
  theme(plot.margin = margin(3, 6, 3, 4, "pt"))

# Assemble 5×2 grid
fig_s9_grid <- (p_s9_bio_1 | p_s9_den_1) /
               (p_s9_bio_2 | p_s9_den_2) /
               (p_s9_bio_3 | p_s9_den_3) /
               (p_s9_bio_4 | p_s9_den_4) /
               (p_s9_bio_5 | p_s9_spacer) +
  plot_layout(heights = c(1, 1, 1, 1, 1))

fig_s9_main <- cowplot::ggdraw() +
  cowplot::draw_plot(fig_s9_grid, x = 0.07, y = 0, width = 0.93, height = 0.97) +
  # Column headers
  cowplot::draw_label("Biomass", x = 0.30, y = 0.99,
                      size = 9, fontface = "bold", fontfamily = "Helvetica") +
  cowplot::draw_label("Density", x = 0.77, y = 0.99,
                      size = 9, fontface = "bold", fontfamily = "Helvetica") +
  # Shared y-axis label
  cowplot::draw_label("MPA / Reference", x = 0.01, y = 0.50,
                      angle = 90, size = 9, fontfamily = "Helvetica") +
  # Trophic row labels
  cowplot::draw_label("Predators", x = 0.045, y = 0.87,
                      angle = 90, size = 8, fontface = "bold",
                      fontfamily = "Helvetica", hjust = 0.5) +
  cowplot::draw_label("Herbivores", x = 0.045, y = 0.49,
                      angle = 90, size = 8, fontface = "bold",
                      fontfamily = "Helvetica", hjust = 0.5) +
  cowplot::draw_label("Producers", x = 0.045, y = 0.10,
                      angle = 90, size = 8, fontface = "bold",
                      fontfamily = "Helvetica", hjust = 0.5)

# Reuse the same legend as Fig 4
fig_s9_legend <- cowplot::ggdraw() +
  cowplot::draw_line(x = c(0.04, 0.12), y = c(0.5, 0.5),
                     color = "grey30", linewidth = 1.2) +
  cowplot::draw_label("Population mean (lmer \u00b1 95% CI)", x = 0.13, y = 0.5,
                      size = 7.5, fontfamily = "Helvetica", hjust = 0) +
  cowplot::draw_line(x = c(0.44, 0.52), y = c(0.5, 0.5),
                     color = "grey55", linewidth = 0.4) +
  cowplot::draw_label("Individual MPA", x = 0.53, y = 0.5,
                      size = 7.5, fontfamily = "Helvetica", hjust = 0) +
  cowplot::draw_line(x = c(0.70, 0.70), y = c(0.2, 0.8),
                     color = "grey50", linewidth = 0.5, linetype = "dotted") +
  cowplot::draw_label("t = 11 years", x = 0.72, y = 0.5,
                      size = 7.5, fontfamily = "Helvetica", hjust = 0)

fig_s9_final <- cowplot::plot_grid(
  fig_s9_main, fig_s9_legend, ncol = 1, rel_heights = c(1, 0.04)
)

save_fig(fig_s9_final, "fig_s10_recovery_bio_den", FIG_S10_DIMS["w"], FIG_S10_DIMS["h"])
cat("  Figure S10 saved: plots/fig_s10_recovery_bio_den.pdf\n")

} # end fig_s10 inner block

} # end shared data prep for fig04/fig_s10


# =============================================================================
# Figure S3 (Supplemental): All taxa log response ratios at example MPAs
# Not in current manuscript - kept as supplemental
# =============================================================================

if (should_render("fig_s03")) {
cat("Building Figure S3 (Supplemental): All taxa time series at example MPAs...\n")

fig_s2_mpas <- c("Naples SMCA", "Scorpion SMR", "Anacapa Island SMR 2003")

# Look up start years from Site table with fallback defaults
fig_s2_starts <- tibble::tibble(
  CA_MPA_Name_Short = fig_s2_mpas,
  MPA_Start = c(2012, 2005, 2005)
) %>%
  dplyr::rows_update(
    Site %>%
      dplyr::filter(CA_MPA_Name_Short %in% fig_s2_mpas) %>%
      dplyr::select(CA_MPA_Name_Short, MPA_Start) %>%
      dplyr::filter(!is.na(MPA_Start)),
    by = "CA_MPA_Name_Short"
  ) %>%
  tibble::deframe()

exclude_species <- c("legal", "sublegal", "all urchins")

# Taxa levels for S3: exclude M. pyrifera (biomass only, no density data)
taxa_levels_s3 <- taxa_levels[taxa_levels != "M. pyrifera"]

fig_s2_data <- All.RR.sub.trans %>%
  dplyr::filter(
    CA_MPA_Name_Short %in% fig_s2_mpas,
    resp == "Den",
    !stringr::str_to_lower(y) %in% exclude_species,
    # M. pyrifera has no density data — exclude to avoid empty facets
    y != "Macrocystis pyrifera"
  ) %>%
  dplyr::mutate(
    source = factor(source, levels = source_levels),
    CA_MPA_Name_Short = factor(CA_MPA_Name_Short, levels = fig_s2_mpas)
  )

# Build species_short column
fig_s2_data <- fig_s2_data %>%
  dplyr::mutate(
    species_short = forcats::fct_recode(
      factor(y),
      "S. purpuratus"   = "Strongylocentrotus purpuratus",
      "M. franciscanus" = "Mesocentrotus franciscanus",
      "P. interruptus"  = "Panulirus interruptus",
      "S. pulcher"      = "Semicossyphus pulcher"
    ) %>%
      forcats::fct_relevel(taxa_levels_s3)
  )

# Short MPA display names
fig_s2_data <- fig_s2_data %>%
  dplyr::mutate(
    site_display = shorten_mpa_name(as.character(CA_MPA_Name_Short)),
    site_display = factor(site_display, levels = shorten_mpa_name(fig_s2_mpas))
  )

# MPA start year lookup
fig_s2_vline_data <- fig_s2_data %>%
  dplyr::distinct(site_display, CA_MPA_Name_Short) %>%
  dplyr::mutate(mpa_start_yr = fig_s2_starts[as.character(CA_MPA_Name_Short)])

# Join start years for before/after shading and trend filtering
fig_s2_data <- fig_s2_data %>%
  dplyr::left_join(
    tibble::enframe(fig_s2_starts,
                    name = "CA_MPA_Name_Short", value = "mpa_start_yr"),
    by = "CA_MPA_Name_Short"
  ) %>%
  dplyr::mutate(period = ifelse(year > mpa_start_yr, "After", "Before"))

# Before-period shading rectangles (one per site column)
fig_s2_shade <- fig_s2_data %>%
  dplyr::group_by(site_display, CA_MPA_Name_Short) %>%
  dplyr::summarise(xmin = min(year, na.rm = TRUE) - 0.5,
                   xmax = fig_s2_starts[as.character(CA_MPA_Name_Short[1])],
                   .groups = "drop")

# Panel tag labels (uses taxa_levels_s3 — M. pyrifera excluded)
s2_tagged_levels <- paste0("(", letters[seq_along(taxa_levels_s3)], ")  ", taxa_levels_s3)
fig_s2_data$species_short <- factor(fig_s2_data$species_short,
  levels = taxa_levels_s3, labels = s2_tagged_levels)
col_taxa_s3 <- col_taxa[taxa_levels_s3]
col_taxa_s2 <- setNames(unname(col_taxa_s3), s2_tagged_levels)

fig_s2 <- ggplot(fig_s2_data, aes(x = year, y = lnDiff,
                                   color = species_short)) +
  # Before-period shading
  geom_rect(data = fig_s2_shade,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "grey92", alpha = 0.6) +
  # Reference line at RR = 1
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30",
             linewidth = 0.5) +
  # MPA implementation vertical line
  geom_vline(data = fig_s2_vline_data,
             aes(xintercept = mpa_start_yr),
             linetype = MPA_LINE_TYPE,
             color = MPA_LINE_COLOR,
             linewidth = MPA_LINE_WIDTH,
             inherit.aes = FALSE) +
  # Post-MPA linear trend
  geom_smooth(data = fig_s2_data %>% dplyr::filter(period == "After"),
              method = "lm", formula = y ~ x, se = TRUE,
              linewidth = 0.8, alpha = 0.15) +
  # Data points — slightly transparent, before period more muted
  geom_point(aes(alpha = period), size = 1.1, shape = 16, na.rm = TRUE) +
  scale_alpha_manual(values = c("Before" = 0.35, "After" = 0.65), guide = "none") +
  scale_color_manual(values = col_taxa_s2, name = "Species") +
  # Free y per row so each species gets its natural range
  facet_grid(species_short ~ site_display, scales = "free") +
  scale_y_rr_free() +
  labs(x = "Year") +
  theme_mpa(base_size = 9) +
  theme(
    legend.position = "none",
    strip.text.x = element_text(size = 9, face = "bold",
                                margin = margin(b = 3, t = 3)),
    strip.text.y = element_text(size = 8.5, face = "italic",
                                margin = margin(l = 3, r = 3)),
    axis.title = element_text(size = 8.5),
    axis.text = element_text(size = 7.5),
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(fill = NA, colour = "grey80", linewidth = 0.4),
    panel.spacing = unit(3, "mm"),
    plot.margin = margin(6, 6, 6, 6)
  )

save_fig(fig_s2, "fig_s03_all_taxa_timeseries", FIG_S3_DIMS["w"], FIG_S3_DIMS["h"])
} # end fig_s03

# (Old Figures S3 and S4 removed — superseded by species-level versions
#  in 10_temporal_analysis.R: fig_s03 through fig_s06)


# =============================================================================
# Figure S7: Model Selection & Heterogeneity
# =============================================================================
# Statistical transparency: variance components and model selection

if (should_render("fig_s07")) {
cat("\n--- Figure S7: Statistical Transparency ---\n")

FIG_S7_DIMS <- c(w = 17, h = 10)  # Conservation Letters max width

# Panel A: Model selection distribution by taxa
model_dist <- SumStats.Final %>%
  group_by(Taxa, Model) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Taxa) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

model_dist$Taxa <- factor(model_dist$Taxa, levels = taxa_levels)

# Model colors — from centralized palette (00b_color_palette.R)

# Pre-compute cumulative positions for Sigmoid labels
sigmoid_labels <- model_dist %>%
  arrange(Taxa, desc(Model)) %>%
  group_by(Taxa) %>%
  mutate(ymax = cumsum(prop), ymin = ymax - prop) %>%
  ungroup() %>%
  filter(Model == "Sigmoid", prop >= 0.03) %>%
  mutate(ymid = (ymin + ymax) / 2,
         label = paste0(round(prop * 100), "%"))

panel_S7A <- ggplot(model_dist, aes(x = Taxa, y = prop, fill = Model)) +
  geom_col(position = "stack", width = 0.7, color = "white", linewidth = 0.3) +
  # Label small Sigmoid slivers so readers know they are intentional

  {if (nrow(sigmoid_labels) > 0) {
    geom_text(
      data = sigmoid_labels,
      aes(x = Taxa, y = ymax, label = label),
      inherit.aes = FALSE,
      size = 2.5, vjust = -0.4, color = "grey25"
    )
  }} +
  scale_fill_manual(
    values = col_model,
    name = "Effect size method",
    guide = guide_legend(nrow = 1)
  ) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  labs(
    x = NULL,
    y = "Proportion of MPAs"
  ) +
  theme_mpa(base_size = 10) +
  theme(
    axis.text.x = element_text(face = "italic", angle = 30, hjust = 1),
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
      values = col_variance,
      name = "Variance component",
      guide = guide_legend(nrow = 1, byrow = TRUE)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    scale_y_discrete(expand = expansion(add = 0.5)) +
    labs(
      x = expression(tau^2),
      y = NULL
    ) +
    theme_mpa(base_size = 10) +
    theme(
      axis.title.x = element_text(size = 10),
      legend.position = "bottom",
      legend.direction = "horizontal",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  fig_s6 <- panel_S7A / panel_S7B +
    # Put the collected guides into an explicit guide row to avoid the large
    # whitespace patchwork can leave when legends are collected implicitly.
    patchwork::guide_area() +
    plot_layout(heights = c(1, 0.65, 0.08), guides = "collect") +
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
    theme(
      plot.tag = element_text(size = 9, face = "bold"),
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
      axis.line.x.bottom = element_line(colour = "black", linewidth = 0.5),
      axis.line.y.left = element_line(colour = "black", linewidth = 0.5)
    )

} else {
  fig_s6 <- panel_S7A +
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
    theme(plot.tag = element_text(size = 9, face = "bold"))
}

save_fig(fig_s6, "fig_s07_statistical_transparency", FIG_S7_DIMS["w"], FIG_S7_DIMS["h"])
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
    resp_label <- ifelse(alt_resp == "Bio",
                         "Biomass MPA / Reference", "Density MPA / Reference")
  } else {
    resp_label <- ifelse(primary_resp == "Bio",
                         "Biomass MPA / Reference", "Density MPA / Reference")
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
    abs(plot_dat$mean_lnRR - 1.96 * plot_dat$se_lnRR),
    abs(plot_dat$mean_lnRR + 1.96 * plot_dat$se_lnRR),
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
    # 95% CI (Conservation Letters requirement)
    geom_ribbon(aes(ymin = mean_lnRR - 1.96 * se_lnRR,
                    ymax = mean_lnRR + 1.96 * se_lnRR),
                fill = taxa_color, alpha = 0.22) +
    geom_line(color = taxa_color, linewidth = 0.6, alpha = 0.8) +
    geom_point(color = taxa_color, size = 1.2, alpha = 0.7) +
    scale_x_continuous(limits = c(year_min, year_max), breaks = year_breaks) +
    coord_cartesian(ylim = c(-y_lim, y_lim)) +
    scale_y_rr(-y_lim, y_lim, name = resp_label) +
    facet_wrap(~ CA_MPA_Name_Short, ncol = n_cols, scales = "fixed",
               labeller = label_wrap_gen(18)) +
    labs(
      title = taxa_i,
      x = "Year",
      y = resp_label
    ) +
    theme_mpa(base_size = 9) +
    theme(
      plot.title = element_text(size = 10, face = "italic",
                                 margin = margin(b = 6)),
      strip.text = element_text(size = 8.5, face = "plain"),
      strip.background = element_blank(),
      axis.text = element_text(size = 8),
      panel.grid.major = element_blank(),
      panel.spacing = unit(0.6, "lines"),
      # Reinforce L-shaped axes
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.line.x.bottom = element_line(colour = "black", linewidth = 0.5),
      axis.line.y.left = element_line(colour = "black", linewidth = 0.5)
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
    # Wrap long MPA names and add dagger to excluded sites
    wrap_fn <- function(x) {
      mapped <- new_labels[x]
      stringr::str_wrap(mapped, width = 18)
    }
    fig_appendix <- fig_appendix +
      facet_wrap(~ CA_MPA_Name_Short, ncol = n_cols, scales = "fixed",
                 labeller = labeller(CA_MPA_Name_Short = wrap_fn))
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


# =============================================================================
# FIGURE 2: Trophic Cascade Case Studies — 3×3 trophic grid
# =============================================================================
# Three Channel Islands MPAs with long time series (19-20 years) that show
# the full trophic cascade: predators ↑, urchins ↓, kelp ↑.
# Layout: 3 rows (Predators / Urchins / Kelp) × 3 columns (sites)
# Each panel shows 1-2 species within the same trophic group.

FIG2_DIMS <- c(w = 17, h = 16)  # 170mm × 160mm, double-column (increased from 130mm to reduce vertical cramping in 3×3 grid)

if (should_render("fig02")) {
cat("\n--- Figure 2 (Main Text): Trophic cascade case studies (3×3 grid) ---\n")

# Sites selected based on cascade consistency analysis (table_s_cascade_consistency.csv)
# and per-MPA slope analysis (table_s_per_mpa_slopes.csv):
#   Scorpion SMR:       4/5 cascade score, strongest kelp recovery slope (+0.39/yr), 20 yr
#   Gull Island SMR:    5/5 cascade score, dramatic kelp effect (lnRR=5.2), 19-20 yr
#   Santa Barbara Island SMR: 5/5 cascade score, all 5 species consistent, 19 yr
cascade_sites <- c("Scorpion SMR", "Gull Island SMR", "Santa Barbara Island SMR")
cascade_site_labels <- c("Scorpion SMR", "Gull Island SMR", "SB Island SMR")

# Map full scientific names to abbreviations
taxa_name_map_fig4 <- c(
  "Strongylocentrotus purpuratus" = "S. purpuratus",
  "Mesocentrotus franciscanus" = "M. franciscanus",
  "Macrocystis pyrifera" = "M. pyrifera",
  "Panulirus interruptus" = "P. interruptus",
  "Semicossyphus pulcher" = "S. pulcher"
)

# Preferred response type per species: density for animals, biomass for kelp
preferred_resp <- c(
  "P. interruptus"  = "Den",
  "S. pulcher"      = "Den",
  "S. purpuratus"   = "Den",
  "M. franciscanus" = "Den",
  "M. pyrifera"     = "Bio"
)

# Trophic groupings for the 3 rows
fig4_trophic_groups <- list(
  Predators = c("P. interruptus", "S. pulcher"),
  Urchins   = c("S. purpuratus", "M. franciscanus"),
  Producers = c("M. pyrifera")
)

# Get MPA implementation years from Site metadata
if (exists("Site")) {
  fig4_impl_df <- Site %>%
    filter(CA_MPA_Name_Short %in% cascade_sites) %>%
    dplyr::select(CA_MPA_Name_Short, MPA_Start) %>%
    distinct()
} else {
  fig4_impl_df <- data.frame(CA_MPA_Name_Short = cascade_sites, MPA_Start = 2003)
}

# Prepare data: filter to cascade sites, compute time since MPA, select preferred response
fig4_data <- All.RR.sub.trans %>%
  mutate(
    Taxa_Short = ifelse(y %in% names(taxa_name_map_fig4),
                        taxa_name_map_fig4[y], y),
    year_num = as.numeric(year)
  ) %>%
  filter(
    CA_MPA_Name_Short %in% cascade_sites,
    Taxa_Short %in% names(preferred_resp)
  ) %>%
  # Keep only the preferred response type per species
  filter(mapply(function(tx, rp) preferred_resp[tx] == rp, Taxa_Short, resp)) %>%
  # Join MPA start year and compute years since implementation
  left_join(fig4_impl_df, by = "CA_MPA_Name_Short") %>%
  mutate(
    time_since = year_num - MPA_Start,
    period = ifelse(time_since <= 0, "Before", "After")
  )

# Compute annual means per site x species x time_since
fig4_annual <- fig4_data %>%
  group_by(CA_MPA_Name_Short, Taxa_Short, time_since, period) %>%
  summarise(
    mean_lnRR = mean(lnDiff, na.rm = TRUE),
    se_lnRR = sd(lnDiff, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(is.finite(mean_lnRR)) %>%
  mutate(se_lnRR = ifelse(is.na(se_lnRR), 0, se_lnRR))

# Factor ordering: trophic order for legend
fig4_annual$Taxa_Short <- factor(fig4_annual$Taxa_Short, levels = taxa_levels)

# Assign trophic group for row-based y-limits
fig4_annual$trophic_group <- trophic_assignment[as.character(fig4_annual$Taxa_Short)]
if (any(is.na(fig4_annual$trophic_group))) {
  warning("Fig 2: Unmapped taxa found: ",
          paste(unique(fig4_annual$Taxa_Short[is.na(fig4_annual$trophic_group)]), collapse = ", "))
}

# Compute y-limits shared per trophic row (not per panel)
fig4_row_ylims <- fig4_annual %>%
  group_by(trophic_group) %>%
  summarise(
    q_lo = min(mean_lnRR - 1.96 * se_lnRR, na.rm = TRUE),
    q_hi = max(mean_lnRR + 1.96 * se_lnRR, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pad = pmax((q_hi - q_lo) * 0.10, 0.3),
    ymin = q_lo - pad,
    ymax = q_hi + pad
  )

# --- Helper: build one panel (site × trophic group) ---
make_fig4_panel <- function(site, species_vec, tag_label,
                             show_xaxis = FALSE, show_yaxis = TRUE) {
  panel_data <- fig4_annual %>%
    filter(CA_MPA_Name_Short == site, Taxa_Short %in% species_vec)

  # Get trophic group y-limits
  tg <- unique(trophic_assignment[species_vec])
  ylims <- fig4_row_ylims[fig4_row_ylims$trophic_group == tg, ]
  if (nrow(ylims) == 0) stop("No y-limits found for trophic group: ", tg)

  # Split before/after for layering
  panel_before <- panel_data %>% filter(period == "Before")
  panel_after  <- panel_data %>% filter(period == "After")

  p <- ggplot(panel_data,
              aes(x = time_since, y = mean_lnRR,
                  color = Taxa_Short, fill = Taxa_Short)) +
    # Reference line at lnRR = 0 (no MPA effect)
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55",
               linewidth = 0.35) +
    # Vertical line at MPA implementation (time = 0)
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40",
               linewidth = 0.45) +
    # Before-MPA period: lighter alpha, 95% CI (Conservation Letters requirement)
    geom_ribbon(data = panel_before,
                aes(ymin = mean_lnRR - 1.96 * se_lnRR,
                    ymax = mean_lnRR + 1.96 * se_lnRR),
                alpha = 0.06, colour = NA) +
    geom_line(data = panel_before,
              aes(group = Taxa_Short), alpha = 0.25, linewidth = 0.4) +
    geom_point(data = panel_before,
               size = 0.7, alpha = 0.30, shape = 16) +
    # After-MPA period: full visual weight, 95% CI (Conservation Letters requirement)
    geom_ribbon(data = panel_after,
                aes(ymin = mean_lnRR - 1.96 * se_lnRR,
                    ymax = mean_lnRR + 1.96 * se_lnRR),
                alpha = 0.12, colour = NA) +
    geom_line(data = panel_after,
              aes(group = Taxa_Short), alpha = 0.55, linewidth = 0.5) +
    geom_point(data = panel_after,
               size = 1.0, alpha = 0.7, shape = 16) +
    # Linear trend overlay: after-period only (the cascade signal)
    # 95% CI shown per Conservation Letters requirement
    geom_smooth(data = panel_after,
                method = "lm", formula = y ~ x, se = TRUE, level = 0.95,
                linewidth = 1.1, alpha = 0.2) +
    # Scales
    scale_color_taxa(name = NULL) +
    scale_fill_taxa(name = NULL) +
    scale_x_continuous(breaks = seq(-20, 20, by = 10),
                       expand = c(0.02, 0)) +
    coord_cartesian(ylim = c(ylims$ymin, ylims$ymax), clip = "on") +
    # RR-scaled y-axis: data on lnRR scale, labels show back-transformed RR.
    # RR = 1 (lnRR = 0) = no MPA effect. Shared y-axis label added via cowplot.
    scale_y_rr(ylims$ymin, ylims$ymax, name = NULL, sparse = TRUE) +
    labs(x = NULL, y = NULL) +
    theme_mpa(base_size = 9) +
    theme(
      legend.position = "none",
      plot.tag = element_text(size = 8, face = "bold"),
      plot.tag.position = c(0, 1),
      plot.margin = margin(3, 6, 3, 4, "pt")
    ) +
    labs(tag = paste0("(", tag_label, ")"))

  # Strip y-axis for non-leftmost columns
 if (!show_yaxis) {
    p <- p + theme(axis.text.y = element_blank(),
                   axis.ticks.y = element_blank())
  }

  # Strip x-axis text for non-bottom rows
  if (!show_xaxis) {
    p <- p + theme(axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())
  }

  p
}

# --- Build 9 panels: iterate rows (trophic) × columns (sites) ---
tag_letters <- letters[1:9]
tag_idx <- 1
fig4_panels <- list()

for (row_i in seq_along(fig4_trophic_groups)) {
  group_name <- names(fig4_trophic_groups)[row_i]
  species_vec <- fig4_trophic_groups[[row_i]]
  is_bottom <- (row_i == length(fig4_trophic_groups))

  for (col_j in seq_along(cascade_sites)) {
    site <- cascade_sites[col_j]
    is_left <- (col_j == 1)

    panel_key <- paste0("r", row_i, "c", col_j)
    fig4_panels[[panel_key]] <- make_fig4_panel(
      site        = site,
      species_vec = species_vec,
      tag_label   = tag_letters[tag_idx],
      show_xaxis  = is_bottom,
      show_yaxis  = is_left
    )
    tag_idx <- tag_idx + 1
  }
}

# --- Assemble with patchwork ---
fig4_grid <- (fig4_panels[["r1c1"]] | fig4_panels[["r1c2"]] | fig4_panels[["r1c3"]]) /
             (fig4_panels[["r2c1"]] | fig4_panels[["r2c2"]] | fig4_panels[["r2c3"]]) /
             (fig4_panels[["r3c1"]] | fig4_panels[["r3c2"]] | fig4_panels[["r3c3"]]) +
  plot_layout(heights = c(1, 1, 1))

# --- Add shared labels, column headers, and row labels via cowplot ---
fig4_main <- cowplot::ggdraw() +
  cowplot::draw_plot(fig4_grid, x = 0.07, y = 0.04, width = 0.93, height = 0.90) +
  # Shared y-axis label
  cowplot::draw_label("MPA / Reference", x = 0.025, y = 0.50,
                      angle = 90, size = 9, fontfamily = "Helvetica") +
  # Trophic row labels (left margin)
  cowplot::draw_label("Predators", x = 0.050, y = 0.79,
                      angle = 90, size = 8, fontface = "bold",
                      fontfamily = "Helvetica", hjust = 0.5) +
  cowplot::draw_label("Urchins", x = 0.050, y = 0.49,
                      angle = 90, size = 8, fontface = "bold",
                      fontfamily = "Helvetica", hjust = 0.5) +
  cowplot::draw_label("Producers", x = 0.050, y = 0.19,
                      angle = 90, size = 8, fontface = "bold",
                      fontfamily = "Helvetica", hjust = 0.5) +
  # Column headers (site names, bold, above panels)
  cowplot::draw_label(cascade_site_labels[1], x = 0.23, y = 0.96,
                      size = 9, fontface = "bold", fontfamily = "Helvetica") +
  cowplot::draw_label(cascade_site_labels[2], x = 0.54, y = 0.96,
                      size = 9, fontface = "bold", fontfamily = "Helvetica") +
  cowplot::draw_label(cascade_site_labels[3], x = 0.84, y = 0.96,
                      size = 9, fontface = "bold", fontfamily = "Helvetica") +
  # Shared x-axis label
  cowplot::draw_label("Years since MPA implementation", x = 0.54, y = 0.015,
                      size = 9, fontfamily = "Helvetica")

# --- Species legend at bottom ---
fig4_legend_data <- data.frame(
  x = 1:5,
  y = rep(0, 5),
  Taxa_Short = factor(taxa_levels, levels = taxa_levels)
)
fig4_legend_plot <- ggplot(fig4_legend_data,
                            aes(x = x, y = y, color = Taxa_Short)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1.2) +
  scale_color_taxa(name = NULL) +
  guides(color = guide_legend(
    nrow = 1,
    override.aes = list(linewidth = 1.5, alpha = 1, size = 2)
  )) +
  theme_mpa(base_size = 9) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(face = "italic", size = 8),
    legend.key.width = unit(10, "mm"),
    legend.key.height = unit(3, "mm"),
    legend.spacing.x = unit(2, "mm"),
    legend.box.margin = margin(t = 0)
  )
fig4_legend_grob <- cowplot::get_legend(fig4_legend_plot)

fig4_final <- cowplot::plot_grid(
  fig4_main, fig4_legend_grob,
  ncol = 1, rel_heights = c(1, 0.06)
)

save_fig(fig4_final, "fig_02_cascade_case_studies", FIG2_DIMS["w"], FIG2_DIMS["h"])
cat("  Figure 2 saved: plots/fig_02_cascade_case_studies.pdf\n")
} # end fig02


# =============================================================================
# FIGURE S11: DHARMa Model Diagnostics (4-panel)
# =============================================================================
# Reads data/model_diagnostics.csv produced by 08_effect_sizes.R and creates
# a 2×2 panel figure summarizing residual diagnostics for all effect-size models (pBACIPS + NLS).

if (should_render("fig_s11")) {
cat("\n--- Figure S11: DHARMa model diagnostics ---\n")

diag_path <- here::here("data", "model_diagnostics.csv")
if (!file.exists(diag_path)) {
  warning("model_diagnostics.csv not found — skipping fig_s11")
} else {

diag <- read.csv(diag_path, stringsAsFactors = FALSE)
diag$Taxa <- factor(diag$Taxa, levels = taxa_levels)

# --- Classify failure type for panel (b) ---
diag$fail_type <- dplyr::case_when(
  grepl("non-normal", diag$Notes) & grepl("heteroscedastic", diag$Notes) ~ "Both",
  grepl("non-normal", diag$Notes)       ~ "Non-normal",
  grepl("heteroscedastic", diag$Notes)  ~ "Heteroscedastic",
  TRUE                                   ~ "OK"
)
diag$fail_type <- factor(diag$fail_type,
                          levels = c("OK", "Non-normal", "Heteroscedastic", "Both"))

n_models <- nrow(diag)

# -------------------------------------------------------------------------
# Panel (a): Pass rate by diagnostic test — Grouped bar chart
# -------------------------------------------------------------------------
# Pivot the four p-value columns to long form, classify pass/fail at α=0.05
# Heteroscedasticity (|cor| > 0.5) is appended separately below.
diag_long <- tidyr::pivot_longer(
  diag,
  cols = c(Uniformity_p, Dispersion_p, Outlier_p, Shapiro_p),
  names_to = "test_raw", values_to = "p_value"
)
# Heteroscedasticity uses |cor| > 0.5 threshold, not a p-value test
diag_hetero <- diag %>%
  dplyr::mutate(test_raw = "Hetero_Cor",
                p_value  = NA_real_,
                pass     = abs(Hetero_Cor) <= 0.5)
diag_long$pass <- diag_long$p_value > 0.05
# Count NAs (tests that couldn't be computed) — these are excluded from pass rates
na_pass <- sum(is.na(diag_long$pass)) + sum(is.na(diag_hetero$pass))
if (na_pass > 0) {
  cat("  Note:", na_pass, "diagnostic tests returned NA (excluded from pass rates)\n")
}
diag_long <- dplyr::bind_rows(
  diag_long %>% dplyr::select(Taxa, test_raw, pass),
  diag_hetero %>% dplyr::select(Taxa, test_raw, pass)
)

# Clean test names
test_labels <- c(
  "Uniformity_p" = "Uniformity",
  "Dispersion_p" = "Dispersion",
  "Outlier_p"    = "Outlier",
  "Shapiro_p"    = "Normality",
  "Hetero_Cor"   = "Homoscedasticity"
)
diag_long$Test <- test_labels[diag_long$test_raw]

# Order tests by overall pass rate (highest first)
test_pass_rate <- diag_long %>%
  dplyr::group_by(Test) %>%
  dplyr::summarise(rate = mean(pass, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(rate))
diag_long$Test <- factor(diag_long$Test, levels = test_pass_rate$Test)

# Summarise per taxa × test
pass_summary <- diag_long %>%
  dplyr::group_by(Taxa, Test) %>%
  dplyr::summarise(pass_rate = mean(pass, na.rm = TRUE) * 100,
                   .groups = "drop")

pa <- ggplot(pass_summary, aes(x = Test, y = pass_rate, fill = Taxa)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_hline(yintercept = 95, linetype = "dashed", color = "grey40",
             linewidth = 0.4) +
  scale_fill_taxa(name = NULL) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Proportion passing (%)") +
  theme_mpa(base_size = 9) +
  theme(
    axis.text.x  = element_text(size = 7.5, angle = 30, hjust = 1),
    legend.position = "none"
  )

# -------------------------------------------------------------------------
# Panel (b): Failure type breakdown — Stacked count bar
# -------------------------------------------------------------------------
# Reverse taxa levels so coord_flip() puts predators at top
diag_b <- diag
diag_b$Taxa <- factor(diag_b$Taxa, levels = rev(taxa_levels))
fail_counts <- diag_b %>%
  dplyr::count(Taxa, fail_type) %>%
  dplyr::group_by(Taxa) %>%
  dplyr::mutate(total = sum(n)) %>%
  dplyr::ungroup()

pb <- ggplot(fail_counts, aes(x = Taxa, y = n, fill = fail_type)) +
  geom_col(position = "stack", width = 0.65) +
  # Count labels inside segments (only for segments with room)
  geom_text(aes(label = ifelse(n >= 2, n, "")),
            position = position_stack(vjust = 0.5),
            size = 2.5, color = "grey20") +
  # Total n at right margin
  geom_text(aes(x = Taxa, y = total, label = paste0("n=", total)),
            inherit.aes = FALSE,
            hjust = -0.15, size = 2.5, color = "grey30",
            data = fail_counts %>%
              dplyr::distinct(Taxa, total)) +
  scale_fill_manual(values = col_diag, name = "Diagnostic result") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  coord_flip() +
  labs(x = NULL, y = "Number of models") +
  theme_mpa(base_size = 9) +
  theme(
    axis.text.y   = element_text(face = "italic", size = 8),
    legend.position = "none"
  )

# -------------------------------------------------------------------------
# Panel (c): Shapiro-Wilk p-value distribution — Violin + jitter
# -------------------------------------------------------------------------
diag$shapiro_pass <- dplyr::case_when(
  is.na(diag$Shapiro_p) ~ "Not tested",
  diag$Shapiro_p > 0.05 ~ "Pass",
  TRUE                   ~ "Fail"
)

pc <- ggplot(diag, aes(x = Taxa, y = Shapiro_p)) +
  geom_violin(aes(fill = Taxa), alpha = 0.18, color = NA,
              scale = "width", width = 0.7) +
  geom_jitter(aes(color = Taxa, shape = shapiro_pass),
              width = 0.15, size = 1.3, alpha = 0.7) +
  geom_hline(yintercept = 0.05, linetype = "dashed",
             color = "#D55E00", linewidth = 0.4) +
  scale_fill_taxa(name = NULL) +
  scale_color_taxa(name = NULL) +
  scale_shape_manual(values = c("Pass" = 16, "Fail" = 17, "Not tested" = 4),
                     name = "Shapiro-Wilk") +
  scale_y_sqrt(breaks = c(0, 0.05, 0.25, 0.5, 0.75, 1.0),
               limits = c(0, 1.05)) +
  labs(x = NULL, y = "Shapiro-Wilk p-value") +
  theme_mpa(base_size = 9) +
  theme(
    axis.text.x  = element_text(face = "italic", size = 7, angle = 30,
                                hjust = 1),
    legend.position = "none"
  )

# -------------------------------------------------------------------------
# Panel (d): R² distribution — Boxplot + jitter
# -------------------------------------------------------------------------
if (!"Pass_All" %in% names(diag)) {
  warning("Fig S11: 'Pass_All' column not found in model_diagnostics.csv")
  diag$Pass_All <- NA
}
diag$Pass_All <- as.logical(diag$Pass_All)
diag$overall_result <- ifelse(diag$Pass_All, "Pass", "Fail")

pd <- ggplot(diag, aes(x = Taxa, y = R_Squared)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, fill = "grey92",
               color = "grey50", linewidth = 0.3) +
  geom_jitter(aes(color = Taxa, shape = overall_result),
              width = 0.15, size = 1.3, alpha = 0.7) +
  scale_color_taxa(name = NULL) +
  scale_shape_manual(values = c("Pass" = 16, "Fail" = 1),
                     name = "All diagnostics") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = NULL, y = expression(R^2)) +
  theme_mpa(base_size = 9) +
  theme(
    axis.text.x  = element_text(face = "italic", size = 7, angle = 30,
                                hjust = 1),
    legend.position = "none"
  )

# -------------------------------------------------------------------------
# Assemble 2×2 with patchwork
# -------------------------------------------------------------------------
fig_s11 <- (pa + pb) / (pc + pd) +
  patchwork::plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

# --- Shared legends at bottom (3 rows for clarity) ---
leg_theme <- theme_mpa(base_size = 9) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7.5, face = "bold"),
        legend.key.size = unit(3, "mm"),
        legend.spacing.x = unit(1, "mm"),
        legend.margin = margin(0, 0, 0, 0))

# Row 1: Taxa colors (full width)
taxa_leg <- cowplot::get_legend(
  ggplot(diag, aes(x = Taxa, y = Shapiro_p, color = Taxa)) +
    geom_point(size = 2) +
    scale_color_taxa(name = NULL) +
    guides(color = guide_legend(nrow = 1)) +
    leg_theme +
    theme(legend.text = element_text(face = "italic", size = 7))
)

# Row 2: Diagnostic result fills (full width)
diag_leg <- cowplot::get_legend(
  ggplot(fail_counts, aes(x = Taxa, y = n, fill = fail_type)) +
    geom_col() +
    scale_fill_manual(values = col_diag, name = "Diagnostic result") +
    guides(fill = guide_legend(nrow = 1)) +
    leg_theme
)

# Row 3: Shape legends side by side
shape_leg_sw <- cowplot::get_legend(
  ggplot(diag, aes(x = Taxa, y = Shapiro_p, shape = shapiro_pass)) +
    geom_point(size = 2, color = "grey30") +
    scale_shape_manual(values = c("Pass" = 16, "Fail" = 17),
                       name = "Normality (panel c)") +
    guides(shape = guide_legend(nrow = 1)) +
    leg_theme
)
shape_leg_all <- cowplot::get_legend(
  ggplot(diag, aes(x = Taxa, y = R_Squared, shape = overall_result)) +
    geom_point(size = 2, color = "grey30") +
    scale_shape_manual(values = c("Pass" = 16, "Fail" = 1),
                       name = "All diagnostics (panel d)") +
    guides(shape = guide_legend(nrow = 1)) +
    leg_theme
)
legend_row3 <- cowplot::plot_grid(shape_leg_sw, shape_leg_all,
                                   nrow = 1, rel_widths = c(1, 1))

legend_block <- cowplot::plot_grid(taxa_leg, diag_leg, legend_row3,
                                    ncol = 1, rel_heights = c(1, 1, 1))

fig_s11_final <- cowplot::plot_grid(fig_s11, legend_block,
                                     ncol = 1, rel_heights = c(1, 0.15))

save_fig(fig_s11_final, "fig_s11_dharma_diagnostics", w = 17, h = 15)
cat("  Figure S11 saved: plots/fig_s11_dharma_diagnostics.pdf\n")
cat("  Models:", n_models, "| Pass all:", sum(diag$Pass_All),
    "(", round(100 * mean(diag$Pass_All)), "%)\n")

} # end file exists check
} # end fig_s11


# =============================================================================
# FIGURE S12: Meta-Analysis Publication Bias — Funnel Plots
# =============================================================================
# Assesses publication bias of the multilevel meta-analysis.
# Panel a/b: Funnel plots (biomass/density) with Egger's regression test.

if (should_render("fig_s12")) {
cat("\n--- Figure S12: Meta-analysis funnel plots & leave-one-out ---\n")

if (!exists("meta_biomass_full") || !exists("meta_density_full")) {
  warning("meta_biomass_full / meta_density_full not found — skipping fig_s12")
} else {

  # --- Helper: build funnel data from an rma.mv model ---
  build_funnel_df <- function(model) {
    yi <- model$yi
    sei <- sqrt(model$vi)
    # Extract taxa from the moderator matrix (cell-means: columns are TaxaXXX)
    mod_mat <- model.matrix(model)
    taxa_names <- sub("^Taxa", "", colnames(mod_mat))
    taxa_idx <- apply(mod_mat, 1, which.max)
    taxa <- taxa_names[taxa_idx]
    unrecognized <- setdiff(unique(taxa), taxa_levels)
    if (length(unrecognized) > 0) {
      warning("Fig S12: Unrecognized taxa in funnel data: ",
              paste(unrecognized, collapse = ", "))
    }
    data.frame(yi = yi, sei = sei, taxa = taxa, stringsAsFactors = FALSE)
  }

  # --- Panel (a): Funnel plot — Biomass ---
  funnel_bio <- build_funnel_df(meta_biomass_full)
  funnel_bio$taxa <- factor(funnel_bio$taxa, levels = taxa_levels)

  # Egger-style regression: regress residuals on SE. Done manually because
  # metafor::regtest() does not support rma.mv (multilevel) models.
  # Significant slope suggests small-study bias.
  egger_bio <- tryCatch({
    resid_bio <- residuals(meta_biomass_full)
    sei_bio   <- sqrt(meta_biomass_full$vi)
    eg_fit    <- lm(resid_bio ~ sei_bio)
    eg_summ   <- summary(eg_fit)
    list(t = eg_summ$coefficients["sei_bio", "t value"],
         p = eg_summ$coefficients["sei_bio", "Pr(>|t|)"])
  }, error = function(e) { message("  Egger test (biomass) failed: ", e$message); NULL })
  egger_bio_text <- if (!is.null(egger_bio)) {
    p_fmt <- if (egger_bio$p < 0.001) "p < 0.001" else sprintf("p = %.3f", egger_bio$p)
    sprintf("Egger's test: t = %.2f, %s", egger_bio$t, p_fmt)
  } else "Egger's test: N/A"

  p_funnel_bio <- ggplot(funnel_bio, aes(x = yi, y = sei)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_point(aes(color = taxa), size = 1.5, alpha = 0.7) +
    scale_color_taxa(name = "Taxa") +
    scale_y_reverse(name = "Standard error") +
    labs(x = "Effect size (lnRR)", subtitle = "Biomass") +
    annotate("text", x = Inf, y = Inf, label = egger_bio_text,
             hjust = 1.05, vjust = -0.5, size = 2.5, color = "grey30") +
    theme_mpa() +
    theme(legend.position = "none")

  # --- Panel (b): Funnel plot — Density ---
  funnel_den <- build_funnel_df(meta_density_full)
  funnel_den$taxa <- factor(funnel_den$taxa, levels = taxa_levels)

  egger_den <- tryCatch({
    resid_den <- residuals(meta_density_full)
    sei_den   <- sqrt(meta_density_full$vi)
    eg_fit    <- lm(resid_den ~ sei_den)
    eg_summ   <- summary(eg_fit)
    list(t = eg_summ$coefficients["sei_den", "t value"],
         p = eg_summ$coefficients["sei_den", "Pr(>|t|)"])
  }, error = function(e) { message("  Egger test (density) failed: ", e$message); NULL })
  egger_den_text <- if (!is.null(egger_den)) {
    p_fmt <- if (egger_den$p < 0.001) "p < 0.001" else sprintf("p = %.3f", egger_den$p)
    sprintf("Egger's test: t = %.2f, %s", egger_den$t, p_fmt)
  } else "Egger's test: N/A"

  p_funnel_den <- ggplot(funnel_den, aes(x = yi, y = sei)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_point(aes(color = taxa), size = 1.5, alpha = 0.7) +
    scale_color_taxa(name = "Taxa") +
    scale_y_reverse(name = "Standard error") +
    labs(x = "Effect size (lnRR)", subtitle = "Density") +
    annotate("text", x = Inf, y = Inf, label = egger_den_text,
             hjust = 1.05, vjust = -0.5, size = 2.5, color = "grey30") +
    theme_mpa() +
    theme(legend.position = "none")

  # --- Assemble side-by-side funnel plots with shared legend ---
  legend_plot <- ggplot(funnel_bio, aes(x = yi, y = sei, color = taxa)) +
    geom_point() + scale_color_taxa(name = "Taxa") +
    theme_mpa() +
    theme(legend.position = "bottom",
          legend.text = element_text(face = "italic", size = 7))
  shared_legend <- cowplot::get_legend(legend_plot)

  fig_s12 <- (p_funnel_bio + p_funnel_den) +
    patchwork::plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold", size = 10))

  fig_s12_final <- cowplot::plot_grid(
    fig_s12, shared_legend, ncol = 1, rel_heights = c(1, 0.08)
  )

  save_fig(fig_s12_final, "fig_s12_funnel_plots", w = 17, h = 8)
  cat("  Figure S12 saved.\n")

} # end meta model check
} # end fig_s12


# =============================================================================
# FIGURE S13: lmer Temporal Model Residual Diagnostics
# =============================================================================
# 4-panel residual diagnostics for the pooled lmer temporal meta-regression.
# Data read from outputs/lmer_residual_diagnostics.csv (generated by
# 10_temporal_analysis.R).

if (should_render("fig_s13")) {
cat("\n--- Figure S13: lmer residual diagnostics ---\n")

resid_path  <- here::here("outputs", "lmer_residual_diagnostics.csv")
ranef_path  <- here::here("outputs", "lmer_ranef_diagnostics.csv")

if (!file.exists(resid_path)) {
  warning("lmer_residual_diagnostics.csv not found — skipping fig_s13. ",
          "Run 10_temporal_analysis.R first.")
} else {

  resid_df <- read.csv(resid_path, stringsAsFactors = FALSE)
  # Map full species names to abbreviated names matching col_taxa keys
  sp_map <- c(
    "Panulirus interruptus"          = "P. interruptus",
    "Semicossyphus pulcher"          = "S. pulcher",
    "Strongylocentrotus purpuratus"  = "S. purpuratus",
    "Mesocentrotus franciscanus"     = "M. franciscanus",
    "Macrocystis pyrifera"           = "M. pyrifera"
  )
  unmapped <- sum(!resid_df$species %in% names(sp_map))
  if (unmapped > 0) {
    warning("Fig S13: ", unmapped, " observations have unrecognized species names: ",
            paste(unique(resid_df$species[!resid_df$species %in% names(sp_map)]), collapse = ", "))
  }
  resid_df$species <- sp_map[resid_df$species]
  # Use pooled model for main diagnostics
  pooled <- resid_df[resid_df$model == "pooled", ]
  if (nrow(pooled) == 0) {
    warning("Fig S13: No 'pooled' model data found. Available models: ",
            paste(unique(resid_df$model), collapse = ", "), ". Skipping figure.")
  } else {
  pooled$species <- factor(pooled$species, levels = names(col_taxa))

  # --- Panel (a): Residuals vs. Fitted ---
  p_resid_fit <- ggplot(pooled, aes(x = fitted, y = residual)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_point(aes(color = species), size = 0.8, alpha = 0.4) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
                color = "#B85A4C", linewidth = 0.6) +
    scale_color_taxa(name = "Taxa") +
    labs(x = "Fitted values", y = "Residuals", subtitle = "Residuals vs. fitted") +
    theme_mpa() +
    theme(legend.position = "none")

  # --- Panel (b): QQ Plot ---
  p_qq <- ggplot(pooled, aes(sample = std_residual)) +
    stat_qq(aes(color = species), size = 0.8, alpha = 0.4) +
    stat_qq_line(color = "black", linewidth = 0.4) +
    scale_color_taxa(name = "Taxa") +
    labs(x = "Theoretical quantiles", y = "Standardized residuals",
         subtitle = "Normal Q-Q") +
    theme_mpa() +
    theme(legend.position = "none")

  # --- Panel (c): Scale-Location ---
  pooled$sqrt_abs_std <- sqrt(abs(pooled$std_residual))
  p_scale_loc <- ggplot(pooled, aes(x = fitted, y = sqrt_abs_std)) +
    geom_point(aes(color = species), size = 0.8, alpha = 0.4) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
                color = "#B85A4C", linewidth = 0.6) +
    scale_color_taxa(name = "Taxa") +
    labs(x = "Fitted values", y = expression(sqrt("|Standardized residuals|")),
         subtitle = "Scale-location") +
    theme_mpa() +
    theme(legend.position = "none")

  # --- Panel (d): Random Effects QQ ---
  if (file.exists(ranef_path)) {
    ranef_df <- read.csv(ranef_path, stringsAsFactors = FALSE)
    ranef_pooled <- ranef_df[ranef_df$model == "pooled", ]
    # The intercept column may be named "(Intercept)" or "X.Intercept."
    int_col <- grep("Intercept", names(ranef_pooled), value = TRUE)[1]
    if (!is.na(int_col)) {
      ranef_pooled$ranef_int <- ranef_pooled[[int_col]]
      # Label extreme MPAs (beyond 1.5 IQR)
      q <- quantile(ranef_pooled$ranef_int, c(0.25, 0.75), na.rm = TRUE)
      iqr_val <- diff(q)
      ranef_pooled$label <- ifelse(
        ranef_pooled$ranef_int < q[1] - 1.5 * iqr_val |
        ranef_pooled$ranef_int > q[2] + 1.5 * iqr_val,
        shorten_mpa_name(ranef_pooled$MPA), ""
      )

      p_ranef <- ggplot(ranef_pooled, aes(sample = ranef_int)) +
        stat_qq(size = 1.5, color = "#2A7B8E") +
        stat_qq_line(color = "black", linewidth = 0.4) +
        ggrepel::geom_text_repel(
          aes(x = qnorm(ppoints(nrow(ranef_pooled)))[order(order(ranef_int))],
              y = ranef_int, label = label),
          size = 2, max.overlaps = 10, segment.color = "grey60"
        ) +
        labs(x = "Theoretical quantiles", y = "MPA random intercept",
             subtitle = "Random effects Q-Q") +
        theme_mpa()
    } else {
      warning("Fig S13 panel (d): No intercept column found in ranef CSV. ",
              "Available columns: ", paste(names(ranef_pooled), collapse = ", "))
      p_ranef <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "Random effects\nnot available") +
        theme_void()
    }
  } else {
    warning("Fig S13 panel (d): ranef CSV not found at ", ranef_path)
    p_ranef <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "Random effects\nnot available") +
      theme_void()
  }

  # --- Assemble 2×2 with shared legend ---
  legend_plot <- ggplot(pooled, aes(x = fitted, y = residual, color = species)) +
    geom_point() + scale_color_taxa(name = "Taxa") +
    theme_mpa() +
    theme(legend.position = "bottom",
          legend.text = element_text(face = "italic", size = 7))
  shared_legend <- cowplot::get_legend(legend_plot)

  fig_s13 <- (p_resid_fit + p_qq) / (p_scale_loc + p_ranef) +
    patchwork::plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold", size = 10))

  fig_s13_final <- cowplot::plot_grid(
    fig_s13, shared_legend, ncol = 1, rel_heights = c(1, 0.06)
  )

  save_fig(fig_s13_final, "fig_s13_lmer_residuals", w = 17, h = 15)
  cat("  Figure S13 saved.\n")

  } # end pooled nrow check
} # end file check
} # end fig_s13


# =============================================================================
# FIGURE S14: NLS Model Selection & Durbin-Watson Autocorrelation
# =============================================================================
# Panel a: Stacked bar of model selection frequencies by taxa.
# Panel b: Dot plot of Durbin-Watson statistics by taxa.
# Data from tables/table_model_selection.csv.

if (should_render("fig_s14")) {
cat("\n--- Figure S14: NLS model selection & DW autocorrelation ---\n")

ms_path <- here::here("tables", "table_model_selection.csv")
if (!file.exists(ms_path)) {
  warning("table_model_selection.csv not found — skipping fig_s14")
} else {

  ms <- read.csv(ms_path, stringsAsFactors = FALSE)
  ms$Taxa <- factor(ms$Taxa, levels = rev(taxa_levels))
  expected_models <- c("Linear", "Mean", "Sigmoid")
  unexpected <- setdiff(unique(ms$Model), expected_models)
  if (length(unexpected) > 0) {
    warning("Fig S14: Unexpected model types will be excluded: ",
            paste(unexpected, collapse = ", "))
  }
  ms$Model <- factor(ms$Model, levels = expected_models)

  # --- Panel (a): Model selection frequencies ---
  ms_counts <- ms %>%
    dplyr::count(Taxa, Model)

  p_model_bar <- ggplot(ms_counts, aes(x = Taxa, y = n, fill = Model)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = n), position = position_stack(vjust = 0.5),
              size = 2.3, color = "grey10") +
    coord_flip() +
    scale_fill_manual(values = col_model[c("Linear", "Mean", "Sigmoid")], name = "NLS model") +
    labs(x = NULL, y = "Number of effect sizes", subtitle = "Model selection") +
    theme_mpa() +
    theme(axis.text.y = element_text(face = "italic"),
          legend.position = "bottom")

  # --- Panel (b): Durbin-Watson statistic distribution ---
  ms_dw <- ms[!is.na(ms$DW_stat), ]
  ms_dw$Taxa <- factor(ms_dw$Taxa, levels = rev(taxa_levels))

  p_dw <- ggplot(ms_dw, aes(x = DW_stat, y = Taxa)) +
    annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
             fill = "grey90", alpha = 0.5) +
    geom_vline(xintercept = 2, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_jitter(aes(color = Taxa, shape = Model), height = 0.15, size = 1.5, alpha = 0.7) +
    scale_color_taxa(name = "Taxa", guide = "none") +
    scale_shape_manual(values = c("Linear" = 16, "Mean" = 15, "Sigmoid" = 17),
                       name = "Model type (shape)") +
    labs(x = "Durbin-Watson statistic", y = NULL,
         subtitle = "Temporal autocorrelation") +
    theme_mpa() +
    theme(axis.text.y = element_text(face = "italic"),
          legend.position = "bottom")

  # --- Assemble side-by-side ---
  fig_s14 <- p_model_bar + p_dw +
    patchwork::plot_layout(widths = c(1, 1.2), guides = "collect") +
    patchwork::plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold", size = 10),
          legend.position = "bottom")

  save_fig(fig_s14, "fig_s14_model_selection", w = 17, h = 10)
  cat("  Figure S14 saved.\n")

} # end file check
} # end fig_s14


# =============================================================================
# FIGURE S15: Sensitivity & Robustness — Cook's Distance + Outlier Methods
# =============================================================================
# Panel a/b: Cook's distance for each observation (biomass/density) with
#   threshold line at 4/n.
# Panel c: Outlier sensitivity comparison across 4 methods.
# Data from outputs/filter_audit_meta_analysis.csv and
# tables/table_s_outlier_sensitivity.csv.

if (should_render("fig_s15")) {
cat("\n--- Figure S15: Sensitivity & robustness summary ---\n")

audit_path <- here::here("outputs", "filter_audit_meta_analysis.csv")
outlier_path <- here::here("tables", "table_s_outlier_sensitivity.csv")

if (!file.exists(audit_path) || !file.exists(outlier_path)) {
  warning("Required CSV files not found — skipping fig_s15. ",
          "Need: outputs/filter_audit_meta_analysis.csv, tables/table_s_outlier_sensitivity.csv")
} else {

  audit <- read.csv(audit_path, stringsAsFactors = FALSE)
  audit$Taxa <- factor(audit$Taxa, levels = taxa_levels)

  # Split by response
  audit_bio <- audit[audit$Response == "Biomass", ]
  audit_den <- audit[audit$Response == "Density", ]

  # --- Panel (a): Cook's distance — Biomass ---
  if (nrow(audit_bio) > 0 && "Cooks_Distance" %in% names(audit_bio)) {
    audit_bio$obs <- seq_len(nrow(audit_bio))
    threshold_bio <- 4 / nrow(audit_bio)
    audit_bio$above <- audit_bio$Cooks_Distance > threshold_bio
    audit_bio$label <- ifelse(audit_bio$above,
                              shorten_mpa_name(audit_bio$MPA), "")

    p_cooks_bio <- ggplot(audit_bio, aes(x = obs, y = Cooks_Distance)) +
      geom_hline(yintercept = threshold_bio, linetype = "dashed",
                 color = "#D55E00", linewidth = 0.3) +
      geom_point(aes(color = Taxa), size = 1.3, alpha = 0.7) +
      ggrepel::geom_text_repel(aes(label = label), size = 2,
                                max.overlaps = 8, segment.color = "grey60") +
      scale_color_taxa(name = "Taxa") +
      labs(x = "Observation", y = "Cook's distance", subtitle = "Biomass") +
      theme_mpa() +
      theme(legend.position = "none")
  } else {
    p_cooks_bio <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No Cook's D data") +
      theme_void()
  }

  # --- Panel (b): Cook's distance — Density ---
  if (nrow(audit_den) > 0 && "Cooks_Distance" %in% names(audit_den)) {
    audit_den$obs <- seq_len(nrow(audit_den))
    threshold_den <- 4 / nrow(audit_den)
    audit_den$above <- audit_den$Cooks_Distance > threshold_den
    audit_den$label <- ifelse(audit_den$above,
                              shorten_mpa_name(audit_den$MPA), "")

    p_cooks_den <- ggplot(audit_den, aes(x = obs, y = Cooks_Distance)) +
      geom_hline(yintercept = threshold_den, linetype = "dashed",
                 color = "#D55E00", linewidth = 0.3) +
      geom_point(aes(color = Taxa), size = 1.3, alpha = 0.7) +
      ggrepel::geom_text_repel(aes(label = label), size = 2,
                                max.overlaps = 8, segment.color = "grey60") +
      scale_color_taxa(name = "Taxa") +
      labs(x = "Observation", y = "Cook's distance", subtitle = "Density") +
      theme_mpa() +
      theme(legend.position = "none")
  } else {
    p_cooks_den <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No Cook's D data") +
      theme_void()
  }

  # --- Panel (c): Outlier sensitivity comparison ---
  outlier <- read.csv(outlier_path, stringsAsFactors = FALSE)
  outlier$Taxa <- factor(outlier$Taxa, levels = taxa_levels)
  # Use actual method names from the CSV
  method_lvls <- sort(unique(outlier$Method))
  outlier$Method <- factor(outlier$Method, levels = method_lvls)

  # Color by method (use col_diag for methodological distinctions)
  available_colors <- unname(col_diag[c("OK", "Non-normal", "Heteroscedastic")])
  if (length(method_lvls) > length(available_colors)) {
    warning("Fig S15: More outlier methods (", length(method_lvls),
            ") than available colors (", length(available_colors), ").")
  }
  method_colors <- setNames(
    available_colors[seq_along(method_lvls)],
    method_lvls
  )

  # Use t-distribution critical value for proper CIs with small sample sizes
  outlier$t_crit <- qt(0.975, df = pmax(outlier$k - 1, 1))

  p_outlier <- ggplot(outlier, aes(x = Estimate, y = Taxa, color = Method)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_errorbar(aes(xmin = Estimate - t_crit * SE, xmax = Estimate + t_crit * SE),
                  width = 0.15, linewidth = 0.3, orientation = "y",
                  position = position_dodge(width = 0.5)) +
    geom_point(size = 1.5, position = position_dodge(width = 0.5)) +
    scale_color_manual(values = method_colors, name = "Outlier method") +
    labs(x = "Effect size estimate (lnRR)", y = NULL,
         subtitle = "Outlier sensitivity") +
    theme_mpa() +
    theme(axis.text.y = element_text(face = "italic"),
          legend.position = "bottom")

  # --- Shared taxa legend ---
  # Build legend from whichever audit subset has data (prefer bio, fall back to full)
  legend_src <- if (nrow(audit_bio) > 0 && "Cooks_Distance" %in% names(audit_bio)) {
    if (!"obs" %in% names(audit_bio)) audit_bio$obs <- seq_len(nrow(audit_bio))
    audit_bio
  } else {
    data.frame(obs = 1, Cooks_Distance = 0, Taxa = names(col_taxa)[1])
  }
  legend_plot <- ggplot(legend_src, aes(x = obs, y = Cooks_Distance, color = Taxa)) +
    geom_point() + scale_color_taxa(name = "Taxa") +
    theme_mpa() +
    theme(legend.position = "bottom",
          legend.text = element_text(face = "italic", size = 7))
  taxa_legend <- cowplot::get_legend(legend_plot)

  # --- Assemble: top row = Cook's panels, bottom = outlier comparison ---
  top_row <- p_cooks_bio + p_cooks_den +
    patchwork::plot_layout(widths = c(1, 1))

  fig_s15 <- top_row / p_outlier +
    patchwork::plot_layout(heights = c(1, 0.8)) +
    patchwork::plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold", size = 10))

  fig_s15_final <- cowplot::plot_grid(
    fig_s15, taxa_legend, ncol = 1, rel_heights = c(1, 0.05)
  )

  save_fig(fig_s15_final, "fig_s15_sensitivity_summary", w = 17, h = 13)
  cat("  Figure S15 saved.\n")

} # end file check
} # end fig_s15

####################################################################################################
## Shared constants for outlier figures (S16–S18) #################################################
####################################################################################################
# These mappings are used by fig_s17 and fig_s18 which work with full species names.
# Defined once here to avoid duplication.

# Abbreviated → full species name mapping
outlier_abbrev_to_full <- c(
  "M. pyrifera"      = "Macrocystis pyrifera",
  "M. franciscanus"  = "Mesocentrotus franciscanus",
  "S. purpuratus"    = "Strongylocentrotus purpuratus",
  "P. interruptus"   = "Panulirus interruptus",
  "S. pulcher"       = "Semicossyphus pulcher"
)

# Response label mapping
outlier_resp_map <- c("Biomass" = "Bio", "Density" = "Den")

# Canonical species order (trophic: predators → urchins → kelp)
outlier_sp_order_full <- c(
  "Panulirus interruptus", "Semicossyphus pulcher",
  "Strongylocentrotus purpuratus", "Mesocentrotus franciscanus",
  "Macrocystis pyrifera"
)

# Full-name → abbreviated for color lookup
outlier_full_to_abbrev <- c(
  "Macrocystis pyrifera"          = "M. pyrifera",
  "Mesocentrotus franciscanus"    = "M. franciscanus",
  "Strongylocentrotus purpuratus" = "S. purpuratus",
  "Panulirus interruptus"         = "P. interruptus",
  "Semicossyphus pulcher"         = "S. pulcher"
)

# Species colors keyed by full name
outlier_sp_colors <- setNames(col_taxa[outlier_full_to_abbrev],
                               names(outlier_full_to_abbrev))

# Facet ordering: "Species — Response" (excluding M. pyrifera Density)
outlier_facet_order <- as.vector(t(outer(
  outlier_sp_order_full, c("Biomass", "Density"),
  function(sp, r) paste0(sp, " \u2014 ", r)
)))
outlier_facet_order <- outlier_facet_order[
  outlier_facet_order != "Macrocystis pyrifera \u2014 Density"
]

####################################################################################################
## Figure S16: Outlier removal rationale ##########################################################
####################################################################################################
# WHY we chose not to remove Cook's D outliers as the primary analysis.
# Four panels: (a) effect sizes by taxon with outlier highlighting, (b) between- vs
# within-taxa distance decomposition, (c) % flagged bar chart, (d) paired forest plot.

if (should_render("fig_s16")) {
cat("\n--- Figure S16: Outlier removal rationale ---\n")

audit_global_path <- here::here("outputs", "filter_audit_meta_analysis.csv")
audit_pertaxa_path <- here::here("outputs", "filter_audit_pertaxa_meta.csv")
sensitivity_path <- here::here("tables", "table_s_outlier_sensitivity.csv")

if (!file.exists(audit_global_path) || !file.exists(audit_pertaxa_path) ||
    !file.exists(sensitivity_path)) {
  warning("Required CSV files not found -- skipping fig_s16.")
} else {

  # --- Load data ---
  audit_g <- read.csv(audit_global_path, stringsAsFactors = FALSE)
  audit_pt <- read.csv(audit_pertaxa_path, stringsAsFactors = FALSE)
  sensitivity <- read.csv(sensitivity_path, stringsAsFactors = FALSE)

  # Standardize factor levels
  audit_g$Taxa <- factor(audit_g$Taxa, levels = taxa_levels)
  audit_pt$Taxa <- factor(audit_pt$Taxa, levels = taxa_levels)
  sensitivity$Taxa <- factor(sensitivity$Taxa, levels = taxa_levels)

  # Compute derived columns
  grand_mean_val <- mean(audit_g$Mean, na.rm = TRUE)
  taxon_resp_means <- aggregate(Mean ~ Taxa + Response, data = audit_g,
                                 FUN = mean, na.rm = TRUE)
  names(taxon_resp_means)[3] <- "taxon_mean"
  audit_g <- merge(audit_g, taxon_resp_means, by = c("Taxa", "Response"))
  audit_g$dist_grand <- abs(audit_g$Mean - grand_mean_val)
  audit_g$dist_taxon <- abs(audit_g$Mean - audit_g$taxon_mean)
  audit_g$outlier_label <- ifelse(audit_g$Is_Outlier, "Flagged by global 4/n", "Retained")

  # Trophic group for annotation
  audit_g$trophic <- trophic_assignment[as.character(audit_g$Taxa)]

  # ---- Panel (a): Effect sizes by taxon with outlier status ----
  # Y range for annotation placement
  y_max_a <- max(audit_g$Mean, na.rm = TRUE) * 1.05
  y_annot <- y_max_a + diff(range(audit_g$Mean, na.rm = TRUE)) * 0.08

  # Build panel (a) base plot, then add trophic annotations as separate grobs
  p_a <- ggplot(audit_g, aes(x = Taxa, y = Mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_hline(yintercept = grand_mean_val, linetype = "dotted",
               color = "grey40", linewidth = 0.3) +
    geom_jitter(aes(color = Taxa,
                    alpha = ifelse(Is_Outlier, 1.0, 0.25),
                    shape = ifelse(Is_Outlier, "Flagged", "Retained")),
                width = 0.2, height = 0, size = 1.8) +
    scale_alpha_identity() +
    scale_shape_manual(values = c("Flagged" = 16, "Retained" = 1),
                       name = "Global 4/n status") +
    scale_color_taxa(name = "Taxa") +
    facet_wrap(~ Response, ncol = 2) +
    labs(x = NULL, y = "Effect size (lnRR)",
         subtitle = expression(
           bold(phantom("Predators") * "\u2191") ~
           bold(phantom("Urchins") * "\u2193") ~
           bold(phantom("Kelp") * "\u2191")
         )) +
    theme_mpa(base_size = 8) +
    theme(axis.text.x = element_text(face = "italic", angle = 35, hjust = 1, size = 6.5),
          legend.position = "none",
          strip.text = element_text(size = 8, face = "bold"),
          plot.subtitle = element_blank(),
          plot.margin = margin(5, 5, 5, 5, "mm"))

  # Add colored trophic annotations using three separate geom_text layers (avoids scale conflict)
  trophic_bio <- data.frame(Response = "Biomass", stringsAsFactors = FALSE)
  p_a <- p_a +
    geom_text(data = cbind(trophic_bio, x = 1.5, y = y_annot),
              inherit.aes = FALSE, aes(x = x, y = y),
              label = "Predators \u2191", color = col_effect["positive"],
              size = 2.3, fontface = "bold", show.legend = FALSE) +
    geom_text(data = cbind(trophic_bio, x = 3.5, y = y_annot),
              inherit.aes = FALSE, aes(x = x, y = y),
              label = "Urchins \u2193", color = col_effect["negative"],
              size = 2.3, fontface = "bold", show.legend = FALSE) +
    geom_text(data = cbind(trophic_bio, x = 5, y = y_annot),
              inherit.aes = FALSE, aes(x = x, y = y),
              label = "Kelp \u2191", color = col_effect["positive"],
              size = 2.3, fontface = "bold", show.legend = FALSE)

  # ---- Panel (b): Between-taxa vs within-taxa distance ----
  p_b <- ggplot(audit_g, aes(x = dist_grand, y = dist_taxon)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(aes(color = Taxa,
                   shape = ifelse(Is_Outlier, "Flagged", "Retained")),
               size = 1.8, alpha = 0.7) +
    scale_color_taxa(name = "Taxa") +
    scale_shape_manual(values = c("Flagged" = 17, "Retained" = 1),
                       name = "Global 4/n status") +
    # Annotation in the lower-right region (below diagonal)
    annotate("text", x = max(audit_g$dist_grand, na.rm = TRUE) * 0.65,
             y = max(audit_g$dist_grand, na.rm = TRUE) * 0.15,
             label = "Flagged due to\ntaxon, not anomaly",
             size = 2.2, color = "grey40", fontface = "italic") +
    # Annotation in the upper-left region (above diagonal)
    annotate("text", x = max(audit_g$dist_grand, na.rm = TRUE) * 0.15,
             y = max(audit_g$dist_grand, na.rm = TRUE) * 0.75,
             label = "True within-taxon\noutliers",
             size = 2.2, color = "grey40", fontface = "italic") +
    labs(x = "|lnRR \u2212 grand mean|",
         y = "|lnRR \u2212 taxon mean|") +
    theme_mpa(base_size = 8) +
    theme(legend.position = "none")

  # ---- Panel (c): % flagged — global vs per-taxon ----
  # Global flagging rate per taxa (combining bio + den)
  flag_global <- aggregate(Is_Outlier ~ Taxa, data = audit_g,
                           FUN = function(x) round(100 * mean(x), 1))
  names(flag_global)[2] <- "pct_flagged"
  flag_global$method <- "Global (4/n)"

  # Per-taxon flagging rate
  flag_pertaxa <- aggregate(Is_Outlier ~ Taxa, data = audit_pt,
                            FUN = function(x) round(100 * mean(x), 1))
  names(flag_pertaxa)[2] <- "pct_flagged"
  flag_pertaxa$method <- "Per-taxon (4/k)"

  # Sample sizes per taxon
  k_per_taxa <- aggregate(Mean ~ Taxa, data = audit_g, FUN = length)
  names(k_per_taxa)[2] <- "k_total"

  flag_df <- rbind(flag_global, flag_pertaxa)
  flag_df$Taxa <- factor(flag_df$Taxa, levels = taxa_levels)
  flag_df$method <- factor(flag_df$method,
                           levels = c("Global (4/n)", "Per-taxon (4/k)"))

  # Merge sample sizes for annotation
  flag_df <- merge(flag_df, k_per_taxa, by = "Taxa")

  p_c <- ggplot(flag_df, aes(x = Taxa, y = pct_flagged, fill = method)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.3) +
    geom_text(aes(label = paste0(pct_flagged, "%")),
              position = position_dodge(width = 0.7),
              vjust = -0.5, size = 2.2) +
    geom_hline(yintercept = 50, linetype = "dotted", color = "grey50", linewidth = 0.3) +
    # k labels below bars
    geom_text(data = flag_df[flag_df$method == "Global (4/n)", ],
              aes(y = -4, label = paste0("k=", k_total)),
              size = 2, color = "grey40") +
    scale_fill_manual(values = c("Global (4/n)" = "#D55E00",
                                 "Per-taxon (4/k)" = "#56B4E9"),
                      name = "Threshold") +
    scale_y_continuous(limits = c(-6, 105), breaks = seq(0, 100, 25)) +
    labs(x = NULL, y = "Observations flagged (%)") +
    theme_mpa(base_size = 8) +
    theme(axis.text.x = element_text(face = "italic", angle = 35, hjust = 1, size = 6.5),
          legend.position = "bottom",
          legend.key.size = unit(3, "mm"),
          legend.text = element_text(size = 7))

  # ---- Panel (d): Meta-analytic estimates — full vs removed ----
  sens_sub <- sensitivity[sensitivity$Method %in%
    c("Joint model, no removal (primary)", "Joint Cook's D (4/n) (legacy)"), ]
  sens_sub$Method_short <- ifelse(
    grepl("no removal", sens_sub$Method),
    "Full data (primary)", "After removal (legacy)"
  )
  sens_sub$Method_short <- factor(sens_sub$Method_short,
    levels = c("Full data (primary)", "After removal (legacy)"))

  # CIs using t-distribution
  sens_sub$t_crit <- qt(0.975, df = pmax(sens_sub$k - 1, 1))
  sens_sub$ci_lo <- sens_sub$Estimate - sens_sub$t_crit * sens_sub$SE
  sens_sub$ci_hi <- sens_sub$Estimate + sens_sub$t_crit * sens_sub$SE

  # Significance label
  sens_sub$sig_label <- ifelse(sens_sub$pval < 0.05, "", "n.s.")
  # k label
  sens_sub$k_label <- paste0("k=", sens_sub$k)

  # Combine k and significance into one label to reduce clutter
  sens_sub$point_label <- paste0(
    "k=", sens_sub$k,
    ifelse(sens_sub$pval >= 0.05, "  n.s.", "")
  )

  p_d <- ggplot(sens_sub, aes(x = Estimate, y = Taxa, color = Method_short)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                   height = 0.15, linewidth = 0.4,
                   position = position_dodge(width = 0.6)) +
    geom_point(size = 2, position = position_dodge(width = 0.6)) +
    geom_text(aes(label = point_label),
              position = position_dodge(width = 0.6),
              size = 1.7, hjust = -0.1, vjust = -0.9, show.legend = FALSE) +
    scale_x_continuous(limits = c(-4, 4), oob = scales::squish) +
    scale_color_manual(values = c("Full data (primary)" = col_effect["positive"],
                                  "After removal (legacy)" = col_effect["negative"]),
                       name = NULL) +
    facet_wrap(~ Response, ncol = 2) +
    labs(x = "Effect size (lnRR)", y = NULL) +
    theme_mpa(base_size = 8) +
    theme(axis.text.y = element_text(face = "italic", size = 7),
          legend.position = "bottom",
          legend.key.size = unit(3, "mm"),
          legend.text = element_text(size = 7),
          strip.text = element_text(size = 8, face = "bold"))

  # ---- Assemble 2x2 ----
  fig_s16 <- (p_a + p_b) / (p_c + p_d) +
    patchwork::plot_layout(heights = c(1, 0.9)) +
    patchwork::plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold", size = 10))

  save_fig(fig_s16, "fig_s16_outlier_rationale",
           w = FIG_S16_DIMS["w"], h = FIG_S16_DIMS["h"])
  cat("  Figure S16 saved.\n")

} # end file check
} # end fig_s16

####################################################################################################
## Figure S17: Temporal trajectories with outlier status ##########################################
####################################################################################################
# Companion to fig_s16. Shows that "outlier" MPA-taxa combinations follow coherent
# recovery trajectories over time — not erratic noise. If these were genuine outliers,
# their time series would be random; instead, they track the cascade prediction.

if (should_render("fig_s17")) {
cat("\n--- Figure S17: Temporal trajectories by outlier status ---\n")

audit_global_path <- here::here("outputs", "filter_audit_meta_analysis.csv")

if (!file.exists(audit_global_path)) {
  warning("Audit CSV not found -- skipping fig_s17.")
} else if (!exists("All.RR.sub.trans")) {
  warning("All.RR.sub.trans not available -- skipping fig_s17.")
} else {

  # --- Load outlier flags (use shared constants) ---
  audit_g <- read.csv(audit_global_path, stringsAsFactors = FALSE)
  audit_g$resp_short <- outlier_resp_map[audit_g$Response]
  audit_g$Taxa_full <- outlier_abbrev_to_full[audit_g$Taxa]

  # Create lookup key: MPA + full taxa name + resp
  audit_g$join_key <- paste(audit_g$MPA, audit_g$Taxa_full, audit_g$resp_short, sep = "|")
  outlier_keys <- audit_g$join_key[audit_g$Is_Outlier]
  retained_keys <- audit_g$join_key[!audit_g$Is_Outlier]

  # --- Detect taxa column name ---
  taxa_col_ts <- if ("y" %in% names(All.RR.sub.trans)) "y" else "Taxa"

  # --- Build temporal dataset (After period, Bio + Den, t <= 15) ---
  ts_data <- All.RR.sub.trans %>%
    dplyr::filter(BA == "After", time >= 0, time <= 15,
                  !is.na(resp), resp %in% c("Bio", "Den")) %>%
    dplyr::mutate(
      time = as.numeric(time),
      join_key = paste(CA_MPA_Name_Short, .data[[taxa_col_ts]], resp, sep = "|"),
      outlier_status = dplyr::case_when(
        join_key %in% outlier_keys ~ "Flagged by global 4/n",
        join_key %in% retained_keys ~ "Retained",
        TRUE ~ "Not in meta-analysis"
      ),
      resp_label = dplyr::case_when(
        resp == "Bio" ~ "Biomass",
        resp == "Den" ~ "Density",
        TRUE ~ resp
      )
    ) %>%
    dplyr::filter(outlier_status != "Not in meta-analysis")

  # Drop M. pyrifera Density (no density data for kelp)
  ts_data <- ts_data %>%
    dplyr::filter(!(.data[[taxa_col_ts]] == "Macrocystis pyrifera" & resp == "Den"))

  # Use shared constants for species ordering and colors
  ts_data[[taxa_col_ts]] <- factor(ts_data[[taxa_col_ts]], levels = outlier_sp_order_full)

  # Combined facet label: "Species — Response"
  ts_data <- ts_data %>%
    dplyr::mutate(
      facet_label = paste0(.data[[taxa_col_ts]], " \u2014 ", resp_label)
    )
  ts_data$facet_label <- factor(ts_data$facet_label, levels = outlier_facet_order)

  # Status factor
  ts_data$outlier_status <- factor(ts_data$outlier_status,
    levels = c("Flagged by global 4/n", "Retained"))

  cat(sprintf("  Temporal data: %d obs, %d flagged trajectories, %d retained\n",
              nrow(ts_data),
              length(unique(ts_data$join_key[ts_data$outlier_status == "Flagged by global 4/n"])),
              length(unique(ts_data$join_key[ts_data$outlier_status == "Retained"]))))

  # --- Build figure ---
  fig_s17 <- ggplot(ts_data, aes(x = time, y = lnDiff)) +
    # Reference line at lnRR = 0 (no MPA effect)
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey30",
               linewidth = 0.4) +
    # Retained MPA trajectories (darker grey, medium weight)
    geom_line(data = ts_data[ts_data$outlier_status == "Retained", ],
              aes(group = interaction(CA_MPA_Name_Short, resp)),
              color = "grey45", alpha = 0.6, linewidth = 0.5) +
    # Flagged MPA trajectories (colored, lighter weight to reduce visual dominance)
    geom_line(data = ts_data[ts_data$outlier_status == "Flagged by global 4/n", ],
              aes(group = interaction(CA_MPA_Name_Short, resp),
                  color = .data[[taxa_col_ts]]),
              alpha = 0.45, linewidth = 0.5) +
    # GAM smooth across ALL data (both flagged + retained)
    geom_smooth(aes(color = .data[[taxa_col_ts]],
                    fill = .data[[taxa_col_ts]]),
                method = "gam", formula = y ~ s(x, k = 5),
                linewidth = 1.2, alpha = 0.15,
                se = TRUE, level = 0.95) +
    facet_wrap(~ facet_label, ncol = 2, scales = "free_y") +
    scale_color_manual(values = outlier_sp_colors, guide = "none") +
    scale_fill_manual(values = outlier_sp_colors, guide = "none") +
    scale_x_continuous(breaks = seq(0, 15, by = 5), limits = c(0, 16),
                       expand = c(0.02, 0)) +
    scale_y_rr_free() +
    labs(
      x = "Years since MPA implementation",
      y = "MPA / Reference"
    ) +
    theme_mpa(base_size = 9) +
    theme(
      strip.text = element_text(face = "italic", size = 8.5,
                                margin = margin(b = 3, t = 3)),
      axis.text  = element_text(size = 7),
      axis.title = element_text(size = 8),
      panel.grid.major = element_blank(),
      plot.margin = margin(6, 6, 6, 6, "pt")
    )

  # Manual legend for line types
  legend_grob <- cowplot::ggdraw() +
    cowplot::draw_line(x = c(0.08, 0.18), y = c(0.5, 0.5),
                       color = "#D55E00", linewidth = 1.0) +
    cowplot::draw_label("Would be removed by global Cook's D (4/n)",
                        x = 0.20, y = 0.5, hjust = 0, size = 7) +
    cowplot::draw_line(x = c(0.70, 0.80), y = c(0.5, 0.5),
                       color = "grey70", linewidth = 0.5) +
    cowplot::draw_label("Retained",
                        x = 0.82, y = 0.5, hjust = 0, size = 7)

  fig_s17_final <- cowplot::plot_grid(
    fig_s17, legend_grob,
    ncol = 1, rel_heights = c(1, 0.03)
  )

  save_fig(fig_s17_final, "fig_s17_temporal_outlier_trajectories",
           w = FIG_S17_DIMS["w"], h = FIG_S17_DIMS["h"])
  cat("  Figure S17 saved.\n")

} # end data checks
} # end fig_s17


# =============================================================================
# Figure S18: Raw data trajectories by outlier status
# =============================================================================
# Shows raw observed values (not response ratios) for MPA and reference sites
# through time, highlighting MPA-taxa combos flagged by global Cook's D.
# Complements fig_s17 (lnRR trajectories) by showing the underlying data
# are well-behaved — flagged combos don't have erratic raw values.
# =============================================================================

if (should_render("fig_s18")) {
cat("\n--- Figure S18: Raw data trajectories by outlier status ---\n")

audit_global_path <- here::here("outputs", "filter_audit_meta_analysis.csv")

if (!file.exists(audit_global_path)) {
  warning("Audit CSV not found -- skipping fig_s18.")
} else if (!exists("All.Resp.sub")) {
  warning("All.Resp.sub not available -- skipping fig_s18.")
} else {

  # --- Load outlier flags ---
  audit_g <- read.csv(audit_global_path, stringsAsFactors = FALSE)

  # Use shared constants for name mapping
  audit_g$resp_short <- outlier_resp_map[audit_g$Response]
  audit_g$Taxa_full <- outlier_abbrev_to_full[audit_g$Taxa]

  # Outlier lookup: keyed by MPA + taxa + resp
  audit_g$join_key <- paste(audit_g$MPA, audit_g$Taxa_full, audit_g$resp_short, sep = "|")
  outlier_keys <- audit_g$join_key[audit_g$Is_Outlier]
  retained_keys <- audit_g$join_key[!audit_g$Is_Outlier]

  # --- Build raw data time series ---
  raw_ts <- All.Resp.sub %>%
    dplyr::filter(!is.na(resp), resp %in% c("Bio", "Den"),
                  !is.na(value), value >= 0) %>%
    dplyr::mutate(
      year = as.numeric(year),
      join_key = paste(CA_MPA_Name_Short, taxon_name, resp, sep = "|"),
      outlier_status = dplyr::case_when(
        join_key %in% outlier_keys ~ "Flagged by global 4/n",
        join_key %in% retained_keys ~ "Retained",
        TRUE ~ "Not in meta-analysis"
      ),
      resp_label = dplyr::case_when(
        resp == "Bio" ~ "Biomass",
        resp == "Den" ~ "Density",
        TRUE ~ resp
      )
    ) %>%
    dplyr::filter(outlier_status != "Not in meta-analysis")

  # Drop M. pyrifera Density (doesn't exist)
  raw_ts <- raw_ts %>%
    dplyr::filter(!(taxon_name == "Macrocystis pyrifera" & resp == "Den"))

  # Use shared constants
  raw_ts$taxon_name <- factor(raw_ts$taxon_name, levels = outlier_sp_order_full)

  # Combined facet label
  raw_ts <- raw_ts %>%
    dplyr::mutate(
      facet_label = paste0(taxon_name, " \u2014 ", resp_label)
    )
  raw_ts$facet_label <- factor(raw_ts$facet_label, levels = outlier_facet_order)

  # Status factor
  raw_ts$outlier_status <- factor(raw_ts$outlier_status,
    levels = c("Flagged by global 4/n", "Retained"))

  # Prettify status labels for legend
  raw_ts$site_type <- ifelse(raw_ts$status == "mpa", "MPA", "Reference")

  cat(sprintf("  Raw data: %d obs, %d flagged MPA-taxa combos, %d retained\n",
              nrow(raw_ts),
              length(unique(raw_ts$join_key[raw_ts$outlier_status == "Flagged by global 4/n"])),
              length(unique(raw_ts$join_key[raw_ts$outlier_status == "Retained"]))))

  # --- Separate MPA and reference into paired lines per MPA ---
  # Average MPA + reference per MPA×year for cleaner spaghetti
  raw_wide <- raw_ts %>%
    dplyr::select(CA_MPA_Name_Short, year, taxon_name, resp, status,
                  value, outlier_status, facet_label) %>%
    tidyr::pivot_wider(names_from = status, values_from = value,
                       values_fn = mean) %>%
    dplyr::filter(!is.na(mpa) | !is.na(reference))

  # --- Build figure ---
  # MPA lines (solid) and reference lines (dashed) per MPA,
  # colored by species for flagged combos, grey for retained.
  # sqrt y-axis to compress extreme values while preserving relative patterns.
  fig_s18 <- ggplot(raw_wide, aes(x = year)) +
    # --- Retained: grey ---
    # MPA (solid)
    geom_line(data = raw_wide[raw_wide$outlier_status == "Retained", ],
              aes(y = mpa, group = CA_MPA_Name_Short),
              color = "grey55", alpha = 0.45, linewidth = 0.35) +
    # Reference (dashed)
    geom_line(data = raw_wide[raw_wide$outlier_status == "Retained", ],
              aes(y = reference, group = CA_MPA_Name_Short),
              color = "grey55", alpha = 0.35, linewidth = 0.3,
              linetype = "dashed") +
    # --- Flagged: colored ---
    # MPA (solid)
    geom_line(data = raw_wide[raw_wide$outlier_status == "Flagged by global 4/n", ],
              aes(y = mpa, group = CA_MPA_Name_Short,
                  color = taxon_name),
              alpha = 0.55, linewidth = 0.5) +
    # Reference (dashed, lighter)
    geom_line(data = raw_wide[raw_wide$outlier_status == "Flagged by global 4/n", ],
              aes(y = reference, group = CA_MPA_Name_Short,
                  color = taxon_name),
              alpha = 0.35, linewidth = 0.4, linetype = "dashed") +
    facet_wrap(~ facet_label, ncol = 2, scales = "free_y") +
    scale_color_manual(values = outlier_sp_colors, guide = "none") +
    scale_y_sqrt(expand = expansion(mult = c(0, 0.05))) +
    labs(
      x = "Year",
      y = "Raw value (sqrt scale)"
    ) +
    theme_mpa(base_size = 9) +
    theme(
      strip.text = element_text(face = "italic", size = 8.5,
                                margin = margin(b = 3, t = 3)),
      axis.text  = element_text(size = 7),
      axis.title = element_text(size = 8),
      panel.grid.major = element_blank(),
      plot.margin = margin(6, 6, 6, 6, "pt"),
      legend.position = "none"
    )

  # Manual legend
  legend_grob <- cowplot::ggdraw() +
    cowplot::draw_line(x = c(0.06, 0.14), y = c(0.5, 0.5),
                       color = "#D55E00", linewidth = 1.0) +
    cowplot::draw_label("Would be removed by global Cook's D (4/n)",
                        x = 0.15, y = 0.5, hjust = 0, size = 7) +
    cowplot::draw_line(x = c(0.52, 0.60), y = c(0.5, 0.5),
                       color = "grey55", linewidth = 0.5) +
    cowplot::draw_label("Retained",
                        x = 0.61, y = 0.5, hjust = 0, size = 7) +
    cowplot::draw_label("(solid = MPA, dashed = reference)",
                        x = 0.77, y = 0.5, hjust = 0, size = 6.5,
                        fontface = "italic", colour = "grey40")

  fig_s18_final <- cowplot::plot_grid(
    fig_s18, legend_grob,
    ncol = 1, rel_heights = c(1, 0.03)
  )

  save_fig(fig_s18_final, "fig_s18_raw_trajectories_outlier_status",
           w = FIG_S18_DIMS["w"], h = FIG_S18_DIMS["h"])
  cat("  Figure S18 saved.\n")

} # end data checks
} # end fig_s18


cat("\n")
cat("========================================================================\n")
if (identical(RENDER_FIGURES, "all")) {
  cat("  ALL MANUSCRIPT FIGURES GENERATED SUCCESSFULLY\n")
  cat("  Main text: Fig 1-4\n")
  cat("  Supplemental: Fig S1, S2, S7a-e (appendix), S8-S12 (diagnostics)\n")
} else {
  cat("  SELECTED FIGURES GENERATED:", paste(RENDER_FIGURES, collapse = ", "), "\n")
}
cat("========================================================================\n")
cat("\n=== Figures saved to:", here::here("plots"), "===\n")
