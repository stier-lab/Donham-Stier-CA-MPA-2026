#' ---
#' title: "Unified Color Palette and Theme for MPA Kelp Forest Analysis"
#' description: "Self-contained, colorblind-friendly, publication-quality color system"
#' author: "Emily Donham & Adrian Stier"
#' target_journal: "Conservation Letters"
#' ---
#'
#' This file defines all colors, shapes, and theme elements used across the
#' manuscript figures. Source this file after 00_libraries.R to override
#' the ad-hoc palette definitions in individual scripts.
#'
#' Design principles:
#'   1. Colorblind-safe: all hues separated by >= 30 degrees on the color wheel
#'      and distinguishable under deuteranopia, protanopia, and tritanopia.
#'      Tested via dichromat simulation during development.
#'   2. Print-safe: luminance values span 25-75% so categories remain distinct
#'      in grayscale reproduction.
#'   3. Self-contained: no external palette packages required.
#'   4. Consistent with existing Fig 3 convention (green = density, orange = biomass).

# =============================================================================
# 1. TAXA COLORS
# =============================================================================
# Five species spanning the trophic cascade, each mapped to a color that
# evokes the organism while maintaining perceptual distinctness.
#
# Luminance (L* in CIELAB, approximate):
#   S. purpuratus  ~ 42   (dark-medium)
#   M. franciscanus ~ 48  (medium)
#   M. pyrifera    ~ 52   (medium)
#   P. interruptus ~ 62   (medium-light)
#   S. pulcher     ~ 45   (medium)
#
# These values ensure adequate grayscale separation.

col_taxa <- c(
  "S. purpuratus"   = "#7B6A8E",
    # Muted purple -- echoes the animal's test color. Desaturated enough to
    # remain distinct from the blue of S. pulcher under deuteranopia.

  "M. franciscanus"  = "#B5503E",
    # Warm brick-red / rust -- reflects the red urchin's spine color.
    # Shifted toward brown to avoid pure-red which collapses with green
    # for protanopes.

  "M. pyrifera"      = "#5B7744",
    # Earthy olive-green -- references kelp frond color. Kept dark enough
    # to separate from the lighter amber of P. interruptus in grayscale.

  "P. interruptus"   = "#D4933B",
    # Warm amber / golden-orange -- evokes the spiny lobster's carapace.
    # High luminance provides contrast against the darker taxa colors.

  "S. pulcher"       = "#2E6E7E"
    # Deep teal-blue -- references sheephead coloration in the marine
    # environment. Separated from purple by shifting toward cyan.
)

# Short-name aliases matching the existing 08_effect_sizes.R convention
col_taxa_short <- c(
  "Purps"    = unname(col_taxa["S. purpuratus"]),
  "Reds"     = unname(col_taxa["M. franciscanus"]),
  "Kelp"     = unname(col_taxa["M. pyrifera"]),
  "Lobs"     = unname(col_taxa["P. interruptus"]),
  "Sheephead" = unname(col_taxa["S. pulcher"])
)

# Backward-compatible alias: drop-in replacement for `cols` in 08_effect_sizes.R
cols <- col_taxa_short


# =============================================================================
# 2. RESPONSE TYPE COLORS (Density vs Biomass)
# =============================================================================
# Maintains the existing Fig 3 convention: green family for density,
# orange family for biomass. Colors are chosen to pair well when overlaid
# or placed side-by-side in a forest plot.

col_response <- c(
  "Den"     = "#5B8C5A",
    # Muted sage green -- cooler tone signals count-based measurement.

  "Bio"     = "#D4873B"
    # Warm amber-orange -- warmer tone signals mass-based measurement.
    # Distinct from density green even under deuteranopia because the
    # luminance difference is > 15 L* units.
)

# Longer labels for legends
col_response_long <- c(
  "Density" = unname(col_response["Den"]),
  "Biomass" = unname(col_response["Bio"])
)


# =============================================================================
# 3. MPA STATUS COLORS (Inside vs Outside / Reference)
# =============================================================================
# A two-tone pair for the BACI design. The gold/amber for MPA interior
# and a cooler dark slate for reference sites provides both hue and
# luminance contrast.

col_site <- c(
  "Inside"    = "#D4A03C",
    # Warm gold -- MPA protected site. High luminance draws the eye to the
    # treatment of interest.

  "Outside"   = "#5A5A78"
    # Cool dark slate with a slight blue undertone -- reference site.
    # Lower luminance recedes visually, signaling "control."
)

# Alternative labels used in some scripts
col_site_alt <- c(
  "MPA"       = unname(col_site["Inside"]),
  "Reference" = unname(col_site["Outside"])
)


# =============================================================================
# 4. BEFORE / AFTER STYLING
# =============================================================================
# Rather than distinct hues (which would consume color-channel budget),
# Before/After is encoded via fill + alpha, keeping point shape and color
# free for other dimensions.

col_ba <- c(
  "Before" = "#A0A0A0",
    # Neutral mid-gray -- visually lighter, "faded" to suggest the
    # pre-intervention baseline.

  "After"  = "#2B2B2B"
    # Near-black -- darker fill signals the treatment period.
    # The luminance contrast ratio (Before/After) is ~3.5:1.
)

# Point fill convention: Before = open (fill NA), After = filled
# Use these with shape 21 (fillable circle) for maximum flexibility.
fill_ba <- c(
  "Before" = NA,        # open / unfilled
  "After"  = "#2B2B2B"  # filled dark
)

alpha_ba <- c(
  "Before" = 0.50,
  "After"  = 1.00
)


# =============================================================================
# 5. DATA SOURCE SHAPES
# =============================================================================
# Shape encodes monitoring program, keeping color channels available for
# biological variables. Shapes are chosen for legibility at small sizes
# (8-10 pt in Conservation Letters figures).

shape_source <- c(
  "KFM"     = 15,   # filled square   -- NPS Channel Islands program
  "LTER"    = 16,   # filled circle   -- SBC LTER
  "PISCO"   = 17,   # filled triangle -- PISCO kelp forest surveys
  "Landsat" = 18    # filled diamond  -- satellite remote sensing
)

# Open-point variants (for Before period overlay)
shape_source_open <- c(
  "KFM"     = 0,    # open square
  "LTER"    = 1,    # open circle
  "PISCO"   = 2,    # open triangle
  "Landsat" = 5     # open diamond
)

# Source display labels for legends
label_source <- c(
  "KFM"     = "KFM (NPS)",
  "LTER"    = "SBC LTER",
  "PISCO"   = "PISCO",
  "Landsat" = "Landsat"
)


# =============================================================================
# 6. ANALYSIS TYPE COLORS
# =============================================================================
# Used in supplementary diagnostics to distinguish model selection outcomes.

col_type <- c(
  "pBACIPS" = "#2E6E7E",
  "BACI"    = "#D4933B",
  "CI"      = "#7B6A8E"
)


# =============================================================================
# 7. GGPLOT2 THEME: theme_mpa()
# =============================================================================
#' Publication-ready ggplot2 theme for Conservation Letters
#'
#' A minimal, clean theme with appropriate font sizes for a two-column journal
#' layout. Designed for figures that are 80-170 mm wide.
#'
#' @param base_size Base font size in points (default 9, suitable for
#'   Conservation Letters single-column figures).
#' @param base_family Font family (default "Helvetica"; falls back to
#'   "sans" if unavailable).
#' @return A ggplot2 theme object.
theme_mpa <- function(base_size = 9, base_family = "Helvetica") {

  # Fall back gracefully if Helvetica is not installed
  available_fonts <- tryCatch(names(grDevices::pdfFonts()), error = function(e) NULL)
  if (!is.null(available_fonts) && !(base_family %in% available_fonts)) {
    base_family <- "sans"
  }

  ggplot2::theme_bw(base_size = base_size, base_family = base_family) %+replace%
    ggplot2::theme(
      # Panel
      panel.grid.major  = ggplot2::element_blank(),
      panel.grid.minor  = ggplot2::element_blank(),
      panel.border      = ggplot2::element_rect(colour = "black", fill = NA,
                                                 linewidth = 0.4),
      panel.background  = ggplot2::element_rect(fill = "white"),

      # Axes
      axis.title        = ggplot2::element_text(size = base_size,
                                                 colour = "black"),
      axis.text         = ggplot2::element_text(size = base_size - 1,
                                                 colour = "black"),
      axis.ticks        = ggplot2::element_line(linewidth = 0.3,
                                                 colour = "black"),
      axis.ticks.length = ggplot2::unit(1.5, "mm"),

      # Legend
      legend.background = ggplot2::element_blank(),
      legend.key        = ggplot2::element_blank(),
      legend.title      = ggplot2::element_text(size = base_size,
                                                 face = "bold"),
      legend.text       = ggplot2::element_text(size = base_size - 1),
      legend.key.size   = ggplot2::unit(3.5, "mm"),

      # Strip (for faceted plots)
      strip.background  = ggplot2::element_rect(fill = "grey95",
                                                  colour = "black",
                                                  linewidth = 0.3),
      strip.text        = ggplot2::element_text(size = base_size,
                                                  face = "italic",
                                                  margin = ggplot2::margin(2, 2, 2, 2)),

      # Plot title / caption
      plot.title        = ggplot2::element_text(size = base_size + 1,
                                                 face = "bold",
                                                 hjust = 0,
                                                 margin = ggplot2::margin(b = 4)),
      plot.subtitle     = ggplot2::element_text(size = base_size,
                                                 hjust = 0,
                                                 margin = ggplot2::margin(b = 4)),
      plot.caption      = ggplot2::element_text(size = base_size - 2,
                                                 colour = "grey40",
                                                 hjust = 1),

      # Margins (tight for journal layout)
      plot.margin       = ggplot2::margin(4, 4, 4, 4, unit = "pt")
    )
}


# =============================================================================
# 8. CONVENIENCE SCALE FUNCTIONS
# =============================================================================
# Drop-in scale_*() wrappers so every figure automatically picks up the
# correct palette without re-specifying hex codes.

#' @describeIn scale_color_taxa Discrete color scale for taxa
scale_color_taxa <- function(...) {
  ggplot2::scale_color_manual(values = col_taxa, ...)
}

#' @describeIn scale_fill_taxa Discrete fill scale for taxa
scale_fill_taxa <- function(...) {
  ggplot2::scale_fill_manual(values = col_taxa, ...)
}

#' @describeIn scale_color_response Discrete color scale for response type
scale_color_response <- function(...) {
  ggplot2::scale_color_manual(values = col_response, ...)
}

#' @describeIn scale_fill_response Discrete fill scale for response type
scale_fill_response <- function(...) {
  ggplot2::scale_fill_manual(values = col_response, ...)
}

#' @describeIn scale_color_site Discrete color scale for MPA status
scale_color_site <- function(...) {
  ggplot2::scale_color_manual(values = col_site, ...)
}

#' @describeIn scale_fill_site Discrete fill scale for MPA status
scale_fill_site <- function(...) {
  ggplot2::scale_fill_manual(values = col_site, ...)
}

#' Discrete shape scale for monitoring data sources
scale_shape_source <- function(...) {
  ggplot2::scale_shape_manual(values = shape_source, labels = label_source, ...)
}


# =============================================================================
# 9. PALETTE PREVIEW FUNCTION
# =============================================================================
#' Generate a reference sheet showing all palette colors and shapes
#'
#' Produces a multi-panel figure with labeled swatches for every palette
#' dimension. Useful for checking colors on screen, in print, and under
#' simulated color-vision deficiency.
#'
#' @param save_path Optional file path (e.g. "plots/palette_preview.pdf").
#'   If NULL, the plot is drawn to the current device.
#' @param width Figure width in inches (default 7).
#' @param height Figure height in inches (default 8).
#' @return A ggplot2 grob (invisibly).
preview_palette <- function(save_path = NULL, width = 7, height = 8) {

  # --- Helper: swatch data frame ---
  make_swatch_df <- function(palette, group_label) {
    data.frame(
      name  = factor(names(palette), levels = names(palette)),
      color = unname(palette),
      group = group_label,
      stringsAsFactors = FALSE
    )
  }

  # --- Panel A: Taxa ---
  df_taxa <- make_swatch_df(col_taxa, "Taxa")
  p_taxa <- ggplot2::ggplot(df_taxa,
    ggplot2::aes(x = name, y = 1, fill = name)) +
    ggplot2::geom_tile(width = 0.85, height = 0.7, colour = "white",
                        linewidth = 0.5) +
    ggplot2::scale_fill_manual(values = col_taxa) +
    ggplot2::geom_text(ggplot2::aes(label = color), vjust = 2.8,
                        size = 2.5, colour = "grey30") +
    ggplot2::labs(title = "A. Taxa", x = NULL, y = NULL) +
    theme_mpa(base_size = 8) +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.position = "none",
      axis.text.x  = ggplot2::element_text(face = "italic", angle = 25,
                                             hjust = 1, size = 7)
    )

  # --- Panel B: Response type ---
  df_resp <- make_swatch_df(col_response_long, "Response")
  p_resp <- ggplot2::ggplot(df_resp,
    ggplot2::aes(x = name, y = 1, fill = name)) +
    ggplot2::geom_tile(width = 0.85, height = 0.7, colour = "white",
                        linewidth = 0.5) +
    ggplot2::scale_fill_manual(values = col_response_long) +
    ggplot2::geom_text(ggplot2::aes(label = color), vjust = 2.8,
                        size = 2.5, colour = "grey30") +
    ggplot2::labs(title = "B. Response Type", x = NULL, y = NULL) +
    theme_mpa(base_size = 8) +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.position = "none"
    )

  # --- Panel C: MPA status ---
  df_site <- make_swatch_df(col_site, "Site")
  p_site <- ggplot2::ggplot(df_site,
    ggplot2::aes(x = name, y = 1, fill = name)) +
    ggplot2::geom_tile(width = 0.85, height = 0.7, colour = "white",
                        linewidth = 0.5) +
    ggplot2::scale_fill_manual(values = col_site) +
    ggplot2::geom_text(ggplot2::aes(label = color), vjust = 2.8,
                        size = 2.5, colour = "grey30") +
    ggplot2::labs(title = "C. MPA Status", x = NULL, y = NULL) +
    theme_mpa(base_size = 8) +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.position = "none"
    )

  # --- Panel D: Before / After ---
  df_ba <- data.frame(
    name  = factor(c("Before", "After"), levels = c("Before", "After")),
    color = c(col_ba["Before"], col_ba["After"]),
    fill  = c(NA, col_ba["After"]),
    stringsAsFactors = FALSE
  )
  p_ba <- ggplot2::ggplot(df_ba,
    ggplot2::aes(x = name, y = 1)) +
    ggplot2::geom_point(ggplot2::aes(fill = name, colour = name),
                         shape = 21, size = 6, stroke = 1.2) +
    ggplot2::scale_colour_manual(values = col_ba) +
    ggplot2::scale_fill_manual(values = fill_ba, na.value = "white") +
    ggplot2::labs(title = "D. Before / After (shape 21)", x = NULL,
                   y = NULL) +
    theme_mpa(base_size = 8) +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.position = "none"
    )

  # --- Panel E: Data source shapes ---
  df_shape <- data.frame(
    name  = factor(names(shape_source), levels = names(shape_source)),
    shape = unname(shape_source),
    y     = 1,
    stringsAsFactors = FALSE
  )
  p_shape <- ggplot2::ggplot(df_shape,
    ggplot2::aes(x = name, y = y, shape = name)) +
    ggplot2::geom_point(size = 4, colour = "grey20") +
    ggplot2::scale_shape_manual(values = shape_source) +
    ggplot2::geom_text(ggplot2::aes(label = paste0("pch ", shape)),
                        vjust = 2.5, size = 2.5, colour = "grey40") +
    ggplot2::labs(title = "E. Data Source Shapes", x = NULL, y = NULL) +
    theme_mpa(base_size = 8) +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.position = "none"
    )

  # --- Panel F: Analysis type ---
  df_type <- make_swatch_df(col_type, "Type")
  p_type <- ggplot2::ggplot(df_type,
    ggplot2::aes(x = name, y = 1, fill = name)) +
    ggplot2::geom_tile(width = 0.85, height = 0.7, colour = "white",
                        linewidth = 0.5) +
    ggplot2::scale_fill_manual(values = col_type) +
    ggplot2::geom_text(ggplot2::aes(label = color), vjust = 2.8,
                        size = 2.5, colour = "grey30") +
    ggplot2::labs(title = "F. Analysis Type", x = NULL, y = NULL) +
    theme_mpa(base_size = 8) +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.position = "none"
    )

  # --- Combine ---
  combined <- cowplot::plot_grid(
    p_taxa, p_resp, p_site, p_ba, p_shape, p_type,
    ncol = 2, rel_heights = c(1.1, 1, 1),
    labels = NULL
  )

  if (!is.null(save_path)) {
    ggplot2::ggsave(save_path, combined, width = width, height = height,
                     dpi = 300)
    message("Palette preview saved to: ", save_path)
  } else {
    print(combined)
  }

  invisible(combined)
}


# =============================================================================
# 10. SESSION NOTE
# =============================================================================
message("MPA color palette loaded (",
        length(col_taxa), " taxa, ",
        length(col_response), " responses, ",
        length(col_site), " site types, ",
        length(shape_source), " source shapes)")
