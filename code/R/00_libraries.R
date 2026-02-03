# =============================================================================
# 00_libraries.R
# =============================================================================
#
# PURPOSE:
#   Load all R packages required for the MPA kelp forest pBACIPS analysis.
#   This script should be sourced first before any other analysis scripts.
#
# WHAT THIS SCRIPT DOES:
#   1. Loads core data manipulation packages (tidyverse, dplyr, data.table)
#   2. Loads statistical modeling packages (lme4, metafor, mgcv)
#   3. Loads plotting packages (ggplot2, cowplot, ggpubr)
#   4. Loads nonlinear curve fitting packages (minpack.lm, nls2)
#   5. Sets up a basic color palette
#
# DEPENDENCIES:
#   All packages must be installed before running this script.
#   Install missing packages with: install.packages("package_name")
#
# AUTHORS: Emily Donham & Adrian Stier
# PROJECT: CA MPA Kelp Forest pBACIPS Analysis
# =============================================================================

# -----------------------------------------------------------------------------
# IMPORTANT: Load plyr BEFORE dplyr/tidyverse
# -----------------------------------------------------------------------------
# The plyr package has functions (like summarise) that conflict with dplyr.
# Loading plyr first ensures dplyr's versions take precedence, which is what
# we want for this analysis. You'll see masking warnings - this is expected.

library(plyr)

# -----------------------------------------------------------------------------
# Core Data Manipulation
# -----------------------------------------------------------------------------
# tidyverse: A collection of packages for data science (includes ggplot2, dplyr,
#            tidyr, readr, purrr, tibble, stringr, forcats)
# data.table: Fast data manipulation for large datasets
# here: Simplifies file path management - always use here::here() for paths

library(tidyverse)
library(data.table)
library(here)

# Load dplyr again to ensure its functions mask plyr's
library(dplyr)

# -----------------------------------------------------------------------------
# Statistical Modeling
# -----------------------------------------------------------------------------
# lme4: Linear mixed-effects models (lmer, glmer)
# lmerTest: Adds p-values to lme4 model summaries
# car: Companion to Applied Regression - provides Anova() and other functions
# MuMIn: Multi-Model Inference - for AICc calculations and model averaging
# mgcv: Generalized Additive Models (GAMs)
# emmeans: Estimated Marginal Means - for post-hoc comparisons
# stats: Base R statistics (loaded by default, but explicit is good practice)
# lmtest: Diagnostic tests for linear models

library(lme4)
library(lmerTest)
library(car)
library(MuMIn)
library(mgcv)
library(emmeans)
library(stats)
library(lmtest)

# -----------------------------------------------------------------------------
# Meta-Analysis
# -----------------------------------------------------------------------------
# metafor: The primary package for meta-analysis in R
#          Provides rma(), rma.mv() for multilevel meta-analysis
#          Used in 09_meta_analysis.R

library(metafor)

# -----------------------------------------------------------------------------
# Nonlinear Curve Fitting (for pBACIPS models)
# -----------------------------------------------------------------------------
# minpack.lm: Levenberg-Marquardt algorithm for nonlinear least squares
#             More robust than base R nls() for difficult optimization
# nls2: Enhanced nls with brute-force grid search for starting values
# AICcmodavg: AICc calculation and model comparison
# nlstools: Tools for nonlinear regression (confidence intervals, diagnostics)
# investr: Inverse estimation and calibration (predFit for prediction intervals)

require(minpack.lm)
require(nls2)
require(AICcmodavg)
library(nlstools)
library(investr)

# -----------------------------------------------------------------------------
# Plotting and Visualization
# -----------------------------------------------------------------------------
# cowplot: Publication-quality ggplot2 themes and plot arrangement
# ggpubr: ggplot2-based publication ready plots
# wesanderson: Color palettes inspired by Wes Anderson films
# yarrr: Pirate plots and color utilities

library(cowplot)
library(ggpubr)
library(wesanderson)
library(yarrr)

# -----------------------------------------------------------------------------
# Tables and Reporting
# -----------------------------------------------------------------------------
# broom: Tidy model outputs into data frames
# gt: Grammar of Tables - create publication-quality tables
# gtsummary: Summary tables for regression models
# flextable: Flexible tables for Word/PowerPoint export
# formattable: Conditional formatting for tables

library(broom)
library(gt)
library(gtsummary)
library(flextable)
library(formattable)

# -----------------------------------------------------------------------------
# Miscellaneous Utilities
# -----------------------------------------------------------------------------
# Rmisc: Miscellaneous functions including summarySE() for summary statistics
# olsrr: OLS regression diagnostics
# pacman: Package management (p_load for install + load)
# sf: Simple Features for spatial data (used for site coordinates)
# rempsyc: Psychology-focused reporting utilities
# remef: Remove random effects from mixed models (optional)

library(Rmisc)
library(olsrr)
library(pacman)
library(sf)
library(rempsyc)

# Conditionally load remef if available (not critical for core analysis)
if (requireNamespace("remef", quietly = TRUE)) library(remef)

# -----------------------------------------------------------------------------
# Vegan for Community Ecology (optional)
# -----------------------------------------------------------------------------
# vegan: Community ecology package - provides diversity indices, ordination, etc.
#        Used for multivariate analysis if needed

library(vegan)

# -----------------------------------------------------------------------------
# Global Color Palette (basic version)
# -----------------------------------------------------------------------------
# This creates a continuous color palette from the Wes Anderson "Zissou1" theme.
# A more comprehensive color system is defined in 00b_color_palette.R

pal <- wes_palette("Zissou1", 25, type = "continuous")

# -----------------------------------------------------------------------------
# Confirmation message
# -----------------------------------------------------------------------------
cat("Libraries loaded successfully.\n")
cat("Next: source 00b_color_palette.R for the unified color system.\n")
