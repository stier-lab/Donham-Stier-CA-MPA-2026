#' ---
#' title: "Resilience analysis suite"
#' description: "Runs the full MPA-resilience module (scripts 14-23) in one call"
#' author: "Emily Donham & Adrian Stier"
#' ---
#'
#' WHAT THIS RUNS
#'   The complete set of analyses asking whether, how, and how robustly California
#'   MPAs confer resilience to climate/biotic disturbance, and how our approach
#'   compares to Kumagai (2024) and Eisaguirre (2020). Organized by resilience facet:
#'
#'   CORE RESPONSE (does the MPA cascade respond to disturbance?)
#'     14_heatwave_analysis.R        - MPA cascade x the 2014-16 marine heatwave (MHW)
#'   REPEATABILITY (does the resilience recur?)
#'     19_heatwave_replication.R     - same response across TWO heatwaves (2014-16 & 2018-20)
#'   GENERALITY ACROSS STRESSOR TYPES (does it hold for non-thermal disturbance?)
#'     20_compound_disturbance.R     - sea-star wasting (SSWD): patchy sunflower star;
#'                                     the cascade response is heat+protection, not SSWD
#'   RESISTANCE & RECOVERY (state-variable decomposition, Kumagai-comparable)
#'     22_resistance_recovery.R      - kelp grew inside / declined outside; recent urchin suppression inside
#'   MEMORY (do the same reserves repeat?)
#'     23_ecological_memory.R        - reserve MHW1 response predicts its MHW2 response (consistency)
#'   STABILITY (does protection damp variability?)
#'     21_resilience_stability.R     - temporal CV inside vs outside (a distinct facet)
#'   ROBUSTNESS / ATTRIBUTION (how much do method, data, covariates change the story?)
#'     15_methods_comparison.R       - analytical multiverse vs Kumagai (method vs data)
#'     16_environmental_moderators.R - env covariates Kumagai only mapped (null after FDR)
#'     18_mpa_effectiveness_predictors.R - can we predict which MPAs are effective? (no)
#'   CROSS-STUDY REPRODUCTION
#'     17_eisaguirre_reproduction.R  - rebuild Eisaguirre 2020 from raw PISCO
#'
#'   Synthesis of all of the above: docs/RESILIENCE_SYNTHESIS.md
#'
#' DATA NEEDS
#'   14, 19, 21 use only the tracked harmonized CSVs (always run; also in run_all.R).
#'   16, 18 additionally read data/sumstats_final.csv (tracked) and the Kumagai
#'     mirror (~/kumagai2024-comparison) for cold-spell/cross-substrate parts.
#'   17, 20 read raw PISCO from the sibling data repo (~/Donham-Stier-CA-MPA-Data-2026).
#'   Scripts skip gracefully (with a message) if their external inputs are absent.
#'
#' USAGE
#'   source(here::here("code", "R", "run_resilience.R"))

rm(list = ls())
modules <- c("14_heatwave_analysis.R", "19_heatwave_replication.R", "20_compound_disturbance.R",
             "22_resistance_recovery.R", "23_ecological_memory.R", "21_resilience_stability.R",
             "15_methods_comparison.R", "16_environmental_moderators.R",
             "18_mpa_effectiveness_predictors.R", "17_eisaguirre_reproduction.R")
cat("==============================================================\n")
cat("  MPA RESILIENCE ANALYSIS SUITE\n  ", length(modules), "modules\n")
cat("==============================================================\n\n")
for (m in modules) {
  cat("---- Running", m, "----\n")
  t0 <- Sys.time()
  ok <- tryCatch({ source(here::here("code", "R", m), local = new.env()); TRUE },
                 error = function(e) { cat("  ****", m, "FAILED:", conditionMessage(e), "\n"); FALSE })
  if (ok) cat("   ", m, "done (", round(difftime(Sys.time(), t0, units = "secs"), 1), "s)\n\n")
}
cat("==============================================================\n")
cat("  Resilience suite complete. See docs/RESILIENCE_SYNTHESIS.md\n")
cat("==============================================================\n")
