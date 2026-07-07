# =============================================================================
# resilience_modules.R
# =============================================================================
#
# PURPOSE:
#   Single source of truth for which scripts make up the resilience module
#   (scripts 14-23). Both orchestrators source this file instead of hard-coding
#   the list, so adding/removing a resilience analysis is a one-line edit here.
#
#   - run_all.R sources RESILIENCE_MODULES_IN_PIPELINE (the subset that is part of
#     the main pBACIPS pipeline; external-data scripts skip gracefully if their
#     inputs are absent).
#   - run_resilience.R sources RESILIENCE_MODULES_FULL_SUITE (all ten, grouped by
#     resilience facet; see docs/RESILIENCE_SYNTHESIS.md).
#
#   Data needs (documented in run_resilience.R): 14/19/21/22/23 use only the
#   tracked harmonized CSVs; 15/16/18 also read data/sumstats_final.csv + the
#   Kumagai mirror; 17/20 read raw PISCO from the sibling data repo.
#
# AUTHORS: Emily Donham & Adrian Stier
# =============================================================================

# Subset run inside the main pipeline (run_all.R), in pipeline order.
RESILIENCE_MODULES_IN_PIPELINE <- c(
  "14_heatwave_analysis.R",        # core response: cascade x 2014-16 MHW (AR1)
  "15_methods_comparison.R",       # robustness: analytical multiverse vs Kumagai
  "16_environmental_moderators.R", # robustness: Kumagai PCA covariates as moderators
  "19_heatwave_replication.R",     # repeatability: resilience across TWO heatwaves
  "21_resilience_stability.R",     # stability: temporal CV inside vs outside
  "22_resistance_recovery.R",      # resistance/recovery on the state variable
  "23_ecological_memory.R"         # memory: do the same reserves repeat?
)

# The full suite (run_resilience.R), grouped by resilience facet.
RESILIENCE_MODULES_FULL_SUITE <- c(
  "14_heatwave_analysis.R",            # CORE RESPONSE
  "19_heatwave_replication.R",         # REPEATABILITY
  "20_compound_disturbance.R",         # GENERALITY ACROSS STRESSORS (SSWD)
  "22_resistance_recovery.R",          # RESISTANCE & RECOVERY
  "23_ecological_memory.R",            # MEMORY / reserve consistency
  "21_resilience_stability.R",         # STABILITY (temporal CV)
  "15_methods_comparison.R",           # ROBUSTNESS: method vs data (multiverse)
  "16_environmental_moderators.R",     # ROBUSTNESS: environmental moderators
  "18_mpa_effectiveness_predictors.R", # ROBUSTNESS: predictability of effectiveness
  "17_eisaguirre_reproduction.R"       # CROSS-STUDY REPRODUCTION
)
