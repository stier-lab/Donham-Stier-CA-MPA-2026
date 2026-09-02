# =============================================================================
# 24_resilience_pipeline_check.R
# =============================================================================
#
# PURPOSE:
#   Fail fast if the resilience analyses used by the manuscript and Supporting
#   Information drift out of the main analysis pipeline. This is a concordance
#   gate, not a new analysis: it checks that run_all.R generated the expected
#   resilience tables/figures and that the headline Figure 5 values still match
#   the manuscript text.
#
# AUTHORS: Emily Donham & Adrian Stier
# =============================================================================

cat("Checking resilience pipeline integration (script 24)...\n")

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required for path checks.", call. = FALSE)
}

repo_file <- function(...) here::here(...)

fail <- function(...) stop(paste0(...), call. = FALSE)

check_files <- function(paths, label) {
  full <- vapply(paths, function(p) repo_file(p), character(1))
  missing <- paths[!file.exists(full)]
  empty <- paths[file.exists(full) & file.info(full)$size <= 0]
  if (length(missing) > 0 || length(empty) > 0) {
    bits <- c(
      if (length(missing) > 0) paste0("missing: ", paste(missing, collapse = ", ")),
      if (length(empty) > 0) paste0("empty: ", paste(empty, collapse = ", "))
    )
    fail(label, " artifact check failed (", paste(bits, collapse = "; "), ").")
  }
  invisible(TRUE)
}

read_required_csv <- function(path) {
  full <- repo_file(path)
  if (!file.exists(full)) fail("Required table missing: ", path)
  read.csv(full, stringsAsFactors = FALSE, check.names = FALSE)
}

expect_one <- function(df, predicate, label) {
  ix <- which(predicate)
  if (length(ix) != 1) fail("Expected one row for ", label, "; found ", length(ix), ".")
  df[ix, , drop = FALSE]
}

expect_equal <- function(actual, expected, label) {
  if (!identical(actual, expected)) {
    fail(label, " mismatch. Expected ", paste(expected, collapse = ", "),
         "; got ", paste(actual, collapse = ", "), ".")
  }
}

expect_true <- function(ok, label) {
  if (!isTRUE(ok)) fail(label)
}

# The manuscript/SI resilience products that are expected from run_all.R.
expected_tables <- c(
  "tables/table_heatwave_period_effects.csv",
  "tables/table_heatwave_contrasts.csv",
  "tables/table_heatwave_diagnostics.csv",
  "tables/table_heatwave_cascade_regression.csv",
  "tables/table_heatwave_per_mpa_exposure.csv",
  "tables/table_heatwave_pbacips_integrated.csv",
  "tables/table_heatwave_sensitivity.csv",
  "tables/table_s_methods_crosswalk.csv",
  "tables/table_s_methods_multiverse.csv",
  "tables/table_s_methods_decomposition.csv",
  "tables/table_s_mpa_env_covariates.csv",
  "tables/table_s_env_moderators.csv",
  "tables/table_heatwave_replication.csv",
  "tables/table_heatwave_replication_continuous.csv",
  "tables/table_resilience_stability.csv",
  "tables/table_resistance_recovery.csv",
  "tables/table_resistance_recovery_sensitivity.csv",
  "tables/table_resistance_recovery_pairs.csv",
  "tables/table_resistance_recovery_pair_counts.csv",
  "tables/table_ecological_memory.csv"
)

expected_figures <- c(
  "plots/fig_heatwave_cascade.pdf",
  "plots/fig_heatwave_cascade.png",
  "plots/fig_heatwave_cascade_regression.pdf",
  "plots/fig_heatwave_cascade_regression.png",
  "plots/fig_s_methods_multiverse.pdf",
  "plots/fig_s_methods_multiverse.png",
  "plots/fig_s_env_moderators.pdf",
  "plots/fig_s_env_moderators.png",
  "plots/fig_heatwave_replication.pdf",
  "plots/fig_heatwave_replication.png",
  "plots/fig_resilience_stability.pdf",
  "plots/fig_resilience_stability.png",
  "plots/fig_resistance_recovery.pdf",
  "plots/fig_resistance_recovery.png",
  "plots/fig_kelp_resilience.pdf",
  "plots/fig_kelp_resilience.png",
  "plots/fig_kelp_resilience_paired.pdf",
  "plots/fig_kelp_resilience_paired.png",
  "plots/fig_ecological_memory.pdf",
  "plots/fig_ecological_memory.png"
)

expected_docs <- c(
  "docs/RESILIENCE_SYNTHESIS.md",
  "docs/RESILIENCE_PIPELINE_INTEGRATION.md",
  "docs/kelp_resilience_figure_text.md",
  "docs/methods_comparison_supplement.md",
  "docs/heatwave_replication.md"
)

check_files(expected_tables, "Resilience table")
check_files(expected_figures, "Resilience figure")
check_files(expected_docs, "Resilience documentation")

# run_all.R must source the same in-pipeline resilience set the manuscript relies on.
source(repo_file("code", "R", "resilience_modules.R"))
expected_modules <- c(
  "14_heatwave_analysis.R",
  "15_methods_comparison.R",
  "16_environmental_moderators.R",
  "19_heatwave_replication.R",
  "21_resilience_stability.R",
  "22_resistance_recovery.R",
  "23_ecological_memory.R"
)
expect_equal(RESILIENCE_MODULES_IN_PIPELINE, expected_modules,
             "RESILIENCE_MODULES_IN_PIPELINE")

# Figure 5 / main-text resistance-recovery claim.
rr <- read_required_csv("tables/table_resistance_recovery.csv")
kelp_metrics <- c("resistance", "recovery", "recovery_recent")
kelp <- rr[rr$taxon == "Macrocystis pyrifera" & rr$metric %in% kelp_metrics, ]
kelp <- kelp[match(kelp_metrics, kelp$metric), ]
expect_equal(kelp$metric, kelp_metrics, "Kelp resistance/recovery metric order")
expect_equal(kelp$n_mpa, c(10L, 10L, 10L), "Kelp resistance/recovery n")
expect_equal(round(kelp$inside, 1), c(2.2, 2.4, 2.1), "Kelp inside baseline ratios")
expect_equal(round(kelp$outside, 1), c(0.9, 0.8, 0.6), "Kelp reference baseline ratios")
expect_equal(round(kelp$paired_p, 3), c(0.010, 0.014, 0.014),
             "Kelp paired Wilcoxon p-values")

# Pair-count evidence used by the Results sentence.
pair_counts <- read_required_csv("tables/table_resistance_recovery_pair_counts.csv")
res_count <- expect_one(pair_counts, pair_counts$metric == "resistance",
                        "resistance pair count")
recent_count <- expect_one(pair_counts, pair_counts$metric == "recovery_recent",
                           "recent recovery pair count")
expect_equal(res_count$n_inside_gt_reference, 9L,
             "Kelp resistance reserves inside > reference")
expect_equal(recent_count$n_inside_gt_reference, 8L,
             "Kelp recent-recovery reserves inside > reference")

# Sensitivity paragraph: primary baseline, alternative baselines, tests, and leave-one-out.
sens <- read_required_csv("tables/table_resistance_recovery_sensitivity.csv")
primary <- sens[sens$baseline == "2010-2013*" & sens$metric %in% kelp_metrics, ]
primary <- primary[match(kelp_metrics, primary$metric), ]
expect_equal(primary$metric, kelp_metrics, "Primary sensitivity metric order")
expect_equal(primary$n_mpa, c(10L, 10L, 10L), "Primary sensitivity n")
expect_equal(round(primary$ratio_in_over_out, 1), c(4.2, 5.2, 4.9),
             "Primary inside/reference ratios")
expect_true(all(primary$ratio_CI_low > 1 & primary$ratio_CI_high > primary$ratio_CI_low),
            "Primary sensitivity CIs must exclude parity.")
expect_true(max(primary$loo_p_max, na.rm = TRUE) <= 0.028,
            "Leave-one-reserve-out p-value drifted above the reported bound.")
expect_true(min(primary$loo_ratio_min, na.rm = TRUE) >= 2.8,
            "Leave-one-reserve-out ratio drifted below the reported bound.")
expect_true(max(sens$p_wilcoxon, na.rm = TRUE) <= 0.047,
            "Alternative-baseline Wilcoxon p-values drifted above the reported bound.")
expect_true(min(round(sens$ratio_in_over_out, 1), na.rm = TRUE) >= 2.5,
            "Alternative-baseline inside/reference ratio drifted below the reported range.")

# SI/discussion support: two-heatwave repeatability, method-vs-data bridge, environmental null,
# temporal-stability distinction, and reserve-level memory.
rep <- read_required_csv("tables/table_heatwave_replication.csv")
rep_kelp <- expect_one(rep, rep$taxon == "Macrocystis pyrifera", "kelp heatwave replication")
expect_true(rep_kelp$MHW1_dlnRR > 0 && rep_kelp$MHW2_dlnRR > 0,
            "Kelp replication contrasts must remain positive.")
expect_true(rep_kelp$MHW1_p < 0.001 && rep_kelp$MHW2_p < 0.001,
            "Kelp replication p-values no longer support both heatwaves.")
expect_true(grepl("REPEATS", rep_kelp$replication),
            "Kelp replication verdict no longer says REPEATS.")

decomp <- read_required_csv("tables/table_s_methods_decomposition.csv")
expect_true(any(grepl("temporal autocorrelation", decomp$contribution)),
            "Methods decomposition is missing the AR1/autocorrelation row.")
kumagai_subtidal <- path.expand("~/kumagai2024-comparison/repo/Processed_data/MLPA_data_summarized_wo_siteblocks.csv")
expect_true(file.exists(kumagai_subtidal),
            paste0("Kumagai processed subtidal data missing at ", kumagai_subtidal,
                   "; cannot verify the paper-grade method-vs-data comparison."))
expect_true(any(grepl("^DATA:", decomp$contribution)),
            "Methods decomposition is missing the Kumagai data-effect row.")

kumagai_cold_spell <- path.expand("~/kumagai2024-comparison/repo/Processed_data/SST/CS_cummulative_intensity_1km.rds")
expect_true(file.exists(kumagai_cold_spell),
            paste0("Kumagai cold-spell grid missing at ", kumagai_cold_spell,
                   "; cannot verify the full environmental-moderator set."))

env <- read_required_csv("tables/table_s_env_moderators.csv")
expect_true(nrow(env) == 35, "Environmental-moderator table is missing expected tests.")
expect_true(all(env$FDR >= 0.05), "At least one environmental moderator now survives FDR.")
env_kelp <- env[env$taxon == "M. pyrifera", ]
expect_true(nrow(env_kelp) >= 7, "Kelp environmental-moderator rows are incomplete.")
expect_true(all(env_kelp$p >= 0.17), "Kelp environmental-moderator p-values drifted below the reported range.")

stab <- read_required_csv("tables/table_resilience_stability.csv")
stab_kelp <- expect_one(stab, stab$taxon == "Macrocystis pyrifera" &
                          stab$era == "disturbance era (>=2013)",
                        "kelp disturbance-era stability")
expect_true(stab_kelp$CV_inside < stab_kelp$CV_outside && stab_kelp$paired_p > 0.05,
            "Kelp stability should remain a nonsignificant trend, not the headline resilience mechanism.")

mem <- read_required_csv("tables/table_ecological_memory.csv")
mem_kelp <- expect_one(mem, mem$taxon == "Macrocystis pyrifera", "kelp ecological memory")
expect_true(mem_kelp$spearman_rho > 0.6 && mem_kelp$spearman_p < 0.05,
            "Kelp reserve-level memory no longer supports the SI/discussion claim.")

cat("  Resilience pipeline integration check passed: main-text Figure 5 and SI-supporting outputs agree.\n")
