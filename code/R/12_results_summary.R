# =============================================================================
# 12_results_summary.R
# =============================================================================
#
# PURPOSE:
#   Generate comprehensive results summaries for easy review by collaborators.
#   Creates both CSV files and a formatted markdown appendix.
#
# OUTPUTS:
#   1. outputs/model_results_summary.csv - All meta-analysis test statistics
#   2. outputs/replicate_effects.csv - Effect sizes for each MPA-taxa combination
#   3. docs/RESULTS_SUMMARY.md - Formatted markdown with all results
#
# DEPENDENCIES:
#   Requires 00-09 scripts to be sourced first (needs SumStats.Final, meta models, Table2)
#
# AUTHORS: Emily Donham & Adrian Stier
# PROJECT: CA MPA Kelp Forest pBACIPS Analysis
# =============================================================================

cat("\n")
cat("====================================================================\n")
cat("  Generating Results Summary\n")
cat("====================================================================\n")

# Ensure outputs directory exists
if (!dir.exists(here::here("outputs"))) {
  dir.create(here::here("outputs"))
}

# =============================================================================
# 1. META-ANALYSIS SUMMARY (test statistics, p-values, DFs)
# =============================================================================

cat("\n--- Building Meta-Analysis Summary ---\n")

# Build model_results directly from joint-model Table2 (full dataset, no outlier removal).
# This ensures RESULTS_SUMMARY.md matches table_02_meta_analysis.csv exactly.
model_results <- data.frame()

if (exists("Table2") && nrow(Table2) > 0) {
  model_results <- data.frame(
    model = Table2$model_type,
    response = Table2$Response,
    taxa = Table2$Taxa,
    n_effects = Table2$k,
    estimate = round(Table2$Estimate, 4),
    std_error = round(Table2$SE, 4),
    test_statistic = round(Table2$tval, 4),
    df = NA_real_,
    p_value = round(Table2$pval, 6),
    ci_lower = round(Table2$CI_lower, 4),
    ci_upper = round(Table2$CI_upper, 4),
    RR = round(Table2$RR, 3),
    RR_lower = round(exp(Table2$CI_lower), 3),
    RR_upper = round(exp(Table2$CI_upper), 3),
    pct_change = round(Table2$Pct_Change, 1),
    significant = Table2$pval < SIGNIFICANCE_ALPHA,
    direction = ifelse(Table2$Estimate > 0, "positive",
                       ifelse(Table2$Estimate < 0, "negative", "neutral")),
    qualitative_effect = ifelse(
      Table2$pval >= SIGNIFICANCE_ALPHA, "No significant effect",
      ifelse(Table2$Estimate > 0, "Significant increase inside MPA",
             "Significant decrease inside MPA")
    ),
    tau2 = round(Table2$tau2, 4),
    I2 = round(Table2$I2, 1),
    pval_fdr = Table2$pval_fdr,
    flag = if ("Flag" %in% names(Table2)) Table2$Flag else "",
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  cat("  Built model_results from joint-model Table2 (", nrow(model_results), " rows)\n")
} else {
  cat("  WARNING: Table2 not found — model_results will be empty\n")
}

# Write model results CSV
model_results_path <- here::here("outputs", "model_results_summary.csv")
write.csv(model_results, model_results_path, row.names = FALSE)
cat("  Wrote:", model_results_path, "\n")

# =============================================================================
# 2. REPLICATE-LEVEL EFFECTS (every MPA-taxa combination)
# =============================================================================

cat("\n--- Building Replicate Effects Table ---\n")

if (exists("SumStats.Final")) {
  # Create comprehensive replicate effects table
  # Note: CI column is the half-width of the confidence interval
  replicate_effects <- SumStats.Final %>%
    dplyr::mutate(
      CI_Lower = as.numeric(Mean) - as.numeric(CI),
      CI_Upper = as.numeric(Mean) + as.numeric(CI)
    ) %>%
    dplyr::select(
      MPA = MPA,
      Taxa = Taxa,
      Source = Source,
      Response = Resp,
      Effect_Size = Mean,
      SE = SE,
      CI_Lower = CI_Lower,
      CI_Upper = CI_Upper,
      Model_Type = Model
    ) %>%
    dplyr::mutate(
      Effect_Size = round(as.numeric(Effect_Size), 4),
      SE = round(as.numeric(SE), 4),
      CI_Lower = round(CI_Lower, 4),
      CI_Upper = round(CI_Upper, 4),
      Significant = (CI_Lower > 0 | CI_Upper < 0),
      Direction = ifelse(Effect_Size > 0, "Increase inside MPA",
                        ifelse(Effect_Size < 0, "Decrease inside MPA", "No change")),
      Qualitative = ifelse(!Significant, "Not significant",
                          ifelse(Effect_Size > 0, "Significant increase", "Significant decrease"))
    ) %>%
    dplyr::arrange(Taxa, Response, MPA)

  # Write replicate effects CSV
  replicate_path <- here::here("outputs", "replicate_effects.csv")
  write.csv(replicate_effects, replicate_path, row.names = FALSE)
  cat("  Wrote:", replicate_path, "\n")
  cat("  Total replicate effects:", nrow(replicate_effects), "\n")
} else {
  cat("  WARNING: SumStats.Final not found - skipping replicate effects\n")
  replicate_effects <- NULL
}

# =============================================================================
# 3. MARKDOWN SUMMARY DOCUMENT
# =============================================================================

cat("\n--- Building Markdown Summary ---\n")

md_lines <- c(
  "# Results Summary: CA MPA Kelp Forest Analysis",
  "",
  paste0("**Generated:** ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("**Pipeline version:** Modular pBACIPS v2.0"),
  "",
  "---",
  "",
  "## Overview",
  "",
  "This document summarizes the key results from the California MPA kelp forest",
  "pBACIPS analysis. Effect sizes are estimated as log response ratios (lnRR = ln[MPA/Reference])",
  "at t=11 years post-MPA implementation, and back-transformed to response ratios (RR = exp(lnRR)).",
  "",
  "**Interpretation:**",
  "- **lnRR > 0 / RR > 1** = higher values inside MPA vs reference",
  "- **lnRR < 0 / RR < 1** = lower values inside MPA vs reference",
  "- **lnRR = 0 / RR = 1** = no difference",
  "- **% Change** = (RR - 1) x 100 (e.g., RR = 1.72 means +72% higher inside MPA)",
  paste0("- Significance threshold: p < ", SIGNIFICANCE_ALPHA),
  "",
  "---",
  ""
)

# Helper for formatting meta-analysis table rows (no Flag column)
format_meta_row <- function(row) {
  ci_rr_str <- paste0("[", row$RR_lower, ", ", row$RR_upper, "]")
  sig_str <- ifelse(row$p_value < 0.001, "***",
                   ifelse(row$p_value < 0.01, "**",
                         ifelse(row$p_value < 0.05, "*", "")))
  fdr_str <- ifelse(is.na(row$pval_fdr), "\u2014",
                   formatC(row$pval_fdr, format = "f", digits = 4))
  pct_str <- paste0(ifelse(row$pct_change > 0, "+", ""), row$pct_change, "%")
  paste0("| ", row$taxa, " | ", row$n_effects, " | ", row$estimate, " | ",
         row$RR, " | ", pct_str, " | ",
         row$std_error, " | ", row$test_statistic, " | ",
         formatC(row$p_value, format = "f", digits = 4), sig_str, " | ",
         fdr_str, " | ", ci_rr_str, " | ", row$direction, " |")
}

# Add meta-analysis summary table
if (nrow(model_results) > 0) {
  if (!"pval_fdr" %in% names(model_results)) model_results$pval_fdr <- NA

  md_lines <- c(md_lines,
    "## Meta-Analysis Results (Table 2)",
    "",
    "**Estimate** = log response ratio (lnRR). **RR** = back-transformed response ratio = exp(lnRR).",
    "RR > 1 means higher inside MPA; RR < 1 means lower inside MPA.",
    "",
    "### Biomass",
    "",
    "| Taxa | k | Estimate | RR | % Change | SE | t | p-value | p (FDR) | 95% CI (RR) | Effect |",
    "|------|---|----------|----|----------|----|----|---------|---------|-------------|--------|"
  )

  bio_rows <- model_results[model_results$response == "Biomass", ]
  for (i in seq_len(nrow(bio_rows))) md_lines <- c(md_lines, format_meta_row(bio_rows[i, ]))

  md_lines <- c(md_lines, "", "### Density", "",
    "| Taxa | k | Estimate | RR | % Change | SE | t | p-value | p (FDR) | 95% CI (RR) | Effect |",
    "|------|---|----------|----|----------|----|----|---------|---------|-------------|--------|"
  )

  den_rows <- model_results[model_results$response == "Density", ]
  for (i in seq_len(nrow(den_rows))) md_lines <- c(md_lines, format_meta_row(den_rows[i, ]))

  md_lines <- c(md_lines, "",
    "*Significance: \\*p<0.05, \\*\\*p<0.01, \\*\\*\\*p<0.001 (uncorrected). p (FDR) = Benjamini-Hochberg adjusted p-values.*",
    ""
  )

  # Heterogeneity statistics
  md_lines <- c(md_lines,
    "### Heterogeneity (Joint Model)",
    "",
    "tau2 and pseudo-I2 are shared across all taxa within each response type.",
    "",
    "| Response | tau2 | I2 (%) |",
    "|----------|------|--------|"
  )
  for (resp in unique(model_results$response)) {
    row <- model_results[model_results$response == resp, ][1, ]
    md_lines <- c(md_lines, paste0("| ", resp, " | ", row$tau2, " | ", row$I2, " |"))
  }
  md_lines <- c(md_lines, "")
}

# ---- Summary statistics (replicate-level) ----
if (!is.null(replicate_effects)) {
  sig_positive <- sum(replicate_effects$Qualitative == "Significant increase", na.rm = TRUE)
  sig_negative <- sum(replicate_effects$Qualitative == "Significant decrease", na.rm = TRUE)
  not_sig <- sum(replicate_effects$Qualitative == "Not significant", na.rm = TRUE)

  md_lines <- c(md_lines,
    "---", "",
    "## Replicate-Level Summary",
    "",
    paste0("- **Significant increases inside MPA:** ", sig_positive),
    paste0("- **Significant decreases inside MPA:** ", sig_negative),
    paste0("- **Not significant:** ", not_sig),
    paste0("- **Total MPA-level effect sizes:** ", nrow(replicate_effects)),
    "",
    "### Per-Taxa Breakdown",
    "",
    "| Taxa | n | Significant | Raw Mean lnRR |",
    "|------|---|-------------|---------------|"
  )

  taxa_summary <- replicate_effects %>%
    dplyr::group_by(Taxa) %>%
    dplyr::summarise(
      n = dplyr::n(),
      n_significant = sum(Significant, na.rm = TRUE),
      mean_effect = round(mean(Effect_Size, na.rm = TRUE), 3),
      .groups = "drop"
    )
  for (i in seq_len(nrow(taxa_summary))) {
    row <- taxa_summary[i, ]
    md_lines <- c(md_lines, paste0(
      "| ", row$Taxa, " | ", row$n, " | ", row$n_significant,
      " (", round(100*row$n_significant/row$n), "%) | ", row$mean_effect, " |"
    ))
  }

  md_lines <- c(md_lines, "",
    "*Raw Mean lnRR is the arithmetic mean of individual MPA effect sizes (bio + den combined), not the meta-analytic estimate. See Table 2 above for meta-analytic means.*",
    "",
    "Full replicate-level effect sizes are in `outputs/replicate_effects.csv`.",
    ""
  )
}

# ---- Cross-taxa meta-regression (Table 3) ----
table3_path <- here::here("tables", "table_03_cross_taxa_meta_regression.csv")
if (file.exists(table3_path)) {
  table3_data <- read.csv(table3_path, stringsAsFactors = FALSE)
  md_lines <- c(md_lines,
    "---", "",
    "## Cross-Taxa Meta-Regressions (Table 3)",
    "",
    "Knapp-Hartung meta-regressions testing whether predator MPA effects predict prey effects across sites.",
    "",
    "| Relationship | k | Slope | SE | p | Significant? |",
    "|-------------|---|-------|----|---|-------------|"
  )
  for (i in seq_len(nrow(table3_data))) {
    row <- table3_data[i, ]
    sig <- ifelse(row$Slope_pval < 0.05, "**yes**", "no")
    md_lines <- c(md_lines, paste0(
      "| ", row$Relationship, " | ", row$k, " | ",
      round(row$Slope_Est, 3), " | ", round(row$Slope_SE, 3), " | ",
      formatC(row$Slope_pval, format = "f", digits = 3), " | ", sig, " |"
    ))
  }
  md_lines <- c(md_lines, "")
}

# ---- Variance components ----
vc_path <- here::here("tables", "table_s_variance_components.csv")
if (file.exists(vc_path)) {
  vc_data <- read.csv(vc_path, stringsAsFactors = FALSE)
  md_lines <- c(md_lines,
    "---", "",
    "## Variance Components",
    "",
    "| Response | Component | Estimate | 95% CI |",
    "|----------|-----------|----------|--------|"
  )
  for (i in seq_len(nrow(vc_data))) {
    row <- vc_data[i, ]
    ci_lo <- if ("CI_Lower" %in% names(row)) round(row$CI_Lower, 4) else NA
    ci_hi <- if ("CI_Upper" %in% names(row)) round(row$CI_Upper, 4) else NA
    ci_str <- if (!is.na(ci_lo) && !is.na(ci_hi)) paste0("[", ci_lo, ", ", ci_hi, "]") else "\u2014"
    md_lines <- c(md_lines, paste0(
      "| ", row$Response, " | ", row$Component, " | ", round(row$Estimate, 4), " | ", ci_str, " |"
    ))
  }
  md_lines <- c(md_lines, "")
}

# ---- Temporal recovery slopes ----
slopes_path <- here::here("tables", "table_s_species_slopes.csv")
if (file.exists(slopes_path)) {
  slopes_data <- read.csv(slopes_path, stringsAsFactors = FALSE)
  pooled <- slopes_data[slopes_data$model == "pooled", ]
  if (nrow(pooled) > 0) {
    md_lines <- c(md_lines,
      "---", "",
      "## Temporal Recovery Slopes (Pooled lmer Model)",
      "",
      "Species-specific slopes from `lmer(lnRR ~ time*species + (1+time|MPA) + (1|source))`.",
      "",
      "| Species | Slope (lnRR/yr) | SE | p |",
      "|---------|-----------------|-----|---|"
    )
    for (i in seq_len(nrow(pooled))) {
      row <- pooled[i, ]
      p_str <- if (row$p_value < 0.001) "< 0.001" else formatC(row$p_value, format = "f", digits = 3)
      md_lines <- c(md_lines, paste0(
        "| ", row$Species_abbrev, " | ", round(row$slope, 4), " | ",
        round(row$SE, 4), " | ", p_str, " |"
      ))
    }
    md_lines <- c(md_lines, "")
  }
}

md_lines <- c(md_lines,
  "---", "",
  "For known limitations, see `docs/ANALYSIS_REVISIONS.md`.",
  "",
  "*This summary was auto-generated by the analysis pipeline.*",
  "*For questions, contact Emily Donham or Adrian Stier.*"
)

# Write markdown file
md_path <- here::here("docs", "RESULTS_SUMMARY.md")
writeLines(md_lines, md_path)
cat("  Wrote:", md_path, "\n")

# =============================================================================
# 4. COMBINED DATA FLOW SUMMARY
# =============================================================================

cat("\n--- Building Combined Data Flow Summary ---\n")

# Read filter audit files if they exist
filter_audit_file <- here::here("outputs", "filter_audit_effect_sizes.csv")
meta_audit_file <- here::here("outputs", "filter_audit_meta_analysis.csv")

if (file.exists(filter_audit_file) && file.exists(meta_audit_file)) {

  # Read audits
  effect_audit <- read.csv(filter_audit_file, stringsAsFactors = FALSE)
  meta_audit <- read.csv(meta_audit_file, stringsAsFactors = FALSE)

  # Create comprehensive data flow table
  cat("\n")
  cat("============================================\n")
  cat("COMPLETE DATA FLOW: Effect Size Generation\n")
  cat("============================================\n")

  # Stage 1: Effect sizes generated
  stage1 <- effect_audit %>%
    dplyr::group_by(Taxa, Resp) %>%
    dplyr::summarise(
      Stage1_Generated = dplyr::n(),
      .groups = "drop"
    )

  # Stage 2: Pass site merge
  stage2 <- effect_audit %>%
    dplyr::filter(Step1_SiteMerge) %>%
    dplyr::group_by(Taxa, Resp) %>%
    dplyr::summarise(
      Stage2_SiteMerge = dplyr::n(),
      .groups = "drop"
    )

  # Stage 3: Pass type filter (pBACIPS or CI+Primary)
  stage3 <- effect_audit %>%
    dplyr::filter(Step2_TypeFilter) %>%
    dplyr::group_by(Taxa, Resp) %>%
    dplyr::summarise(
      Stage3_TypeFilter = dplyr::n(),
      .groups = "drop"
    )

  # Stage 4: Pass MPA exclusion filter
  stage4 <- effect_audit %>%
    dplyr::filter(Passes_All) %>%
    dplyr::group_by(Taxa, Resp) %>%
    dplyr::summarise(
      Stage4_SumStatsFinal = dplyr::n(),
      .groups = "drop"
    )

  # Stage 5: Enter meta-analysis (same as Stage 4)
  # Stage 6: Pass outlier filter (in meta_audit)
  stage6 <- meta_audit %>%
    dplyr::filter(In_Final_Analysis) %>%
    dplyr::group_by(Taxa, Response) %>%
    dplyr::summarise(
      Stage6_FinalK = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::rename(Resp = Response) %>%
    dplyr::mutate(Resp = ifelse(Resp == "Biomass", "Bio", "Den"))

  # Merge all stages
  data_flow <- stage1 %>%
    dplyr::left_join(stage2, by = c("Taxa", "Resp")) %>%
    dplyr::left_join(stage3, by = c("Taxa", "Resp")) %>%
    dplyr::left_join(stage4, by = c("Taxa", "Resp")) %>%
    dplyr::left_join(stage6, by = c("Taxa", "Resp")) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~replace(., is.na(.), 0))) %>%
    dplyr::arrange(Taxa, Resp)

  # Print data flow table
  cat("\nData Flow by Taxa and Response:\n")
  cat("  Stage 1: Effect sizes generated (08_effect_sizes.R)\n")
  cat("  Stage 2: Pass Site metadata merge\n")
  cat("  Stage 3: Pass Type filter (pBACIPS or CI+Primary)\n")
  cat("  Stage 4: Pass MPA exclusion filter -> SumStats.Final\n")
  cat("  Stage 5: Enter meta-analysis\n")
  cat("  Stage 6: Pass outlier filter -> Final k\n")
  cat("\n")

  # Print as formatted table
  cat(sprintf("%-20s %-5s %8s %8s %8s %8s %8s\n",
              "Taxa", "Resp", "Stage1", "Stage2", "Stage3", "Stage4", "Final_k"))
  cat(paste(rep("-", 75), collapse = ""), "\n")

  for (i in seq_len(nrow(data_flow))) {
    row <- data_flow[i, ]
    cat(sprintf("%-20s %-5s %8d %8d %8d %8d %8d\n",
                row$Taxa, row$Resp,
                row$Stage1_Generated, row$Stage2_SiteMerge,
                row$Stage3_TypeFilter, row$Stage4_SumStatsFinal,
                row$Stage6_FinalK))
  }

  # Write data flow summary
  data_flow_file <- here::here("outputs", "data_flow_summary.csv")
  write.csv(data_flow, data_flow_file, row.names = FALSE)
  cat("\nData flow summary written to:", data_flow_file, "\n")

  # Identify attrition points
  cat("\n--- Major Attrition Points ---\n")

  for (i in seq_len(nrow(data_flow))) {
    row <- data_flow[i, ]
    attrition <- c()

    if (row$Stage2_SiteMerge < row$Stage1_Generated) {
      lost <- row$Stage1_Generated - row$Stage2_SiteMerge
      attrition <- c(attrition, sprintf("Site merge: -%d", lost))
    }
    if (row$Stage3_TypeFilter < row$Stage2_SiteMerge) {
      lost <- row$Stage2_SiteMerge - row$Stage3_TypeFilter
      attrition <- c(attrition, sprintf("Type filter: -%d", lost))
    }
    if (row$Stage4_SumStatsFinal < row$Stage3_TypeFilter) {
      lost <- row$Stage3_TypeFilter - row$Stage4_SumStatsFinal
      attrition <- c(attrition, sprintf("MPA exclusion: -%d", lost))
    }
    if (row$Stage6_FinalK < row$Stage4_SumStatsFinal) {
      lost <- row$Stage4_SumStatsFinal - row$Stage6_FinalK
      attrition <- c(attrition, sprintf("Outlier removal: -%d", lost))
    }

    if (length(attrition) > 0) {
      cat(sprintf("%s (%s): %s\n", row$Taxa, row$Resp, paste(attrition, collapse = ", ")))
    }
  }

  # Detailed exclusion reasons
  cat("\n--- Exclusion Reason Summary ---\n")
  reason_counts <- table(effect_audit$Exclusion_Reason)
  for (reason in names(sort(reason_counts, decreasing = TRUE))) {
    cat(sprintf("  %s: %d\n", reason, reason_counts[reason]))
  }

} else {
  cat("  Filter audit files not found - run 08_effect_sizes.R and 09_meta_analysis.R first\n")
}

# =============================================================================
# 5. SAMPLE SIZE & DATA PROVENANCE TABLES (Conservation Letters requirements)
# =============================================================================
# A Conservation Letters manuscript must fully document:
#   - Total number of MPAs analyzed
#   - Number of MPAs per data source
#   - Number of reference sites per MPA
#   - Years of data per MPA (before/after implementation)
#   - Number of observations per species x response type
#   - Any data exclusions and reasons
#
# These tables complement the existing filter_audit and data_flow CSVs by
# providing concise, publication-ready summaries.

cat("\n--- Building Sample Size & Data Provenance Tables ---\n")

# -------------------------------------------------------------------------
# 5a. MPA-level data availability matrix (Table S1a)
# -------------------------------------------------------------------------
# For each MPA in the final analysis: which taxa have data, data source,
# MPA type, implementation year, location, N time-series observations,
# years of data before/after implementation.

if (exists("SumStats.Final") && exists("All.RR.sub.trans") && exists("Site")) {

  # Build per-MPA x taxa availability from SumStats.Final
  mpa_avail <- SumStats.Final %>%
    dplyr::group_by(MPA) %>%
    dplyr::summarise(
      Data_Sources   = paste(sort(unique(Source)), collapse = ", "),
      Taxa_Biomass   = paste(sort(unique(Taxa[Resp == "Bio"])), collapse = "; "),
      Taxa_Density   = paste(sort(unique(Taxa[Resp == "Den"])), collapse = "; "),
      n_taxa_bio     = dplyr::n_distinct(Taxa[Resp == "Bio"]),
      n_taxa_den     = dplyr::n_distinct(Taxa[Resp == "Den"]),
      n_effects_bio  = sum(Resp == "Bio"),
      n_effects_den  = sum(Resp == "Den"),
      n_effects_total = dplyr::n(),
      .groups = "drop"
    )

  # Merge MPA metadata from Site
  # Detect which column in Site holds the MPA name
  site_mpa_col <- if ("CA_MPA_Name_Short" %in% names(Site)) "CA_MPA_Name_Short" else "Site"
  site_meta <- Site %>%
    dplyr::distinct(.data[[site_mpa_col]], .keep_all = TRUE) %>%
    dplyr::select(dplyr::any_of(c(site_mpa_col, "type", "Location", "MPA_Start",
                                    "Lat", "Lon", "Hectares")))
  names(site_meta)[names(site_meta) == site_mpa_col] <- "MPA"

  mpa_avail <- dplyr::left_join(mpa_avail, site_meta, by = "MPA")

  # Add years before/after from response ratio data
  # Compute from All.RR.sub.trans using BA column
  ba_summary <- All.RR.sub.trans %>%
    dplyr::filter(CA_MPA_Name_Short %in% mpa_avail$MPA) %>%
    dplyr::group_by(CA_MPA_Name_Short) %>%
    dplyr::summarise(
      year_min        = min(year, na.rm = TRUE),
      year_max        = max(year, na.rm = TRUE),
      n_years_total   = dplyr::n_distinct(year),
      n_years_before  = dplyr::n_distinct(year[BA == "Before"]),
      n_years_after   = dplyr::n_distinct(year[BA == "After"]),
      n_obs_rr        = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::rename(MPA = CA_MPA_Name_Short)

  mpa_avail <- dplyr::left_join(mpa_avail, ba_summary, by = "MPA")

  # Add reference site counts from All.Resp.sub (if available)
  if (exists("All.Resp.sub")) {
    status_col <- if ("status" %in% names(All.Resp.sub)) "status"
                  else if ("Status" %in% names(All.Resp.sub)) "Status"
                  else if ("site_status" %in% names(All.Resp.sub)) "site_status"
                  else NULL
    if (!is.null(status_col)) {
      ref_counts <- All.Resp.sub %>%
        dplyr::filter(CA_MPA_Name_Short %in% mpa_avail$MPA,
                      tolower(.data[[status_col]]) %in% c("outside", "reference", "control", "ref", "o", "r")) %>%
        dplyr::group_by(CA_MPA_Name_Short) %>%
        dplyr::summarise(
          n_reference_site_years = dplyr::n(),
          .groups = "drop"
        ) %>%
        dplyr::rename(MPA = CA_MPA_Name_Short)
      mpa_avail <- dplyr::left_join(mpa_avail, ref_counts, by = "MPA")
    }
  }

  # Sort by location then MPA name
  mpa_avail <- mpa_avail %>%
    dplyr::arrange(Location, MPA)

  mpa_avail_path <- here::here("tables", "table_s_mpa_data_availability.csv")
  write.csv(mpa_avail, mpa_avail_path, row.names = FALSE)
  cat("  Wrote:", mpa_avail_path, "\n")
  cat("  Total MPAs in final analysis:", nrow(mpa_avail), "\n")

} else {
  cat("  WARNING: Required objects not found for MPA availability table\n")
}

# -------------------------------------------------------------------------
# 5b. Sample sizes by taxa x response, broken down by data source (Table S1b)
# -------------------------------------------------------------------------
# For each taxa x response type: total k (number of MPA-level effect sizes),
# number of unique MPAs, broken down by data source.

if (exists("SumStats.Final")) {

  # Per taxa x response totals
  taxa_k_total <- SumStats.Final %>%
    dplyr::group_by(Taxa, Resp) %>%
    dplyr::summarise(
      k_total     = dplyr::n(),
      n_MPAs      = dplyr::n_distinct(MPA),
      Sources     = paste(sort(unique(Source)), collapse = ", "),
      .groups = "drop"
    )

  # Per taxa x response x source breakdown
  taxa_k_by_source <- SumStats.Final %>%
    dplyr::group_by(Taxa, Resp, Source) %>%
    dplyr::summarise(
      k = dplyr::n(),
      n_MPAs = dplyr::n_distinct(MPA),
      MPAs = paste(sort(unique(MPA)), collapse = "; "),
      .groups = "drop"
    )

  # Pivot to wide format: one row per taxa x response, columns for each source's k
  taxa_k_wide <- taxa_k_by_source %>%
    dplyr::select(Taxa, Resp, Source, k) %>%
    tidyr::pivot_wider(
      names_from  = Source,
      values_from = k,
      values_fill = 0,
      names_prefix = "k_"
    )

  # Merge with totals
  sample_sizes <- dplyr::left_join(taxa_k_total, taxa_k_wide, by = c("Taxa", "Resp"))

  # Also add meta-analysis final k (after outlier removal) from data_flow if available
  # Use in-memory data_flow object (created in section 4 above) if available;
  # fall back to reading from disk if the section was skipped.
  if (exists("data_flow") && is.data.frame(data_flow) && "Stage6_FinalK" %in% names(data_flow)) {
    data_flow_k <- data_flow %>%
      dplyr::select(Taxa, Resp, k_final_meta = Stage6_FinalK)
    sample_sizes <- dplyr::left_join(sample_sizes, data_flow_k, by = c("Taxa", "Resp"))
  } else {
    data_flow_file <- here::here("outputs", "data_flow_summary.csv")
    if (file.exists(data_flow_file)) {
      data_flow_disk <- read.csv(data_flow_file, stringsAsFactors = FALSE)
      data_flow_k <- data_flow_disk %>%
        dplyr::select(Taxa, Resp, k_final_meta = Stage6_FinalK)
      sample_sizes <- dplyr::left_join(sample_sizes, data_flow_k, by = c("Taxa", "Resp"))
    }
  }

  # Order by canonical trophic order
  trophic_order <- c("P. interruptus", "S. pulcher", "S. purpuratus",
                     "M. franciscanus", "M. pyrifera")
  sample_sizes <- sample_sizes %>%
    dplyr::mutate(Taxa = factor(Taxa, levels = trophic_order)) %>%
    dplyr::arrange(Taxa, Resp)

  sample_sizes_path <- here::here("tables", "table_s_sample_sizes.csv")
  write.csv(sample_sizes, sample_sizes_path, row.names = FALSE)
  cat("  Wrote:", sample_sizes_path, "\n")
  cat("  Taxa x Response combinations:", nrow(sample_sizes), "\n")

} else {
  cat("  WARNING: SumStats.Final not found for sample sizes table\n")
}

# -------------------------------------------------------------------------
# 5c. High-level analysis summary (Table S1c / Methods paragraph)
# -------------------------------------------------------------------------
# Single-row summary with key numbers for the manuscript methods section.

if (exists("SumStats.Final") && exists("All.RR.sub.trans")) {

  analysis_summary <- data.frame(
    Metric = c(
      "Total unique MPAs analyzed",
      "Total MPA x taxa x response combinations (effect sizes)",
      "Effect sizes - Biomass",
      "Effect sizes - Density",
      "Unique taxa",
      "Data sources",
      "MPAs from PISCO",
      "MPAs from KFM",
      "MPAs from LTER",
      "MPAs from Landsat",
      "Effect sizes from before-after analysis (pBACIPS)",
      "Effect sizes from control-impact only (CI)",
      "Channel Islands MPAs",
      "Mainland MPAs",
      "SMR (State Marine Reserve) MPAs",
      "SMCA (State Marine Conservation Area) MPAs",
      "Earliest year in dataset",
      "Latest year in dataset",
      "MPA implementation years"
    ),
    Value = c(
      length(unique(SumStats.Final$MPA)),
      nrow(SumStats.Final),
      sum(SumStats.Final$Resp == "Bio"),
      sum(SumStats.Final$Resp == "Den"),
      length(unique(SumStats.Final$Taxa)),
      paste(sort(unique(SumStats.Final$Source)), collapse = ", "),
      length(unique(SumStats.Final$MPA[SumStats.Final$Source == "PISCO"])),
      length(unique(SumStats.Final$MPA[SumStats.Final$Source == "KFM"])),
      length(unique(SumStats.Final$MPA[SumStats.Final$Source == "LTER"])),
      length(unique(SumStats.Final$MPA[SumStats.Final$Source == "Landsat"])),
      sum(SumStats.Final$AnalysisType == "pBACIPS"),
      sum(SumStats.Final$AnalysisType == "CI"),
      tryCatch(
        length(unique(SumStats.Final$MPA[SumStats.Final$Location == "C"])),
        error = function(e) NA
      ),
      tryCatch(
        length(unique(SumStats.Final$MPA[SumStats.Final$Location == "ML"])),
        error = function(e) NA
      ),
      tryCatch(
        length(unique(SumStats.Final$MPA[grepl("SMR", SumStats.Final$type)])),
        error = function(e) NA
      ),
      tryCatch(
        length(unique(SumStats.Final$MPA[grepl("SMCA", SumStats.Final$type)])),
        error = function(e) NA
      ),
      min(All.RR.sub.trans$year, na.rm = TRUE),
      max(All.RR.sub.trans$year, na.rm = TRUE),
      paste(sort(unique(SumStats.Final$MPA_Start)), collapse = ", ")
    ),
    stringsAsFactors = FALSE
  )

  # Add data exclusion summary
  exclusion_rows <- data.frame(
    Metric = c(
      "--- DATA EXCLUSIONS ---",
      "Excluded KFM sites (short/non-paired time series)",
      "Excluded reference sites (habitat/sampling issues)",
      "Sheephead-only MPAs (excluded from cross-taxa, included for S. pulcher)",
      "Excluded MPAs (fishing allowed or data quality)",
      "BACI-type analyses excluded (linear before period)",
      "Non-primary model selections excluded"
    ),
    Value = c(
      "",
      paste(EXCLUDED_KFM_SITES, collapse = ", "),
      paste(EXCLUDED_REFERENCE_SITES, collapse = ", "),
      paste(SHEEPHEAD_ONLY_MPAS, collapse = ", "),
      paste(EXCLUDED_MPAS[!EXCLUDED_MPAS %in% c("N/A")], collapse = ", "),
      "See filter_audit_effect_sizes.csv for details",
      "See filter_audit_effect_sizes.csv for details"
    ),
    stringsAsFactors = FALSE
  )

  analysis_summary <- rbind(analysis_summary, exclusion_rows)

  analysis_summary_path <- here::here("tables", "table_s_analysis_summary.csv")
  write.csv(analysis_summary, analysis_summary_path, row.names = FALSE)
  cat("  Wrote:", analysis_summary_path, "\n")

} else {
  cat("  WARNING: Required objects not found for analysis summary table\n")
}

# -------------------------------------------------------------------------
# 5d. Per-MPA detailed effect size table with source breakdown
# -------------------------------------------------------------------------
# Cross-tabulation: for each MPA, show which taxa x response combinations
# are available and from which source. This is the "data availability matrix"
# that Conservation Letters reviewers expect in supplementary materials.

if (exists("SumStats.Final")) {

  # Define canonical trophic order for column ordering
  trophic_taxa <- c("P. interruptus", "S. pulcher", "S. purpuratus",
                     "M. franciscanus", "M. pyrifera")
  trophic_resp <- c("Bio", "Den")
  canonical_col_order <- as.vector(outer(trophic_taxa, trophic_resp,
                                          function(t, r) paste0(t, " (", r, ")")))

  mpa_taxa_matrix <- SumStats.Final %>%
    dplyr::mutate(
      taxa_resp = factor(paste0(Taxa, " (", Resp, ")"), levels = canonical_col_order),
      source_model = paste0(Source, "/", Model)
    ) %>%
    dplyr::select(MPA, taxa_resp, source_model) %>%
    tidyr::pivot_wider(
      names_from  = taxa_resp,
      values_from = source_model,
      values_fill = "",
      values_fn   = function(x) paste(x, collapse = "; ")
    ) %>%
    dplyr::arrange(MPA)

  # Reorder columns: MPA first, then canonical trophic order
  col_order <- c("MPA", intersect(canonical_col_order, names(mpa_taxa_matrix)))
  mpa_taxa_matrix <- mpa_taxa_matrix[, col_order]

  mpa_matrix_path <- here::here("tables", "table_s_mpa_taxa_matrix.csv")
  write.csv(mpa_taxa_matrix, mpa_matrix_path, row.names = FALSE)
  cat("  Wrote:", mpa_matrix_path, "\n")
  cat("  MPAs x taxa combinations documented\n")

} else {
  cat("  WARNING: SumStats.Final not found for MPA-taxa matrix\n")
}

cat("\n--- Sample Size & Provenance Tables Complete ---\n")

# =============================================================================
# 6. SUMMARY
# =============================================================================

cat("\n")
cat("====================================================================\n")
cat("  Results Summary Complete\n")
cat("====================================================================\n")
cat("\nPipeline outputs (generated across scripts 08-12):\n")
cat("  1. outputs/model_results_summary.csv          - Meta-analysis test statistics (12)\n")
cat("  2. outputs/replicate_effects.csv              - Effect sizes by MPA-taxa (12)\n")
cat("  3. outputs/filter_audit_effect_sizes.csv      - Detailed filtering trace (08)\n")
cat("  4. outputs/filter_audit_meta_analysis.csv     - Meta-analysis outlier trace (09)\n")
cat("  5. outputs/filter_summary_by_taxa.csv         - Taxa-level filtering summary (08)\n")
cat("  6. outputs/data_flow_summary.csv              - Complete data flow by stage (12)\n")
cat("  7. docs/RESULTS_SUMMARY.md                    - Formatted markdown appendix (12)\n")
cat("  8. tables/table_s_mpa_data_availability.csv   - Per-MPA data availability (12)\n")
cat("  9. tables/table_s_sample_sizes.csv            - Sample sizes by taxa x source (12)\n")
cat(" 10. tables/table_s_analysis_summary.csv        - High-level analysis summary (12)\n")
cat(" 11. tables/table_s_mpa_taxa_matrix.csv         - MPA x taxa cross-tabulation (12)\n")
cat(" 12. outputs/lmer_residual_diagnostics.csv      - lmer residual diagnostics (10)\n")
cat(" 13. outputs/lmer_ranef_diagnostics.csv         - lmer random effect BLUPs (10)\n")
cat("\n")
