# =============================================================================
# 11_results_summary.R
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

# Function to extract comprehensive model summary
extract_model_summary <- function(model, data, response_type, model_name) {
  coef_table <- coef(summary(model))
  stat_col <- if ("tval" %in% colnames(coef_table)) "tval" else "zval"

  # Extract taxa names
  taxa_names <- gsub("^Taxa", "", rownames(coef_table))

  # Sample sizes per taxa
  k_per_taxa <- vapply(taxa_names, function(taxon) {
    sum(data$Taxa == taxon, na.rm = TRUE)
  }, integer(1))

  # Build comprehensive summary
  data.frame(
    model = model_name,
    response = response_type,
    taxa = taxa_names,
    n_effects = k_per_taxa,
    estimate = round(coef_table[, "estimate"], 4),
    std_error = round(coef_table[, "se"], 4),
    test_statistic = round(coef_table[, stat_col], 4),
    df = if ("ddf" %in% colnames(coef_table)) round(coef_table[, "ddf"], 1) else NA,
    p_value = round(coef_table[, "pval"], 6),
    ci_lower = round(coef_table[, "ci.lb"], 4),
    ci_upper = round(coef_table[, "ci.ub"], 4),
    significant = coef_table[, "pval"] < SIGNIFICANCE_ALPHA,
    direction = ifelse(coef_table[, "estimate"] > 0, "positive",
                       ifelse(coef_table[, "estimate"] < 0, "negative", "neutral")),
    qualitative_effect = ifelse(
      coef_table[, "pval"] >= SIGNIFICANCE_ALPHA, "No significant effect",
      ifelse(coef_table[, "estimate"] > 0, "Significant increase inside MPA",
             "Significant decrease inside MPA")
    ),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

# Extract summaries from meta-analysis models (if they exist)
model_results <- data.frame()

if (exists("meta_biomass") && exists("biomass_clean")) {
  bio_summary <- extract_model_summary(meta_biomass, biomass_clean, "Biomass", "Meta-analysis (REML)")
  model_results <- rbind(model_results, bio_summary)
  cat("  Added biomass meta-analysis results\n")
}

if (exists("meta_density") && exists("density_clean")) {
  den_summary <- extract_model_summary(meta_density, density_clean, "Density", "Meta-analysis (REML)")
  model_results <- rbind(model_results, den_summary)
  cat("  Added density meta-analysis results\n")
}

# Add heterogeneity statistics
if (exists("meta_biomass")) {
  model_results$tau2_MPA[model_results$response == "Biomass"] <- round(meta_biomass$sigma2[1], 4)
  model_results$tau2_Source[model_results$response == "Biomass"] <- round(meta_biomass$sigma2[2], 4)
}
if (exists("meta_density")) {
  model_results$tau2_MPA[model_results$response == "Density"] <- round(meta_density$sigma2[1], 4)
  model_results$tau2_Source[model_results$response == "Density"] <- round(meta_density$sigma2[2], 4)
}

# Write model results CSV
model_results_path <- here::here("outputs", "model_results_summary.csv")
write_csv(model_results, model_results_path)
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
  write_csv(replicate_effects, replicate_path)
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
  "pBACIPS analysis. Effect sizes represent the log response ratio (ln[MPA/Reference])",
  "at t=11 years post-MPA implementation.",
  "",
  "**Interpretation:**",
  "- Positive effect size = higher values inside MPA vs reference",
  "- Negative effect size = lower values inside MPA vs reference",
  paste0("- Significance threshold: p < ", SIGNIFICANCE_ALPHA),
  "",
  "---",
  ""
)

# Add meta-analysis summary table
if (nrow(model_results) > 0) {
  md_lines <- c(md_lines,
    "## Meta-Analysis Results (Table 2)",
    "",
    "### Biomass",
    "",
    "| Taxa | k | Estimate | SE | t | p-value | 95% CI | Effect |",
    "|------|---|----------|----|----|---------|--------|--------|"
  )

  bio_rows <- model_results[model_results$response == "Biomass", ]
  for (i in seq_len(nrow(bio_rows))) {
    row <- bio_rows[i, ]
    ci_str <- paste0("[", row$ci_lower, ", ", row$ci_upper, "]")
    sig_str <- ifelse(row$p_value < 0.001, "***",
                     ifelse(row$p_value < 0.01, "**",
                           ifelse(row$p_value < 0.05, "*", "")))
    md_lines <- c(md_lines, paste0(
      "| ", row$taxa, " | ", row$n_effects, " | ", row$estimate, " | ",
      row$std_error, " | ", row$test_statistic, " | ",
      formatC(row$p_value, format = "f", digits = 4), sig_str, " | ",
      ci_str, " | ", row$direction, " |"
    ))
  }

  md_lines <- c(md_lines, "", "### Density", "",
    "| Taxa | k | Estimate | SE | t | p-value | 95% CI | Effect |",
    "|------|---|----------|----|----|---------|--------|--------|"
  )

  den_rows <- model_results[model_results$response == "Density", ]
  for (i in seq_len(nrow(den_rows))) {
    row <- den_rows[i, ]
    ci_str <- paste0("[", row$ci_lower, ", ", row$ci_upper, "]")
    sig_str <- ifelse(row$p_value < 0.001, "***",
                     ifelse(row$p_value < 0.01, "**",
                           ifelse(row$p_value < 0.05, "*", "")))
    md_lines <- c(md_lines, paste0(
      "| ", row$taxa, " | ", row$n_effects, " | ", row$estimate, " | ",
      row$std_error, " | ", row$test_statistic, " | ",
      formatC(row$p_value, format = "f", digits = 4), sig_str, " | ",
      ci_str, " | ", row$direction, " |"
    ))
  }

  md_lines <- c(md_lines, "",
    "*Significance: \\*p<0.05, \\*\\*p<0.01, \\*\\*\\*p<0.001*",
    ""
  )

  # Add heterogeneity statistics
  md_lines <- c(md_lines,
    "### Heterogeneity Statistics",
    "",
    "| Response | tau2 (MPA) | tau2 (Source) |",
    "|----------|------------|---------------|"
  )

  for (resp in c("Biomass", "Density")) {
    row <- model_results[model_results$response == resp, ][1, ]
    md_lines <- c(md_lines, paste0(
      "| ", resp, " | ", row$tau2_MPA, " | ", row$tau2_Source, " |"
    ))
  }
  md_lines <- c(md_lines, "")
}

# Add replicate-level appendix
if (!is.null(replicate_effects) && nrow(replicate_effects) > 0) {
  md_lines <- c(md_lines,
    "---",
    "",
    "## Appendix: Replicate-Level Effect Sizes",
    "",
    "Effect sizes for each MPA-taxa-response combination.",
    ""
  )

  # Group by taxa
  for (taxon in unique(replicate_effects$Taxa)) {
    md_lines <- c(md_lines, paste0("### ", taxon), "")

    taxon_data <- replicate_effects[replicate_effects$Taxa == taxon, ]

    md_lines <- c(md_lines,
      "| MPA | Response | Source | Effect Size | SE | 95% CI | Result |",
      "|-----|----------|--------|-------------|-------|--------|--------|"
    )

    for (i in seq_len(nrow(taxon_data))) {
      row <- taxon_data[i, ]
      ci_str <- paste0("[", row$CI_Lower, ", ", row$CI_Upper, "]")
      md_lines <- c(md_lines, paste0(
        "| ", row$MPA, " | ", row$Response, " | ", row$Source, " | ",
        row$Effect_Size, " | ", row$SE, " | ", ci_str, " | ", row$Qualitative, " |"
      ))
    }
    md_lines <- c(md_lines, "")
  }
}

# Add summary statistics
md_lines <- c(md_lines,
  "---",
  "",
  "## Summary Statistics",
  ""
)

if (!is.null(replicate_effects)) {
  # Count significant effects by direction
  sig_positive <- sum(replicate_effects$Qualitative == "Significant increase", na.rm = TRUE)
  sig_negative <- sum(replicate_effects$Qualitative == "Significant decrease", na.rm = TRUE)
  not_sig <- sum(replicate_effects$Qualitative == "Not significant", na.rm = TRUE)

  md_lines <- c(md_lines,
    "### Effect Direction Summary",
    "",
    paste0("- **Significant increases inside MPA:** ", sig_positive),
    paste0("- **Significant decreases inside MPA:** ", sig_negative),
    paste0("- **Not significant:** ", not_sig),
    paste0("- **Total effects analyzed:** ", nrow(replicate_effects)),
    ""
  )

  # Summary by taxa
  md_lines <- c(md_lines, "### Summary by Taxa", "")

  taxa_summary <- replicate_effects %>%
    dplyr::group_by(Taxa) %>%
    dplyr::summarise(
      n = dplyr::n(),
      n_significant = sum(Significant, na.rm = TRUE),
      mean_effect = round(mean(Effect_Size, na.rm = TRUE), 3),
      .groups = "drop"
    )

  md_lines <- c(md_lines,
    "| Taxa | n | Significant | Mean Effect |",
    "|------|---|-------------|-------------|"
  )

  for (i in seq_len(nrow(taxa_summary))) {
    row <- taxa_summary[i, ]
    md_lines <- c(md_lines, paste0(
      "| ", row$Taxa, " | ", row$n, " | ", row$n_significant,
      " (", round(100*row$n_significant/row$n), "%) | ", row$mean_effect, " |"
    ))
  }
}

md_lines <- c(md_lines, "",
  "---",
  "",
  "*This summary was auto-generated by the analysis pipeline.*",
  "*For questions, contact Emily Donham or Adrian Stier.*"
)

# Write markdown file
md_path <- here::here("docs", "RESULTS_SUMMARY.md")
writeLines(md_lines, md_path)
cat("  Wrote:", md_path, "\n")

# =============================================================================
# 4. SUMMARY
# =============================================================================

cat("\n")
cat("====================================================================\n")
cat("  Results Summary Complete\n")
cat("====================================================================\n")
cat("\nOutputs generated:\n")
cat("  1. outputs/model_results_summary.csv - Meta-analysis test statistics\n")
cat("  2. outputs/replicate_effects.csv    - Effect sizes by MPA-taxa\n")
cat("  3. docs/RESULTS_SUMMARY.md          - Formatted markdown appendix\n")
cat("\n")
