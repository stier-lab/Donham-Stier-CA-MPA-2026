# =============================================================================
# 23_ecological_memory.R
# =============================================================================
#
# PURPOSE:
#   Script 19 showed the MPA kelp-resilience response repeats across two heatwaves
#   ON AVERAGE. This asks a finer question: do the SAME reserves repeat? I.e., is a
#   reserve's response to the first heatwave (2014-16) predictive of its response to
#   the second (2018-20)? A positive across-reserve correlation means reserve
#   performance is a CONSISTENT property (reliable resilience / ecological memory),
#   not a different set of reserves responding each time. We also test whether the
#   second-event response exceeds the first (priming / acclimation at the reserve
#   level) vs simply ongoing recovery.
#
# METRIC (per MPA): mean response ratio lnDiff = ln(MPA / Reference) during
#   MHW1 (2014-2016) and during MHW2 (2018-2020), for giant kelp and purple urchin.
#   Across reserves: Pearson + Spearman correlation MHW1 vs MHW2 (consistency), and
#   a paired test of MHW2 - MHW1 (priming). Small k -> exploratory.
#
#   SCOPE: Southern California Bight.
#
# OUTPUTS (tables/, plots/):
#   table_ecological_memory.csv       - per-taxon MHW1<->MHW2 correlation + paired delta
#   fig_ecological_memory.{pdf,png}   - per-MPA MHW1 vs MHW2 response scatter
#
# AUTHORS: Emily Donham & Adrian Stier
# =============================================================================

cat("Running ecological-memory analysis (do the same reserves repeat?) (script 23)...\n")
source(here::here("code", "R", "00_libraries.R"))
source(here::here("code", "R", "00b_color_palette.R"))
source(here::here("code", "R", "01_utils.R"))
suppressMessages(library(data.table))

MHW1 <- 2014:2016; MHW2 <- 2018:2020
taxa <- c("Macrocystis pyrifera", "Strongylocentrotus purpuratus")
short <- c("Macrocystis pyrifera" = "M. pyrifera", "Strongylocentrotus purpuratus" = "S. purpuratus")

rr <- as.data.table(read.csv(here::here("data", "harmonized", "harmonized_response_ratios.csv"), stringsAsFactors = FALSE))
rr <- rr[is.finite(lnDiff)]

mem_rows <- list(); pts <- list()
for (tx in taxa) {
  s <- rr[y == tx]
  e1 <- s[year %in% MHW1, .(mhw1 = mean(lnDiff)), by = CA_MPA_Name_Short]
  e2 <- s[year %in% MHW2, .(mhw2 = mean(lnDiff)), by = CA_MPA_Name_Short]
  w <- merge(e1, e2, by = "CA_MPA_Name_Short")
  w <- w[is.finite(mhw1) & is.finite(mhw2)]
  if (nrow(w) < 5) next
  pe <- suppressWarnings(cor.test(w$mhw1, w$mhw2))
  sp <- suppressWarnings(cor.test(w$mhw1, w$mhw2, method = "spearman"))
  pt <- suppressWarnings(wilcox.test(w$mhw2, w$mhw1, paired = TRUE))
  mem_rows[[tx]] <- data.frame(taxon = tx, n_mpa = nrow(w),
    pearson_r = round(unname(pe$estimate), 3), pearson_p = signif(pe$p.value, 3),
    spearman_rho = round(unname(sp$estimate), 3), spearman_p = signif(sp$p.value, 3),
    median_MHW1 = round(median(w$mhw1), 3), median_MHW2 = round(median(w$mhw2), 3),
    delta_MHW2_minus_MHW1 = round(median(w$mhw2 - w$mhw1), 3), priming_p = signif(pt$p.value, 3))
  pts[[tx]] <- transform(w, taxon = tx)
}
mem_tab <- do.call(rbind, mem_rows)
write.csv(mem_tab, here::here("tables", "table_ecological_memory.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# Figure: per-MPA MHW1 vs MHW2 response (consistency)
# ---------------------------------------------------------------------------
pd <- do.call(rbind, pts)
if (!is.null(pd)) {
  pd$lab <- factor(short[pd$taxon], levels = short[taxa])
  rng <- range(c(pd$mhw1, pd$mhw2), na.rm = TRUE)
  p <- ggplot(pd, aes(mhw1, mhw2)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey50") +
    geom_hline(yintercept = 0, linewidth = 0.2, color = "grey80") +
    geom_vline(xintercept = 0, linewidth = 0.2, color = "grey80") +
    geom_point(aes(color = lab), size = 1.6, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.6) +
    facet_wrap(~ lab) + scale_color_manual(values = setNames(col_taxa[c(5, 3)], short[taxa]), guide = "none") +
    coord_equal(xlim = rng, ylim = rng) +
    labs(x = "Reserve response in MHW1 (2014-16), lnRR", y = "Reserve response in MHW2 (2018-20), lnRR",
         title = "Ecological memory: do the same reserves respond in both heatwaves?") +
    theme_mpa(base_size = 9) + theme(plot.title = element_text(size = 9, face = "bold"),
                                     strip.text = element_text(face = "italic"))
  ggsave(here::here("plots", "fig_ecological_memory.pdf"), p, width = 150, height = 85, units = "mm", device = cairo_pdf)
  ggsave(here::here("plots", "fig_ecological_memory.png"), p, width = 150, height = 85, units = "mm", dpi = 600)
}

cat("\n=== Ecological memory: reserve-level MHW1 <-> MHW2 consistency ===\n")
print(mem_tab, row.names = FALSE)
cat("\n  Tables + figure written.\n")
