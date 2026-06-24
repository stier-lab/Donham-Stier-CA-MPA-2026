# =============================================================================
# 21_resilience_stability.R
# =============================================================================
#
# PURPOSE:
#   Temporal STABILITY as a resilience facet. Scripts 14/19/20 test the MEAN MPA
#   effect (how much more kelp / fewer urchins inside reserves) and its response to
#   disturbances. Stability is a distinct dimension of resilience: do reserves damp
#   year-to-year VARIABILITY, holding the system in a more constant state? Here we
#   ask whether kelp and the cascade taxa are temporally LESS variable inside MPAs
#   than at paired reference sites (lower coefficient of variation), over the full
#   record and through the disturbance era.
#
# METRIC:
#   Per MPA x status (mpa / reference) x taxon, the temporal coefficient of variation
#   CV = sd / mean of annual abundance (>= 5 years required). CV normalizes for the
#   mean, so a "more stable inside" result is not just "more abundant inside". We
#   compare CV inside vs the paired reference across MPAs (paired Wilcoxon signed-rank
#   on log CV), full series (>=2002) and the disturbance era (>=2013).
#
#   SCOPE: Southern California Bight. Density for animals, biomass for giant kelp.
#
# OUTPUTS (tables/, plots/):
#   table_resilience_stability.csv  - per taxon: median CV inside/outside, paired p
#   fig_resilience_stability.{pdf,png} - paired CV inside vs outside, by taxon
#
# AUTHORS: Emily Donham & Adrian Stier
# =============================================================================

cat("Running temporal-stability resilience analysis (script 21)...\n")
source(here::here("code", "R", "00_libraries.R"))
source(here::here("code", "R", "00b_color_palette.R"))
source(here::here("code", "R", "01_utils.R"))
suppressMessages(library(data.table))

taxa <- c("Panulirus interruptus", "Semicossyphus pulcher", "Strongylocentrotus purpuratus",
          "Mesocentrotus franciscanus", "Macrocystis pyrifera")
short <- c("Panulirus interruptus" = "P. interruptus", "Semicossyphus pulcher" = "S. pulcher",
           "Strongylocentrotus purpuratus" = "S. purpuratus",
           "Mesocentrotus franciscanus" = "M. franciscanus", "Macrocystis pyrifera" = "M. pyrifera")
resp_of <- c("Panulirus interruptus" = "Den", "Semicossyphus pulcher" = "Den",
             "Strongylocentrotus purpuratus" = "Den", "Mesocentrotus franciscanus" = "Den",
             "Macrocystis pyrifera" = "Bio")
MINYRS <- 5

d <- as.data.table(read.csv(here::here("data", "harmonized", "harmonized_raw_responses.csv"), stringsAsFactors = FALSE))

cv_table <- function(year_min, era_label) {
  rows <- list()
  for (tx in taxa) {
    s <- d[taxon_name == tx & resp == resp_of[tx] & year >= year_min]
    # collapse to ONE value per year (mean over monitoring programs) so the temporal CV
    # is not contaminated by within-year between-program variance, and .N counts years
    s <- s[, .(value = mean(value)), by = .(CA_MPA_Name_Short, status, year)]
    # CV per MPA x status (require >= MINYRS distinct years and positive mean)
    cv <- s[, .(cv = if (.N >= MINYRS && mean(value) > 0) sd(value) / mean(value) else NA_real_,
                ny = .N), by = .(CA_MPA_Name_Short, status)]
    w <- dcast(cv[!is.na(cv)], CA_MPA_Name_Short ~ status, value.var = "cv")
    if (!all(c("mpa", "reference") %in% names(w))) next
    w <- w[is.finite(mpa) & is.finite(reference)]
    if (nrow(w) < 5) next
    pt <- suppressWarnings(wilcox.test(w$mpa, w$reference, paired = TRUE))
    rows[[tx]] <- data.frame(era = era_label, taxon = tx, n_mpa = nrow(w),
      CV_inside = round(median(w$mpa), 3), CV_outside = round(median(w$reference), 3),
      pct_lower_inside = round(100 * mean(w$mpa < w$reference)),
      paired_p = signif(pt$p.value, 3),
      stabilizing = ifelse(median(w$mpa) < median(w$reference), "more stable inside", "not"))
    rows[[tx]]$.pairs <- list(w)   # keep for figure
  }
  rows
}
full <- cv_table(2002, "full series (>=2002)")
dist <- cv_table(2013, "disturbance era (>=2013)")

mk <- function(rows) do.call(rbind, lapply(rows, function(r) r[, setdiff(names(r), ".pairs")]))
stab_tab <- rbind(mk(full), mk(dist))
write.csv(stab_tab, here::here("tables", "table_resilience_stability.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# Figure: paired CV inside vs outside (full series), by taxon
# ---------------------------------------------------------------------------
pairs_df <- do.call(rbind, lapply(names(full), function(tx) {
  w <- full[[tx]]$.pairs[[1]]; if (is.null(w)) return(NULL)
  data.frame(taxon = tx, MPA = w$CA_MPA_Name_Short, inside = w$mpa, outside = w$reference)
}))
if (!is.null(pairs_df)) {
  pl <- rbind(data.frame(taxon = pairs_df$taxon, MPA = pairs_df$MPA, site = "Reference", cv = pairs_df$outside),
              data.frame(taxon = pairs_df$taxon, MPA = pairs_df$MPA, site = "MPA", cv = pairs_df$inside))
  pl$lab <- factor(short[pl$taxon], levels = short[taxa])
  pl$site <- factor(pl$site, levels = c("Reference", "MPA"))
  p <- ggplot(pl, aes(site, cv, group = MPA)) +
    geom_line(color = "grey75", linewidth = 0.3) +
    geom_point(aes(color = site), size = 1.3) +
    facet_wrap(~ lab, nrow = 1, scales = "free_y") +
    scale_color_manual(values = c(Reference = "#D55E00", MPA = "#0072B2"), guide = "none") +
    labs(x = NULL, y = "Temporal CV of annual abundance (lower = more stable)",
         title = "Does MPA protection stabilize the kelp forest? (paired CV inside vs outside)") +
    theme_mpa(base_size = 9) +
    theme(axis.text.x = element_text(size = 7), plot.title = element_text(size = 9, face = "bold"))
  ggsave(here::here("plots", "fig_resilience_stability.pdf"), p, width = 180, height = 80, units = "mm", device = cairo_pdf)
  ggsave(here::here("plots", "fig_resilience_stability.png"), p, width = 180, height = 80, units = "mm", dpi = 600)
}

cat("\n=== Temporal stability (CV inside vs outside) ===\n")
print(stab_tab, row.names = FALSE)
cat("\n  Tables + figure written.\n")
