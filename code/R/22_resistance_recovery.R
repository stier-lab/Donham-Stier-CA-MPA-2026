# =============================================================================
# 22_resistance_recovery.R
# =============================================================================
#
# PURPOSE:
#   Decompose resilience into RESISTANCE and RECOVERY on the STATE variable
#   (abundance), inside vs outside MPAs -- the framework Kumagai et al. (2024) used
#   (resistance/recovery permutation tests). Our scripts 14/19 work on the lnRR
#   (a ratio), which hides whether a rising MPA effect means "kelp held inside" or
#   "kelp crashed outside". Working on the absolute trajectory distinguishes these
#   and yields communicable quantities ("kelp retained X% inside vs Y% outside").
#
# METRICS (per MPA x status, relative to a pre-heatwave baseline 2010-2013):
#   resistance      = during(2014-16) / baseline   (1 = fully resisted; <1 = declined)
#   recovery        = after(2017-19)  / baseline   (1 = returned to baseline)
#   recovery_recent = recent(2020-23) / baseline   (completeness; <1 = incomplete / possible state shift)
#   We compare each inside (MPA) vs the paired reference (paired Wilcoxon on log-ratio)
#   and report the absolute trajectory.
#
#   SCOPE: Southern California Bight. Giant kelp (biomass), purple/red urchin,
#   lobster, sheephead (density).
#
# OUTPUTS (tables/, plots/):
#   table_resistance_recovery.csv   - per taxon: resistance/recovery inside vs outside, paired p
#   fig_resistance_recovery.{pdf,png} - abundance trajectories inside vs outside (kelp + urchins)
#
# AUTHORS: Emily Donham & Adrian Stier
# =============================================================================

cat("Running resistance/recovery decomposition (script 22)...\n")
source(here::here("code", "R", "00_libraries.R"))
source(here::here("code", "R", "00b_color_palette.R"))
source(here::here("code", "R", "00c_analysis_constants.R"))
source(here::here("code", "R", "01_utils.R"))
suppressMessages(library(data.table))

# Shared resilience constants/helpers live in 00c_analysis_constants.R + 01_utils.R
taxa <- RESILIENCE_TAXA
short <- RESILIENCE_TAXA_SHORT
resp_of <- RESILIENCE_RESP_OF
BASE <- RESILIENCE_BASELINE_YEARS; DUR <- MHW1_YEARS; AFT <- RESILIENCE_AFTER_YEARS; REC <- RESILIENCE_RECENT_YEARS

d <- as.data.table(load_harmonized_raw())
# one abundance per MPA x year x taxon x status (mean over source)
d <- d[, .(value = mean(value)), by = .(CA_MPA_Name_Short, year, taxon_name, status, resp)]

rr_rows <- list()
for (tx in taxa) {
  st <- d[taxon_name == tx & resp == resp_of[tx]]
  # per MPA x status windows (window_mean from 01_utils replaces the old win() helper)
  agg <- st[, .(base = window_mean(value, year, BASE), dur = window_mean(value, year, DUR),
                aft = window_mean(value, year, AFT), rec = window_mean(value, year, REC)),
            by = .(CA_MPA_Name_Short, status)]
  agg <- agg[is.finite(base) & base > 0]
  agg[, `:=`(resistance = dur / base, recovery = aft / base, recovery_recent = rec / base)]
  w <- dcast(agg, CA_MPA_Name_Short ~ status, value.var = c("resistance", "recovery", "recovery_recent"))
  paircmp <- function(metric) {
    a <- w[[paste0(metric, "_mpa")]]; b <- w[[paste0(metric, "_reference")]]
    ok <- is.finite(a) & is.finite(b) & a > 0 & b > 0
    if (sum(ok) < 5) return(c(inside = NA, outside = NA, p = NA, n = sum(ok)))
    p <- suppressWarnings(wilcox.test(log(a[ok]), log(b[ok]), paired = TRUE)$p.value)
    c(inside = median(a[ok]), outside = median(b[ok]), p = p, n = sum(ok))
  }
  for (mt in c("resistance", "recovery", "recovery_recent")) {
    r <- paircmp(mt)
    rr_rows[[paste(tx, mt)]] <- data.frame(taxon = tx, metric = mt, n_mpa = unname(r["n"]),
      inside = round(unname(r["inside"]), 2), outside = round(unname(r["outside"]), 2),
      paired_p = signif(unname(r["p"]), 3))
  }
}
rr_tab <- do.call(rbind, rr_rows)
write.csv(rr_tab, here::here("tables", "table_resistance_recovery.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# Figure: abundance trajectories inside vs outside (kelp + purple urchin)
# ---------------------------------------------------------------------------
traj <- d[taxon_name %in% c("Macrocystis pyrifera", "Strongylocentrotus purpuratus")]
traj <- traj[(taxon_name == "Macrocystis pyrifera" & resp == "Bio") |
             (taxon_name == "Strongylocentrotus purpuratus" & resp == "Den")]
traj <- traj[year >= 2008 & year <= 2023]
tr <- traj[, .(abund = mean(value), se = sd(value) / sqrt(.N)), by = .(taxon_name, year, status)]
tr$lab <- factor(short[tr$taxon_name], levels = short[c("Macrocystis pyrifera", "Strongylocentrotus purpuratus")])
tr$Site <- factor(tr$status, levels = c("reference", "mpa"), labels = c("Reference", "MPA"))
p <- ggplot(tr, aes(year, abund, color = Site, fill = Site)) +
  annotate("rect", xmin = 2013.5, xmax = 2016.5, ymin = -Inf, ymax = Inf, alpha = 0.10, fill = "#D55E00") +
  geom_ribbon(aes(ymin = abund - se, ymax = abund + se), alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 1) +
  facet_wrap(~ lab, nrow = 1, scales = "free_y") +
  scale_color_manual(values = c(Reference = "#D55E00", MPA = "#0072B2"), name = NULL, aesthetics = c("color", "fill")) +
  labs(x = NULL, y = "Abundance (kelp biomass / urchin density)",
       title = "Resistance & recovery: kelp and urchin state trajectories inside vs outside MPAs") +
  theme_mpa(base_size = 9) +
  theme(legend.position = "bottom", plot.title = element_text(size = 9, face = "bold"))
ggsave(here::here("plots", "fig_resistance_recovery.pdf"), p, width = 170, height = 85, units = "mm", device = cairo_pdf)
ggsave(here::here("plots", "fig_resistance_recovery.png"), p, width = 170, height = 85, units = "mm", dpi = 600)

# Main-text single-panel: giant-kelp resistance & recovery (Conservation Letters, single column 80 mm)
trk <- tr[taxon_name == "Macrocystis pyrifera"]
lab_pos <- trk[, .SD[which.max(year)], by = Site]
p_kelp <- ggplot(trk, aes(year, abund, color = Site, fill = Site)) +
  annotate("rect", xmin = 2013.5, xmax = 2016.5, ymin = -Inf, ymax = Inf, alpha = 0.10, fill = "#D55E00") +
  annotate("text", x = 2015, y = max(trk$abund + trk$se), label = "2014-16\nMHW",
           size = 2.3, color = "grey35", lineheight = 0.85, vjust = 1) +
  geom_ribbon(aes(ymin = abund - se, ymax = abund + se), alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 0.9) +
  geom_text(data = lab_pos, aes(label = Site), hjust = 0, nudge_x = 0.3, size = 2.7, fontface = "bold") +
  scale_color_manual(values = c(Reference = "#D55E00", MPA = "#0072B2"), aesthetics = c("color", "fill"), guide = "none") +
  scale_x_continuous(limits = c(2008, 2028.5), breaks = seq(2008, 2024, 4)) +
  labs(x = NULL, y = "Giant kelp biomass") +
  theme_mpa(base_size = 8) + theme(plot.margin = margin(6, 10, 4, 4))
ggsave(here::here("plots", "fig_kelp_resilience.pdf"), p_kelp, width = 80, height = 72, units = "mm", device = cairo_pdf)
ggsave(here::here("plots", "fig_kelp_resilience.png"), p_kelp, width = 80, height = 72, units = "mm", dpi = 600)

# Per-reserve PAIRED view (each reserve vs its own control) -- shows the resistance/recovery
# pattern is consistent across pairs, not a pooled-mean artifact. Ratios to 2010-13 baseline.
kb <- d[taxon_name == "Macrocystis pyrifera" & resp == "Bio", .(value = mean(value)),
        by = .(CA_MPA_Name_Short, status, year)]
agb <- kb[, .(base = window_mean(value, year, BASE), dur = window_mean(value, year, DUR),
              rec = window_mean(value, year, REC)),
          by = .(CA_MPA_Name_Short, status)][is.finite(base) & base > 0]
wb <- dcast(agb, CA_MPA_Name_Short ~ status, value.var = c("base", "dur", "rec"))
wb <- wb[is.finite(base_mpa) & is.finite(base_reference)]
wb[, `:=`(Res_in = dur_mpa / base_mpa, Res_out = dur_reference / base_reference,
          Rec_in = rec_mpa / base_mpa, Rec_out = rec_reference / base_reference)]
ordb <- wb[order(Res_in), CA_MPA_Name_Short]
segb <- rbind(wb[, .(reserve = CA_MPA_Name_Short, metric = "Resistance (during)", out = Res_out, ins = Res_in)],
              wb[, .(reserve = CA_MPA_Name_Short, metric = "Recovery (2020-23)", out = Rec_out, ins = Rec_in)])
segb <- segb[is.finite(out) & is.finite(ins) & out > 0 & ins > 0]
segb$reserve <- factor(segb$reserve, levels = ordb)
segb$metric  <- factor(segb$metric, levels = c("Resistance (during)", "Recovery (2020-23)"))
ptsb <- rbind(segb[, .(reserve, metric, Site = "Reference", value = out)],
              segb[, .(reserve, metric, Site = "MPA", value = ins)])
ptsb$Site <- factor(ptsb$Site, levels = c("Reference", "MPA"))
p_pair <- ggplot() +
  geom_vline(xintercept = 1, linetype = "dotted", color = "grey45") +
  geom_segment(data = segb, aes(x = out, xend = ins, y = reserve, yend = reserve), color = "grey75", linewidth = 0.5) +
  geom_point(data = ptsb, aes(value, reserve, color = Site), size = 1.8) +
  facet_wrap(~ metric) + scale_x_log10() +
  scale_color_manual(values = c(Reference = "#D55E00", MPA = "#0072B2"), name = NULL) +
  labs(x = "Kelp biomass relative to 2010-13 baseline (log scale; >1 = above baseline)", y = NULL,
       title = "Giant-kelp resilience in each reserve vs its paired control") +
  theme_mpa(base_size = 8) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 6.5), plot.title = element_text(size = 9, face = "bold"))
ggsave(here::here("plots", "fig_kelp_resilience_paired.pdf"), p_pair, width = 150, height = 82, units = "mm", device = cairo_pdf)
ggsave(here::here("plots", "fig_kelp_resilience_paired.png"), p_pair, width = 150, height = 82, units = "mm", dpi = 600)

cat("\n=== Resistance / recovery (state-based, inside vs outside; ratio to 2010-13 baseline) ===\n")
print(rr_tab, row.names = FALSE)
cat("\n  Tables + figure written.\n")
