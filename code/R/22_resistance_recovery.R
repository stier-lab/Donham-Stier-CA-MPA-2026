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
BASE <- 2010:2013; DUR <- 2014:2016; AFT <- 2017:2019; REC <- 2020:2023

d <- as.data.table(read.csv(here::here("data", "harmonized", "harmonized_raw_responses.csv"), stringsAsFactors = FALSE))
# one abundance per MPA x year x taxon x status (mean over source)
d <- d[, .(value = mean(value)), by = .(CA_MPA_Name_Short, year, taxon_name, status, resp)]

win <- function(s, yrs) { v <- s$value[s$year %in% yrs]; if (length(v) >= 1) mean(v) else NA_real_ }
rr_rows <- list()
for (tx in taxa) {
  st <- d[taxon_name == tx & resp == resp_of[tx]]
  # per MPA x status windows
  agg <- st[, .(base = win(.SD, BASE), dur = win(.SD, DUR), aft = win(.SD, AFT), rec = win(.SD, REC)),
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

cat("\n=== Resistance / recovery (state-based, inside vs outside; ratio to 2010-13 baseline) ===\n")
print(rr_tab, row.names = FALSE)
cat("\n  Tables + figure written.\n")
