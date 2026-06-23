# =============================================================================
# 19_heatwave_replication.R
# =============================================================================
#
# PURPOSE:
#   Does the MPA-conferred resilience repeat across MORE THAN ONE marine heatwave?
#   Kumagai et al. (2024) and our script 14 evaluate a SINGLE thermal event
#   (the 2014-2016 MHW). Our series runs to 2023 and the SBC MHW catalog shows a
#   clear SECOND event in 2018-2020 ("Blob 2.0"; 122/26/83 heatwave-days) that
#   script 14 lumps into a flat "after" baseline. Here we re-code the period
#   structure to separate the two events and test whether the inside-MPA effect
#   (especially giant-kelp resilience) rises in BOTH heatwaves -- a replication
#   test that a single-event design cannot perform, and which also de-contaminates
#   the "after" baseline used in script 14.
#
# DESIGN:
#   period5 = before (<=2013) / MHW1 (2014-2016) / interim (2017, quiet) /
#             MHW2 (2018-2020) / recent (>=2021).
#   Per taxon, glmmTMB AR1(year|MPA) + site + Source random intercepts on the
#   proportion-based response ratio (lnDiff = ln[MPA / Reference]), as in script 14.
#   Key contrasts (emmeans): MHW1-before and MHW2-before (do BOTH events elevate
#   the MPA effect?), and MHW2-interim (cleanest: second heatwave vs the adjacent
#   quiet year, holding reserve maturation ~constant). Replication = the two events
#   agree in sign/significance.
#   A continuous check regresses lnDiff on annual MHW intensity across the full
#   post-establishment record (controlling for time-since-establishment), so the
#   2019 intensity dip is handled naturally and both events contribute.
#
#   SCOPE: Southern California Bight only (asserted, < 34.45 N).
#
# OUTPUTS (tables/, plots/):
#   table_heatwave_replication.csv            - per-taxon MHW1 vs MHW2 contrasts + verdict
#   table_heatwave_replication_continuous.csv - recovery-controlled MHW-intensity slope
#   fig_heatwave_replication.{pdf,png}        - per-taxon RR across the 5 periods
#
# AUTHORS: Emily Donham & Adrian Stier
# =============================================================================

cat("Running heatwave replication analysis (two MHWs; Southern California)...\n")
source(here::here("code", "R", "00_libraries.R"))
source(here::here("code", "R", "00b_color_palette.R"))
source(here::here("code", "R", "01_utils.R"))
suppressMessages({library(glmmTMB); library(emmeans); library(lme4); library(lmerTest)})

taxa <- c("Panulirus interruptus", "Semicossyphus pulcher",
          "Strongylocentrotus purpuratus", "Mesocentrotus franciscanus",
          "Macrocystis pyrifera")
role <- c("Panulirus interruptus" = "Predator", "Semicossyphus pulcher" = "Predator",
          "Strongylocentrotus purpuratus" = "Herbivore",
          "Mesocentrotus franciscanus" = "Herbivore", "Macrocystis pyrifera" = "Producer")
short <- c("Panulirus interruptus" = "P. interruptus", "Semicossyphus pulcher" = "S. pulcher",
           "Strongylocentrotus purpuratus" = "S. purpuratus",
           "Mesocentrotus franciscanus" = "M. franciscanus", "Macrocystis pyrifera" = "M. pyrifera")

# ---------------------------------------------------------------------------
# 1. Load + join + 5-level period
# ---------------------------------------------------------------------------
rr   <- read.csv(here::here("data", "harmonized", "harmonized_response_ratios.csv"), stringsAsFactors = FALSE)
hw   <- read.csv(here::here("data", "heatwave_exposure_SBC_annual.csv"), stringsAsFactors = FALSE)
meta <- read.csv(here::here("data", "harmonized", "harmonized_site_metadata.csv"), stringsAsFactors = FALSE)
stopifnot("Non-Southern-California sites (lat > 34.45 N)" = all(meta$Lat <= 34.45, na.rm = TRUE))

rr <- merge(rr, hw[, c("year", "mhw_days", "mhw_icum")], by = "year", all.x = TRUE)
rr <- subset(rr, !is.na(lnDiff) & !is.na(mhw_icum))
period5 <- function(y) factor(
  ifelse(y <= 2013, "before",
  ifelse(y <= 2016, "MHW1",
  ifelse(y == 2017, "interim",
  ifelse(y <= 2020, "MHW2", "recent")))),
  levels = c("before", "MHW1", "interim", "MHW2", "recent"))
rr$period <- period5(rr$year)
rr$yrf <- factor(rr$year)

fit_ar1 <- function(form_fixed, data) {
  nsrc <- length(unique(data$source))
  re <- if (nsrc > 1) "+ (1|source) + ar1(yrf + 0 | CA_MPA_Name_Short)" else
                      "+ ar1(yrf + 0 | CA_MPA_Name_Short)"
  f <- as.formula(paste(form_fixed, "+ (1|CA_MPA_Name_Short)", re))
  m <- tryCatch(suppressWarnings(glmmTMB(f, data = data)), error = function(e) NULL)
  if (!is.null(m) && isTRUE(m$sdr$pdHess)) m else NULL
}

# ---------------------------------------------------------------------------
# 2. Per-taxon period model + the two-event contrasts
# ---------------------------------------------------------------------------
# contrast coefficient vectors over c(before, MHW1, interim, MHW2, recent)
ctr_defs <- list("MHW1-before"  = c(-1, 1, 0, 0, 0),
                 "MHW2-before"  = c(-1, 0, 0, 1, 0),
                 "MHW2-interim" = c(0, 0, -1, 1, 0),
                 "MHW1-MHW2"    = c(0, 1, 0, -1, 0))
em_list <- list(); rep_rows <- list()
for (tx in taxa) {
  sub <- subset(rr, y == tx)
  # need all five periods present to estimate the full contrast set
  m <- fit_ar1("lnDiff ~ period", sub)
  if (is.null(m)) { m <- suppressWarnings(lmer(lnDiff ~ period + (1 | CA_MPA_Name_Short), data = sub, REML = TRUE))
                    message("  [19] AR1 failed for ", tx, "; naive LMM used.") }
  em <- emmeans(m, ~ period)
  emdf <- as.data.frame(em)
  names(emdf)[names(emdf) %in% c("asymp.LCL", "lower.CL")] <- "lower"
  names(emdf)[names(emdf) %in% c("asymp.UCL", "upper.CL")] <- "upper"
  em_list[[tx]] <- transform(emdf[, c("period", "emmean", "SE", "lower", "upper")], taxon = tx)
  present <- levels(rr$period) %in% as.character(emdf$period)
  defs_ok <- Filter(function(v) all(v[!present] == 0), ctr_defs)  # only contrasts whose cells exist
  cc <- as.data.frame(contrast(em, method = defs_ok))
  get1 <- function(nm, col) { i <- which(cc$contrast == nm); if (length(i)) cc[i, col] else NA }
  d1 <- get1("MHW1-before", "estimate"); p1 <- get1("MHW1-before", ncol(cc))
  d2 <- get1("MHW2-before", "estimate"); p2 <- get1("MHW2-before", ncol(cc))
  d2i <- get1("MHW2-interim", "estimate"); p2i <- get1("MHW2-interim", ncol(cc))
  verdict <- if (any(is.na(c(d1, d2)))) "incomplete" else
    if (sign(d1) == sign(d2) && p1 < 0.05 && p2 < 0.05) "REPEATS (both sig, same sign)" else
    if (sign(d1) == sign(d2)) "same direction, not both sig" else "diverges"
  rep_rows[[tx]] <- data.frame(taxon = tx, role = role[tx],
    MHW1_dlnRR = round(d1, 3), MHW1_p = signif(p1, 3),
    MHW2_dlnRR = round(d2, 3), MHW2_p = signif(p2, 3),
    MHW2_vs_interim_dlnRR = round(d2i, 3), MHW2_vs_interim_p = signif(p2i, 3),
    replication = verdict)
}
rep_tab <- do.call(rbind, rep_rows)
write.csv(rep_tab, here::here("tables", "table_heatwave_replication.csv"), row.names = FALSE)
emm <- do.call(rbind, em_list)

# ---------------------------------------------------------------------------
# 3. Continuous, recovery-controlled MHW-intensity effect across BOTH events
#    (post-establishment data only; controls for time-since-establishment)
# ---------------------------------------------------------------------------
cont_rows <- list()
for (tx in taxa) {
  s <- subset(rr, y == tx & BA == "After" & !is.na(time))
  s$z <- as.numeric(scale(s$mhw_icum))
  m <- fit_ar1("lnDiff ~ time + z", s)
  if (!is.null(m)) { co <- summary(m)$coefficients$cond
    cont_rows[[tx]] <- data.frame(taxon = tx, role = role[tx], n = nrow(s),
      recovery_slope = round(co["time", 1], 3), recovery_p = signif(co["time", 4], 3),
      mhw_intensity_slope = round(co["z", 1], 3), mhw_intensity_p = signif(co["z", 4], 3)) }
}
write.csv(do.call(rbind, cont_rows), here::here("tables", "table_heatwave_replication_continuous.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 4. Figure: per-taxon RR across the 5 periods, both heatwaves highlighted
# ---------------------------------------------------------------------------
emm$lab <- factor(short[emm$taxon], levels = short[taxa])
emm$period <- factor(emm$period, levels = c("before", "MHW1", "interim", "MHW2", "recent"),
                     labels = c("Before", "MHW1\n'14-16", "Interim\n'17", "MHW2\n'18-20", "Recent\n'21-23"))
emm$heat <- ifelse(grepl("MHW", emm$period), "heatwave", "non-heatwave")
yr <- range(c(emm$lower, emm$upper), na.rm = TRUE)
p <- ggplot(emm, aes(period, emmean, group = lab)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  geom_line(linewidth = 0.4, color = "grey65") +
  geom_pointrange(aes(ymin = lower, ymax = upper, color = heat), linewidth = 0.55, size = 0.32) +
  facet_wrap(~ lab, nrow = 1) +
  scale_color_manual(values = c(heatwave = "#D55E00", `non-heatwave` = "#0072B2"), name = NULL) +
  scale_y_rr(yr[1], yr[2], name = "MPA effect on abundance  (RR = MPA / Reference)") +
  labs(x = NULL, title = "Does MPA resilience repeat across two marine heatwaves? (Southern California)") +
  theme_mpa(base_size = 9) +
  theme(legend.position = "bottom", axis.text.x = element_text(size = 6),
        plot.title = element_text(size = 9, face = "bold"))
ggsave(here::here("plots", "fig_heatwave_replication.pdf"), p, width = 180, height = 90, units = "mm", device = cairo_pdf)
ggsave(here::here("plots", "fig_heatwave_replication.png"), p, width = 180, height = 90, units = "mm", dpi = 600)

cat("\n=== Heatwave replication (MHW1 2014-16 vs MHW2 2018-20) ===\n")
print(rep_tab[, c("taxon", "MHW1_dlnRR", "MHW1_p", "MHW2_dlnRR", "MHW2_p", "replication")], row.names = FALSE)
cat("\n  Tables + figure written.\n")
