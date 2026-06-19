# =============================================================================
# 14_heatwave_analysis.R
# =============================================================================
#
# PURPOSE:
#   Test whether the MPA-driven trophic cascade changed strength with the
#   2014-2016 marine heatwave (MHW) in Southern California, and provide a
#   like-for-like comparison with Kumagai et al. (2024, Global Change Biology),
#   who showed Southern California no-take MPAs buffered kelp forests against the
#   same heatwave via a predator-urchin-kelp cascade.
#
# APPROACH (our framing):
#   Protection is already encoded in our proportion-based response ratio
#   (lnRR = ln(MPA / Reference)). We classify each survey year as before
#   (<=2013), during (2014-2016), or after (>=2017) the MHW -- the same windows
#   Kumagai et al. used -- and fit, per taxon, lnRR ~ period + (1|MPA). The
#   period term tests directly whether the inside-vs-reference effect shifted
#   with the heatwave. Marginal means and during/after-vs-before contrasts via
#   emmeans.
#
#   SCOPE: Southern California Bight only. All MPAs in this study lie south of
#   Point Conception (32.7-34.45 N); the script asserts this.
#
#   Marine-heatwave timing/intensity: Hobday (2016) heatwaveR catalog on daily
#   NOAA OISST v2.1 for the Santa Barbara Coastal region (maintained in the
#   sbc-oceanography source-of-truth repo). Annual exposure is provided as a
#   tracked input: data/heatwave_exposure_SBC_annual.csv.
#
# INPUTS:
#   - data/harmonized/harmonized_response_ratios.csv
#   - data/harmonized/harmonized_site_metadata.csv  (latitude check)
#   - data/heatwave_exposure_SBC_annual.csv
#
# OUTPUTS:
#   - tables/table_heatwave_period_effects.csv
#   - tables/table_heatwave_contrasts.csv
#   - plots/fig_heatwave_cascade.pdf / .png
#
# AUTHORS: Emily Donham & Adrian Stier
# =============================================================================

cat("Running heatwave analysis (Southern California)...\n")

source(here::here("code", "R", "00_libraries.R"))
source(here::here("code", "R", "00b_color_palette.R"))
source(here::here("code", "R", "01_utils.R"))
suppressMessages({library(lme4); library(lmerTest); library(emmeans)})

# ---------------------------------------------------------------------------
# 1. Load + join
# ---------------------------------------------------------------------------
rr   <- read.csv(here::here("data", "harmonized", "harmonized_response_ratios.csv"),
                 stringsAsFactors = FALSE)
hw   <- read.csv(here::here("data", "heatwave_exposure_SBC_annual.csv"),
                 stringsAsFactors = FALSE)
meta <- read.csv(here::here("data", "harmonized", "harmonized_site_metadata.csv"),
                 stringsAsFactors = FALSE)

# Southern California assertion (no sites north of Point Conception, 34.45 N)
stopifnot("Non-Southern-California sites present (lat > 34.45 N)" =
            all(meta$Lat <= 34.45, na.rm = TRUE))

rr <- merge(rr, hw[, c("year", "mhw_days", "mhw_icum", "period")], by = "year", all.x = TRUE)
rr <- subset(rr, !is.na(period) & !is.na(lnDiff))
rr$period <- factor(rr$period, levels = c("before", "during", "after"))

# ---------------------------------------------------------------------------
# 2. Per-taxon mixed models: lnRR ~ period + (1|MPA)
# ---------------------------------------------------------------------------
taxa <- c("Panulirus interruptus", "Semicossyphus pulcher",
          "Strongylocentrotus purpuratus", "Mesocentrotus franciscanus",
          "Macrocystis pyrifera")
role <- c("Panulirus interruptus" = "Predator", "Semicossyphus pulcher" = "Predator",
          "Strongylocentrotus purpuratus" = "Herbivore",
          "Mesocentrotus franciscanus" = "Herbivore",
          "Macrocystis pyrifera" = "Producer")

em_list <- list(); ct_list <- list()
for (tx in taxa) {
  sub <- subset(rr, y == tx)
  m <- lmer(lnDiff ~ period + (1 | CA_MPA_Name_Short), data = sub, REML = TRUE)
  em <- emmeans(m, ~ period)
  em_list[[tx]] <- transform(as.data.frame(em), taxon = tx,
                             n = nrow(sub),
                             n_mpa = length(unique(sub$CA_MPA_Name_Short)))
  cc <- as.data.frame(contrast(em, method = list(
    "during-before" = c(-1, 1, 0), "after-before" = c(-1, 0, 1))))
  ct_list[[tx]] <- transform(cc, taxon = tx)
}
emm <- do.call(rbind, em_list)
ctr <- do.call(rbind, ct_list)

# ---------------------------------------------------------------------------
# 3. Tables (RR scale + contrasts)
# ---------------------------------------------------------------------------
tab <- reshape(
  data.frame(taxon = emm$taxon, period = emm$period, RR = round(exp(emm$emmean), 2)),
  idvar = "taxon", timevar = "period", direction = "wide")
names(tab) <- gsub("RR\\.", "RR_", names(tab))
tab$Role <- role[tab$taxon]
tab <- tab[match(taxa, tab$taxon), c("taxon", "Role", "RR_before", "RR_during", "RR_after")]
write.csv(tab, here::here("tables", "table_heatwave_period_effects.csv"), row.names = FALSE)

ctr_out <- data.frame(
  Taxon = ctr$taxon, Role = role[ctr$taxon], contrast = ctr$contrast,
  dlnRR = round(ctr$estimate, 3), SE = round(ctr$SE, 3),
  t = round(ctr$t.ratio, 2), p = signif(ctr$p.value, 3))
ctr_out <- ctr_out[order(match(ctr_out$Taxon, taxa), ctr_out$contrast), ]
write.csv(ctr_out, here::here("tables", "table_heatwave_contrasts.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 4. Figure (project theme + RR-scaled axis)
# ---------------------------------------------------------------------------
short <- c("Panulirus interruptus" = "P. interruptus", "Semicossyphus pulcher" = "S. pulcher",
           "Strongylocentrotus purpuratus" = "S. purpuratus",
           "Mesocentrotus franciscanus" = "M. franciscanus",
           "Macrocystis pyrifera" = "M. pyrifera")
emm$lab <- factor(short[emm$taxon], levels = short[taxa])
emm$period <- factor(emm$period, levels = c("before", "during", "after"),
                     labels = c("Before\n(≤2013)", "During\n(2014-16)", "After\n(≥2017)"))
yr <- range(c(emm$lower.CL, emm$upper.CL))

p <- ggplot(emm, aes(period, emmean, color = lab, group = lab)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  geom_line(linewidth = 0.5, alpha = 0.6) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL), linewidth = 0.6, size = 0.35) +
  facet_wrap(~ lab, nrow = 1) +
  scale_color_taxa() +
  scale_y_rr(yr[1], yr[2], name = "MPA effect on abundance  (RR = MPA / Reference)") +
  labs(x = NULL) +
  theme_mpa(base_size = 9) +
  theme(legend.position = "none", axis.text.x = element_text(size = 7))

ggsave(here::here("plots", "fig_heatwave_cascade.pdf"), p,
       width = 180, height = 80, units = "mm", device = cairo_pdf)
ggsave(here::here("plots", "fig_heatwave_cascade.png"), p,
       width = 180, height = 80, units = "mm", dpi = 600)

cat("  Heatwave analysis complete: tables/table_heatwave_period_effects.csv,",
    "tables/table_heatwave_contrasts.csv, plots/fig_heatwave_cascade.{pdf,png}\n")
