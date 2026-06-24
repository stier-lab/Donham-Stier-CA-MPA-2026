# =============================================================================
# 16_environmental_moderators.R
# =============================================================================
#
# PURPOSE:
#   Test whether the MPA effect on each taxon varies across the ENVIRONMENTAL
#   gradients that Kumagai et al. (2024) mapped in their PCA (Fig. 3) but never
#   entered into a statistical model. This closes the one substantive covariate
#   gap surfaced by the methods comparison (script 15): both studies model
#   protection, heatwave timing, and trophic densities, but neither tests whether
#   the reserve effect scales with thermal stress, thermal regime, biogeography,
#   or reserve size. Here we do, as an inverse-variance meta-regression on our
#   per-MPA effect sizes.
#
# MODERATORS (per-MPA; standardized to per-SD slopes):
#   mhw_during   marine-heatwave cumulative intensity during 2014-2016
#                (Kumagai's 1-km MHW raster sampled at MPA centroids;
#                 data/per_mpa_mhw_exposure.csv) -- THERMAL STRESS
#   cs_mean      cold-spell cumulative intensity, climatological mean
#                (Kumagai's 1-km cold-spell grid, nearest cell) -- UPWELLING / THERMAL REGIME
#   lat          latitude -- BIOGEOGRAPHIC / THERMAL GRADIENT
#   log_size     log MPA area (Hectares) -- RESERVE DESIGN
#
#   NOT INCLUDED (no per-MPA data without fabrication): wave exposure (Kumagai's
#   per-site `hsmax` came from a kelp-canopy NetCDF not in the mirror; only a
#   single REGIONAL SBC wave series exists locally, with no among-MPA variation),
#   nitrate/NPP, depth, human gravity. Documented as a recommended addition.
#
# MODEL: per taxon, a random-effects meta-regression rma(yi = effect size,
#   sei = SE, mods = ~ z(moderator), method = "REML", test = "knha") on the
#   density effect sizes in SumStats.Final. The Knapp-Hartung adjustment
#   (small-sample t-test on k-2 df) is ESSENTIAL here: with k = 9-15 MPAs, the
#   default normal test treats SE^2 as known and is badly anticonservative. Even
#   so these tests are EXPLORATORY -- low k, univariate moderators screened one at
#   a time, and the thermal/biogeographic moderators are mutually collinear (see
#   Known Limitations -> Moderator power in CLAUDE.md). Benjamini-Hochberg FDR
#   across all moderator tests.
#
# OUTPUTS:
#   tables/table_s_mpa_env_covariates.csv  - the per-MPA environmental covariate table
#   tables/table_s_env_moderators.csv      - taxon x moderator: slope/CI/p/FDR/k
#   plots/fig_s_env_moderators.{pdf,png}   - per-SD moderator slopes (forest, by taxon)
#
# AUTHORS: Emily Donham & Adrian Stier
# =============================================================================

cat("Running environmental moderator meta-regression (script 16)...\n")
source(here::here("code", "R", "00_libraries.R"))
source(here::here("code", "R", "00b_color_palette.R"))
source(here::here("code", "R", "01_utils.R"))
suppressMessages(library(metafor))

taxa <- c("P. interruptus", "S. pulcher", "S. purpuratus", "M. franciscanus", "M. pyrifera")
role <- c("P. interruptus" = "Predator", "S. pulcher" = "Predator", "S. purpuratus" = "Herbivore",
          "M. franciscanus" = "Herbivore", "M. pyrifera" = "Producer")
safe <- function(expr) tryCatch(suppressWarnings(suppressMessages(expr)), error = function(e) NULL)

# ---------------------------------------------------------------------------
# 1. Per-MPA environmental covariate table
# ---------------------------------------------------------------------------
# effect sizes (per MPA x taxon) -- carries MPA, Source, Lat, Lon.
# Response: density for the animals (the currency comparable to Kumagai), but BIOMASS
# for giant kelp -- our pBACIPS pipeline computes kelp effect sizes on biomass only
# (SumStats.Final has 0 density effects for M. pyrifera). Each taxon's moderator
# meta-regression is within-taxon, so the response type need only be internally
# consistent; the `response` column records which was used.
ss <- read.csv(here::here("data", "sumstats_final.csv"), stringsAsFactors = FALSE)
ss <- subset(ss, Taxa %in% taxa & is.finite(Mean) & is.finite(SE) & SE > 0)
ss <- subset(ss, (Taxa != "M. pyrifera" & Resp == "Den") | (Taxa == "M. pyrifera" & Resp == "Bio"))
ss <- ss[!duplicated(ss[, c("MPA", "Taxa")]), ]   # one effect per MPA x taxon

# MPA size (Hectares) from response ratios
rr <- read.csv(here::here("data", "harmonized", "harmonized_response_ratios.csv"), stringsAsFactors = FALSE)
size <- unique(rr[, c("CA_MPA_Name_Short", "Hectares")]); names(size)[1] <- "MPA"

# MHW thermal stress during the heatwave (per-MPA centroid sample of Kumagai raster)
mhw <- read.csv(here::here("data", "per_mpa_mhw_exposure.csv"), stringsAsFactors = FALSE)
mhw_during <- aggregate(mhw_icum_mpa ~ CA_MPA_Name_Short,
                        data = subset(mhw, year >= 2014 & year <= 2016), FUN = mean)
names(mhw_during) <- c("MPA", "mhw_during")

# Cold-spell climatology (per-MPA nearest cell of Kumagai's 1-km cold-spell grid)
cs_path <- "~/kumagai2024-comparison/repo/Processed_data/SST/CS_cummulative_intensity_1km.rds"
mpa_ll <- unique(ss[, c("MPA", "Lat", "Lon")])
if (file.exists(path.expand(cs_path))) {
  cs <- readRDS(path.expand(cs_path))
  csagg <- aggregate(CS_cummulative ~ long + lat, data = cs, FUN = mean)  # climatological mean per cell
  nearest_cs <- function(la, lo) {
    d <- (csagg$lat - la)^2 + (csagg$long - lo)^2
    csagg$CS_cummulative[which.min(d)]
  }
  mpa_ll$cs_mean <- mapply(nearest_cs, mpa_ll$Lat, mpa_ll$Lon)
} else {
  message("  [16] cold-spell grid not found; cs_mean omitted.")
  mpa_ll$cs_mean <- NA_real_
}

# wave exposure: per-MPA mean HSMAX (max significant wave height) from Bell (2023)
# island-inclusive kelp-canopy-env product (EDI knb-lter-sbc.162; CDIP MOP v1.1 waves),
# derived by code/R/extract_wave_exposure.R -> data/per_mpa_wave_exposure.csv (all 34 MPAs)
wave_path <- here::here("data", "per_mpa_wave_exposure.csv")
wave_dt <- if (file.exists(wave_path)) {
  w <- read.csv(wave_path, stringsAsFactors = FALSE); w[, c("MPA", "wave_hs")]
} else NULL

env <- Reduce(function(a, b) merge(a, b, by = "MPA", all.x = TRUE),
              c(list(mpa_ll, mhw_during, size), if (!is.null(wave_dt)) list(wave_dt)))
env$log_size <- log(env$Hectares)
keep_cols <- c("MPA", "Lat", "Lon", "mhw_during", "cs_mean", "log_size", "Hectares")
if (!is.null(wave_dt)) keep_cols <- c(keep_cols, "wave_hs")
env <- env[, keep_cols]
names(env)[names(env) == "Lat"] <- "lat"
write.csv(env, here::here("tables", "table_s_mpa_env_covariates.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 2. Meta-regression: per taxon x moderator (inverse-variance, Source RE)
# ---------------------------------------------------------------------------
dat <- merge(ss[, c("MPA", "Taxa", "Resp", "Mean", "SE", "Source")], env, by = "MPA")
moderators <- c(mhw_during = "MHW intensity (during)", cs_mean = "Cold-spell intensity",
                lat = "Latitude", log_size = "log MPA size")
if (!is.null(wave_dt)) moderators <- c(moderators, wave_hs = "Wave exposure (HSMAX)")

rows <- list()
for (tx in taxa) {
  d0 <- subset(dat, Taxa == tx)
  for (mk in names(moderators)) {
    d <- d0[is.finite(d0[[mk]]), ]
    if (nrow(d) < 5 || length(unique(d[[mk]])) < 4) next
    d$z <- as.numeric(scale(d[[mk]]))
    # random-effects meta-regression with Knapp-Hartung small-sample t-test
    m <- safe(rma(yi = Mean, sei = SE, mods = ~ z, data = d, method = "REML", test = "knha"))
    if (is.null(m)) next
    rows[[paste(tx, mk)]] <- data.frame(
      taxon = tx, role = role[tx], response = d$Resp[1], moderator = moderators[mk], k = nrow(d),
      slope_perSD = as.numeric(m$b[2]), se = m$se[2],
      ci_lo = m$ci.lb[2], ci_hi = m$ci.ub[2], p = m$pval[2])
  }
}
mod_tab <- do.call(rbind, rows)
mod_tab$FDR <- p.adjust(mod_tab$p, method = "BH")
num <- c("slope_perSD", "se", "ci_lo", "ci_hi", "p", "FDR")
mod_out <- mod_tab; mod_out[num] <- lapply(mod_out[num], function(x) signif(x, 3))
write.csv(mod_out, here::here("tables", "table_s_env_moderators.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 3. Figure: per-SD moderator slopes (lnRR scale), forest by taxon, faceted by moderator
# ---------------------------------------------------------------------------
mt <- mod_tab
mt$taxon <- factor(mt$taxon, levels = rev(taxa))                 # trophic order top-to-bottom
mt$moderator <- factor(mt$moderator, levels = moderators)
mt$sig <- ifelse(mt$FDR < 0.05, "FDR < 0.05", ifelse(mt$p < 0.05, "p < 0.05", "ns"))
p <- ggplot(mt, aes(slope_perSD, taxon, colour = taxon)) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey40") +
  geom_pointrange(aes(xmin = ci_lo, xmax = ci_hi, shape = sig), linewidth = 0.5, size = 0.35) +
  facet_wrap(~ moderator, nrow = 1) +
  scale_colour_manual(values = setNames(col_taxa, taxa), guide = "none") +
  scale_shape_manual(values = c("FDR < 0.05" = 16, "p < 0.05" = 17, "ns" = 1), name = NULL) +
  labs(x = "Change in MPA effect (lnRR) per SD of moderator", y = NULL,
       title = "MPA effect vs environmental gradients (Kumagai PCA covariates, here modeled)") +
  theme_mpa(base_size = 9) +
  theme(legend.position = "bottom", plot.title = element_text(size = 9, face = "bold"),
        axis.text.y = element_text(face = "italic", size = 7),
        strip.text = element_text(size = 7))
ggsave(here::here("plots", "fig_s_env_moderators.pdf"), p, width = 180, height = 80, units = "mm", device = cairo_pdf)
ggsave(here::here("plots", "fig_s_env_moderators.png"), p, width = 180, height = 80, units = "mm", dpi = 600)

cat("  Environmental moderators complete:", nrow(mod_tab), "taxon x moderator tests;",
    sum(mod_tab$p < 0.05), "with p<0.05,", sum(mod_tab$FDR < 0.05), "surviving FDR.\n")
cat("  Tables: table_s_mpa_env_covariates.csv, table_s_env_moderators.csv; figure: fig_s_env_moderators.\n")
