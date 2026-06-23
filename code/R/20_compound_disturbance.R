# =============================================================================
# 20_compound_disturbance.R
# =============================================================================
#
# PURPOSE:
#   Evaluate MPA resilience to a COMPOUND disturbance -- the 2013-14 sea-star
#   wasting disease (SSWD) immediately followed by the 2014-16 marine heatwave
#   (MHW) -- as a complement to the heat-only framing of Kumagai (2024) / script 14
#   and the disease-focused Eisaguirre (2020). This is a DIFFERENT STRESSOR TYPE
#   (biotic disease / keystone-predator removal) layered on the thermal event.
#
# THE HONEST FRAMING:
#   In the Southern California Bight, SSWD (2013-14) and the MHW (2014-16) are
#   temporally inseparable, so this is a *compound* stressor test, NOT an
#   SSWD-specific one. The new empirical content our 5-taxon harmonized data lacks
#   is the SEA-STAR LOSS itself, which we document from the raw PISCO swath.
#
#   Logical key (established here): the sea star crashed EQUALLY inside and outside
#   reserves (disease is not stopped by protection), so any MPA buffering of the
#   urchin->kelp cascade must come from the PROTECTED fished predators (lobster,
#   sheephead), not from protecting the sea star -- the compensation hypothesis.
#
# ANALYSES:
#   (A) Document the disturbance: sea-star (Pycnopodia + Pisaster) density crash in
#       SoCal, inside vs outside, 2008-2022 (raw PISCO swath).
#   (B) Cascade reorganization: did the inside-outside gap for urchins and kelp
#       reorganize across the compound event? log(raw) ~ status*period + (1|MPA),
#       period = before(<=2013) / compound(2014-16) / after(>=2017).
#   (C) Resilience (lnRR): MPA effect on urchin suppression & kelp across the
#       compound event (glmmTMB AR1, as script 14).
#   (D) Cross-study synthesis: Eisaguirre (SSWD) / Kumagai (MHW) / ours (compound).
#
#   SCOPE: Southern California Bight (lat <= 34.45 N).
#
# OUTPUTS (tables/, plots/):
#   table_compound_seastar_crash.csv       - sea-star density by status x year
#   table_compound_cascade_reorg.csv       - status x period interactions (urchin, kelp)
#   table_compound_resilience.csv          - lnRR MPA effect by period (AR1)
#   table_compound_crossstudy.csv          - 3-study synthesis
#   fig_compound_disturbance.{pdf,png}     - sea-star crash + cascade response
#
# NOTE: standalone (reads raw PISCO swath from the sibling data repo for sea stars).
#
# AUTHORS: Emily Donham & Adrian Stier
# =============================================================================

cat("Running compound-disturbance analysis (SSWD + MHW; Southern California)...\n")
source(here::here("code", "R", "00_libraries.R"))
source(here::here("code", "R", "00b_color_palette.R"))
source(here::here("code", "R", "01_utils.R"))
suppressMessages({library(glmmTMB); library(emmeans); library(lme4); library(lmerTest); library(data.table)})
safe <- function(e) tryCatch(suppressWarnings(suppressMessages(e)), error = function(x) NULL)
period3 <- function(y) factor(ifelse(y <= 2013, "before", ifelse(y <= 2016, "compound", "after")),
                              levels = c("before", "compound", "after"))

# ---------------------------------------------------------------------------
# (A) Document the SSWD: sea-star crash in SoCal, inside vs outside
# ---------------------------------------------------------------------------
DR <- path.expand("~/Donham-Stier-CA-MPA-Data-2026/data/PISCO")
have_swath <- dir.exists(DR)
if (have_swath) {
  st <- fread(file.path(DR, "master_site_table_Emilyedit.csv"))
  soc <- unique(st[latitude <= 34.45 & site_status %in% c("mpa", "reference"),
                   .(site, status = site_status)])[!duplicated(site)]
  sw <- fread(file.path(DR, "MLPA_kelpforest_swath_2024.csv"), nThread = 2,
              select = c("site", "year", "zone", "transect", "classcode", "count"))
  sw <- merge(sw, soc, by = "site")[year %in% 2008:2022]
  sw$tid <- paste(sw$site, sw$year, sw$zone, sw$transect)
  star <- sw[classcode %in% c("PYCHEL", "PISGIG"), .(star = sum(count)), by = .(tid, status, year)]
  tr <- unique(sw[, .(tid, status, year)])
  tr <- merge(tr, star, by = c("tid", "status", "year"), all.x = TRUE); tr$star[is.na(tr$star)] <- 0
  crash <- tr[, .(star_density = round(mean(star), 2), n = .N), by = .(year, status)]
  crash <- crash[order(year, status)]
  write.csv(crash, here::here("tables", "table_compound_seastar_crash.csv"), row.names = FALSE)
  # drop magnitude pre vs post, by status
  pre  <- tr[year <= 2013, .(pre = mean(star)), by = status]
  post <- tr[year >= 2014, .(post = mean(star)), by = status]
  drop <- merge(pre, post, by = "status"); drop$pct_drop <- round(100 * (1 - drop$post / drop$pre), 1)
} else message("  [20] raw PISCO swath not found; sea-star section skipped.")

# ---------------------------------------------------------------------------
# (B) Cascade reorganization across the compound event (harmonized raw responses)
#     log(raw) ~ status*period + (1|MPA): does the inside-outside gap reorganize?
# ---------------------------------------------------------------------------
oraw <- as.data.table(read.csv(here::here("data", "harmonized", "harmonized_raw_responses.csv"), stringsAsFactors = FALSE))
oraw <- oraw[year >= 2002]
oraw$period <- period3(oraw$year)
oraw$status <- factor(oraw$status, levels = c("reference", "mpa"))
log_eps <- function(x) { nz <- x[x > 0 & is.finite(x)]; log(x + (if (length(nz)) min(nz) / 2 else 1e-6)) }
reorg_one <- function(taxon, rsp) {
  d <- oraw[taxon_name == taxon & resp == rsp]
  if (nrow(d) < 30) return(NULL)
  d$lv <- log_eps(d$value)
  m <- safe(lmer(lv ~ status * period + (1 | CA_MPA_Name_Short), data = d, REML = TRUE))
  if (is.null(m)) return(NULL)
  cc <- coef(summary(m))
  g <- function(r) if (r %in% rownames(cc)) cc[r, c("Estimate", "Pr(>|t|)")] else c(NA, NA)
  data.frame(taxon = taxon, resp = rsp,
             inside_before = round(g("statusmpa")[1], 3),       p_before = signif(g("statusmpa")[2], 3),
             reorg_compound = round(g("statusmpa:periodcompound")[1], 3), p_compound = signif(g("statusmpa:periodcompound")[2], 3),
             reorg_after = round(g("statusmpa:periodafter")[1], 3),       p_after = signif(g("statusmpa:periodafter")[2], 3))
}
reorg <- do.call(rbind, list(
  reorg_one("Strongylocentrotus purpuratus", "Den"),
  reorg_one("Mesocentrotus franciscanus", "Den"),
  reorg_one("Macrocystis pyrifera", "Bio"),
  reorg_one("Panulirus interruptus", "Den"),
  reorg_one("Semicossyphus pulcher", "Den")))
write.csv(reorg, here::here("tables", "table_compound_cascade_reorg.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# (C) Resilience: MPA effect (lnRR) by period, AR1 (as script 14)
# ---------------------------------------------------------------------------
rr <- read.csv(here::here("data", "harmonized", "harmonized_response_ratios.csv"), stringsAsFactors = FALSE)
rr <- subset(rr, !is.na(lnDiff) & year >= 2002); rr$period <- period3(rr$year); rr$yrf <- factor(rr$year)
fit_ar1 <- function(form_fixed, data) {
  nsrc <- length(unique(data$source))
  re <- if (nsrc > 1) "+ (1|source) + ar1(yrf + 0 | CA_MPA_Name_Short)" else "+ ar1(yrf + 0 | CA_MPA_Name_Short)"
  m <- safe(glmmTMB(as.formula(paste(form_fixed, "+ (1|CA_MPA_Name_Short)", re)), data = data))
  if (!is.null(m) && isTRUE(m$sdr$pdHess)) m else NULL
}
taxa <- c("Panulirus interruptus", "Semicossyphus pulcher", "Strongylocentrotus purpuratus",
          "Mesocentrotus franciscanus", "Macrocystis pyrifera")
short <- c("Panulirus interruptus" = "P. interruptus", "Semicossyphus pulcher" = "S. pulcher",
           "Strongylocentrotus purpuratus" = "S. purpuratus",
           "Mesocentrotus franciscanus" = "M. franciscanus", "Macrocystis pyrifera" = "M. pyrifera")
res_rows <- list(); em_list <- list()
for (tx in taxa) {
  sub <- subset(rr, y == tx)
  m <- fit_ar1("lnDiff ~ period", sub)
  if (is.null(m)) m <- suppressWarnings(lmer(lnDiff ~ period + (1 | CA_MPA_Name_Short), data = sub, REML = TRUE))
  em <- emmeans(m, ~ period); emdf <- as.data.frame(em)
  names(emdf)[names(emdf) %in% c("asymp.LCL", "lower.CL")] <- "lower"
  names(emdf)[names(emdf) %in% c("asymp.UCL", "upper.CL")] <- "upper"
  em_list[[tx]] <- transform(emdf[, c("period", "emmean", "SE", "lower", "upper")], taxon = tx)
  cc <- as.data.frame(contrast(em, method = list("compound-before" = c(-1, 1, 0), "after-before" = c(-1, 0, 1))))
  res_rows[[tx]] <- data.frame(taxon = tx,
    compound_before = round(cc$estimate[1], 3), p_compound = signif(cc[[ncol(cc)]][1], 3),
    after_before = round(cc$estimate[2], 3), p_after = signif(cc[[ncol(cc)]][2], 3))
}
write.csv(do.call(rbind, res_rows), here::here("tables", "table_compound_resilience.csv"), row.names = FALSE)
emm <- do.call(rbind, em_list)

# ---------------------------------------------------------------------------
# (D) Cross-study synthesis
# ---------------------------------------------------------------------------
cross <- data.frame(
  study = c("Eisaguirre et al. 2020", "Kumagai et al. 2024", "This study (compound)"),
  region = c("N. Channel Islands", "S. + Central California", "S. California Bight"),
  disturbance_framing = c("SSWD (sea-star loss)", "2014-16 marine heatwave",
                          "compound SSWD + MHW (inseparable in SoCal)"),
  design = c("paired inside/outside, pre/post-SSWD", "pooled protected-vs-unprotected, heatwave windows",
             "pBACIPS paired lnRR + AR1, before/compound/after"),
  kelp_conclusion = c("MPA predators buffer kelp after sea-star loss",
                      "MPAs preserve cascade -> kelp resilience to MHW",
                      "kelp resilient to the compound stressor (repeats across two MHWs, script 19)"),
  sea_star = c("Pycnopodia, key NCI urchin predator", "cited as SSWD context (not modeled)",
               "Pisaster + Pycnopodia crash documented, equal inside/outside (script 20)"))
write.csv(cross, here::here("tables", "table_compound_crossstudy.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# Figure: (a) sea-star crash inside/outside; (b) cascade MPA effect by period
# ---------------------------------------------------------------------------
plots <- list()
if (have_swath) {
  cr <- crash; cr$Status <- factor(cr$status, levels = c("reference", "mpa"), labels = c("Reference", "MPA"))
  plots$a <- ggplot(cr, aes(year, star_density, color = Status)) +
    annotate("rect", xmin = 2013.5, xmax = 2016.5, ymin = -Inf, ymax = Inf, alpha = 0.12, fill = "#D55E00") +
    geom_line(linewidth = 0.7) + geom_point(size = 1.2) +
    scale_color_manual(values = c(Reference = "#D55E00", MPA = "#0072B2"), name = NULL) +
    labs(x = NULL, y = expression("Sea-star density (per 60 m"^2*")"),
         title = "a  Sea-star wasting crash (Pycnopodia + Pisaster), equal inside & outside") +
    theme_mpa(base_size = 9) + theme(legend.position = "bottom", plot.title = element_text(size = 8.5, face = "bold"))
}
emm$lab <- factor(short[emm$taxon], levels = short[taxa])
emm$period <- factor(emm$period, levels = c("before", "compound", "after"),
                     labels = c("Before\n<=2013", "Compound\n'14-16", "After\n>=2017"))
yr <- range(c(emm$lower, emm$upper), na.rm = TRUE)
plots$b <- ggplot(emm, aes(period, emmean, group = lab, color = lab)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  geom_line(linewidth = 0.5, alpha = 0.6) +
  geom_pointrange(aes(ymin = lower, ymax = upper), linewidth = 0.5, size = 0.3) +
  facet_wrap(~ lab, nrow = 1) + scale_color_taxa() +
  scale_y_rr(yr[1], yr[2], name = "MPA effect (RR = MPA / Reference)") +
  labs(x = NULL, title = "b  Cascade response to the compound disturbance") +
  theme_mpa(base_size = 9) +
  theme(legend.position = "none", axis.text.x = element_text(size = 6),
        plot.title = element_text(size = 8.5, face = "bold"))
fig <- if (!is.null(plots$a)) patchwork::wrap_plots(plots$a, plots$b, ncol = 1, heights = c(1, 1.1)) else plots$b
ggsave(here::here("plots", "fig_compound_disturbance.pdf"), fig,
       width = 180, height = if (!is.null(plots$a)) 175 else 90, units = "mm", device = cairo_pdf)
ggsave(here::here("plots", "fig_compound_disturbance.png"), fig,
       width = 180, height = if (!is.null(plots$a)) 175 else 90, units = "mm", dpi = 600)

cat("\n=== Sea-star crash (% drop pre->post 2014, by status) ===\n"); if (have_swath) print(drop)
cat("\n=== Cascade resilience (lnRR MPA effect by period, AR1) ===\n")
print(do.call(rbind, res_rows), row.names = FALSE)
cat("\n  Tables + figure written.\n")
