# Heatwave analysis — draft manuscript text (Southern California)

*Draft for inclusion in Donham & Stier. All numbers verified from
`code/R/14_heatwave_analysis.R`. Models are glmmTMB with site + Source random
intercepts and an AR1(year|MPA) residual structure (the PRIMARY specification;
naive lme4 results are reported only to show autocorrelation inflation). Scope:
Southern California Bight only (all 15 MPAs south of Point Conception,
32.7–34.45 °N). Needs a `stier-writing-voice` polish.*

---

## Methods — Marine heatwave analysis

We tested whether the MPA effect on each focal taxon changed with the 2014–2016
marine heatwave (MHW). Each survey year was classified as **before** (≤2013),
**during** (2014–2016), or **after** (≥2017), matching Kumagai et al. (2024).
Heatwave timing/intensity came from a Hobday et al. (2016) `heatwaveR` event
catalog on daily NOAA OISST v2.1 for the Santa Barbara Coastal region (90th-
percentile threshold, ≥5-day duration, 30-yr climatology; the lab's source-of-truth
product). In this catalog the 2014–2016 MHW is the most extreme of the 1982–2026
record (two category-III "Severe" events; 189/277/72 heatwave-days in 2014/15/16
vs a pre-2014 median near zero).

Protection is encoded in our proportion-based response ratio (lnRR = ln[MPA /
Reference]), so the period term tests directly whether the MPA effect shifted with
the heatwave. Diagnostics on naive mixed models revealed strong residual temporal
autocorrelation in annual lnRR (within-MPA lag-1 ACF 0.2–0.6, largest for kelp and
purple urchin) and non-normal residuals (DHARMa KS p < 0.01); ignoring the
autocorrelation inflated significance by many orders of magnitude. We therefore
fit, per taxon, glmmTMB models with site and Source (PISCO/KFM/LTER) random
intercepts and a first-order autoregressive AR1(year | MPA) residual structure
(`lnRR ~ period + (1|MPA) + (1|Source) + AR1`), with during/after-vs-before
contrasts via `emmeans`. We tested the cascade mechanism with the same model class
(density lnRR regressed on density lnRR; kelp~urchin, urchin~predator, kelp~predator;
cf. Kumagai Fig. 7), and the spatial gradient by adding per-MPA annual MHW
cumulative intensity sampled from Kumagai et al.'s 1-km raster. Finally, because
our reserves were established over 2003–2012 while the heatwave is fixed in
calendar time, we separated the ongoing pBACIPS *recovery* trajectory from any
*heatwave-specific* effect by fitting, on post-establishment data,
`lnRR ~ time-since-establishment + exposure + (1|MPA) + (1|Source) + AR1` — a
separation Kumagai et al. could not make without a pre-establishment baseline.

## Results

**Giant kelp.** The robust result is for the foundation species. Giant kelp shifted
from roughly half as abundant inside reserves before the heatwave (RR = 0.43) to
2.5- and 3.4-fold more abundant inside during and after it (RR = 2.52, 3.38;
during- and after-vs-before ΔlnRR = +1.78 and +2.07, p < 10⁻⁴ and < 10⁻⁵ under the
AR1 model). This held when biomass and density were modelled separately and across
sensitivities — MPAs sampled in all three periods, no-take SMRs only, and dropping
the sheephead-only MPAs (after-vs-before +1.93 to +2.56, all p < 10⁻⁴). Critically,
the kelp effect was not merely ongoing recovery: controlling for time since
establishment (recovery slope +0.16 lnRR yr⁻¹, p = 0.009), kelp's inside-MPA
advantage still rose with contemporaneous heatwave exposure (+0.28 per SD,
p = 0.03), and the per-MPA spatial gradient pointed the same way (+0.23 per SD,
p = 0.06). MPAs thus delivered an acute, heatwave-specific buffer for giant kelp,
beyond the slow recovery as reserves matured.

**Predators.** Spiny lobster was more abundant inside reserves after the heatwave
(RR 0.88 → 2.05; after-vs-before +0.85, p = 0.005), but this reflected MPA recovery
rather than the heatwave per se: controlling for time since establishment, lobster's
trend was recovery-driven (recovery p = 0.003) with at most a marginal heatwave-
specific component (+0.15 per SD, p = 0.06). California sheephead were modestly more
abundant inside reserves throughout (RR ≈ 1.8) and did not change with the heatwave
(p ≥ 0.7).

**Urchins — weaker than naive models imply.** Once temporal autocorrelation was
accounted for, the urchin signal was much weaker than an autocorrelation-naive
analysis suggests. Purple urchins were marginally fewer inside reserves during the
heatwave (RR 0.97 → 0.55, during-vs-before p = 0.055) but not significantly so
afterward (p = 0.06), and red urchins showed no significant period shift (p ≥ 0.15);
neither showed a heatwave-specific effect beyond recovery (p ≥ 0.48). (The same
contrasts under a naive no-AR1 model returned p < 10⁻¹⁰, illustrating the
autocorrelation inflation.)

**Cascade mechanism — partial support.** Across MPA-years, giant-kelp lnRR was
positively associated with spiny-lobster lnRR (+0.31, p = 0.003) and negatively with
urchin lnRR (purple −0.21, p = 0.045; red −0.26, p = 0.046), consistent with a
lobster→kelp pathway. However, the intermediate link — urchin lnRR vs lobster lnRR —
was not supported under the AR1 model (−0.01, p = 0.91), and sheephead predicted
neither urchins nor kelp (p ≥ 0.45). The data therefore robustly support an
association between lobster and kelp resilience, but provide only weak support for
urchin grazing as the mediating step once autocorrelation is accounted for.

## Discussion

Our analysis corroborates the central message of Kumagai et al. (2024) — that
Southern California's no-take MPAs conferred resilience of **giant kelp** to the
2014–2016 marine heatwave — using an independent, broader dataset (PISCO + NPS KFM +
SBC LTER) and a before/after-establishment design. We add two things their
protected-vs-unprotected, heatwave-window comparison could not. First, by exploiting
staggered MPA establishment we separate slow MPA *recovery* from an acute
*heatwave-specific* effect, and show the kelp benefit is genuinely heatwave-specific,
not merely the coincidence of maturing reserves with a warm period. Second, by
explicitly modelling the strong temporal autocorrelation in annual response ratios,
we show that the predator/urchin components of the cascade are weaker and less
certain than autocorrelation-naive models (including period-mean comparisons) imply:
the robust, mechanistically interpretable signals are the lobster–kelp association
and the kelp resilience itself, whereas urchin suppression and the lobster→urchin
step are not robustly resolved in our data.

This tempers, rather than overturns, the trophic-cascade narrative. The lobster–kelp
link and the marginal urchin signals are directionally consistent with a cascade,
and our power to resolve the urchin step is limited by ~11 MPAs and high
year-to-year variability. But it cautions against strong claims of urchin-mediated
control from observational reserve data without accounting for temporal dependence.
That the cascade is carried by lobster rather than sheephead matches Kumagai (kelp
tracked lobster, not sheephead; the 2015–2016 sheephead recruitment pulse likely
decoupled sheephead abundance from grazing pressure), and the weaker red- than
purple-urchin response is consistent with partial release of red urchins from
fishing inside reserves (Malakhoff & Miller 2021).

**Caveats.** The recovery-vs-heatwave decomposition uses *contemporaneous* exposure,
so it isolates the acute effect and attributes sustained post-2016 kelp elevation to
recovery — making the kelp heatwave benefit a conservative estimate; the full
NLS-pBACIPS extraction-time machinery is not yet integrated. Residuals remain
non-normal (heavy-tailed from extreme kelp lnRR where reference kelp is near zero),
which the main analysis handles by inverse-variance weighting; a robust/heavy-tailed
likelihood is a sensible refinement. Per-MPA exposure is sampled at MPA centroids
(centroid + buffer) from Kumagai et al.'s 1-km raster; polygon-area weighting would
be a modest improvement. Analysis is restricted to the Southern California Bight; we
do not engage the Central/Northern California (sea-otter / sunflower-star) dynamics.

---

## Direct data cross-check vs Kumagai et al. (2024)

Computing the inside-vs-outside response *descriptively* on each dataset
(`table_heatwave_crosscheck_kumagai.csv`; raw Full-vs-Reference density ratios for
theirs, exponentiated paired lnRR for ours) shows the same period pattern — kelp
rising sharply (after-period inside-advantage 3.1–3.5× in both), lobster high, urchins
lower, sheephead flat — with the main difference being the pre-heatwave baseline for
lobster and kelp (higher inside in their raw-density framing, near parity in our
proportion-based, paired-reference framing). This descriptive convergence is
reassuring; the inferential conclusions above rest on the autocorrelation- and
recovery-corrected models, which are more conservative than period-mean comparisons.
