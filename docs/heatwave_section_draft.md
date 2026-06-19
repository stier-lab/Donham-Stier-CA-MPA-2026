# Heatwave analysis — draft manuscript text (Southern California)

*Draft for inclusion in Donham & Stier. All numbers verified from
`code/R/14_heatwave_analysis.R` on harmonized response ratios joined to the SBC
LTER marine-heatwave catalog. Scope: Southern California Bight only (all 15 MPAs
south of Point Conception, 32.7–34.45 °N). Needs a `stier-writing-voice` polish.*

---

## Methods — Marine heatwave analysis

We tested whether the strength of the MPA-driven trophic cascade changed with the
2014–2016 marine heatwave (MHW). Each survey year was classified as **before**
(≤2013), **during** (2014–2016), or **after** (≥2017) the MHW, the windows used by
Kumagai et al. (2024) for the same region. Heatwave timing and intensity came from
a Hobday et al. (2016) event catalog built with the `heatwaveR` package on the
daily NOAA OISST v2.1 record for the Santa Barbara Coastal region (90th-percentile
threshold, ≥5-day duration, 30-year climatology), the lab's source-of-truth
physical-environment product. In this catalog the 2014–2016 MHW is the most
extreme of the 1982–2026 record (two category-III "Severe" events, Oct 2014–Apr
2015 [206 d] and Jul–Nov 2015 [129 d]; 189/277/72 heatwave-days in 2014/2015/2016
vs a pre-2014 median near zero).

Because protection is already encoded in our proportion-based response ratio
(lnRR = ln[MPA / Reference]), the period term tests directly whether the MPA effect
shifted with the heatwave. For each focal taxon we fit `lnRR ~ period + (1|MPA)`
(linear mixed model, `lme4`), estimating marginal means and during/after-vs-before
contrasts with `emmeans`. We confirmed results held when biomass and density
response ratios were modelled separately. To test the cascade mechanism (cf.
Kumagai et al. 2024, their Fig. 7), we regressed density lnRR on density lnRR across
all MPA-years with a site random intercept (kelp~urchin, urchin~predator,
kelp~predator). Finally, to ask whether the MPA effect scales with *local*
heatwave intensity (not only with calendar period), we extracted per-MPA annual
MHW cumulative intensity (centroid + 3-km buffer) from the 1-km marine-heatwave
raster of Kumagai et al. (2024) — the same product, ensuring spatial
comparability — and fit `lnRR ~ exposure(standardized) + (1|MPA)` per taxon.
Finally, and most importantly, because MPAs in our network were established over
2003–2012 while the heatwave is fixed in calendar time, we can separate the
ongoing pBACIPS *recovery* trajectory from any *heatwave-specific* effect. For
post-establishment data we fit `lnRR ~ time-since-establishment + exposure + (1|MPA)`
and used a likelihood-ratio test against a recovery-only model to ask whether the
heatwave shifted the cascade beyond the recovery trend — a separation Kumagai
et al. (2024) could not make without a pre-establishment baseline.

## Results — The MPA cascade strengthens with the heatwave

Before the heatwave, Southern California MPAs showed the expected predator signal
but little herbivore or producer response: sheephead and (weakly) lobster were
more abundant inside reserves, while urchins and giant kelp did not differ from
reference sites (Table X; Fig. X). The heatwave sharpened every link of the
cascade. The inside-MPA advantage for spiny lobster rose from no detectable
difference before (RR = 0.93) to two- and three-fold during and after the MHW
(RR = 2.05, 3.03; after-vs-before ΔlnRR = +1.18, p < 10⁻¹³). Urchins shifted from
neutral to strongly suppressed inside reserves — purple urchins RR 1.08 → 0.37
(ΔlnRR = −1.07, p < 10⁻¹⁶) and red urchins 1.22 → 0.63 (ΔlnRR = −0.65, p < 10⁻⁹).
Most strikingly, giant kelp reversed from roughly half as abundant inside reserves
before (RR = 0.47) to three- to four-fold more abundant inside during and after
the heatwave (RR = 2.59, 3.51; ΔlnRR = +2.01, p < 10⁻²⁶) — the signature of MPAs
conferring resistance and recovery. California sheephead were the exception:
consistently more abundant inside reserves (RR ≈ 1.8–2.0) but unchanged by the
heatwave (ΔlnRR, p ≥ 0.23). Every shift held when biomass and density were modelled
separately (after-vs-before p < 0.01 in both for lobster, both urchins, and kelp),
and were robust to restricting the panel to MPAs sampled in all three periods, to
no-take SMRs only, and to dropping the sheephead-only MPAs
(`table_heatwave_sensitivity.csv`; lobster after-vs-before ΔlnRR +0.97 to +1.28,
purple urchin −1.02 to −1.25, kelp +1.93 to +2.36, all p < 10⁻⁶).

The cascade regressions corroborate the mechanism. Across MPA-years, kelp lnRR
declined with urchin lnRR (purple: slope −0.51, p = 2×10⁻⁹; red: −0.52, p = 9×10⁻⁶)
and urchin lnRR declined with lobster lnRR (−0.33, p = 2×10⁻⁶), while kelp lnRR rose
with lobster lnRR (+0.71, p = 1×10⁻¹¹). Sheephead lnRR predicted neither urchins
(p = 0.20) nor kelp (p = 0.28), indicating a lobster-led cascade.

The same pattern emerges spatially: where local heatwave exposure was more intense,
the inside-MPA advantage was larger for kelp (+0.41 lnRR per SD of exposure,
p = 2×10⁻⁶) and lobster (+0.23, p = 0.002) and the disadvantage stronger for purple
urchins (−0.19, p = 0.001), while sheephead and red urchins did not respond to local
exposure (p ≥ 0.20).

Separating recovery from heatwave, however, shows these signals have different
origins. Controlling for time since MPA establishment (the pBACIPS recovery
trajectory), the inside-MPA shifts in lobster, sheephead, and red urchins are
explained by recovery and add no detectable heatwave-specific effect (heatwave
term p ≥ 0.34; recovery term p < 10⁻³ for all). For purple urchins the heatwave
adds a marginal further suppression beyond recovery (−0.11 per SD exposure,
p = 0.08). Only giant kelp shows a clear heatwave-specific effect: beyond a strong
recovery trend (+0.12 lnRR yr⁻¹, p < 10⁻⁷), kelp's inside-MPA advantage rose
further with contemporaneous heatwave exposure (+0.31 per SD, p = 0.003;
likelihood-ratio p = 0.003). In other words, the predator and urchin components of
the cascade strengthened steadily as reserves matured, whereas MPAs delivered an
additional, acute buffering of giant kelp specifically during the heatwave.

## Discussion — Convergence with, and extension of, Kumagai et al. (2024)

These results independently corroborate Kumagai et al. (2024), who found that
Southern California's no-take MPAs buffered kelp forests against the 2014–2016
heatwave through a predator–urchin–kelp cascade. We reach the same conclusion from
a different design and a broader dataset: where Kumagai contrasted protected and
unprotected sites during the heatwave, we tracked the inside-vs-reference response
ratio across before/during/after periods within a network monitored by PISCO, the
NPS Kelp Forest Monitoring program (Channel Islands), and SBC LTER. The agreement
across data sources and analytical framings strengthens the inference that the
cascade — not a confound of MPA placement — drives the resilience signal.

Our before/after-establishment design also lets us go a step further than studies
that compare protected and unprotected sites only during a heatwave. Because our
reserves were established over a decade (2003–2012), we can separate the steady
recovery of the cascade as reserves mature from any acute, heatwave-specific
effect. Doing so reveals that most of the apparent heatwave "strengthening" — in
lobster, sheephead, and red urchins — reflects ongoing recovery rather than a
distinct response to the 2014–2016 event. The clear exception is giant kelp, for
which MPAs delivered an acute resilience benefit during the heatwave over and above
the recovery trend. This refines the prevailing narrative: predator recovery and
urchin control build slowly with protection, but their payoff for the foundation
species becomes most visible precisely when a climate shock arrives.

Two finer-grained results also converge. The cascade is carried by spiny lobster,
not sheephead: lobster's inside-MPA advantage tripled across the heatwave and
predicted both urchin suppression and kelp gain, whereas sheephead's advantage was
large but static and predicted neither — matching Kumagai's finding that kelp
tracked lobster (not sheephead), plausibly because the 2015–2016 sheephead
recruitment pulse decoupled sheephead numbers from their grazing effect. Second,
red urchins were suppressed less than purple urchins, consistent with the partial
release of red urchins from fishing inside reserves (Malakhoff & Miller 2021)
acting against the cascade — a complication we document elsewhere in this study.

**Caveats.** The recovery-vs-heatwave decomposition controls for time since
establishment but uses *contemporaneous* heatwave exposure, so it isolates the
acute effect and attributes the sustained post-2016 kelp elevation to recovery;
the kelp heatwave benefit is therefore a conservative (lower-bound) estimate, and
the full NLS-pBACIPS extraction-time machinery is not yet integrated. Per-MPA
heatwave exposure is derived by sampling Kumagai et al.'s
1-km raster at MPA centroids (centroid + buffer); polygon-area-weighted exposure
would be a modest refinement. Cross-MPA tests of the cascade are underpowered (≈11 MPAs with paired data), so we
rest the mechanistic inference on the well-replicated MPA-year regressions rather
than on per-MPA before→after change correlations (which were non-significant).
Analysis is restricted to the Southern California Bight; we do not engage the
Central/Northern California (sea-otter / sunflower-star) dynamics that Kumagai and
others address.

---

## Direct data cross-check vs Kumagai et al. (2024)

Computing the inside-vs-outside response independently on each dataset
(`table_heatwave_crosscheck_kumagai.csv`) confirms convergence. From Kumagai's
processed Southern California subtidal data we computed the Full-vs-Reference raw
density ratio per period; from ours, the proportion-based paired lnRR (exponentiated).
Both show the same cascade through the heatwave — lobster high/rising, sheephead flat,
purple urchins collapsing (nearly identical: 1.1 → 0.45 → 0.4 in both), red urchins
declining, and kelp rising sharply (after-period inside-advantage 3.1–3.5× in both).
The principal difference is the pre-heatwave baseline for lobster and kelp, which sit
higher inside reserves in their raw-density framing but near or below parity in our
proportion-based, paired-reference, KFM-inclusive framing — the baseline difference
the pBACIPS design is intended to absorb. The agreement across two datasets and two
analytical traditions is strong evidence that the predator–urchin–kelp cascade, not a
methodological artifact, underlies MPA-conferred kelp resilience to marine heatwaves
in Southern California.
