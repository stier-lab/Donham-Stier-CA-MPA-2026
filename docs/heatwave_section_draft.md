# Heatwave analysis — draft manuscript text (Southern California)

*Draft for inclusion in Donham & Stier. Numbers verified from `01_formal_models.R`
(emmeans/contrasts) on harmonized response ratios joined to the SBC LTER marine
heatwave catalog. Scope: Southern California Bight only (all 15 MPAs lie south of
Point Conception, 32.7–34.45°N). Needs a `stier-writing-voice` polish pass.*

---

## Methods — Marine heatwave analysis

To test whether the strength of the MPA-driven trophic cascade changed with the
2014–2016 marine heatwave (MHW), we classified each survey year as **before**
(≤2013), **during** (2014–2016), or **after** (≥2017) the MHW, following the
windows used by Kumagai et al. (2024) for the same region. Heatwave timing and
intensity were taken from a Hobday et al. (2016) event catalog built with the
`heatwaveR` package on the daily NOAA OISST v2.1 record for the Santa Barbara
Coastal region (90th-percentile threshold, ≥5-day minimum duration, 30-year
climatology), maintained as the lab's single source-of-truth physical-environment
product. In this catalog the 2014–2016 MHW is the most extreme of the 1982–2026
record: two "Severe" (Hobday category III) events spanning Oct 2014–Apr 2015
(206 d) and Jul–Nov 2015 (129 d), with 189, 277, and 72 heatwave-days in 2014,
2015, and 2016 versus a median of ~0 in pre-2014 years.

For each focal taxon we modelled the log response ratio of inside-MPA to
reference abundance (lnRR; the same proportion-based response ratio used
throughout) as a function of heatwave period with a site random intercept:
`lnRR ~ period + (1 | MPA)` (linear mixed model, `lme4`). Because protection is
already encoded in the response ratio, the period term tests directly whether the
MPA effect on each taxon shifted with the heatwave. We estimated marginal means
per period and contrasts of during-vs-before and after-vs-before with `emmeans`.

## Results — The MPA cascade strengthens with the heatwave

Before the heatwave, MPAs in Southern California showed the expected predator
signal but little herbivore or producer response: California sheephead and (more
weakly) spiny lobster were more abundant inside reserves, while urchins and giant
kelp did not differ from reference sites (Table X; Figure X). The heatwave
sharpened every link of the cascade. The inside-MPA advantage for spiny lobster
rose from no detectable difference before (RR = 0.93) to roughly two- and
three-fold during and after the MHW (RR = 2.05 and 3.03; ΔlnRR after = +1.18,
p < 10⁻¹³). Concurrently, urchins shifted from neutral or slightly higher inside
reserves to markedly suppressed: purple urchins fell from RR = 1.08 to 0.37
(ΔlnRR after = −1.07, p < 10⁻¹⁶) and red urchins from 1.22 to 0.63
(Δlnrr after = −0.65, p < 10⁻⁹). Most strikingly, giant kelp reversed from being
roughly half as abundant inside reserves before the heatwave (RR = 0.47) to three-
to four-fold more abundant inside during and after it (RR = 2.59 and 3.51;
Δlnrr after = +2.01, p < 10⁻²⁶) — the signature of MPAs conferring resistance to
and recovery from the climate shock. California sheephead were the exception: they
were consistently more abundant inside reserves (RR ≈ 1.8–2.0) but this advantage
did not change with the heatwave (Δlnrr, p ≥ 0.23).

## Discussion — Convergence with, and extension of, Kumagai et al. (2024)

These results independently corroborate Kumagai et al. (2024), who found that
Southern California's no-take MPAs buffered kelp forests against the 2014–2016
heatwave through a predator–urchin–kelp cascade. We reach the same conclusion from
a different design and a broader dataset: where Kumagai contrasted protected and
unprotected sites during the heatwave, we tracked the inside-vs-reference response
ratio across before/during/after periods within a network monitored by PISCO, the
NPS Kelp Forest Monitoring program (Channel Islands), and SBC LTER. The agreement
across data sources and analytical framings strengthens the inference that the
cascade — not some confound of MPA placement — drives the resilience signal.

Two finer-grained results also converge. First, the cascade is carried by spiny
lobster rather than sheephead: lobster's inside-MPA advantage tripled across the
heatwave while sheephead's was large but static, matching Kumagai's finding that
giant kelp tracked lobster (not sheephead) abundance, plausibly because the
2015–2016 sheephead recruitment pulse decoupled sheephead numbers from their
grazing effect. Second, red urchins show a weaker suppression than purple urchins,
consistent with the partial release of red urchins from fishing inside reserves
(Malakhoff & Miller 2021) acting against the cascade — a complication we document
elsewhere in this study.

[Caveats to fold in: this is a period-level test on the response ratio and does
not yet use the full pBACIPS baseline correction or per-MPA heatwave exposure;
analysis is restricted to the Southern California Bight; the before-period kelp
deficit inside reserves (RR = 0.47) warrants the pBACIPS treatment that motivates
the main analysis.]
