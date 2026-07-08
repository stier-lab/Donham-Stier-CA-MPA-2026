# Supplement — How analytical method shapes inferred MPA × heatwave effect sizes

*Draft for the Supporting Information of Donham & Stier. All numbers verified from
`code/R/15_methods_comparison.R` and `tables/table_s_methods_*.csv`. Companion to
the heatwave analysis (`docs/heatwave_section_draft.md`) and the direct
data cross-check against Kumagai et al. (2024). Needs a `stier-writing-voice`
polish. Scope: Southern California Bight, density responses, heatwave-window
contrasts (before ≤2013 / during 2014–16 / after ≥2017).*

---

## Why this supplement

Our study and Kumagai et al. (2024, *Global Change Biology*) ask nearly the same
question of an overlapping system — did Southern California's no-take MPAs buffer
the kelp-forest trophic cascade against the 2014–2016 marine heatwave? — and reach
a convergent headline (yes, carried by giant kelp). Yet the two studies differ on
almost every analytical choice at once: we use a paired-reference, proportion-based
log response ratio synthesised in a multilevel inverse-variance meta-analysis,
whereas they pool protected against unprotected sites within a region and fit raw
density in one-stage GLMMs. When two pipelines this different agree on some taxa
and diverge on others, the natural question is whether the disagreements are
*biological*, *data-driven*, or merely *artifacts of method*. This supplement
answers that question directly, by re-running the analysis many ways and measuring
how much each choice moves the result.

## Approach — an analytical multiverse

We decomposed the difference between the two pipelines into two parts.

**A method bridge (holds the data constant; isolates method).** Starting from a
faithful reimplementation of the Kumagai-style endpoint and walking to our own, we
flipped exactly one analytical choice at a time, on a single common dataset (our
harmonized PISCO + NPS-KFM + SBC-LTER data), and recorded how the per-taxon effect
size moved at each step:

| Waypoint | Reference | Scale | Estimator | = endpoint |
|---|---|---|---|---|
| WP0 | pooled | raw density | Tweedie GLMM (log link) | **Kumagai-style** |
| WP1 | pooled | proportion | Tweedie GLMM (log link) | *flip: response scale* |
| WP2 | pooled | proportion | Gaussian LMM on log | *flip: distribution / link* |
| WP3 | paired | proportion | LMM on lnRR | *flip: reference / pairing* |
| WP4 | paired | proportion | LMM on lnRR + AR1(year\|MPA) + Source RE | *flip: temporal autocorrelation* — **our heatwave endpoint** |
| WP5 | paired | proportion | inverse-variance meta-analysis (`rma.mv`) | *flip: weighting* — **our main-paper synthesis** |

Each model is coded so that the interaction (pooled designs) or period (paired
designs) coefficients are directly the *during-versus-before* and
*after-versus-before* changes in ln(MPA / Reference), with standard errors — a
common currency across all six waypoints. The change at each step is attributed to
the single choice flipped at that step.

**A data effect (holds method constant; isolates data).** The Kumagai dataset is
structured for pooled, not paired, analysis (it carries no MPA-to-reference
mapping), so it can support only the three pooled waypoints (WP0–WP2). We ran those
three on *both* datasets — ours and theirs (Kumagai's
`MLPA_data_summarized_wo_siteblocks.csv`, Southern Coast, Full vs Reference) — and
read the gap between substrates at a matched waypoint as the contribution of the
data itself.

The full grid (44 cells: 6 waypoints × 5 taxa on our data + 3 pooled waypoints × 5
taxa on theirs) is in `table_s_methods_multiverse.csv`; the per-axis decomposition
is in `table_s_methods_decomposition.csv`; a qualitative crosswalk of all the
analytical choices is in `table_s_methods_crosswalk.csv`.

## Results

**Giant kelp resilience is method-invariant.** The after-versus-before increase in
giant kelp's inside-MPA advantage is large and highly significant at *every*
waypoint — from the Kumagai-style endpoint (Δln = +1.06, *p* < 10⁻⁵) through our
meta-analytic endpoint (Δln = +1.44, *p* < 10⁻⁹) — and on both datasets. No single
analytical choice, and no combination, makes the kelp signal appear or disappear.
The foundation-species result that anchors both papers does not depend on how it is
analysed.

**The lobster signal depends on pairing.** Spiny lobster's post-heatwave increase
inside reserves is invisible under the Kumagai-style pooled, raw-density model
(Δln = +0.04, *p* = 0.86) and emerges only once the response is standardised and,
especially, once reserves are paired to their own references (Δln = +0.49 to +0.86,
*p* ≤ 0.008 from WP1 onward). Pooling protected against unprotected sites buries
the lobster effect in among-site abundance differences that pairing removes.

**The urchin inference — not its magnitude — is set by the autocorrelation
choice.** This is the sharpest result. Across waypoints the urchin point estimates
are stable (purple urchin after-versus-before Δln ≈ −0.6 to −0.9 throughout), but
the *inference* swings entirely on whether temporal autocorrelation is modelled.
Without an AR1 structure, purple-urchin suppression is "highly significant"
(WP3, *p* = 1.4 × 10⁻⁶); adding AR1(year | MPA) — the spec our heatwave analysis
adopts as primary — renders the same effect non-significant (WP4, *p* = 0.12); and
an autocorrelation-naive inverse-variance meta-analysis returns it to spurious high
significance (WP5, *p* ~ 10⁻¹⁰). Red urchin behaves identically (WP3 *p* = 4 × 10⁻⁶
→ WP4 *p* = 0.10). The autocorrelation flip moves the point estimate by only 0.14
ln units on average yet flips the significance of two of five taxa — it is an
*inferential* lever, not a magnitude lever.

**California sheephead is a robust null.** Sheephead shows no heatwave-related shift
in its inside-MPA advantage at any waypoint, on either dataset (|Δln| ≤ 0.20,
*p* > 0.05 throughout).

**The data matter more than any single method choice.** Of all the axes we varied,
the largest mover of effect sizes is the dataset itself: switching from our data to
Kumagai's at a matched pooled waypoint moves the after-versus-before estimate by
0.63 ln units on average (up to 1.83 for giant kelp) and flips the significance of
all five taxa. The biggest single methodological lever on the point estimate is the
distribution/link choice (mean |Δln| = 0.29), followed by the response scale and the
weighting (≈ 0.20 each); pairing and autocorrelation move the point estimate least
(0.09 and 0.14) even though autocorrelation is decisive for inference.

| Axis (one analytical flip) | mean \|Δln\| | max \|Δln\| | taxa with significance flip |
|---|---|---|---|
| Response scale (raw → proportion) | 0.20 | 0.44 | 1 |
| Distribution / link (Tweedie GLMM → log-LMM) | 0.29 | 0.51 | 1 |
| Reference (pooled → paired) | 0.09 | 0.21 | 0 |
| Temporal autocorrelation (no AR1 → AR1 + Source RE) | 0.14 | 0.31 | **2** |
| Weighting (AR1 LMM → IV meta-analysis) | 0.20 | 0.65 | 2 |
| **Data (ours vs Kumagai), matched pooled waypoints** | **0.63** | **1.83** | **5** |

## What this means

The disagreement between the two studies is not arbitrary. Where they agree — giant
kelp's resilience and sheephead's flat response — the result is robust to every
analytical choice and to the dataset, and can be reported with confidence. Where
they appear to diverge — the strength of urchin suppression — the divergence is
largely an artifact of one choice: whether the strong year-to-year autocorrelation
in annual response ratios is modelled. An autocorrelation-naive analysis (whether a
simple mixed model or an inverse-variance meta-analysis, and including the
period-mean comparisons both papers lean on) returns urchin significance that does
not survive an AR1 correction. This is the methodological basis for the tempered
urchin conclusion in our main heatwave analysis, and it cautions against strong
claims of urchin-mediated control from observational reserve data without accounting
for temporal dependence.

Two broader points follow. First, the single most influential choice is not a
modelling option at all but the data: our broader, three-program dataset and
Kumagai's PISCO/MLPA core can yield effect sizes differing by more than an order of
magnitude at the same waypoint (most strikingly for giant kelp, whose pooled
after-period advantage is far larger in their raw-density framing). Convergent
conclusions from independent datasets are therefore more reassuring than any single
analysis. Second, the pairing and standardisation choices that distinguish our
pBACIPS design — pairing each reserve to its own reference and standardising by site
composition — chiefly *reshape the baseline* (they lower the apparent pre-heatwave
inside-MPA advantage for lobster and kelp) rather than the heatwave response itself,
which is why our before-period response ratios sit below theirs while the
during/after estimates converge.

---

## Environmental covariates: testing what Kumagai mapped but did not model

Kumagai et al. extracted a rich environmental layer — sea-surface temperature,
significant wave height, nitrate, marine-heatwave and cold-spell cumulative
intensity, depth, and a human-gravity index — but used it only descriptively, in a
principal-components ordination (their Fig. 3); none of these covariates entered any
of their inferential models. Because this is the one substantive predictor set that
neither study had formally tested, we asked directly whether the MPA effect varies
across these gradients (`code/R/16_environmental_moderators.R`).

We assembled, per MPA, the covariates obtainable at adequate spatial resolution
without fabrication: marine-heatwave cumulative intensity during 2014–2016
(thermal stress; Kumagai's 1-km MHW raster sampled at MPA centroids), cold-spell
cumulative intensity (a thermal-regime / upwelling proxy; Kumagai's matching 1-km
cold-spell grid, nearest cell), latitude (a biogeographic / thermal gradient), and
log reserve area (reserve design). We then ran, per taxon, a random-effects
meta-regression of the per-MPA effect size on each standardised moderator
(`metafor::rma`, REML, Knapp–Hartung small-sample *t*-test on *k* − 2 df; density
effect sizes for the four animals, biomass for giant kelp, which has no density
effect sizes in our pipeline), with Benjamini–Hochberg FDR across all tests.
We added **wave exposure** and **nitrate (nutrients)** for all MPAs (including the
Channel Islands) from Bell's (2023) island-inclusive kelp-canopy-environment product —
the per-pixel NetCDF (EDI `knb-lter-sbc.162`; 332,640 Landsat-pixel stations ×
quarterly 1984–2021; CDIP MOP v1.1 waves, SST-derived nitrate) that underlies Kumagai
et al.'s `hsmax` and Wanner et al. (2024)'s Channel-Islands kelp-synchrony analysis.
Per-MPA values = mean climatological covariate over kelp-pixel stations within 3 km of
each MPA (`code/R/extract_kelp_env_covariates.R` → `data/per_mpa_kelp_env.csv`, which
also holds temperature/npp/depth; provenance in `~/sbc-kelp-env/`); all 34 MPAs are
covered (median 0.03 km to the nearest station) and the values recover the expected
gradients (wave: San Miguel Island 4.9 m → Catalina-lee Blue Cavern 0.65 m; nitrate:
cold upwelling 3.0 µM at Pt Conception / NW islands → warm 0.39 µM at San Diego).
We also added **human gravity** (a Cinner et al. 2018-style population-pressure /
accessibility index — population within 50 km weighted by inverse-square distance)
from Kumagai et al.'s per-kelp-patch grid, matched to the nearest point; it spans
0 at the remote NW islands (San Miguel, Santa Rosa, Santa Barbara Is.) to ~22,000 at
urban-mainland reserves (Campus Point, Palos Verdes, San Diego), correctly capturing
the human-pressure gradient. With this, the **entire set of environmental covariates
Kumagai mapped in their PCA is now obtained per MPA** (island-inclusive) and tested.
(An earlier attempt used the mainland-only tabular product `knb-lter-sbc.144`, whose
500-m coastline segments left the island MPAs 32–105 km from data; it is superseded.)

The result is a clean null. **No environmental moderator survives FDR correction**
(`table_s_env_moderators.csv`). The strongest signal is a decline in the red-urchin
(*M. franciscanus*) effect with marine-heatwave intensity (−1.84 lnRR per SD,
*p* = 0.003, FDR = 0.21) — suggestive but not robust, and directionally sensible
(red urchins fared relatively worse inside the most heat-exposed reserves).
**Wave exposure**, **nitrate**, and **human gravity** (all 34 MPAs incl. islands,
*k* = 10–17) do not significantly modulate any taxon's effect either (none survive
FDR; nitrate's only nominal signal, red urchin *p* = 0.04, mirrors its MHW signal
because nitrate is SST-derived and tracks the cold-upwelling gradient; human gravity's
nominal lobster signal, *p* = 0.07, is in the expected direction — weaker reserve
benefit near cities — but does not survive correction).
Critically, **giant kelp's effect is not detectably modulated by any gradient**
(all *p* ≥ 0.17): its heatwave resilience is not an artifact of where reserves
happen to sit thermally, biogeographically, by size, or by wave exposure. The Knapp–Hartung
adjustment is essential here — without it the same *k* = 9–15 regressions return
implausibly tiny *p*-values (down to 10⁻¹⁸) because inverse-variance weighting with
a handful of high-precision MPAs is overconfident; these tests are exploratory
and the thermal/biogeographic moderators are mutually collinear. The practical
upshot for the main analysis is reassuring: incorporating the environmental
covariates Kumagai only visualised does not reveal hidden modulation of the reserve
effect, so their omission from the primary models does not appear to bias the
headline conclusions.

---

### Reproducibility

`code/R/15_methods_comparison.R` and `code/R/16_environmental_moderators.R` (both
wired into `run_all.R` after the heatwave module). Script 15 requires the Kumagai
mirror at `~/kumagai2024-comparison/repo/` for the cross-substrate panel (if absent,
the method bridge on our data still runs and the data-effect panel is skipped);
script 16 requires the cold-spell grid in that mirror for the `cs_mean` covariate
(otherwise it is omitted). Outputs — script 15: `tables/table_s_methods_crosswalk.csv`,
`tables/table_s_methods_multiverse.csv`, `tables/table_s_methods_decomposition.csv`,
`plots/fig_s_methods_multiverse.{pdf,png}`; script 16:
`tables/table_s_mpa_env_covariates.csv`, `tables/table_s_env_moderators.csv`,
`plots/fig_s_env_moderators.{pdf,png}`.

### Figure legends

**Figure S — Analytical multiverse: how method and data shape inferred MPA ×
heatwave effect sizes (Southern California).** Each point is the after-versus-before
change in ln(MPA / Reference) abundance for a taxon (back-transformed to an x-fold
change on the RR-scaled axis), ± 95% CI. **(a)** The method bridge on a single
common dataset (our harmonized data), walking from the Kumagai-style endpoint (WP0,
orange) through one analytical flip per step to our endpoints (WP4–WP5, blue);
intermediate waypoints in grey. Giant kelp (*M. pyrifera*) sits well above the
no-effect line at every step; urchins (*S. purpuratus*, *M. franciscanus*) drift
toward the line at the AR1 step (WP4). **(b)** The data effect: the three pooled
waypoints (WP0–WP2) run on our data (blue) versus Kumagai's MLPA data (green),
showing that the dataset moves estimates as much as or more than any single method
choice. Waypoint definitions and the per-axis decomposition are in
`table_s_methods_decomposition.csv`.

**Figure S — The MPA effect across environmental gradients.** Random-effects
meta-regression slopes (change in MPA effect, lnRR, per SD of each standardised
moderator; ± 95% CI, Knapp–Hartung) for each taxon, faceted by moderator: marine-
heatwave intensity during 2014–2016, cold-spell intensity, latitude, and log reserve
area. Open circles are non-significant; filled triangles mark *p* < 0.05 (none
survive FDR). The covariates are those Kumagai et al. (2024) mapped in their PCA but
did not model. Effect sizes are density for the animals and biomass for giant kelp;
full statistics in `table_s_env_moderators.csv`, per-MPA covariates in
`table_s_mpa_env_covariates.csv`.
