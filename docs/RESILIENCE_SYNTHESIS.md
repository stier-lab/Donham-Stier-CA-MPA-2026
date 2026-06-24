# MPA resilience: synthesis of the analysis suite (scripts 14–21)

*Integrating the resilience analyses for Donham & Stier. Each section links to the
script that produces it and its standalone doc. Run the whole suite with
`source(here::here("code","R","run_resilience.R"))`.*

---

## The one-line theme

**California MPAs confer a robust, repeatable increase in giant-kelp resilience to
climate shocks — but the finer claims (the urchin-mediated mechanism, predator
diversity/compensation, and which reserves will be most effective) are
specification-sensitive, stressor-specific, or unpredictable.** The strongest,
most defensible statement is about the *foundation species* and the *average* effect;
the weakest links are *mechanism attribution* and *among-reserve prediction*.

## How the analyses map to facets of resilience

| Facet of resilience | Question | Script | Verdict |
|---|---|---|---|
| **Core response** | Does the MPA cascade respond to a heatwave? | `14` | Giant kelp far more abundant inside reserves during/after the 2014–16 MHW; predator/urchin steps weaker under AR1 |
| **Repeatability** | Does the resilience recur for a second event? | `19` | **Yes** — kelp elevated inside in *both* the 2014–16 and 2018–20 heatwaves (only taxon significant in both); tracks heat beyond recovery |
| **Resistance & recovery** | Did kelp hold/return on the *state* variable? | `22` | **Yes** — kelp *grew* inside (resistance 2.18× baseline) while declining outside (0.94×), recovered to 2.38× inside vs 0.83× outside, and by 2020–23 inside held 2.10× vs outside 0.64× (all p ≤ 0.014); recent urchin suppression inside (p ≈ 0.02) |
| **Memory / consistency** | Do the *same* reserves repeat? | `23` | **Yes** — a reserve's MHW1 response predicts its MHW2 response (kelp ρ=0.71, p=0.028; purple urchin r=0.87, p=0.001); resilience is a consistent reserve property, with no priming |
| **Generality across stressor types** | Does it hold for disease, not just heat? | `20` | The cascade response is **heat + fishing protection, not sea-star wasting** — the sunflower star was patchy (present at ~1/3 of reserves) and the response is the same where it was absent |
| **Stability** | Does protection damp variability? | `21` | **No clear effect** — kelp trends more stable inside (ns); urchins more variable inside. Resilience is in the *mean*, not the *variance* |
| **Robustness / attribution** | How much do method & data change the story? | `15` | Kelp resilience is method-invariant; urchin inference hinges on the AR1 choice; the **dataset is the single largest lever** |
| **Moderators** | What environmental gradients modulate the effect? | `16` | **None survive FDR** (30 tests); kelp effect not modulated by MHW/cold-spell intensity, latitude, reserve size, wave exposure, or nutrients (nitrate) — wave & nitrate now island-inclusive from Bell 2023 (EDI sbc.162) |
| **Predictability** | Can we predict which MPAs are most effective? | `18` | **No** — leave-one-out CV R² is negative; among-reserve variation is real but idiosyncratic |
| **Cross-study reproduction** | Does the paired-design literature reproduce? | `17` | Eisaguirre (2020) design-level results reproduce; the partial size-structure/diversity mechanism does **not** robustly reproduce |

## What is robust (report with confidence)

1. **Giant-kelp resilience to marine heatwaves.** Kelp is strongly more abundant
   inside reserves through the 2014–16 MHW and again through the 2018–20 MHW (script
   19: MHW1 +1.65, MHW2 +2.21 lnRR, both p < 10⁻³), tracks heatwave intensity beyond
   ongoing recovery, and the result is invariant to nearly every analytical choice
   (script 15) and to the environmental gradients reserves sit on (script 16). This is
   the headline, and it is as solid as observational reserve data allow.
2. **It is driven by heat + protection, not by sea-star wasting** (script 20): the
   kelp recovery is identical at reserves that never had a sunflower star, so the
   keystone-removal mechanism Eisaguirre documented at the cold Channel Islands does
   not generalize to most of Southern California.
3. **Convergence with independent work.** Our paired pBACIPS design reaches the same
   kelp-resilience conclusion as Kumagai's (2024) unpaired GLMM and Eisaguirre's
   (2020) Channel-Islands study — three datasets/designs, one conclusion (script 15
   cross-substrate; script 17 reproduction).
4. **Resistance, recovery, and reserve-level consistency.** On the state variable
   (script 22), kelp *grew* inside reserves during the heatwave while it *declined*
   outside, recovered above baseline inside while staying depressed outside, and the
   gap widened through 2023 — significant resistance *and* recovery. And the *same*
   reserves repeat across both heatwaves (script 23) — resilience is a consistent
   reserve property, not a different set of reserves each time.

## What is fragile or unresolved (state honestly, do not overclaim)

1. **The urchin-mediated cascade step.** Under the autocorrelation-robust AR1 model
   the predator→urchin→kelp mediation is weak (lobster→urchin dissolves; urchin
   suppression marginal) — and whether it appears at all depends on the AR1 choice
   (script 15) and on site/island selection (scripts 17, the Eisaguirre robustness
   sweep).
2. **Predator diversity / redundancy.** No independent diversity signal survives
   beyond total predation pressure; it is confounded with protection and reverses by
   subset. The "compensation" pattern is consistent with protection/abundance, not
   with diversity per se.
3. **Among-reserve effectiveness is consistent but unexplained.** Reserve responses
   are *repeatable* across the two heatwaves (script 23) — yet *not predictable* from
   environmental, design, or trophic covariates (script 18, negative LOO-CV R²). So
   reserves differ reliably and persistently, but we cannot yet say *why*: report the
   average effect and its heterogeneity (τ²/I²), and note that the heterogeneity is a
   stable reserve property awaiting better (e.g. wave, governance, connectivity)
   covariates — not a modelling target with current data.

## Limitations that cut across the suite

- **Southern California only** (lat ≤ 34.45 °N); the colder Channel-Islands
  sunflower-star dynamics are out of scope for our harmonized data.
- **Confounded stressors / time.** SSWD, the 2014–16 MHW, and the El Niño overlap in
  2013–2016; reserve maturation runs alongside. We separate heat-specific from
  recovery effects where possible (scripts 14, 19) but cannot fully isolate each
  driver.
- **Small k for moderators/prediction** (≈10–17 MPAs with full covariates) — genuine
  power limits, honestly reported via cross-validation and Knapp–Hartung tests.
- **Observational design** — partial-mechanism attribution (which predator, diversity
  vs abundance) is the recurring weak point, as the multiverse and reproductions show.

## Other analyses (status)

- [x] **Resistance vs recovery decomposition (Kumagai's framework) on the state
  variable** — done, script 22.
- [x] **Recovery completeness / return** — done, script 22 (`recovery_recent`: outside
  kelp stays at 0.64× baseline through 2023, incomplete recovery; inside fully recovered).
- [x] **Ecological memory / consistency** — done, script 23 (same reserves repeat).
- [ ] **Recovery *rate* / return time** — per-MPA recovery slopes (e.g. 2016→2020) to
  quantify *how fast*, not just whether, kelp recovered inside vs outside. *Feasible now.*
- [ ] **Alternative-stable-state / hysteresis** — formal test for persistent urchin
  barrens (a reserve's cascade failing to return); script 22 already hints outside
  sites did not recover. *Feasible with the time series.*
- [x] **Spatial wave-exposure (storm) resilience** — wave exposure (HSMAX) obtained
  for **all 34 MPAs including the Channel Islands** from Bell's (2023) island-inclusive
  per-pixel kelp-canopy-env NetCDF (EDI knb-lter-sbc.162, CDIP MOP v1.1 waves;
  `code/R/extract_wave_exposure.R` → `data/per_mpa_wave_exposure.csv`; recovers the
  San Miguel-exposed → Catalina-lee-sheltered gradient) and added to the moderator
  screen (script 16): no significant modulation of any taxon (k=10–17, all p≥0.35).
  Wave exposure does not explain the MPA effect or its among-reserve variation.

## Outputs and reproducibility

Run `code/R/run_resilience.R` for the whole suite. Harmonized-data scripts (14, 19,
21) also run inside `run_all.R`; scripts that read the Kumagai mirror (15, 16, 18) or
raw PISCO (17, 20) skip gracefully if that external data is absent. Per-analysis docs:
`heatwave_section_draft.md`, `heatwave_replication.md`, `compound_disturbance.md`,
`methods_comparison_supplement.md`, `eisaguirre_reproduction.md`,
`mpa_effectiveness_predictors.md`; tables in `tables/`, figures in `plots/`.
