# Proposed main-text addition: giant-kelp resilience to the marine heatwave

*Conservative, robust addition for Donham & Stier (Journal of Applied Ecology). Figure:
`plots/fig_kelp_resilience.{pdf,png}` (single column, 80 mm). Analysis:
`code/R/22_resistance_recovery.R`; statistics from `tables/table_resistance_recovery.csv`.
All numbers are from the giant-kelp (biomass) row, paired Wilcoxon signed-rank across
MPAs, relative to a 2010–2013 pre-heatwave baseline.*

---

## Results paragraph (draft)

Giant kelp — the foundation species of these forests — was resilient to the
2014–2016 marine heatwave inside reserves but degraded outside them. During the
heatwave, kelp biomass rose to 2.2× its pre-heatwave (2010–2013) baseline inside MPAs
while falling to 0.9× at paired reference sites (paired Wilcoxon across reserve–control
pairs, *p* = 0.010), and the divergence widened rather than closed: by 2020–2023 kelp
remained at 2.1× baseline inside reserves but had dropped to 0.6× outside (*p* = 0.014).
This was a within-pair effect, not a pooled-mean artifact: kelp held or grew more
inside than at its **own** paired control in 9 of 10 reserves during the heatwave and
8 of 10 reserves afterward (the lone resistance exception, Anacapa Island SMR,
recovered more slowly at the reserve than at its reference). Protection thus buffered
the foundation species against the most severe marine heatwave in the satellite
record, while unprotected reefs lost kelp and did not recover.

## Two figures

- **Main text (`fig_kelp_resilience`, single column):** the inside-vs-outside biomass
  trajectory through the heatwave era — the intuitive, communicable headline.
- **Supporting (`fig_kelp_resilience_paired`):** each reserve against its own control
  for resistance and recovery — shows the effect is consistent across the individual
  pairs (9/10, 8/10), with the two exceptions visible rather than hidden.

## Figure caption (draft)

**Figure X. Giant-kelp biomass inside (blue) and outside (orange) Channel Islands
marine protected areas, 2008–2023.** Points are annual mean biomass across MPAs
(± 1 SE); the 2014–2016 marine heatwave is shaded. Kelp increased and remained
elevated inside reserves through and after the heatwave while declining and failing to
recover at paired reference sites — MPA-conferred resistance and recovery of the
foundation species (during heatwave: 2.2× the 2010–2013 baseline inside vs 0.9×
outside, paired Wilcoxon *p* = 0.010; 2020–2023: 2.1× inside vs 0.6× outside,
*p* = 0.014).

## Why this is the conservative choice

This result is the only one that survives every robustness check applied to the
resilience analyses: the analytical multiverse (method-invariant), Benjamini–Hochberg
correction, the AR1 autocorrelation-robust model, the sea-star-wasting control (the
kelp response is the same where the sunflower star was absent), and all environmental
moderators (none significant). It is reported on the absolute state variable rather
than a derived ratio, uses a simple paired nonparametric test, and concerns the
foundation species — so it can be stated without the caveats that attach to the
urchin-mediation, predator-diversity, and among-reserve-prediction analyses.

*(Optional companion for the SI: the two-heatwave replication, `fig_heatwave_replication`,
showing the same kelp elevation inside reserves through both the 2014–16 and 2018–20
heatwaves — evidence the resilience is not a one-off.)*

## Robustness of the headline (SI paragraph, ready to paste)

*Source: `table_resistance_recovery_sensitivity.csv` (giant kelp, `22_resistance_recovery.R`).
Reports the inside/outside geometric-mean ratio with 95% CI, three paired tests, and
leave-one-reserve-out, across six pre-heatwave baseline windows. Pair-level reserve
counts used in the Results sentence are exported by the same script to
`table_resistance_recovery_pairs.csv` and `table_resistance_recovery_pair_counts.csv`;
`24_resilience_pipeline_check.R` checks both the values and the counts.*

The MPA advantage in giant kelp was robust to the choices underlying the paired test.
Across six pre-heatwave baseline windows (2008–2013 through 2013 alone), giant kelp was
2.5–10× more abundant inside reserves than at paired references during the heatwave,
after it, and through 2020–2023, and the paired Wilcoxon signed-rank test on
reserve-level log ratios remained significant in every window (all *p* ≤ 0.047; primary
2010–2013 baseline: *p* = 0.010, 0.014, 0.014 for resistance, recovery, and the
2020–2023 window). Reported on the primary baseline, the inside/outside ratio was 4.2×
(95% CI 1.3–13.0) during the heatwave, 5.2× (1.3–20.0) in the recovery window, and 4.9×
(1.5–16.0) by 2020–2023 — every confidence interval excluding parity. The result did not
depend on any single reserve: leaving each reserve out in turn, the worst-case *p*
stayed ≤ 0.027 and the inside/outside ratio never fell below 2.9× (direction never
reversed). The one qualification is that the 2020–2023 elevation is carried by
magnitude rather than by a unanimous count — the ratio and both the Wilcoxon and *t*
tests are significant, but only 8 of 10 reserves individually exceed their reference, so
the sign test is not significant (*p* = 0.11); we therefore describe the sustained
2020–2023 advantage in terms of magnitude rather than as a per-reserve tally.
