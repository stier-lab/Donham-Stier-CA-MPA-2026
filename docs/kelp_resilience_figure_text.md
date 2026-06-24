# Proposed main-text addition: giant-kelp resilience to the marine heatwave

*Conservative, robust addition for Donham & Stier (Conservation Letters). Figure:
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
