# Results Summary: CA MPA Kelp Forest Analysis

**Generated:** 2026-06-17 06:42:29
**Pipeline version:** Modular pBACIPS v2.0

---

## Overview

This document summarizes the key results from the California MPA kelp forest
pBACIPS analysis. Effect sizes are estimated as log response ratios (lnRR = ln[MPA/Reference])
at t=11 years post-MPA implementation, and back-transformed to response ratios (RR = exp(lnRR)).

**Interpretation:**
- **lnRR > 0 / RR > 1** = higher values inside MPA vs reference
- **lnRR < 0 / RR < 1** = lower values inside MPA vs reference
- **lnRR = 0 / RR = 1** = no difference
- **% Change** = (RR - 1) x 100 (e.g., RR = 1.72 means +72% higher inside MPA)
- Significance threshold: p < 0.05

---

## Meta-Analysis Results (Table 2)

**Estimate** = log response ratio (lnRR). **RR** = back-transformed response ratio = exp(lnRR).
RR > 1 means higher inside MPA; RR < 1 means lower inside MPA.

### Biomass

| Taxa | k | Estimate | RR | % Change | SE | t | p-value | p (FDR) | 95% CI (RR) | Effect |
|------|---|----------|----|----------|----|----|---------|---------|-------------|--------|
| S. purpuratus | 10 | -0.4461 | 0.64 | -36% | 0.2368 | -1.884 | 0.0637 | 0.0955 | [0.399, 1.026] | negative |
| M. franciscanus | 14 | 0.5978 | 1.818 | +81.8% | 0.2214 | 2.7008 | 0.0086** | 0.0156 | [1.169, 2.827] | positive |
| M. pyrifera | 31 | 0.6221 | 1.863 | +86.3% | 0.2302 | 2.7026 | 0.0086** | 0.0156 | [1.177, 2.948] | positive |
| P. interruptus | 5 | 0.2223 | 1.249 | +24.9% | 0.305 | 0.7289 | 0.4685 | 0.5270 | [0.68, 2.294] | positive |
| S. pulcher | 16 | 1.0961 | 2.992 | +199.2% | 0.2261 | 4.8477 | 0.0000*** | 0.0001 | [1.906, 4.697] | positive |

### Density

| Taxa | k | Estimate | RR | % Change | SE | t | p-value | p (FDR) | 95% CI (RR) | Effect |
|------|---|----------|----|----------|----|----|---------|---------|-------------|--------|
| S. purpuratus | 13 | -1.4692 | 0.23 | -77% | 0.3827 | -3.839 | 0.0003*** | 0.0013 | [0.107, 0.495] | negative |
| M. franciscanus | 16 | -0.3991 | 0.671 | -32.9% | 0.3744 | -1.0659 | 0.2906 | 0.3736 | [0.317, 1.418] | negative |
| P. interruptus | 17 | 1.0836 | 2.955 | +195.5% | 0.3756 | 2.8851 | 0.0054** | 0.0156 | [1.395, 6.261] | positive |
| S. pulcher | 20 | -0.1728 | 0.841 | -15.9% | 0.3736 | -0.4625 | 0.6454 | 0.6454 | [0.399, 1.775] | negative |

*Significance: \*p<0.05, \*\*p<0.01, \*\*\*p<0.001 (uncorrected). p (FDR) = Benjamini-Hochberg adjusted p-values.*

### Heterogeneity (Joint Model)

tau2 and pseudo-I2 are shared across all taxa within each response type.

| Response | tau2 | I2 (%) |
|----------|------|--------|
| Biomass | 0.2262 | 77.7 |
| Density | 0.839 | 95.5 |

---

## Replicate-Level Summary

- **Significant increases inside MPA:** 59
- **Significant decreases inside MPA:** 34
- **Not significant:** 49
- **Total MPA-level effect sizes:** 142

### Per-Taxa Breakdown

| Taxa | n | Significant | Raw Mean lnRR |
|------|---|-------------|---------------|
| P. interruptus | 22 | 13 (59%) | 1.036 |
| S. pulcher | 36 | 25 (69%) | 0.558 |
| S. purpuratus | 23 | 18 (78%) | -1.22 |
| M. franciscanus | 30 | 21 (70%) | -0.378 |
| M. pyrifera | 31 | 16 (52%) | 1.376 |

*Raw Mean lnRR is the arithmetic mean of individual MPA effect sizes (bio + den combined), not the meta-analytic estimate. See Table 2 above for meta-analytic means.*

Full replicate-level effect sizes are in `outputs/replicate_effects.csv`.

---

## Cross-Taxa Meta-Regressions (Table 3)

Knapp-Hartung meta-regressions testing whether predator MPA effects predict prey effects across sites.

| Relationship | k | Slope | SE | p | Significant? |
|-------------|---|-------|----|---|-------------|
| S. purpuratus density -> M. pyrifera biomass | 11 | 0.057 | 0.235 | 0.807 | no |
| M. franciscanus density -> M. pyrifera biomass | 10 | -0.126 | 0.134 | 0.349 | no |
| P. interruptus density -> S. purpuratus density | 11 | 0.154 | 0.41 | 0.707 | no |
| S. pulcher density -> S. purpuratus density | 11 | -0.125 | 0.688 | 0.856 | no |
| P. interruptus biomass -> S. purpuratus biomass | 4 | -0.422 | 0.188 | 0.025 | **yes** |
| S. pulcher biomass -> S. purpuratus biomass | 9 | 0.29 | 0.183 | 0.114 | no |

---

## Variance Components

| Response | Component | Estimate | 95% CI |
|----------|-----------|----------|--------|
| Biomass | MPA | 0.1938 | [0.0306, 0.7792] |
| Biomass | Source | 0 | [0, 3.0592] |
| Density | MPA | 0.3996 | [0.0297, 1.9906] |
| Density | Source | 0.7438 | [0, 21.4301] |

---

## Temporal Recovery Slopes (Pooled lmer Model)

Species-specific slopes from `lmer(lnRR ~ time*species + (1+time|MPA) + (1|source))`.

| Species | Slope (lnRR/yr) | SE | p |
|---------|-----------------|-----|---|
| P. interruptus | 0.1115 | 0.0286 | — |
| S. pulcher | -0.0238 | 0.0244 | — |
| S. purpuratus | -0.0781 | 0.0256 | — |
| M. franciscanus | -0.0596 | 0.0256 | — |
| M. pyrifera | 0.1642 | 0.0251 | — |

---

## Outlier Sensitivity

The primary analysis uses the **full dataset** (no outlier removal, k = 142).
A global Cook's D threshold of 4/n flagged **93 of 142** observations (65%) as influential.
These flagged observations span all taxa and arise primarily from between-taxon divergence
(the expected trophic cascade signal) rather than within-taxon anomalies.

A per-taxon Cook's D threshold (4/k) flagged only **14 of 142** observations,
confirming the global threshold's over-flagging.

See Table S9 for 4-method sensitivity comparison and Figures S13-S15 for visualizations.

---

For known limitations, see `docs/ANALYSIS_REVISIONS.md`.

*This summary was auto-generated by the analysis pipeline.*
*For questions, contact Emily Donham or Adrian Stier.*
