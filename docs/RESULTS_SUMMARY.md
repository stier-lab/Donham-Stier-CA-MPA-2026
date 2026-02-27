# Results Summary: CA MPA Kelp Forest Analysis

**Generated:** 2026-02-27 11:50:00
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
| S. purpuratus | 10 | -0.6694 | 0.512 | -48.8% | 0.2636 | -2.5399 | 0.0137* | 0.0261 | [0.302, 0.868] | negative |
| M. franciscanus | 14 | 0.6102 | 1.841 | +84.1% | 0.2422 | 2.5195 | 0.0145* | 0.0261 | [1.134, 2.989] | positive |
| M. pyrifera | 17 | 0.5838 | 1.793 | +79.3% | 0.2653 | 2.2004 | 0.0317* | 0.0476 | [1.054, 3.048] | positive |
| P. interruptus | 5 | 0.2542 | 1.289 | +28.9% | 0.3192 | 0.7961 | 0.4292 | 0.4828 | [0.681, 2.442] | positive |
| S. pulcher | 18 | 1.0475 | 2.851 | +185.1% | 0.2476 | 4.23 | 0.0001*** | 0.0007 | [1.737, 4.679] | positive |

### Density

| Taxa | k | Estimate | RR | % Change | SE | t | p-value | p (FDR) | 95% CI (RR) | Effect |
|------|---|----------|----|----------|----|----|---------|---------|-------------|--------|
| S. purpuratus | 13 | -1.4643 | 0.231 | -76.9% | 0.3731 | -3.9247 | 0.0002*** | 0.0010 | [0.11, 0.487] | negative |
| M. franciscanus | 16 | -0.394 | 0.674 | -32.6% | 0.3645 | -1.0809 | 0.2838 | 0.3649 | [0.326, 1.397] | negative |
| P. interruptus | 17 | 1.0864 | 2.964 | +196.4% | 0.3658 | 2.9699 | 0.0042** | 0.0126 | [1.427, 6.155] | positive |
| S. pulcher | 22 | -0.1657 | 0.847 | -15.3% | 0.3633 | -0.4562 | 0.6498 | 0.6498 | [0.41, 1.751] | negative |

*Significance: \*p<0.05, \*\*p<0.01, \*\*\*p<0.001 (uncorrected). p (FDR) = Benjamini-Hochberg adjusted p-values.*

### Heterogeneity (Joint Model)

tau2 and pseudo-I2 are shared across all taxa within each response type.

| Response | tau2 | I2 (%) |
|----------|------|--------|
| Biomass | 0.2238 | 80.5 |
| Density | 0.7658 | 95.1 |

---

## Replicate-Level Summary

- **Significant increases inside MPA:** 56
- **Significant decreases inside MPA:** 35
- **Not significant:** 41
- **Total MPA-level effect sizes:** 132

### Per-Taxa Breakdown

| Taxa | n | Significant | Raw Mean lnRR |
|------|---|-------------|---------------|
| P. interruptus | 22 | 13 (59%) | 1.036 |
| S. pulcher | 40 | 26 (65%) | 0.524 |
| S. purpuratus | 23 | 19 (83%) | -1.237 |
| M. franciscanus | 30 | 21 (70%) | -0.337 |
| M. pyrifera | 17 | 12 (71%) | 2.111 |

*Raw Mean lnRR is the arithmetic mean of individual MPA effect sizes (bio + den combined), not the meta-analytic estimate. See Table 2 above for meta-analytic means.*

Full replicate-level effect sizes are in `outputs/replicate_effects.csv`.

---

## Cross-Taxa Meta-Regressions (Table 3)

Knapp-Hartung meta-regressions testing whether predator MPA effects predict prey effects across sites.

| Relationship | k | Slope | SE | p | Significant? |
|-------------|---|-------|----|---|-------------|
| S. purpuratus density -> M. pyrifera biomass | 11 | 0.169 | 0.415 | 0.683 | no |
| M. franciscanus density -> M. pyrifera biomass | 10 | -0.035 | 0.254 | 0.891 | no |
| P. interruptus density -> S. purpuratus density | 11 | 0.154 | 0.41 | 0.707 | no |
| S. pulcher density -> S. purpuratus density | 11 | -0.095 | 0.692 | 0.890 | no |
| P. interruptus biomass -> S. purpuratus biomass | 4 | -0.42 | 0.187 | 0.025 | **yes** |
| S. pulcher biomass -> S. purpuratus biomass | 9 | 0.29 | 0.166 | 0.081 | no |

---

## Variance Components

| Response | Component | Estimate | 95% CI |
|----------|-----------|----------|--------|
| Biomass | MPA | 0.1626 | [0.0085, 0.7434] |
| Biomass | Source | 0.19 | [0, 8.1082] |
| Density | MPA | 0.2461 | [0.0166, 1.3064] |
| Density | Source | 0.7086 | [0, 19.4104] |

---

## Temporal Recovery Slopes (Pooled lmer Model)

Species-specific slopes from `lmer(lnRR ~ time*species + (1+time|MPA) + (1|source))`.

| Species | Slope (lnRR/yr) | SE | p |
|---------|-----------------|-----|---|
| P. interruptus | 0.1082 | 0.0281 | < 0.001 |
| S. pulcher | -0.0274 | 0.0236 | 0.246 |
| S. purpuratus | -0.0796 | 0.0252 | 0.002 |
| M. franciscanus | -0.0631 | 0.0252 | 0.012 |
| M. pyrifera | 0.1618 | 0.0245 | < 0.001 |

---

For known limitations, see `docs/ANALYSIS_REVISIONS.md`.

*This summary was auto-generated by the analysis pipeline.*
*For questions, contact Emily Donham or Adrian Stier.*
