# Before/After Results Comparison

**Purpose:** Side-by-side comparison of every number that changed between the V5 manuscript and the current pipeline. Use this to update Results, Abstract, and Discussion.

---

## At a Glance

| | MS V5 | Current Pipeline |
|---|---|---|
| Significant effects | 2 of 9 | **6 of 9** |
| Main text figures | 4 (map, pipeline diagram, mean effects, urchin-kelp scatter) | 4 (map, cascade case studies, mean effects, recovery curves) |
| Supplemental figures | ~6 | 15 (S1-S15) |
| Main tables | 2 | 3 |
| Supplemental tables | ~3 | 9 (S1-S9) |
| Effect sizes (k) | ~144 | 132 |
| MPAs analyzed | 19 | 23 |
| Meta-analysis approach | Joint multi-taxa model | Joint multi-taxa model (same structure, improved random effects) |
| Outlier detection | Joint Cook's D (flagged 62% of data) | No removal (primary); Cook's D as sensitivity (Figs S12-S15) |
| P-value correction | None | Benjamini-Hochberg FDR |

---

## Results That Flipped Significance

### Now SIGNIFICANT (weren't before)

| Taxa | Metric | V5 Estimate | V5 p-value | → | Pipeline Estimate | Pipeline p (FDR) |
|------|--------|------------|-----------|---|------------------|-----------------|
| S. pulcher (sheephead) | Biomass | — | p = 0.56 | → | **+185% (RR = 2.85)** | **p < 0.001 (FDR = 0.001)** |
| P. interruptus (lobster) | Density | — | p = 0.047 (marginal) | → | **+196% (RR = 2.96)** | **p = 0.004 (FDR = 0.013)** |
| M. franciscanus (red urchin) | Biomass | — | non-sig | → | **+84% (RR = 1.84)** | **p = 0.015 (FDR = 0.026)** |

### Now NON-SIGNIFICANT (were before)

| Taxa | Metric | V5 Estimate | V5 p-value | → | Pipeline Estimate | Pipeline p (FDR) |
|------|--------|------------|-----------|---|------------------|-----------------|
| M. franciscanus (red urchin) | Density | ~-51% | p = 0.049 | → | -33% (RR = 0.67) | p = 0.284 (FDR = 0.365) |
| P. interruptus (lobster) | Biomass | — | p = 0.047 | → | +29% (RR = 1.29) | p = 0.429 (FDR = 0.483) |

### Still significant (direction unchanged)

| Taxa | Metric | V5 Estimate | → | Pipeline Estimate | Pipeline p (FDR) |
|------|--------|------------|---|------------------|-----------------|
| S. purpuratus (purple urchin) | Biomass | ~-58% | → | **-49% (RR = 0.51)** | **p = 0.014 (FDR = 0.026)** |
| S. purpuratus (purple urchin) | Density | ~-78% | → | **-77% (RR = 0.23)** | **p < 0.001 (FDR = 0.001)** |
| M. pyrifera (kelp) | Biomass | ~+266% | → | **+79% (RR = 1.79)** | **p = 0.032 (FDR = 0.048)** |

### Still non-significant

| Taxa | Metric | Pipeline Estimate | Pipeline p (FDR) |
|------|--------|------------------|-----------------|
| S. pulcher (sheephead) | Density | -15% (RR = 0.85) | p = 0.650 (FDR = 0.650) |

---

## Full Current Results (Table 2 values)

### Biomass

| Taxa | k | lnRR | SE | RR | % Change | p (FDR) | 95% CI (RR) | Sig? |
|------|---|------|----|----|----------|---------|-------------|------|
| S. purpuratus | 10 | -0.669 | 0.264 | 0.51 | -49% | 0.026 | [0.30, 0.87] | **YES** |
| M. franciscanus | 14 | +0.610 | 0.242 | 1.84 | +84% | 0.026 | [1.13, 2.99] | **YES** |
| M. pyrifera | 17 | +0.584 | 0.265 | 1.79 | +79% | 0.048 | [1.05, 3.05] | **YES** |
| P. interruptus | 5 | +0.254 | 0.319 | 1.29 | +29% | 0.483 | [0.68, 2.44] | no |
| S. pulcher | 18 | +1.048 | 0.248 | 2.85 | +185% | 0.001 | [1.74, 4.68] | **YES** |

### Density

| Taxa | k | lnRR | SE | RR | % Change | p (FDR) | 95% CI (RR) | Sig? |
|------|---|------|----|----|----------|---------|-------------|------|
| S. purpuratus | 13 | -1.464 | 0.373 | 0.23 | -77% | 0.001 | [0.11, 0.49] | **YES** |
| M. franciscanus | 16 | -0.394 | 0.365 | 0.67 | -33% | 0.365 | [0.33, 1.40] | no |
| P. interruptus | 17 | +1.086 | 0.366 | 2.96 | +196% | 0.013 | [1.43, 6.16] | **YES** |
| S. pulcher | 22 | -0.166 | 0.363 | 0.85 | -15% | 0.650 | [0.41, 1.75] | no |

---

## Abstract: Line-by-Line Revision Guide

| V5 Claim | What to Write Instead |
|----------|----------------------|
| "we did not detect effects of MPA implementation on predator density and biomass" | Sheephead biomass +185% (p < 0.001) and lobster density +196% (p = 0.004) — **significant predator responses detected** |
| "sea urchin densities declined" | **Purple** urchin density declined 77% (p < 0.001); red urchin density non-significant (-33%, p = 0.28) |
| "kelp ~266% higher" | Kelp +79% higher (p = 0.032) — same direction, smaller magnitude in joint model |

## Results Section: Key Sentences to Replace

| V5 Statement | Replacement |
|-------------|-------------|
| "We failed to detect a significant effect on P. interruptus or S. pulcher biomass or density" | S. pulcher biomass increased 185% (RR = 2.85, p < 0.001) and P. interruptus density increased 196% (RR = 2.96, p = 0.004) inside MPAs |
| "significant effect on densities of BOTH purple and red sea urchins" | Only purple urchin density declined significantly (-77%, p < 0.001); red urchin density was non-significant (-33%, p = 0.28) |
| "densities were ~78% and 51% lower for purple and red urchins" | Purple urchin density was 77% lower (p < 0.001); red urchin density was 33% lower but non-significant (p = 0.28) |
| (no V5 equivalent) | Red urchin **biomass** increased 84% (RR = 1.84, p = 0.015) — urchins grow larger inside MPAs even as density declines, consistent with reduced harvest pressure |
| "19 MPAs" | 23 MPAs |

## Discussion: Key Claims to Revise

| V5 Claim | Revision |
|----------|---------|
| "we did not detect an increase in density or biomass of key predatory species" | 2 of 4 predator metrics are significant (sheephead biomass, lobster density) — discuss predator recovery |
| Urchin decline claims (both species) | Specify only S. purpuratus density; M. franciscanus effects are mixed (biomass up, density non-sig) |
| (new) | Discuss red urchin biomass increase: +84% inside MPAs suggests fishing release on individual size, even as density effects are non-significant |

---

## Figure Reference Mapping (V5 → Current)

| When V5 says... | Change to... | Notes |
|----------------|-------------|-------|
| Fig 1 | Fig 1 | Same concept (map) but inset time series removed |
| Fig 2 | Fig S1 | Pipeline diagram demoted to supplemental |
| Fig 3 | Fig 3 | Mean effects — same concept, different number path |
| Fig 4 | Table 3 | Urchin-kelp scatter replaced by cross-taxa meta-regressions |
| — | Fig 2 (NEW) | Cascade case studies: 3x3 before/after grid |
| — | Fig 4 (NEW) | Recovery trajectories: lmer predictions with 95% CI |
| Fig S1 | Fig S2 | Forest plot renumbered |

## Table Reference Mapping

| When V5 says... | Change to... |
|----------------|-------------|
| Table 1 (meta-analysis) | Table 2 (Table 1 is now average responses) |
| Table 2 | Table 2 (same) |
| Table 3 "linear models" | Table 3 (cross-taxa meta-regression — different analysis) |

---

## New Table 3: Trophic Cascade Meta-Regressions

This table is new — V5 had scatter plots (old Fig 4) instead.

| Relationship | k | Slope | p | Sig? |
|-------------|---|-------|---|------|
| S. purpuratus density → M. pyrifera biomass | 11 | +0.169 | 0.683 | no |
| M. franciscanus density → M. pyrifera biomass | 10 | -0.035 | 0.891 | no |
| P. interruptus density → S. purpuratus density | 11 | +0.154 | 0.707 | no |
| S. pulcher density → S. purpuratus density | 11 | -0.095 | 0.890 | no |
| **P. interruptus biomass → S. purpuratus biomass** | 4 | **-0.420** | **0.025** | **YES** |
| S. pulcher biomass → S. purpuratus biomass | 9 | +0.290 | 0.081 | no |

Note: Cross-taxa meta-regressions have limited power (k=4-11) and are attenuated by errors-in-variables bias (estimated effect sizes used as moderators). See known limitations in `docs/ANALYSIS_REVISIONS.md` Section 8.

---

*Generated from pipeline output (joint multilevel rma.mv model). All p-values are FDR-corrected (Benjamini-Hochberg).*
