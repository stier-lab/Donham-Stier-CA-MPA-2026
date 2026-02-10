# Results Summary: CA MPA Kelp Forest Analysis

**Generated:** 2026-02-09 14:38:31
**Pipeline version:** Modular pBACIPS v2.0

---

## Overview

This document summarizes the key results from the California MPA kelp forest
pBACIPS analysis. Effect sizes represent the log response ratio (ln[MPA/Reference])
at t=11 years post-MPA implementation.

**Interpretation:**
- Positive effect size = higher values inside MPA vs reference
- Negative effect size = lower values inside MPA vs reference
- Significance threshold: p < 0.05

---

## Meta-Analysis Results (Table 2)

### Biomass

| Taxa | k | Estimate | SE | t | p-value | p (FDR) | 95% CI | Effect | Flag |
|------|---|----------|----|----|---------|---------|--------|--------|------|
| S. purpuratus | 3 | -0.5161 | 0.4426 | -1.1662 | 0.2530 | 0.4112 | [-1.4213, 0.389] | negative | preliminary (k<5) |
| M. franciscanus | 2 | 0.4807 | 0.475 | 1.0121 | 0.3198 | 0.4112 | [-0.4907, 1.4521] | positive | preliminary (k<5) |
| M. pyrifera | 17 | 0.543 | 0.3721 | 1.4592 | 0.1552 | 0.4112 | [-0.2181, 1.3041] | positive |  |
| P. interruptus | 2 | 0.5861 | 0.5633 | 1.0404 | 0.3067 | 0.4112 | [-0.566, 1.7382] | positive | preliminary (k<5) |
| S. pulcher | 10 | 1.1657 | 0.3738 | 3.1182 | 0.0041** | 0.0185 | [0.4011, 1.9302] | positive |  |

### Density

| Taxa | k | Estimate | SE | t | p-value | p (FDR) | 95% CI | Effect | Flag |
|------|---|----------|----|----|---------|---------|--------|--------|------|
| S. purpuratus | 6 | -2.4015 | 0.5307 | -4.5254 | 0.0004*** | 0.0036 | [-3.5327, -1.2704] | negative |  |
| M. franciscanus | 1 | -0.4432 | 0.7477 | -0.5928 | 0.5621 | 0.5621 | [-2.0369, 1.1504] | negative | preliminary (k<5) |
| P. interruptus | 2 | 0.7315 | 0.5839 | 1.2528 | 0.2294 | 0.4112 | [-0.513, 1.9761] | positive | preliminary (k<5) |
| S. pulcher | 10 | -0.3244 | 0.5124 | -0.6331 | 0.5362 | 0.5621 | [-1.4166, 0.7678] | negative |  |

*Significance: \*p<0.05, \*\*p<0.01, \*\*\*p<0.001 (uncorrected). p (FDR) = Benjamini-Hochberg adjusted p-values.*

### Heterogeneity Statistics

| Response | tau2 (MPA) | tau2 (Source) |
|----------|------------|---------------|
| Biomass | 0.1416 | 0.2388 |
| Density | 0.0742 | 0.5131 |

---

## Appendix: Replicate-Level Effect Sizes

Effect sizes for each MPA-taxa-response combination.

### S. purpuratus

| MPA | Response | Source | Effect Size | SE | 95% CI | Result |
|-----|----------|--------|-------------|-------|--------|--------|
| Abalone Cove SMCA | Bio | PISCO | 0.2129 | 0.471 | [-0.8133, 1.2391] | Not significant |
| Anacapa Island SMR 2003 | Bio | KFM | -2.2454 | 0.8248 | [-4.0035, -0.4874] | Significant decrease |
| Campus Point SMCA | Bio | PISCO | -1.7674 | 0.6992 | [-3.3254, -0.2094] | Significant decrease |
| Gull Island SMR | Bio | PISCO | -1.5178 | 0.4907 | [-2.5778, -0.4577] | Significant decrease |
| Harris Point SMR | Bio | PISCO | -1.0384 | 0.2412 | [-1.5557, -0.5212] | Significant decrease |
| Naples SMCA | Bio | PISCO | -0.7267 | 0.1586 | [-1.0723, -0.3812] | Significant decrease |
| Point Vicente SMCA | Bio | PISCO | -0.4846 | 0.2893 | [-1.1096, 0.1403] | Not significant |
| Scorpion SMR | Bio | PISCO | -0.9521 | 0.3822 | [-1.7623, -0.1419] | Significant decrease |
| South Point SMR | Bio | PISCO | -1.3862 | 0.532 | [-2.5454, -0.2271] | Significant decrease |
| South Point SMR | Bio | KFM | -0.9506 | 0.2888 | [-1.5661, -0.3351] | Significant decrease |
| Abalone Cove SMCA | Den | PISCO | -2.8135 | 1.1871 | [-5.4262, -0.2007] | Significant decrease |
| Anacapa Island SMR 2003 | Den | KFM | -2.6173 | 0.8097 | [-4.3432, -0.8914] | Significant decrease |
| Anacapa Island SMR 2003 | Den | PISCO | -2.2179 | 0.403 | [-3.0822, -1.3535] | Significant decrease |
| Campus Point SMCA | Den | PISCO | -2.5823 | 0.5035 | [-3.6793, -1.4852] | Significant decrease |
| Gull Island SMR | Den | PISCO | -1.6996 | 0.4605 | [-2.6873, -0.7119] | Significant decrease |
| Harris Point SMR | Den | PISCO | -0.388 | 0.2193 | [-0.8528, 0.0768] | Not significant |
| Matlahuayl SMR | Den | PISCO | -2.3777 | 0.4292 | [-3.4811, -1.2743] | Significant decrease |
| Naples SMCA | Den | PISCO | -0.576 | 0.1503 | [-0.9006, -0.2514] | Significant decrease |
| Point Vicente SMCA | Den | PISCO | -0.3472 | 0.3731 | [-1.1475, 0.4531] | Not significant |
| Scorpion SMR | Den | PISCO | -1.4648 | 0.3309 | [-2.1599, -0.7697] | Significant decrease |
| South Point SMR | Den | KFM | -1.3069 | 0.2597 | [-1.8604, -0.7533] | Significant decrease |
| South Point SMR | Den | PISCO | -1.4105 | 0.3935 | [-2.2543, -0.5666] | Significant decrease |

### M. franciscanus

| MPA | Response | Source | Effect Size | SE | 95% CI | Result |
|-----|----------|--------|-------------|-------|--------|--------|
| Abalone Cove SMCA | Bio | PISCO | -1.6149 | 0.4914 | [-2.6965, -0.5333] | Significant decrease |
| Anacapa Island SMR 2003 | Bio | KFM | -2.6218 | 0.436 | [-3.5512, -1.6925] | Significant decrease |
| Campus Point SMCA | Bio | PISCO | -0.2588 | 0.3803 | [-1.0959, 0.5783] | Not significant |
| Gull Island SMR | Bio | KFM | 1.146 | 0.3909 | [0.3516, 1.9404] | Significant increase |
| Gull Island SMR | Bio | PISCO | 0.247 | 0.2923 | [-0.3799, 0.874] | Not significant |
| Harris Point SMR | Bio | PISCO | 0.0103 | 0.1497 | [-0.3107, 0.3313] | Not significant |
| Naples SMCA | Bio | LTER | 1.0984 | 0.2057 | [0.6402, 1.5567] | Significant increase |
| Naples SMCA | Bio | PISCO | 1.2056 | 0.2469 | [0.6677, 1.7436] | Significant increase |
| Point Vicente SMCA | Bio | PISCO | 0.2549 | 0.126 | [-0.0174, 0.5272] | Not significant |
| Scorpion SMR | Bio | PISCO | 0.0656 | 0.1477 | [-0.2461, 0.3773] | Not significant |
| South Point SMR | Bio | PISCO | 1.2023 | 0.2143 | [0.7393, 1.6652] | Significant increase |
| South Point SMR | Bio | KFM | 1.2126 | 0.0961 | [1.009, 1.4163] | Significant increase |
| Abalone Cove SMCA | Den | PISCO | -2.4162 | 0.3481 | [-3.1824, -1.65] | Significant decrease |
| Anacapa Island SMR 2003 | Den | KFM | -2.7358 | 0.4302 | [-3.6528, -1.8188] | Significant decrease |
| Anacapa Island SMR 2003 | Den | PISCO | -1.2279 | 0.258 | [-1.7813, -0.6745] | Significant decrease |
| Campus Point SMCA | Den | PISCO | -1.0448 | 0.3686 | [-1.8478, -0.2417] | Significant decrease |
| Gull Island SMR | Den | PISCO | -0.8908 | 0.3627 | [-1.6688, -0.1129] | Significant decrease |
| Harris Point SMR | Den | KFM | -0.5323 | 0.3681 | [-1.2775, 0.2129] | Not significant |
| Harris Point SMR | Den | PISCO | 0.5331 | 0.1157 | [0.2879, 0.7784] | Significant increase |
| Matlahuayl SMR | Den | PISCO | -6.323 | 1.0662 | [-9.0636, -3.5823] | Significant decrease |
| Naples SMCA | Den | PISCO | 1.0694 | 0.1653 | [0.7123, 1.4266] | Significant increase |
| Point Vicente SMCA | Den | PISCO | 0.2004 | 0.1234 | [-0.0642, 0.465] | Not significant |
| Scorpion SMR | Den | PISCO | -0.4569 | 0.1575 | [-0.7877, -0.1261] | Significant decrease |
| South Point SMR | Den | PISCO | 0.5759 | 0.1842 | [0.1833, 0.9684] | Significant increase |
| South Point SMR | Den | KFM | 0.6094 | 0.097 | [0.4037, 0.8152] | Significant increase |

### M. pyrifera

| MPA | Response | Source | Effect Size | SE | 95% CI | Result |
|-----|----------|--------|-------------|-------|--------|--------|
| Abalone Cove SMCA | Bio | PISCO | 0.0929 | 0.4395 | [-0.8646, 1.0504] | Not significant |
| Abalone Cove SMCA | Bio | Landsat | 3.0387 | 0.7799 | [1.4538, 4.6236] | Significant increase |
| Anacapa Island SMR 2003 | Bio | PISCO | 1.536 | 0.4059 | [0.6708, 2.4012] | Significant increase |
| Anacapa Island SMR 2003 | Bio | KFM | 0.1113 | 1.8342 | [-3.777, 3.9995] | Not significant |
| Cabrillo SMR | Bio | Landsat | -0.1213 | 1.8258 | [-3.8405, 3.5978] | Not significant |
| Campus Point SMCA | Bio | LTER | 1.7789 | 0.5273 | [0.6822, 2.8755] | Significant increase |
| Campus Point SMCA | Bio | PISCO | 1.6169 | 0.6343 | [0.2348, 2.9989] | Significant increase |
| Campus Point SMCA | Bio | Landsat | 0.4308 | 0.4877 | [-0.5593, 1.4209] | Not significant |
| Carrington Pt SMR | Bio | Landsat | 2.5496 | 0.6268 | [1.2771, 3.8222] | Significant increase |
| Farnsworth Onshore SMCA | Bio | Landsat | -1.9067 | 0.9737 | [-3.9261, 0.1127] | Not significant |
| Gull Island SMR | Bio | PISCO | 1.6757 | 0.5108 | [0.58, 2.7713] | Significant increase |
| Gull Island SMR | Bio | KFM | 5.2348 | 1.385 | [2.41, 8.0596] | Significant increase |
| Gull Island SMR | Bio | Landsat | 0.8811 | 0.2825 | [0.3076, 1.4546] | Significant increase |
| Harris Point SMR | Bio | PISCO | -0.5664 | 0.3456 | [-1.2991, 0.1663] | Not significant |
| Harris Point SMR | Bio | Landsat | 0.7468 | 0.5278 | [-0.3259, 1.8194] | Not significant |
| Matlahuayl SMR | Bio | PISCO | 1.962 | 0.5776 | [0.5487, 3.3754] | Significant increase |
| Naples SMCA | Bio | LTER | 3.7129 | 1.0927 | [1.4405, 5.9854] | Significant increase |
| Naples SMCA | Bio | Landsat | -0.0922 | 0.7464 | [-1.6074, 1.423] | Not significant |
| Naples SMCA | Bio | PISCO | -0.7152 | 0.2081 | [-1.1647, -0.2656] | Significant decrease |
| Point Dume SMCA | Bio | Landsat | 0.2394 | 0.4933 | [-0.762, 1.2408] | Not significant |
| Point Dume SMR | Bio | Landsat | 1.0933 | 1.0334 | [-1.0046, 3.1913] | Not significant |
| Point Vicente SMCA | Bio | Landsat | 1.8299 | 0.7655 | [0.2725, 3.3873] | Significant increase |
| Point Vicente SMCA | Bio | PISCO | -0.5287 | 0.3316 | [-1.2399, 0.1826] | Not significant |
| Santa Barbara Island SMR | Bio | KFM | 3.9147 | 1.1253 | [1.6302, 6.1991] | Significant increase |
| Santa Barbara Island SMR | Bio | Landsat | 2.0391 | 0.4963 | [1.0316, 3.0466] | Significant increase |
| Scorpion SMR | Bio | PISCO | 2.568 | 0.6118 | [1.2827, 3.8534] | Significant increase |
| Skunk Pt SMR | Bio | Landsat | 0.7773 | 0.611 | [-0.4644, 2.0191] | Not significant |
| South La Jolla SMR | Bio | Landsat | 3.0581 | 1.4819 | [0.0466, 6.0697] | Significant increase |
| South Point SMR | Bio | PISCO | 0.5977 | 0.2879 | [-0.016, 1.2113] | Not significant |
| South Point SMR | Bio | KFM | 3.7711 | 1.475 | [0.6272, 6.9149] | Significant increase |
| South Point SMR | Bio | Landsat | 0.0178 | 0.2574 | [-0.5048, 0.5403] | Not significant |
| Swamis SMCA | Bio | Landsat | -1.8229 | 1.0917 | [-4.0391, 0.3933] | Not significant |

### P. interruptus

| MPA | Response | Source | Effect Size | SE | 95% CI | Result |
|-----|----------|--------|-------------|-------|--------|--------|
| Abalone Cove SMCA | Bio | PISCO | -1.674 | 0.5994 | [-2.9933, -0.3547] | Significant decrease |
| Campus Point SMCA | Bio | PISCO | 2.9683 | 0.9072 | [0.9917, 4.945] | Significant increase |
| Campus Point SMCA | Bio | LTER | 1.6562 | 0.7007 | [0.095, 3.2174] | Significant increase |
| Matlahuayl SMR | Bio | PISCO | 0.2317 | 0.6473 | [-1.3521, 1.8155] | Not significant |
| Naples SMCA | Bio | LTER | 0.8481 | 0.3136 | [0.1578, 1.5383] | Significant increase |
| Naples SMCA | Bio | PISCO | 1.1476 | 0.4302 | [0.2183, 2.0769] | Significant increase |
| Point Vicente SMCA | Bio | PISCO | 0.3924 | 0.4786 | [-0.661, 1.4459] | Not significant |
| Abalone Cove SMCA | Den | PISCO | -0.8885 | 0.4568 | [-1.8838, 0.1068] | Not significant |
| Anacapa Island SMR 2003 | Den | KFM | 3.3458 | 0.3004 | [2.7089, 3.9827] | Significant increase |
| Anacapa Island SMR 2003 | Den | PISCO | 2.0849 | 0.2378 | [1.578, 2.5918] | Significant increase |
| Campus Point SMCA | Den | PISCO | 2.3747 | 0.7748 | [0.6865, 4.0629] | Significant increase |
| Campus Point SMCA | Den | LTER | -0.0553 | 0.2314 | [-0.5645, 0.454] | Not significant |
| Gull Island SMR | Den | PISCO | 2.1034 | 0.3673 | [1.3206, 2.8862] | Significant increase |
| Harris Point SMR | Den | PISCO | 0.6936 | 0.2641 | [0.1307, 1.2565] | Significant increase |
| Matlahuayl SMR | Den | PISCO | -0.1091 | 0.3758 | [-1.0287, 0.8105] | Not significant |
| Naples SMCA | Den | LTER | 0.5586 | 0.2986 | [-0.0986, 1.2158] | Not significant |
| Naples SMCA | Den | PISCO | 0.6104 | 0.2984 | [-0.0342, 1.255] | Not significant |
| Point Vicente SMCA | Den | PISCO | -0.288 | 0.3857 | [-1.1154, 0.5393] | Not significant |
| Santa Barbara Island SMR | Den | KFM | 1.9917 | 0.3406 | [1.3021, 2.6812] | Significant increase |
| Scorpion SMR | Den | PISCO | 2.6353 | 0.256 | [2.0995, 3.1712] | Significant increase |
| Scorpion SMR | Den | KFM | 2.5971 | 0.3072 | [1.9747, 3.2196] | Significant increase |
| South Point SMR | Den | PISCO | 2.5756 | 0.4202 | [1.6745, 3.4768] | Significant increase |
| South Point SMR | Den | KFM | 0.8924 | 0.2611 | [0.3388, 1.446] | Significant increase |

### S. pulcher

| MPA | Response | Source | Effect Size | SE | 95% CI | Result |
|-----|----------|--------|-------------|-------|--------|--------|
| Abalone Cove SMCA | Bio | PISCO | 2.3533 | 0.9139 | [0.362, 4.3446] | Significant increase |
| Anacapa Island SMR 2003 | Bio | PISCO | 0.3783 | 0.1302 | [0.1047, 0.6518] | Significant increase |
| Blue Cavern Onshore SMCA | Bio | PISCO | 1.6715 | 0.5452 | [0.3822, 2.9608] | Significant increase |
| Campus Point SMCA | Bio | LTER | -2.2354 | 1.6333 | [-5.6423, 1.1716] | Not significant |
| Campus Point SMCA | Bio | PISCO | 0.0886 | 0.2769 | [-0.5052, 0.6824] | Not significant |
| Cat Harbor SMCA | Bio | PISCO | 0.663 | 0.1526 | [0.3023, 1.0238] | Significant increase |
| Dana Point SMCA | Bio | PISCO | -0.1 | 0.7991 | [-2.0554, 1.8554] | Not significant |
| Farnsworth Onshore SMCA | Bio | PISCO | 0.3805 | 0.2431 | [-0.2443, 1.0053] | Not significant |
| Gull Island SMR | Bio | PISCO | 0.6297 | 0.1425 | [0.3303, 0.929] | Significant increase |
| Harris Point SMR | Bio | PISCO | 2.5255 | 0.6108 | [1.2307, 3.8203] | Significant increase |
| Long Point SMR | Bio | PISCO | 1.5277 | 0.3739 | [0.6436, 2.4119] | Significant increase |
| Matlahuayl SMR | Bio | PISCO | 0.3191 | 0.4515 | [-0.7856, 1.4239] | Not significant |
| Naples SMCA | Bio | PISCO | 2.0512 | 0.2049 | [1.6118, 2.4906] | Significant increase |
| Point Vicente SMCA | Bio | PISCO | 1.2727 | 0.38 | [0.4673, 2.0782] | Significant increase |
| Santa Barbara Island SMR | Bio | PISCO | 0.566 | 0.1302 | [0.2657, 0.8663] | Significant increase |
| Scorpion SMR | Bio | PISCO | 0.5102 | 0.0681 | [0.3676, 0.6529] | Significant increase |
| South Point SMR | Bio | PISCO | 0.8413 | 0.1074 | [0.6165, 1.0662] | Significant increase |
| Swamis SMCA | Bio | PISCO | 3.0836 | 0.7764 | [1.1839, 4.9834] | Significant increase |
| Abalone Cove SMCA | Den | PISCO | 1.4808 | 0.5818 | [0.2131, 2.7485] | Significant increase |
| Anacapa Island SMR 2003 | Den | KFM | 1.2379 | 0.4273 | [0.327, 2.1487] | Significant increase |
| Anacapa Island SMR 2003 | Den | PISCO | -0.6778 | 0.2673 | [-1.2417, -0.1139] | Significant decrease |
| Blue Cavern Onshore SMCA | Den | PISCO | -1.8543 | 0.7458 | [-3.6793, -0.0293] | Significant decrease |
| Campus Point SMCA | Den | LTER | -0.5359 | 0.6937 | [-1.983, 0.9111] | Not significant |
| Campus Point SMCA | Den | PISCO | -0.0571 | 0.2412 | [-0.5745, 0.4602] | Not significant |
| Cat Harbor SMCA | Den | PISCO | 0.4123 | 0.1889 | [-0.0342, 0.8589] | Not significant |
| Dana Point SMCA | Den | PISCO | -0.3392 | 0.5732 | [-1.7417, 1.0633] | Not significant |
| Farnsworth Onshore SMCA | Den | PISCO | 0.0311 | 0.0948 | [-0.2125, 0.2748] | Not significant |
| Gull Island SMR | Den | PISCO | 0.3184 | 0.0959 | [0.117, 0.5199] | Significant increase |
| Harris Point SMR | Den | PISCO | 1.5312 | 0.2686 | [0.9619, 2.1005] | Significant increase |
| Harris Point SMR | Den | KFM | -1.1723 | 0.5113 | [-2.2103, -0.1343] | Significant decrease |
| Long Point SMR | Den | PISCO | 0.6287 | 0.2231 | [0.1013, 1.1562] | Significant increase |
| Matlahuayl SMR | Den | PISCO | 0.1079 | 0.3409 | [-0.7263, 0.9421] | Not significant |
| Naples SMCA | Den | PISCO | 1.6366 | 0.157 | [1.2999, 1.9734] | Significant increase |
| Point Vicente SMCA | Den | PISCO | 0.3289 | 0.1083 | [0.1003, 0.5575] | Significant increase |
| Santa Barbara Island SMR | Den | PISCO | -0.1035 | 0.1049 | [-0.3409, 0.1339] | Not significant |
| Scorpion SMR | Den | PISCO | 0.0884 | 0.0683 | [-0.0546, 0.2315] | Not significant |
| South Point SMR | Den | KFM | 0.2348 | 0.1109 | [-2e-04, 0.4699] | Not significant |
| South Point SMR | Den | PISCO | 0.5243 | 0.0896 | [0.3368, 0.7119] | Significant increase |
| Swamis SMCA | Den | PISCO | 2.139 | 0.6512 | [0.5457, 3.7323] | Significant increase |

---

## Summary Statistics

### Effect Direction Summary

- **Significant increases inside MPA:** 61
- **Significant decreases inside MPA:** 32
- **Not significant:** 48
- **Total effects analyzed:** 141

### Summary by Taxa

| Taxa | n | Significant | Mean Effect |
|------|---|-------------|-------------|
| S. purpuratus | 22 | 18 (82%) | -1.394 |
| M. franciscanus | 25 | 18 (72%) | -0.428 |
| M. pyrifera | 32 | 17 (53%) | 1.235 |
| P. interruptus | 23 | 15 (65%) | 1.161 |
| S. pulcher | 39 | 25 (64%) | 0.577 |

---

*This summary was auto-generated by the analysis pipeline.*
*For questions, contact Emily Donham or Adrian Stier.*
