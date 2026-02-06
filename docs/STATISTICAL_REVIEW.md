# Statistical Methods Documentation

## CA MPA Kelp Forest pBACIPS Analysis

**Document Version:** 2.0
**Last Updated:** 2026-02-06
**Authors:** Emily Donham & Adrian Stier
**Target Journal:** Conservation Letters

---

## Table of Contents

1. [Methods Summary for Manuscript](#1-methods-summary-for-manuscript)
2. [Glossary of Key Terms](#2-glossary-of-key-terms)
3. [Statistical Equations](#3-statistical-equations)
4. [Methodological Justifications](#4-methodological-justifications)
5. [Diagnostic Summary](#5-diagnostic-summary)
6. [Methods Checklist for Reviewers](#6-methods-checklist-for-reviewers)
7. [Implementation Verification](#7-implementation-verification)
8. [Known Limitations](#8-known-limitations)
9. [References](#9-references)

---

## 1. Methods Summary for Manuscript

### Copy-Paste Ready Text

> **Statistical Analysis**
>
> We used progressive-change Before-After-Control-Impact-Pairs with Symmetry (pBACIPS) analysis to estimate MPA effects on kelp forest taxa (Thiault et al. 2017). For each MPA, we paired observations from protected sites with nearby unprotected reference sites and calculated log response ratios (lnRR) as the difference in log-transformed proportional biomass or density between MPA and reference sites at each sampling time.
>
> We fitted four candidate temporal models to the difference between MPA and reference sites (delta): (1) a step model assuming abrupt change at MPA implementation, (2) a linear model allowing gradual constant-rate change, (3) an asymptotic Michaelis-Menten model representing rapid initial change that saturates, and (4) a sigmoid Hill function model representing S-shaped change. We selected among models using AICc weights (Burnham & Anderson 2002). Effect sizes were calculated as the predicted difference at t = 11 years post-implementation to ensure comparability across MPAs with different protection durations.
>
> We aggregated effect sizes across MPAs using multilevel meta-analysis with restricted maximum likelihood estimation (REML) implemented in the metafor package (Viechtbauer 2010). Models included taxa as a fixed-effect moderator and MPA and data source as crossed random effects. We applied the Knapp-Hartung adjustment for improved inference with small sample sizes (Knapp & Hartung 2003). Influential outliers were identified using Cook's distance (threshold: 4/n) and removed before final model fitting. Heterogeneity was quantified using tau-squared variance components and pseudo-I-squared statistics.
>
> To handle zeros in proportional data before log transformation, we applied an adaptive pseudocount correction using half the minimum non-zero proportion observed in each MPA-taxa combination (Aitchison 1986; Martin-Fernandez et al. 2003). Biomass was estimated from size-frequency data using 1,000 bootstrap resamples. Standard errors for effect sizes were obtained using proper covariance-aware contrasts via the emmeans package (Lenth 2021).

### Condensed Version (for space-limited journals)

> We used pBACIPS analysis to estimate MPA effects (Thiault et al. 2017), fitting step, linear, asymptotic, and sigmoid models to log response ratios and selecting among them using AICc. Effect sizes were aggregated via multilevel meta-analysis (REML with Knapp-Hartung adjustment) using MPA and data source as random effects. Zeros were handled using adaptive pseudocount correction (Aitchison 1986), and influential outliers were identified via Cook's distance.

---

## 2. Glossary of Key Terms

| Term | Definition | Mathematical Notation |
|------|------------|----------------------|
| **pBACIPS** | Progressive-change Before-After-Control-Impact-Pairs with Symmetry. A statistical design that allows MPA effects to develop gradually over time rather than assuming an instantaneous step change. | - |
| **lnRR** | Log response ratio. The natural logarithm of the ratio of a response variable (e.g., biomass) between treatment (MPA) and control (reference) sites. | lnRR = ln(MPA/Reference) |
| **Delta (delta)** | The difference between MPA and reference site values at each time point. In this analysis, delta = lnRR(MPA) - lnRR(Reference). | delta_t = y_MPA,t - y_Ref,t |
| **Effect size** | The magnitude of MPA effect, measured as the change in delta from the before period to a standardized post-implementation time point (t = 11 years). | ES = delta_after - delta_before |
| **Sampling variance (vi)** | The squared standard error of an effect size estimate, used as a weight in meta-analysis. | vi = SE^2 |
| **tau-squared** | Between-study variance in meta-analysis; measures heterogeneity not explained by sampling error. | tau^2 |
| **I-squared (pseudo)** | Percentage of total variability due to true heterogeneity rather than sampling error. | I^2 = 100 * tau^2 / (tau^2 + v) |
| **emmeans** | Estimated marginal means. A method for obtaining predictions and contrasts from fitted models that properly accounts for covariance structure. | - |
| **AICc** | Corrected Akaike Information Criterion. A model selection criterion that penalizes complexity, with a correction for small sample sizes. | AICc = AIC + 2K(K+1)/(n-K-1) |
| **Cook's distance** | A diagnostic measuring the influence of each observation on model estimates; values > 4/n indicate influential outliers. | D_i |
| **REML** | Restricted Maximum Likelihood. An estimation method for mixed models that provides unbiased variance component estimates. | - |
| **Knapp-Hartung** | An adjustment to standard errors in meta-analysis that accounts for uncertainty in heterogeneity estimation, using t-distribution instead of normal. | test = "t" in metafor |

---

## 3. Statistical Equations

### 3.1 Log Response Ratio Calculation

The log response ratio quantifies the proportional difference between MPA and reference sites:

```
lnRR = ln(p_MPA / p_Reference)
     = ln(p_MPA) - ln(p_Reference)
```

where `p` is the proportion of total biomass or density for a given taxon.

**Zero-correction (Adaptive Pseudocount):**

```
p_corrected = p + (min(p_nonzero) / 2)
```

This adaptive method (Aitchison 1986) scales the correction to the data, avoiding inflation of effect sizes for rare species.

### 3.2 pBACIPS Candidate Models

All models predict delta (MPA - Reference difference) as a function of time since MPA implementation:

**Step Model (ANOVA):**
```
delta_t = mu_before      if t = 0 (Before period)
delta_t = mu_after       if t > 0 (After period)
```

**Linear Model:**
```
delta_t = beta_0 + beta_1 * t + epsilon_t
```

**Asymptotic Model (Michaelis-Menten):**
```
delta_t = (M * t) / (L + t) + B + epsilon_t

where:
  M = maximum effect (asymptote)
  L = half-saturation time (time to reach M/2)
  B = baseline (before-period mean)
```

**Sigmoid Model (Hill Function):**
```
delta_t = (M * (t/L)^K) / (1 + (t/L)^K) + B + epsilon_t

where:
  M = maximum effect
  L = inflection point (time to half-maximum)
  K = Hill coefficient (steepness)
  B = baseline
```

### 3.3 Model Selection (AICc Weights)

```
AICc = -2 * ln(L) + 2K + (2K(K+1)) / (n-K-1)

delta_AICc_i = AICc_i - min(AICc)

w_i = exp(-0.5 * delta_AICc_i) / sum(exp(-0.5 * delta_AICc_j))
```

where `L` is likelihood, `K` is number of parameters, and `n` is sample size.

### 3.4 Effect Size Extraction

For step models:
```
ES = mean(delta_after) - mean(delta_before)
SE = sqrt(var(delta_after)/n_after + var(delta_before)/n_before)
```

For trend models (linear, asymptotic, sigmoid):
```
ES = predicted(t = 11) - predicted(t = 0)
SE = obtained from emmeans contrast covariance matrix
```

The choice of t = 11 years ensures comparability across MPAs with different establishment dates (see Section 4.3).

### 3.5 Multilevel Meta-Analysis Model

```
y_i = mu + Taxa_j + u_MPA,k + u_Source,l + epsilon_i

where:
  y_i = effect size for observation i
  mu = overall intercept (removed when using ~ Taxa - 1)
  Taxa_j = fixed effect for taxon j
  u_MPA,k ~ N(0, tau^2_MPA) = random effect for MPA k
  u_Source,l ~ N(0, tau^2_Source) = random effect for source l
  epsilon_i ~ N(0, vi) = sampling error with known variance vi
```

**Implementation in R:**
```r
rma.mv(yi = Mean, V = vi,
       mods = ~ Taxa - 1,
       random = list(~ 1 | MPA, ~ 1 | Source),
       data = data, method = "REML", test = "t")
```

---

## 4. Methodological Justifications

### 4.1 Why pBACIPS Instead of Traditional BACI?

**Choice:** We use pBACIPS rather than assuming an instantaneous step change.

**Justification:** Ecological responses to protection rarely occur instantaneously (Thiault et al. 2017). Recovery trajectories may be:
- **Linear:** Constant rate of change as populations rebuild
- **Asymptotic:** Rapid initial recovery that slows as carrying capacity is approached
- **Sigmoidal:** Delayed response followed by rapid change (e.g., threshold dynamics)

The pBACIPS approach allows the data to select among these ecologically meaningful alternatives, improving effect size estimation accuracy.

**Literature Support:**
- Thiault et al. (2017) Methods in Ecology and Evolution 8:288-296
- Claudet et al. (2008) Ecology Letters 11:481-489

### 4.2 Why Adaptive Pseudocount Instead of Fixed Constant?

**Choice:** We use `min(p_nonzero)/2` rather than a fixed constant (e.g., 0.01).

**Justification:**
1. **Scale-appropriate:** The correction scales with observed data, avoiding artificial inflation for rare species
2. **Bias reduction:** A fixed 0.01 can inflate effect sizes by 5-10x for species with low proportions
3. **Compositional data theory:** Recommended by Aitchison (1986) for log-ratio analysis of proportions

**Literature Support:**
- Aitchison, J. (1986) The Statistical Analysis of Compositional Data. Chapman & Hall
- Martin-Fernandez et al. (2003) Mathematical Geology 35:253-278

### 4.3 Why t = 11 Years for Effect Size Extraction?

**Choice:** Effect sizes are calculated at t = 11 years post-implementation for all MPAs.

**Justification:**
1. **Comparability:** Different MPAs have different protection durations; using maximum observed time would conflate MPA age with effect magnitude
2. **Biological relevance:** 11 years represents the age of our youngest established MPA in the dataset (MLPA South Coast MPAs, implemented 2012, analyzed through 2023)
3. **Conservative estimate:** For older MPAs, this may underestimate full effects, but ensures fair comparison

**Alternative considered:** Using maximum observed time per MPA would yield larger effects for older MPAs, confounding protection duration with effect magnitude.

### 4.4 Why Include Source as a Random Effect?

**Choice:** The meta-analysis model includes both MPA and data source (PISCO, KFM, LTER) as random effects.

**Justification:**
1. **Non-independence:** Effect sizes from the same monitoring program share methodological characteristics (survey protocols, observer training, equipment)
2. **Hierarchical structure:** MPAs are nested within geographic regions that different programs cover
3. **Variance partitioning:** Allows separation of MPA-level and program-level heterogeneity

**Caveat:** With only 3-4 source levels, variance component estimates for the Source effect have high uncertainty. Sensitivity analyses comparing models with and without the Source random effect are performed (see `09_meta_analysis.R`).

**Literature Support:**
- Viechtbauer (2010) Journal of Statistical Software 36:1-48
- Konstantopoulos (2011) Research Synthesis Methods 2:61-76

### 4.5 Why Knapp-Hartung Adjustment?

**Choice:** We use `test = "t"` in metafor for t-distribution based inference.

**Justification:**
1. **Small sample correction:** Standard meta-analysis assumes large samples; Knapp-Hartung uses residual-based variance estimation
2. **Proper uncertainty:** Accounts for the fact that heterogeneity (tau-squared) is estimated, not known
3. **Conservative inference:** Produces wider confidence intervals and higher p-values than normal approximation

**Literature Support:**
- Knapp & Hartung (2003) Statistics in Medicine 22:2693-2710
- IntHout et al. (2014) BMC Medical Research Methodology 14:25

### 4.6 Why Cook's Distance for Outlier Detection?

**Choice:** Outliers are identified using Cook's distance > 4/n threshold.

**Justification:**
1. **Influence-based:** Cook's distance measures how much removing an observation changes model estimates, not just how extreme the value is
2. **Standard practice:** The 4/n threshold is conventional in meta-analysis
3. **Transparency:** Outliers are reported before removal so readers can assess impact

**Literature Support:**
- Viechtbauer & Cheung (2010) Research Synthesis Methods 1:112-125

### 4.7 Why Bootstrap for Biomass Estimation?

**Choice:** Biomass is estimated using 1,000 bootstrap resamples from size-frequency distributions.

**Justification:**
1. **Missing individual sizes:** When only aggregate counts are available, individual sizes must be imputed from contemporaneous size-frequency data
2. **Uncertainty quantification:** Bootstrap SE captures variability due to unknown individual sizes
3. **Non-parametric:** Makes no distributional assumptions about size distributions

**Implementation:** For each observation with count N but no individual sizes:
1. Find size-frequency records from same MPA, year, and site status
2. Resample N individuals from pooled size distribution
3. Convert sizes to biomass using allometric equations
4. Repeat 1,000 times; report mean and SE

---

## 5. Diagnostic Summary

### 5.1 Model Diagnostics Performed

| Diagnostic | Method | Threshold | Purpose |
|------------|--------|-----------|---------|
| Residual normality | Shapiro-Wilk test | p > 0.05 | Verify error distribution |
| Residual uniformity | DHARMa simulation | p > 0.05 | Detect systematic patterns |
| Heteroscedasticity | Correlation(abs(resid), fitted) | r < 0.5 | Check constant variance |
| Outliers | DHARMa outlier test | p > 0.05 | Identify extreme values |
| Overdispersion | DHARMa dispersion test | p > 0.05 | Check variance assumptions |
| Model influence | Cook's distance | D < 4/n | Identify influential points |

### 5.2 Diagnostic Results Summary

The following checks are performed automatically in `08_effect_sizes.R`:

1. **Before-period trend test:** Linear regression on pre-MPA data
   - If significant (p < 0.05): Analysis downgraded to BACI (step only)
   - If non-significant: Full pBACIPS model selection proceeds

2. **NLS convergence monitoring:**
   - Multiple starting value strategies attempted (5 per model type)
   - Fallback models available (piecewise linear, quadratic, GAM)
   - Convergence logged for post-hoc review via `get_model_fit_summary()`

3. **DHARMa diagnostics:**
   - Run on all linear models when DHARMa package is available
   - Custom diagnostics for NLS models (residual checks)
   - Results stored in `ModelDiagnostics` dataframe

### 5.3 Heterogeneity Assessment

Meta-analysis heterogeneity is assessed via:

1. **Tau-squared components:**
   - tau^2_MPA: Between-MPA variance
   - tau^2_Source: Between-source variance

2. **Pseudo-I-squared:**
   - < 25%: Low heterogeneity
   - 25-75%: Moderate heterogeneity
   - > 75%: High heterogeneity

3. **Variance component confidence intervals:**
   - Computed using profile likelihood method
   - Exported to `data/table_s_variance_components.csv`

---

## 6. Methods Checklist for Reviewers

Use this checklist to verify methodological rigor:

### Data Processing
- [ ] Raw data from multiple sources (PISCO, KFM, LTER) standardized to common units
- [ ] Biomass conversions use published allometric equations (see Section 7.3)
- [ ] Size cutoffs applied per source protocol (PISCO: 25mm minimum for urchins)
- [ ] Survey areas normalized to per-m2 densities

### BACI Design
- [ ] Each MPA paired with appropriate reference site(s)
- [ ] Before and after periods defined by MPA implementation year
- [ ] Before-period trends tested; confounded sites flagged
- [ ] Log response ratios calculated with adaptive zero-correction

### Effect Size Estimation
- [ ] Four candidate models fitted (step, linear, asymptotic, sigmoid)
- [ ] Model selection via AICc weights (appropriate for small samples)
- [ ] Effect sizes extracted at standardized time point (t = 11 years)
- [ ] Standard errors obtained from emmeans contrasts (covariance-aware)

### Meta-Analysis
- [ ] Biomass and density analyzed separately (non-independent)
- [ ] MPA and Source included as random effects
- [ ] REML estimation with Knapp-Hartung adjustment
- [ ] Cook's distance outlier detection (4/n threshold)
- [ ] Heterogeneity statistics reported (tau^2, pseudo-I^2)

### Sensitivity Analyses
- [ ] Models with/without Source random effect compared
- [ ] Variance component confidence intervals computed
- [ ] Outlier removal impact assessed

### Reproducibility
- [ ] Random seeds set for bootstrap procedures
- [ ] 1,000 bootstrap iterations used throughout
- [ ] All code available in numbered R scripts

---

## 7. Implementation Verification

### 7.1 Code Locations

| Component | File | Key Lines |
|-----------|------|-----------|
| pBACIPS model fitting | `02_pBACIPS_function.R` | 774-1277 |
| Bootstrap biomass | `01_utils.R` | 252-310 |
| Effect size calculation | `08_effect_sizes.R` | 560-681 |
| Meta-analysis | `09_meta_analysis.R` | 140-285 |
| Zero-correction | `01_utils.R` | `calculate_proportions()` |

### 7.2 Verified Implementations

| Component | Status | Notes |
|-----------|--------|-------|
| Sampling variance (vi = SE^2) | CORRECT | Line 86 of 09_meta_analysis.R |
| Multilevel random effects | CORRECT | `random = list(~1|MPA, ~1|Source)` |
| Knapp-Hartung adjustment | CORRECT | `test = "t"` in rma.mv calls |
| Cook's distance threshold | CORRECT | 4/n threshold |
| Bootstrap iterations | CORRECT | 1,000 in all scripts |
| Adaptive pseudocount | CORRECT | Standardized 2026-02-05 |
| Effect size SE from emmeans | CORRECT | Covariance-aware contrasts |

### 7.3 Biomass Conversion Formulas

| Species | Formula | Code Location |
|---------|---------|---------------|
| S. purpuratus | W = 0.00059 x TD^2.870 | `01_utils.R:189` |
| M. franciscanus | W = 0.00059 x TD^2.917 | `01_utils.R:174` |
| S. pulcher | W = 0.0144 x TL^3.04 | `06_lter_processing.R:640` |
| P. interruptus | W = 0.00135 x CL^2.91 | `01_utils.R:159` |
| M. pyrifera | B = Stipes x 0.085 x 1000 | `01_utils.R:211` |

---

## 8. Known Limitations

### 8.1 Acknowledged Limitations

1. **Source random effect uncertainty:** With only 3-4 source levels (PISCO, KFM, LTER, Landsat), variance component estimates for Source have wide confidence intervals. Sensitivity analyses are performed.

2. **Measurement error in predictors:** Cross-taxa regression (Table 3) treats predictor effect sizes as fixed values, ignoring their uncertainty. This may underestimate standard errors.

3. **Spatial pseudoreplication:** Multiple transects within an MPA are averaged before pBACIPS analysis, but some spatial non-independence may remain.

4. **Missing before-period data:** Some MPAs lack pre-implementation data; these are analyzed as Control-Impact only (CI analysis).

### 8.2 Potential Reviewer Concerns

| Concern | Response |
|---------|----------|
| "Why not use a Bayesian meta-analysis?" | Frequentist REML with Knapp-Hartung provides valid inference for our sample sizes. Bayesian approach would require prior specification and may not substantially change conclusions. |
| "How do you handle spatial autocorrelation?" | Reference sites are selected to be nearby but outside MPA boundaries. Within-MPA spatial structure is averaged over transects. |
| "Why separate biomass and density models?" | These metrics from the same surveys are non-independent; pooling would violate meta-analysis assumptions. |
| "What if the best model varies across MPAs?" | This is expected and appropriate - pBACIPS allows each MPA to have its own best-fit trajectory. Effect sizes are comparable because they're extracted at a common time point. |

---

## 9. References

Aitchison, J. (1986). The Statistical Analysis of Compositional Data. Chapman & Hall, London.

Burnham, K.P. & Anderson, D.R. (2002). Model Selection and Multimodel Inference: A Practical Information-Theoretic Approach, 2nd ed. Springer, New York.

Claudet, J., et al. (2008). Marine reserves: size and age do matter. Ecology Letters, 11, 481-489.

IntHout, J., Ioannidis, J.P.A., & Borm, G.F. (2014). The Hartung-Knapp-Sidik-Jonkman method for random effects meta-analysis is straightforward and considerably outperforms the standard DerSimonian-Laird method. BMC Medical Research Methodology, 14, 25.

Knapp, G. & Hartung, J. (2003). Improved tests for a random effects meta-regression with a single covariate. Statistics in Medicine, 22, 2693-2710.

Konstantopoulos, S. (2011). Fixed effects and variance components estimation in three-level meta-analysis. Research Synthesis Methods, 2, 61-76.

Lenth, R.V. (2021). emmeans: Estimated Marginal Means, aka Least-Squares Means. R package version 1.7.0. https://CRAN.R-project.org/package=emmeans

Martin-Fernandez, J.A., Barcelo-Vidal, C., & Pawlowsky-Glahn, V. (2003). Dealing with zeros and missing values in compositional data sets using nonparametric imputation. Mathematical Geology, 35, 253-278.

Thiault, L., Kernaléguen, L., Osenberg, C.W., & Claudet, J. (2017). Progressive-Change BACIPS: a flexible approach for environmental impact assessment. Methods in Ecology and Evolution, 8, 288-296.

Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. Journal of Statistical Software, 36, 1-48.

Viechtbauer, W. & Cheung, M.W.L. (2010). Outlier and influence diagnostics for meta-analysis. Research Synthesis Methods, 1, 112-125.

---

## Appendix: Change Log

### Version 2.0 (2026-02-05)
- Added methods summary for manuscript (Section 1)
- Added glossary of key terms (Section 2)
- Added statistical equations (Section 3)
- Added methodological justifications with literature support (Section 4)
- Added diagnostic summary (Section 5)
- Added reviewer checklist (Section 6)
- Reorganized document structure for clarity
- Standardized adaptive pseudocount method across all processing scripts

### Version 1.0 (2026-02-03)
- Initial statistical review document
- Identified discrepancies between manuscript and pipeline
- Fixed effect size SE calculation (emmeans contrasts)
- Changed prediction to confidence intervals
- Added heterogeneity reporting
