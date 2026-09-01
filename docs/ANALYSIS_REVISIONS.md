**From:** Adrian  |  **Date:** February 27, 2026

> **Historical / superseded.** This memo predates the later heatwave,
> environmental-moderator, SSWD, and resilience modules now integrated in scripts
> 14-23. Treat old figure numbers, table numbers, action items, and "no climate
> covariates" language below as provenance from that revision stage, not current
> repo guidance.

Emily, this walks through every meaningful difference between our original analysis (V5 manuscript / `pBACIPS_PISCO_V10.R`) and the rebuilt pipeline. I cover the stats, results, figures, tables, and repo layout. It's long but worth reading through before we sit down to rewrite the manuscript.

---

## Executive Summary

The big picture: 6 of 9 taxa-response combinations are now significant. I switched to a joint multilevel model on the full dataset, dropped the outlier removal step, and added FDR correction. The story flips from "MPAs don't detectably affect predators" to "MPAs boost predator biomass and density, knock down urchins, and grow kelp."

Compared to V5: sheephead biomass, lobster density, and red urchin biomass **gained** significance. Red urchin density **lost** it. And here's the weird one: red urchin biomass is *positive* (higher inside MPAs). More on that below.

**What didn't change:** The core pBACIPS approach, the data, the study sites, and the finding that kelp goes up inside MPAs.

**What we need to rewrite:** Pretty much every number in Results and the Abstract. The Discussion needs to flip its conclusion about predators. See Section 9 for the full action list.

### At a Glance

| | MS V5 | Current Pipeline |
|---|---|---|
| Significant effects | 2 of 9 | **6 of 9** |
| Main text figures | 4 (map, pipeline diagram, mean effects, urchin-kelp scatter) | 5 (map, cascade case studies, mean effects, recovery curves, kelp resilience) |
| Supplemental figures | ~6 | 15 (S1-S15) |
| Main tables | 2 | 3 |
| Supplemental tables | ~3 | 9 (S1-S9) |
| Effect sizes (k) | ~144 | 142 |
| MPAs analyzed | 19 | 23 |
| Meta-analysis approach | Joint multi-taxa model | Joint multi-taxa model (same structure, improved random effects) |
| Outlier detection | Joint Cook's D (flagged 62% of data) | No removal (primary); Cook's D as sensitivity (Figs S12-S15) |
| P-value correction | None | Benjamini-Hochberg FDR |

---

## 1. Meta-Analysis Model Structure

### What V5 did

V5 fit two joint `rma.mv()` models, one for density and one for biomass (`pBACIPS_PISCO_V10.R`, lines 6941-6970). Both fit all 5 species at once with Taxa as a moderator. The density model used `random = ~1|MPA`. The biomass model used `random = ~1|MPA/Type.x`. After fitting, V5 ran Cook's D outlier detection with hardcoded global thresholds: `4/62` for density (62 observations across all taxa) and `4/79` for biomass. Anything above the threshold got tossed and the model was refit.

### What I changed

**Primary analysis:** Two joint `rma.mv()` models (one biomass, one density) with `mods = ~Taxa - 1` (cell-means), `random = list(~1|MPA, ~1|Source)`, and **no outlier removal**. This borrows strength across taxa for shared variance components and gives us per-taxon estimates. We keep the full dataset, which is best practice (Viechtbauer & Cheung 2010; Noble et al. 2022; Cochrane Handbook).

**Sensitivity:** I also ran (1) per-taxa separate `rma.mv()` models without outlier removal, (2) per-taxa models with Cook's D at 4/k, and (3) joint Cook's D at 4/n (V5's approach). All compared in Table S9.

### Why the joint model makes more sense

Meta-regression with a categorical moderator is the textbook approach (Borenstein & Higgins 2013; Cochrane Handbook Ch. 10; Nakagawa & Santos 2012). Three reasons it works better here:

1. **Borrowing strength:** Shared τ²(MPA) and τ²(Source) come from all data (~56 bio, ~76 den), not per-taxon k=2 to 22. Way more stable variance estimates.
2. **Small-k taxa actually get real df:** P. interruptus biomass (k=5) and S. purpuratus biomass (k=10) get proper degrees of freedom instead of falling back to z-tests.
3. **One model per response type.** Cleaner than fitting 9 separate models and making ad hoc choices about random effects for each one.

### Why V5's Cook's D was a problem

This is the key issue. V5 used taxa as moderators, which means we *expect* different taxa to have different effect sizes. But then it applied Cook's D across all taxa pooled together. So any observation that differed from the *overall* pattern got flagged as an "outlier," even if it was totally normal for its taxon.

Here's what that looked like. The density model pooled all 62 observations (lobster + sheephead + both urchins + kelp) and flagged any with Cook's D > 4/62 = 0.065. Lobster density has large positive lnRR values (~1.3). Urchins have negative ones (~-1.5). So the lobster observations looked extreme relative to the pooled distribution and got flagged. Same thing in the biomass model with kelp (lnRR up to ~9) vs. urchins (~-0.8).

Those aren't outliers. They're the cascade signal. V5 threw out **62%** of the data:

| Taxa | Response | Total k | Kept by V5 Cook's D | Lost |
|------|----------|---------|---------------------|------|
| P. interruptus | Density | 17 | 2 | **15 (88%)** |
| M. pyrifera | Biomass | 17 | 5 | **12 (71%)** |
| S. purpuratus | Biomass | 10 | 3 | **7 (70%)** |
| S. pulcher | Biomass | 18 | 11 | **7 (39%)** |
| P. interruptus | Biomass | 5 | 2 | **3 (60%)** |

### Why keeping the full dataset is the right call

Beyond the Cook's D problem above, there are strong reasons to report the full dataset as primary.

First, the observations V5 removed are real data from real MPAs. Every one came from a standardized monitoring program (PISCO, KFM, or LTER) using established transect protocols. None of them are data-entry errors, equipment failures, or sites with known confounds. Removing them means ignoring legitimate ecological variation across MPAs.

Second, the meta-analysis literature is clear on this. Viechtbauer and Cheung (2010) argue that influential observations in meta-analysis often reflect meaningful biological heterogeneity and that removing them can bias results toward the null. Noble et al. (2022) make the same point for ecology: dropping "outliers" from meta-analyses risks discarding the very signal you're trying to detect. The Cochrane Handbook (Ch. 10) recommends reporting the full dataset as the primary analysis and using outlier removal only as a sensitivity check.

Third, this is a trophic cascade study. We *expect* large effect sizes for predators and kelp. That's the whole hypothesis. Lobster density going up 3x inside MPAs is exactly what MPA theory predicts. Kelp biomass ratios of 5-10x at some sites make sense when reference-site kelp is near zero after urchin barrens. Calling these "outliers" is basically saying the cascade is too strong to be real.

Fourth, the multilevel model already handles variation the right way. The `~1|MPA` random effect accounts for site-to-site differences in effect magnitude. Sites with unusually large effects get partially shrunk toward the overall taxon mean through the hierarchical structure. Inverse-variance weighting downweights imprecise estimates. These are the standard tools for handling heterogeneity, and they work better than throwing data away.

The four-way sensitivity analysis confirms all of this. Significant results hold under approaches 1 through 3 (joint no removal, per-taxa no removal, per-taxa Cook's D at 4/k). Only V5's aggressive joint Cook's D at 4/n (approach 4) kills the predator signals, because it removes the taxa that diverge most from the grand mean (Table S9). Figures S13–S15 provide visual evidence: Figure S13 shows the flagged observations span the full within-taxon distribution and that the global threshold flags 50–82% of observations per taxon (vs. 4–13% for per-taxon thresholds); Figure S14 shows flagged MPA–taxon combinations follow coherent recovery trajectories, not erratic noise; Figure S15 shows the underlying raw data are well-behaved.

---

## 2. Random Effects Structure

| Component | V5 (density) | V5 (biomass) | Current (both) |
|-----------|-------------|--------------|----------------|
| MPA | `~1\|MPA` | `~1\|MPA/Type.x` | `~1\|MPA` |
| Source | — | nested within MPA | `~1\|Source` (crossed) |

V5 used different random-effect structures for density vs biomass (the biomass model nested Source within MPA). Because the same program (PISCO/KFM/LTER) samples multiple MPAs, Source and MPA are crossed rather than nested, so the current pipeline uses crossed random effects for both responses. This is a small change in practice: a sensitivity check shows estimates move less than 10% with vs. without the Source random effect, and significance is identical either way.

---

## 3. Results Comparison

### Summary table

| Taxa | Metric | V5 Estimate | V5 p-value | Current Estimate | Current p (FDR) | Change |
|------|--------|-------------|------------|------------------|-----------------|--------|
| S. purpuratus | Biomass | -0.77 (RR=0.46, -54%) | p=0.002 | **-0.67 (RR=0.51, -49%)** | **p=0.026** | Still significant |
| S. purpuratus | Density | -1.46 (RR=0.23, -77%) | p<0.001 | **-1.46 (RR=0.23, -77%)** | **p=0.001** | No change |
| M. franciscanus | Biomass | — | — | **+0.61 (RR=1.84, +84%)** | **p=0.026** | **Gained significance** |
| M. franciscanus | Density | RR~0.49, -51% | p=0.049 | **-0.39 (RR=0.67, -33%)** | **p=0.365** | **Lost significance** |
| M. pyrifera | Biomass | "~266% higher" | p<0.05 | **+0.58 (RR=1.79, +79%)** | **p=0.048** | Still significant |
| P. interruptus | Biomass | — | non-sig | **+0.25 (RR=1.29, +29%)** | **p=0.483** | Non-significant in both |
| P. interruptus | Density | — | non-sig† | **+1.09 (RR=2.96, +196%)** | **p=0.013** | **Gained significance** |
| S. pulcher | Biomass | — | p=0.56 | **+1.05 (RR=2.85, +185%)** | **p<0.001** | **Gained significance** |
| S. pulcher | Density | — | — | **-0.17 (RR=0.85, -15%)** | **p=0.650** | Non-significant |

*†V5 text (p.8) says "We failed to detect a significant effect of MPA implementation on P. interruptus or S. pulcher biomass or density." But V5 Fig 3 shows a significance star on P. interruptus. I couldn't resolve this V5 contradiction.*

### How the narrative changes

**V5:** "MPAs reduced urchins and increased kelp, but we failed to detect effects on predators."

**Current:** "MPAs increased sheephead biomass (+185%), lobster density (+196%), reduced purple urchin biomass (-49%) and density (-77%), and increased kelp biomass (+79%). Red urchin biomass was unexpectedly *higher* inside MPAs (+84%, p=0.026), which could reflect larger individual size even as density dropped. Red urchin density was non-significant."

Bottom line: the trophic cascade from predators through urchins to kelp now has statistical support at multiple levels, not just the bottom two.

---

## 4. Taxa-by-Taxa Details

### S. pulcher (sheephead) biomass: non-significant → significant

- **V5:** p=0.56, described as "we did not detect an effect"
- **Current:** RR=2.85, +185%, p<0.001 (FDR=0.0007), k=18 from 17 MPAs
- **Why:** V5's Cook's D tossed sheephead observations because their positive effects looked extreme next to urchin negatives. When you keep everything, there's a strong, consistent signal.
- **Robustness:** Significant under all four sensitivity approaches (joint no removal: p<0.001, per-taxa no removal: p<0.001, per-taxa Cook's D: p<0.001, joint Cook's D: p=0.010). Rock solid.

### P. interruptus (lobster) density: non-significant → significant

- **V5:** Text says "failed to detect a significant effect" (but V5 Fig 3 shows a star, which is confusing)
- **Current:** RR=2.96, +196%, p=0.004 (FDR=0.013), k=17 from 11 MPAs
- **Why:** V5's Cook's D kept only 2 of 17 lobster density observations. That's the whole problem.
- **Robustness:** Significant under joint no removal (p=0.004), per-taxa no removal (p=0.016), per-taxa Cook's D (p=0.002). Only non-significant under joint Cook's D (k=2, p=0.124).

### M. franciscanus (red urchin) biomass: non-significant → significant

- **V5:** Not reported separately
- **Current:** RR=1.84, +84%, p=0.014 (FDR=0.026), k=14 from 9 MPAs
- **Why:** The joint model borrows strength from shared τ²(MPA) and τ²(Source). The per-taxa model shows a near-zero, non-significant effect (-0.11, p=0.79) because k=14 with high heterogeneity just doesn't have the power on its own.
- **The tricky part:** Red urchin *biomass* is higher inside MPAs even as density trends negative. This is actually consistent with prior work. Malakhoff & Miller (2021, *Proc. R. Soc. B* 288: 20203061) found that red urchin biomass nearly quadrupled (+397%) inside Channel Islands MPAs after 15 years of protection, driven by release from fishing pressure rather than top-down cascade effects. Teck et al. (2017, *Biol. Conserv.* 210: 321–330) similarly showed red urchins were larger inside MPAs with greater adult biomass density and reproductive biomass density, consistent with reduced harvest allowing individuals to grow larger. So this isn't anomalous — it's a well-documented direct fishing effect that runs counter to the indirect cascade prediction. We should frame it in the Discussion as the expected outcome of protecting a commercially harvested species, distinct from the predator-mediated urchin suppression seen in purple urchins.
- **Robustness:** Significant only under the joint model (p=0.014). Non-significant under per-taxa no removal (p=0.794), per-taxa Cook's D (p=0.479), and joint Cook's D (p=0.459). **This one is model-dependent.** We should flag that clearly.

### M. franciscanus (red urchin) density: significant → non-significant

- **V5:** p=0.049 (barely made it), described as ~51% lower
- **Current:** RR=0.67, -33%, p=0.284 (FDR=0.365), k=16 from 10 MPAs
- **Why:** Lots of heterogeneity across MPAs. The effect points the right direction but isn't consistent enough.

### P. interruptus (lobster) biomass: non-significant (unchanged)

- **V5:** Non-significant
- **Current:** RR=1.29, +29%, p=0.429 (FDR=0.483), k=5 from 5 MPAs
- **Why:** Just not enough data. Only 5 MPAs have lobster biomass. The density metric captures the lobster signal much better (k=17, strong).

### M. pyrifera (kelp) biomass: still significant, smaller effect

- **V5:** "~266% higher" (I think they meant RR≈2.66, i.e. +166%, but the wording is ambiguous)
- **Current:** RR=1.79, +79%, p=0.032 (FDR=0.048), k=17 from 11 MPAs
- **Why:** The joint model shares variance components across taxa, which shifts estimates vs. per-taxa models. The per-taxa model gives a bigger number (+1.04 lnRR, RR=2.84) with a wider CI.
- **Leave-one-out:** Significant in 16 of 16 iterations. Not going anywhere.

---

## 5. Cross-Taxa Meta-Regression (Table 3)

This is entirely new. V5 only had an informal urchin-kelp scatter plot (V5 Fig 4). I ran Knapp-Hartung meta-regressions testing whether predator MPA effects predict prey effects across sites.

| Relationship | k | Slope | p | Significant? |
|-------------|---|-------|---|-------------|
| S. purpuratus density → M. pyrifera biomass | 11 | +0.169 | 0.683 | no |
| M. franciscanus density → M. pyrifera biomass | 10 | -0.035 | 0.891 | no |
| P. interruptus density → S. purpuratus density | 11 | +0.154 | 0.707 | no |
| S. pulcher density → S. purpuratus density | 11 | -0.095 | 0.890 | no |
| **P. interruptus biomass → S. purpuratus biomass** | **4** | **-0.420** | **0.025** | **yes** |
| S. pulcher biomass → S. purpuratus biomass | 9 | +0.290 | 0.081 | no |

Only 1 of 6 pathways hits: lobster biomass predicts purple urchin biomass (p=0.025). But it's based on just k=4 MPAs, so we should interpret carefully. The other pathways are flat. Using estimated effect sizes as predictors introduces errors-in-variables bias that pulls slopes toward zero, sample sizes are small, and between-MPA heterogeneity is high. The temporal stuff (recovery curves, cascade consistency, phase portraits) tells a more convincing cascade story than these cross-sectional regressions.

---

## 6. Other Methodological Changes

### Effect size estimation

| Change | V5 | Current | Impact |
|--------|----|---------|---------|
| **NLS standard errors** | Assumed independence of NLS parameters | Delta method using full variance-covariance matrix | V5 underestimated uncertainty |
| **NLS fitting** | `nls.lm()` → `nls2(brute-force)` only | Multi-step cascade: `nls(port)` → `nls.lm()` → `nls2(brute-force)` → SSasymp/SSlogis selfStart → 4PL reparameterization | Recovers fits that V5 dropped. Fallbacks flagged in `outputs/model_fallback_audit.csv` |
| **Zero correction** | Fixed +0.01 added to all zero proportions | Half the minimum non-zero proportion (adaptive) | Standard practice. The fixed constant was dominating small proportions |

### Multiple testing

**New:** Benjamini-Hochberg FDR correction across all 9 taxon-response coefficients. All 6 significant results survive.

### Temporal recovery models (entirely new)

V5 had nothing temporal. I added `lmer(lnRR ~ time*species + (1+time|MPA) + (1|source))`.

Key slopes from the pooled model:

| Species | Slope (lnRR/yr) | p |
|---------|-----------------|---|
| M. pyrifera (kelp) | +0.16 | <0.001 |
| P. interruptus (lobster) | +0.11 | <0.001 |
| S. purpuratus (purple urchin) | -0.08 | 0.002 |
| M. franciscanus (red urchin) | -0.06 | 0.014 |
| S. pulcher (sheephead) | +0.02 | 0.261 |

### Publication bias and diagnostics (entirely new)

V5 had Cook's D removal but no real sensitivity analysis, no publication bias tests, and no model diagnostics. I added:

- Funnel plots + Egger's test from the joint models (Fig S9)
- Leave-one-out sensitivity for kelp biomass (16/16 significant)
- DHARMa residual diagnostics for the NLS models (Fig S8)
- lmer residual diagnostics (Fig S10)
- Four-method outlier/model sensitivity comparison (Table S9)

### Cascade consistency scoring and phase portraits (entirely new)

8 of 11 MPAs with data for all 5 species score 4 to 5 out of 5 on cascade consistency (predators up, urchins down, kelp up). Phase portraits show network-wide mean lnRR trajectories over time. These are descriptive, not formal tests. The real evidence comes from the meta-analysis and temporal models.

---

## 7. Figures

### Main text figures

| Figure | V5 | Current |
|--------|----|---------|
| **Fig 1** | Map with 6 inset time series panels | Map only (time series moved to Fig 2) |
| **Fig 2** | Data processing pipeline | **New.** Trophic cascade case studies: 3x3 grid (predators/urchins/kelp x 3 Channel Islands sites) with before/after time series and trends |
| **Fig 3** | Meta-analytic mean effects by taxa | Same idea, now with RR-scaled axis, FDR significance stars, and common names |
| **Fig 4** | Urchin-kelp scatter plot | **New.** Recovery trajectories: 3x2 trophic grid, lmer predictions with 95% CI |

- V5's pipeline figure → now **Fig S1**
- V5's urchin-kelp scatter → **gone** (stats now in Table 3)

### Supplemental figures

| SI Figure | Disk File | Description | New? |
|-----------|-----------|-------------|------|
| S1 | `fig_s01_data_processing` | Pipeline illustration | Relocated from V5 Fig 2 |
| S2 | `fig_s02_forest_plot` | Forest plot by MPA | Revised |
| S3 | `fig_s04_recovery_curves` | GAM + lmer recovery overlay | **New** |
| S4 | `fig_s05_cascade_phase` | Phase portraits | **New** |
| S5 | `fig_s06_slope_comparison` | Slope comparison + cascade consistency | **New** |
| S6 | `fig_s09_moderator_comparisons` | SMR vs SMCA, Islands vs mainland | **New** |
| S7a-e | `fig_s08_appendix_*` | Site-level time series (5 files) | Revised |
| S8 | `fig_s11_dharma_diagnostics` | DHARMa residuals | **New** |
| S9 | `fig_s12_funnel_plots` | Funnel plots + Egger's test | **New** |
| S10 | `fig_s13_lmer_residuals` | lmer residual diagnostics | **New** |
| S11 | `fig_s14_model_selection` | NLS model selection + DW | **New** |
| S12 | `fig_s15_sensitivity_summary` | Cook's D + outlier sensitivity | **New** |

---

## 8. Known Limitations

We should hit these in the Discussion:

1. **Bio/Den non-independence:** Biomass and density come from the same individuals. The 9 tests aren't fully independent. FDR helps but doesn't fully fix this.
2. **Cross-taxa attenuation bias:** Using estimated effect sizes as predictors introduces errors-in-variables bias. Table 3 regressions are conservative (biased toward the null).
3. **Climate/stressor attribution:** The core pBACIPS analysis is paired but observational. Scripts 14-23 now add marine-heatwave timing/intensity, environmental-moderator, repeatability, resistance/recovery, and SSWD/sunflower-star checks; ENSO and other overlapping regional drivers still cannot be fully isolated.
4. **Cross-program biomass bootstrap:** KFM and LTER urchin biomass comes from applying PISCO size-frequency distributions to density counts. That assumes similar size structure across programs.
5. **Proportion-based lnRR:** Effect sizes are on proportions (species % of community), not raw densities. Non-standard. We need a clear justification in Methods.
6. **M. franciscanus biomass-density split:** Red urchin *biomass* is up +84% (p=0.026) inside MPAs but density is down -33% (p=0.365, non-sig). This pattern is consistent with prior work: Malakhoff & Miller (2021) found red urchin biomass nearly quadrupled inside Channel Islands reserves due to release from fishing pressure, and Teck et al. (2017) showed red urchins were larger inside MPAs with greater adult and reproductive biomass density. The most parsimonious explanation is a direct fishing effect — reduced harvest allows individuals to grow larger (more biomass per capita) even as density may decline from predation. This complicates the simple cascade narrative but is ecologically coherent. Only significant under the joint model; we should discuss as a direct fishery effect overlaid on the indirect cascade.
7. **Cascade consistency has no null:** 8 of 11 MPAs scored 4 to 5 out of 5. Looks great but there's no formal test against a null expectation.
8. **Moderator power:** SMR vs SMCA and Islands vs mainland comparisons are exploratory. Small subgroups.
9. **Extreme kelp values:** Some MPAs have lnRR≈9 for kelp because reference-site kelp is near zero. Not measurement error. Inverse-variance weighting handles it.
10. **Short before-period:** Only 3 to 5 years before MPA establishment. Hard to formally test parallel trends.
11. **SSWD:** Sea Star Wasting Disease (2013-2014) hit urchin populations differently inside vs outside MPAs.
12. **Fishing displacement:** Effort may pile up near MPA boundaries. Not modeled.

---

## 9. Action Items for the Manuscript

### Must do before submission

| Priority | Task | Details |
|----------|------|---------|
| **1** | **Rewrite Results** | Almost every number changed. Big ones: sheephead biomass (+185%, p<0.001), lobster density (+196%, p=0.013), red urchin biomass (+84%, p=0.026) are now significant. Red urchin density is now non-significant. Pull values from `tables/table_02_meta_analysis.csv`. |
| **2** | **Rewrite Abstract** | V5 says "we did not detect effects on predator density and biomass." Now sheephead biomass (+185%, p<0.001) and lobster density (+196%, p=0.013) are both significant. Kelp is +79% (V5 said "~266% higher"). Red urchin biomass is unexpectedly positive (+84%). |
| **3** | **Rewrite Discussion** | V5 says "we did not detect an increase in density or biomass of key predatory species." Now 2 of 4 predator metrics are significant. Purple urchin density and biomass both decline, but red urchin biomass *increases* and density is non-sig. We need to discuss the size-density tradeoff idea. |
| **4** | **Update Methods** | Describe joint multilevel meta-analysis, FDR correction, Source RE, delta method SEs. Draft text below. |
| **5** | **Update figure/table refs** | V5 figure numbers are completely different. See mapping above. |
| **6** | **Add Limitations paragraph** | Cover bio/den non-independence, cross-taxa attenuation, climate/stressor attribution, and red urchin heterogeneity. |
| **7** | **Fill reference placeholders** | V5 has "REF", "REFS", "CITE" on pages 5, 6, 7, 12. |
| **8** | **Check SI document** | Browse `docs/supporting_information.html`. Make sure figure/table refs match and update SI refs in the manuscript. |
| **9** | **Upload Dryad data** | Upload `dryad_staging/donham_stier_mpa_kelp_data.zip`, then update the DOI placeholder in `code/R/00_download_data.R` and Data Availability. |
| **10** | **Browse results report** | Open `docs/results_report.html`. Every pipeline output on one page. Good for sanity checking. |

### Line-by-line revision guide

**Abstract — V5 claims to replace:**

| V5 Claim | What to Write Instead |
|----------|----------------------|
| "we did not detect effects of MPA implementation on predator density and biomass" | Sheephead biomass +185% (p < 0.001) and lobster density +196% (p = 0.004) — **significant predator responses detected** |
| "sea urchin densities declined" | **Purple** urchin density declined 77% (p < 0.001); red urchin density non-significant (-33%, p = 0.28) |
| "kelp ~266% higher" | Kelp +79% higher (p = 0.032) — same direction, smaller magnitude in joint model |

**Results — key sentences to replace:**

| V5 Statement | Replacement |
|-------------|-------------|
| "We failed to detect a significant effect on P. interruptus or S. pulcher biomass or density" | S. pulcher biomass increased 185% (RR = 2.85, p < 0.001) and P. interruptus density increased 196% (RR = 2.96, p = 0.004) inside MPAs |
| "significant effect on densities of BOTH purple and red sea urchins" | Only purple urchin density declined significantly (-77%, p < 0.001); red urchin density was non-significant (-33%, p = 0.28) |
| "densities were ~78% and 51% lower for purple and red urchins" | Purple urchin density was 77% lower (p < 0.001); red urchin density was 33% lower but non-significant (p = 0.28) |
| (no V5 equivalent) | Red urchin **biomass** increased 84% (RR = 1.84, p = 0.015) — urchins grow larger inside MPAs even as density declines, consistent with reduced harvest pressure |
| "19 MPAs" | 23 MPAs |

**Discussion — key claims to revise:**

| V5 Claim | Revision |
|----------|---------|
| "we did not detect an increase in density or biomass of key predatory species" | 2 of 4 predator metrics are significant (sheephead biomass, lobster density) — discuss predator recovery |
| Urchin decline claims (both species) | Specify only S. purpuratus density; M. franciscanus effects are mixed (biomass up, density non-sig) |
| (new) | Discuss red urchin biomass increase: +84% inside MPAs suggests fishing release on individual size, even as density effects are non-significant |

### Figure and table reference mappings

**Figures (V5 → Current):**

| When V5 says... | Change to... | Notes |
|----------------|-------------|-------|
| Fig 1 | Fig 1 | Same concept (map) but inset time series removed |
| Fig 2 | Fig S1 | Pipeline diagram demoted to supplemental |
| Fig 3 | Fig 3 | Mean effects — same concept, different number path |
| Fig 4 | Table 3 | Urchin-kelp scatter replaced by cross-taxa meta-regressions |
| — | Fig 2 (NEW) | Cascade case studies: 3x3 before/after grid |
| — | Fig 4 (NEW) | Recovery trajectories: lmer predictions with 95% CI |
| Fig S1 | Fig S2 | Forest plot renumbered |

**Tables:**

| When V5 says... | Change to... |
|----------------|-------------|
| Table 1 (meta-analysis) | Table 2 (Table 1 is now average responses) |
| Table 2 | Table 2 (same) |
| Table 3 "linear models" | Table 3 (cross-taxa meta-regression — different analysis) |

### I need your expertise on two things

1. **Kelp biomass conversion:** The `bio_macro()` function uses stipe-calibrated slopes, but LTER data gives frond density. Can we apply stipe-calibrated allometry to frond counts? (See `Donham-Stier-CA-MPA-Data-2026/code/R/06_lter_processing.R`)
2. **KFM macrocystis column:** Can you confirm the column name `SurveyYear` in `KFM_Macrocystis_RawData_1984-2023.csv`?

---

## 10. Suggested Methods Text

### Joint multilevel meta-analysis

> We estimated the overall effect of MPA protection on each taxon-response combination using joint multilevel meta-analytic models (metafor::rma.mv, Viechtbauer 2010). For each response type (biomass and density separately), we fit a single model with taxon as a categorical moderator in cell-means parameterization (mods = ~ Taxa - 1). MPA identity and data source (PISCO, KFM/MBON, LTER, Landsat) were crossed random effects. The joint model borrows strength across taxa for shared variance components, yielding more stable estimates for small-sample taxa. The primary analysis reports the full dataset without outlier removal (best practice, Viechtbauer & Cheung 2010). We conducted sensitivity analyses comparing four approaches: (1) joint model without outlier removal (primary), (2) per-taxa separate models without outlier removal, (3) per-taxa Cook's distance at 4/k, and (4) joint Cook's distance at 4/n (Table S9). P-values were corrected for multiple comparisons across all 9 taxon-response coefficients using the Benjamini-Hochberg false discovery rate method.

### Cross-taxa meta-regression

> To test for trophic cascade effects across MPAs, we used Knapp-Hartung meta-regression (metafor::rma, method = "REML", test = "knha") relating each prey's log response ratio to the corresponding predator's log response ratio across MPAs. Six predator-prey relationships were tested: lobster and sheephead effects on purple and red urchin (biomass and density), and purple and red urchin density effects on kelp biomass.

---

## 11. Repository Structure

### Two repos, one pipeline

I split the analysis into two repos so the raw monitoring data (~1.3 GB) doesn't bog down the analysis code.

- **[`Donham-Stier-CA-MPA-Data-2026`](https://github.com/stier-lab/Donham-Stier-CA-MPA-Data-2026)** — raw data processing (PISCO, KFM, LTER, Landsat). Takes the original monitoring files and produces 4 harmonized CSVs.
- **[`Donham-Stier-CA-MPA-2026`](https://github.com/stier-lab/Donham-Stier-CA-MPA-2026)** (this repo) — loads those CSVs, runs the full meta-analysis, and produces every figure, table, and summary in the manuscript.

The harmonized CSVs (~1 MB) are tracked in git under `data/harmonized/`. You don't need the raw data to run anything in this repo.

### Running the pipeline

```r
# Full pipeline (~2.3 min in the latest complete run)
source(here::here("code", "R", "run_all.R"))

# Figures only (~17 sec, uses cached snapshot)
source(here::here("code", "R", "run_figures_only.R"))
```

`run_all.R` runs everything in order: loads data, calculates effect sizes, runs meta-analysis, fits temporal models, makes all figures, runs the in-pipeline resilience subset, and writes all tables. `run_figures_only.R` skips the computation and just regenerates plots from a cached snapshot, which is useful when you only need to tweak figure styling.

### Directory layout

Here's what lives where. The most important folders for you are `plots/`, `tables/`, and `docs/`.

| Folder | What's in it | You'll use it for |
|--------|-------------|-------------------|
| `plots/` | All figures as PDF + PNG (600 DPI) | Reviewing figures, grabbing files for manuscript submission |
| `tables/` | All manuscript tables as CSV | Checking numbers, copying into manuscript |
| `docs/` | Manuscript draft (V5), this revisions doc, SI source, interactive reports | Reading and reviewing |
| `code/R/` | All analysis scripts (numbered 00-23 plus orchestrators) | Understanding methodology, modifying analyses |
| `data/harmonized/` | The 4 input CSVs from the data repo | Checking raw input data |
| `data/cache/` | Intermediate results (bootstrap, map data, figure snapshot) | Don't need to touch these |
| `outputs/` | Audit trails, filter logs, replicate-level effect sizes | Checking data flow and individual MPA results |

### Pipeline scripts

Scripts run in numerical order. Here's what each one does and what it produces.

| Script | What it does | Key outputs |
|--------|-------------|-------------|
| `00_libraries.R` | Loads all R packages | — |
| `00b_color_palette.R` | Defines all colors, shapes, and `theme_mpa()` | — |
| `00c_analysis_constants.R` | Named constants (survey areas, size thresholds, excluded sites) | — |
| `01_utils.R` | Utility functions (effect size helpers, figure rendering, RR axis scales) | — |
| `02_pBACIPS_function.R` | Core pBACIPS statistical methodology | — |
| `03_load_harmonized_data.R` | Loads the 4 harmonized CSVs into R objects | `All.RR.sub.trans`, `All.Resp.sub`, `Landsat.RR`, `Site` |
| `08_effect_sizes.R` | Fits NLS models per MPA/taxa, extracts effect sizes at t=11 | `SumStats.Final` (142 effect sizes), Table S1b, model diagnostics |
| `09_meta_analysis.R` | Joint multilevel meta-analysis, per-taxa sensitivity, cross-taxa regressions | **Table 2**, **Table 3**, Tables S2-S3, S8-S9, variance components |
| `10_temporal_analysis.R` | lmer recovery models, GAMs, phase portraits, cascade consistency | SI Figs S3-S5, Tables S4-S5, S7, species slopes |
| `11_figures.R` | Core main text figures + most SI figures | **Figs 1-4**, SI Figs S1-S2, S7a-e, S8-S12 |
| `13_additional_analyses.R` | Moderator analyses (SMR vs SMCA, Islands vs mainland) | SI Fig S6, Table S6 |
| `14/15/16/19/21/22/23` | In-pipeline resilience / heatwave robustness subset | Heatwave, multiverse, environmental-moderator, replication, stability, resistance/recovery including Fig 5, and memory outputs |
| `17/18/20` | Full resilience-suite modules needing external raw-data or comparison inputs | Eisaguirre reproduction, effectiveness prediction, and SSWD/sunflower-star checks |
| `12_results_summary.R` | Writes all summary CSVs and `RESULTS_SUMMARY.md` | Data flow audit, replicate effects, model results |
| `run_all.R` | Runs everything above in order | Pipeline log in `logs/` |

---

## 12. Supporting Information: V5 vs. Current

The current pipeline produces a complete, auto-generated SI document (`docs/supporting_information.html`, source in `docs/supporting_information.Rmd`) with two main sections: Supplementary Methods (S1) and Supplementary Results (S2). Every table and figure is data-driven — the numbers, significance statements, and narrative text all pull directly from the pipeline output CSVs, so they stay in sync with the analysis automatically.

### What's new in the SI

**Supplementary Methods (Section S1)** — entirely new:

| Section | What it covers |
|---------|---------------|
| S1.1 Data Sources | Detailed descriptions of KFM, LTER, PISCO, and Landsat programs; data harmonization procedure (biomass conversions, bootstrap resampling, proportion-of-max standardization); site selection and exclusion criteria |
| S1.2 Effect Size Estimation | Full pBACIPS methodology: 4 candidate models (step-change, linear, asymptotic, sigmoidal), AICc selection, control-impact fallback, delta method SE estimation, NLS convergence strategies |
| S1.3 Meta-Analytic Framework | Joint multilevel `rma.mv` model specification (equations), outlier identification (Cook's D, 4 approaches), FDR correction, cross-taxa meta-regressions for cascade tests |
| S1.4 Temporal Recovery | lmer specification with random slopes, GAM linearity checks, phase portrait construction, cascade consistency scoring |
| S1.5 Moderator Analyses | SMR vs SMCA, Islands vs mainland meta-regression specification |
| S1.6 Sensitivity Analyses | Six formal sensitivity analyses: outlier methods, kelp leave-one-out, variance components, source RE, AR1 autocorrelation, GAM linearity |

V5 had none of this. The original manuscript just said "see S1" for program details and described the statistical approach in the main text in about two paragraphs.

**Supplementary Results (Section S2)** — entirely new:

| Section | Figures | Tables | What it adds |
|---------|---------|--------|-------------|
| S2.1 Model Selection | — | S1 (data availability + taxa matrix), S1b (sample sizes) | V5 referenced "Table S1" for model frequencies but never produced it |
| S2.2 Meta-Analytic Results | S2 (forest plot) | — | Data-driven narrative: pulls exact estimates, RR, p-values from Table 2 |
| S2.3 Variance Components | — | S2 (variance components), S3a-b (source sensitivity) | New. V5 had MPA as the only random effect; no variance decomposition |
| S2.4 Cross-Taxa Regressions | — | — | New. V5 had the urchin-kelp scatter as a main figure; current analysis tests 6 trophic pairs formally |
| S2.5 Temporal Recovery | S3 (GAM + lmer curves), S4 (phase portraits), S5 (slope comparison) | S4 (temporal meta-regression), S5 (cascade consistency) | Entirely new. V5 had no temporal analysis |
| S2.6 Moderator Analyses | S6 (SMR vs SMCA, Islands vs mainland) | S6 (moderator meta-regression + diagnostics) | Entirely new. V5 didn't test moderators |
| S2.7 Sensitivity Analyses | — | S7a-b (GAM + AR1 diagnostics), S8 (kelp leave-one-out), S9 (outlier sensitivity: 4 methods) | Entirely new. V5 used Cook's D removal as primary with no sensitivity comparison |
| S2.8 Site-Level Appendix | S7a-e (5 taxa, per-MPA time series) | — | Revised. V5 had no per-site appendix |
| S2.9 Model Diagnostics | S8 (DHARMa), S9 (funnel plots), S10 (lmer residuals), S11 (NLS model selection + DW), S12 (Cook's D + sensitivity summary) | — | Entirely new. V5 had no formal diagnostics |

### Key methodological differences exposed by the SI

1. **Zero correction.** V5 used a fixed constant (a = 0.01). The current pipeline uses adaptive zero-correction: half the minimum non-zero proportion within each time series. This avoids inflating effect sizes for rare species.

2. **SE estimation.** V5 computed CIs "from the pooled standard deviation between time points." The current pipeline uses the delta method for NLS models and `emmeans::pairs()` for linear/step-change models, both of which correctly account for covariance between predictions at t=0 and t=11 from the same model.

3. **Random effects.** V5 had only `(1|MPA)`. The current model adds `(1|Source)` to account for systematic differences between monitoring programs (PISCO, KFM, LTER, Landsat).

4. **Outlier handling.** V5 removed outliers via Cook's D as the primary analysis. The current pipeline keeps the full dataset as primary and reports Cook's D as one of four sensitivity approaches (Table S9).

5. **Model structure.** V5 ran separate `rma.mv` per taxon. The current primary is a joint model with `~ Taxa - 1` so all taxa share common variance components, with per-taxa models as sensitivity.

6. **Temporal analysis.** V5 had none. The current SI includes lmer recovery models, GAM linearity diagnostics, phase portraits, cascade consistency scoring, and AR1 sensitivity checks.

7. **Cross-taxa tests.** V5 had a single urchin-kelp scatter plot as Fig 4. The current analysis tests 6 trophic-pair meta-regressions formally (Table 3) and reports cascade consistency scores per MPA.

### Totals

|  | V5 | Current SI |
|--|----|-----------:|
| SI figures | 0 | 12 (S1-S12, with S7 split into 5 sub-figures) |
| SI tables | 1 (referenced, never produced) | 9 (S1-S9, all auto-generated) |
| Sensitivity analyses | 0 | 6 formal approaches |
| Methods equations | 2 (inline) | 5 (numbered, with full derivations) |
| Diagnostic figures | 0 | 5 (DHARMa, funnel, lmer residuals, model selection, Cook's D) |
| References in SI | 0 | 12 |

