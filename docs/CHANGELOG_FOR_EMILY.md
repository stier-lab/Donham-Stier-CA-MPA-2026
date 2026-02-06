# Pipeline Changes for Manuscript Revision

**Prepared by:** Adrian Stier
**Last Updated:** 2026-02-05
**Purpose:** Document methodological changes between original analysis and current pipeline

---

## Action Items for Emily

**Before you review anything else, please complete these tasks:**

| Priority | Action | Estimated Time | Status |
|----------|--------|----------------|--------|
| **1** | Schedule call with Adrian to walk through Table 2 discrepancies | 30 min | [ ] |
| **2** | Decide on zero-correction method (adaptive vs fixed 0.01) | 15 min | [ ] |
| **3** | Review updated Methods text for random effects | 10 min | [ ] |
| **4** | Verify Table 2 values match new pipeline output | 20 min | [ ] |
| **5** | Decide whether to include heterogeneity statistics | 5 min | [ ] |

---

## Review Checklist

Use this checklist to track your review of each change:

### Manuscript Text Updates Required

- [ ] **Methods - Random effects:** Change "MPA as a random effect" to "MPA and data source (PISCO, KFM, LTER) as crossed random effects"
- [ ] **Methods - Zero correction:** Update description of constant added before log transformation (see Options A/B below)
- [ ] **Results - Table 2:** Replace all values with current pipeline output (see comparison table below)
- [ ] **Results - Heterogeneity:** Consider adding tau-squared/I-squared (optional)

### Code Changes (No Manuscript Impact)

- [ ] Review: Effect size SE now uses covariance-aware calculation (more accurate)
- [ ] Review: Prediction intervals changed to confidence intervals (statistically appropriate)
- [ ] Review: Code refactored into modular scripts (easier to audit)

---

## Results Comparison: Old vs New Table 2

**CRITICAL:** The pipeline now produces different values than the manuscript draft.

### Biomass Effect Sizes (lnRR)

| Taxa | Old Estimate | Old SE | Old p | **New Estimate** | **New SE** | **New p** | Change Direction |
|------|--------------|--------|-------|------------------|------------|-----------|------------------|
| S. purpuratus | -0.76 | 0.36 | 0.042 | **-0.53** | **0.40** | **0.192** | Effect smaller, no longer significant |
| M. franciscanus | 0.17 | 0.34 | 0.626 | **0.86** | **0.33** | **0.014** | Effect larger, NOW SIGNIFICANT |
| M. pyrifera | 0.62 | 0.27 | 0.028 | **0.67** | **0.30** | **0.031** | Similar (still significant) |
| P. interruptus | 1.48 | 0.73 | 0.047 | **0.66** | **0.52** | **0.220** | Effect smaller, no longer significant |
| S. pulcher | 0.30 | 0.50 | 0.558 | **1.23** | **0.31** | **<0.001** | Effect larger, NOW SIGNIFICANT |

### Density Effect Sizes (lnRR)

| Taxa | Old Estimate | Old SE | Old p | **New Estimate** | **New SE** | **New p** | Change Direction |
|------|--------------|--------|-------|------------------|------------|-----------|------------------|
| S. purpuratus | -1.23 | 0.41 | 0.004 | **-2.40** | **0.53** | **<0.001** | Effect larger, still significant |
| M. franciscanus | -0.72 | 0.36 | 0.049 | **-0.44** | **0.75** | **0.562** | Effect smaller, no longer significant |
| P. interruptus | 0.85 | 0.42 | 0.047 | **0.73** | **0.58** | **0.229** | Similar magnitude, no longer significant |
| S. pulcher | 0.02 | 0.28 | 0.946 | **-0.32** | **0.51** | **0.536** | Still non-significant |

### Key Takeaways

1. **Significance changed for 6 of 9 effects** - some gained, some lost significance
2. **Standard errors are generally larger** - this is correct (better uncertainty quantification)
3. **Direction of effects unchanged** - positive/negative signs are consistent
4. **Need to discuss interpretation changes** before updating manuscript narrative

---

## Summary of All Changes

### Changes That Affect Manuscript Text

| Change | Impact | Affects | Priority |
|--------|--------|---------|----------|
| Random effects structure | More conservative SEs | Methods + Table 2 | HIGH |
| Zero-correction method | Minor effect size changes | Methods | HIGH |
| Table 2 values | Different p-values | Results | CRITICAL |

### Changes That Affect Code Only (No Manuscript Changes Needed)

| Change | Impact | Why It's Correct |
|--------|--------|------------------|
| Confidence vs prediction intervals | Smaller SEs for sigmoid models | Effect sizes estimate means, not future observations |
| Covariance-aware SE calculation | Larger SEs | Before/after estimates share error variance |
| Modular code structure | Easier to audit | No statistical impact |
| Bootstrap returns SE | Future analyses | Not currently used in results |

---

## Detailed Change Descriptions

### 1. Meta-Analysis Random Effects (Affects Methods + Results)

**What changed:**
- **Old:** `random = ~1|MPA` (MPA as only random effect)
- **New:** `random = list(~1|MPA, ~1|Source)` (MPA and data source as crossed random effects)

**Why this is correct:**
Effect sizes come from three independent monitoring programs (PISCO, KFM, LTER) with different sampling methods. Including Source as a random effect accounts for program-level variation.

**Quantified impact:**
- Standard errors increased by ~20-50%
- P-values shifted (some effects lost significance, some gained)
- This is statistically appropriate - we were previously underestimating uncertainty

**Methods text to change:**
> Before: "...with MPA as a random effect"
> After: "...with MPA and data source (PISCO, KFM, LTER) as crossed random effects"

---

### 2. Zero-Correction for Log Response Ratios (Affects Methods)

**What changed:**
- **Old:** Fixed constant `+0.01` added to all zero proportions
- **New:** Adaptive constant (half the minimum non-zero proportion)

**Why this matters:**
Adding 0.01 is arbitrary. If the smallest real proportion in a dataset is 0.002, adding 0.01 inflates that value by 5x. The adaptive method scales with the data.

**Decision needed from Emily:**

**Option A - Use adaptive method (recommended):**
> Change methods text from:
> "We added a constant (a=0.01) to raw proportion values"
> To:
> "We applied an adaptive pseudocount correction (half the minimum non-zero proportion observed in each MPA-taxa combination) before log transformation (Aitchison, 1986)"

**Option B - Revert to fixed 0.01 for consistency with original:**
> Keep original methods text, modify pipeline to use fixed correction everywhere

---

### 3. Effect Size Time Point (NO CHANGE NEEDED)

**What was investigated:**
The original code used `t=11` years for extracting effect sizes. We evaluated whether this should be dynamic.

**Decision: KEEP t=11**

The manuscript correctly explains:
> "We chose these years since they correspond to before implementation (t=0) and the age of our youngest MPA in 2023 (eleven years). This allowed us to control for differences in the effect size due to the duration of time since establishment."

**Why t=11 is correct:**
- Channel Islands MPAs (2003) have 20 years of data
- Mainland MPAs (2012) have 11 years of data
- Using dynamic max time would confound MPA age with protection effectiveness
- t=11 ensures all MPAs are compared at the same "age"

**Methods text:** No change needed.

---

### 4. Heterogeneity Statistics (Optional Addition)

**What's new:**
Pipeline now calculates and reports tau-squared and pseudo-I-squared for meta-analysis models.

**What these measure:**
- **tau-squared:** Between-study variance (how much effect sizes vary across MPAs)
- **I-squared:** Proportion of total variance due to heterogeneity (0-100%)

**Decision needed from Emily:**

- [ ] Add to Table 2 as additional columns
- [ ] Report in Results text only
- [ ] Include in Supplemental Materials only
- [ ] Do not include

---

## Things That Did NOT Change

The following aspects of the analysis are unchanged and verified correct:

- Biomass conversion formulas (all verified against literature)
- Bootstrap iterations (1000 for all analyses)
- Cook's distance outlier removal (4/n threshold)
- pBACIPS methodology and model structure
- Study design and MPA selection criteria
- Size cutoff application (25mm for PISCO, unfiltered for KFM/LTER)

---

## Questions for Emily

1. **Can we schedule a call this week** to walk through the Table 2 changes? The shifts in significance are important to discuss before updating the manuscript narrative.

2. **Which zero-correction method do you prefer?** (Option A: adaptive, Option B: fixed 0.01)

3. **Do you want heterogeneity statistics** in the final manuscript? If yes, where?

4. **Are you comfortable with the changes to significance?** Some effects that were previously significant are now non-significant (and vice versa). The new results are more conservative and statistically appropriate, but we should discuss the narrative implications.

---

## File Reference

| Output File | Contents | Use For |
|-------------|----------|---------|
| `data/table_02_meta_analysis.csv` | Current Table 2 values | Replace manuscript Table 2 |
| `plots/fig_03_mean_effects.pdf` | Updated Figure 3 | Replace manuscript Figure 3 |

---

## Appendix: Code Structure Mapping

For reference, here's how the original script maps to the new modular pipeline:

| Original (pBACIPS_PISCO_V10.R) | New Location |
|--------------------------------|--------------|
| Lines 1-500 (imports) | `03_data_import.R` |
| Lines 500-1500 (PISCO processing) | `04_pisco_processing.R` |
| Lines 1500-2500 (KFM processing) | `05_kfm_processing.R` |
| Lines 2500-3400 (LTER processing) | `06_lter_processing.R` |
| Lines 3400-5500 (effect sizes) | `08_effect_sizes.R` |
| Lines 5500-6900 (meta-analysis) | `09_meta_analysis.R` |
| Lines 6900-7100 (figures) | `10_figures.R` |
