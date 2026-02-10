# Scope C: Data Processing Module Audit Report

**Auditor:** Claude Opus 4.6 (Scope Agent C)
**Date:** 2026-02-09
**Scope:** 6 files, ~3,112 lines
**Files audited:**
- `code/R/03_data_import.R` (301 lines)
- `code/R/04_pisco_processing.R` (900 lines)
- `code/R/05_kfm_processing.R` (763 lines)
- `code/R/06_lter_processing.R` (710 lines)
- `code/R/06b_landsat_processing.R` (192 lines)
- `code/R/07_combine_data.R` (247 lines)

**Supporting files also reviewed:**
- `code/R/01_utils.R` (875 lines) -- utility functions and constants
- `code/R/00c_analysis_constants.R` (148 lines) -- analysis constants

---

## CRITICAL Issues (Breaks pipeline or produces fundamentally wrong results)

### C-01: `SizeFreq.Urch.OG` preserves FILTERED data, not unfiltered [CRITICAL]

**File:** `code/R/04_pisco_processing.R`, lines 446-450
**Impact:** KFM and LTER urchin biomass estimation uses 25mm-filtered size frequency data instead of the full size distribution

The comment on line 442-445 states: "We preserve the original unfiltered size frequency data (SizeFreq.Urch.OG) for use in those scripts." However, the assignment occurs AFTER the 25mm filter:

```r
# Line 446: FILTER applied first
SizeFreq.Urch <- subset(SizeFreq.Urch, size >= 25)

# Line 450: "OG" assigned AFTER filter
SizeFreq.Urch.OG <- SizeFreq.Urch
```

Both `05_kfm_processing.R` (line 192) and `06_lter_processing.R` (line 142) set `SizeFreq.Urch <- SizeFreq.Urch.OG`, believing they are getting unfiltered data. In reality, they receive the 25mm-filtered PISCO subset.

**Consequence:** KFM and LTER urchin biomass bootstrap resamples only from urchins >= 25mm. If KFM/LTER protocols include smaller urchins in their counts, biomass per individual is overestimated (sampling from larger-than-actual size distribution), inflating total biomass estimates.

**Note:** Whether this materially affects results depends on whether KFM/LTER actually record urchins < 25mm. The comment explicitly says "KFM and LTER use different protocols and may include smaller urchins," so the intent was clearly to NOT apply the 25mm filter for these programs. The fix is to move line 450 before line 446.

---

### C-02: `sites.short.edit` not deduplicated before merge -- root cause of 260 duplicates [CRITICAL]

**Files:** `code/R/03_data_import.R` (lines 212-217), `code/R/04_pisco_processing.R` (lines 851-852), and all files that merge with it
**Impact:** Inflated row counts in response ratio datasets fed into meta-analysis

`sites.short.edit` is extracted from the `Site` table which has one row per monitoring site, not one row per MPA. Multiple monitoring sites can belong to the same MPA:

```r
# 03_data_import.R, lines 212-217:
sites.short.edit <- subset(Site, select = c(
  CA_MPA_Name_Short, type, Location, Hectares
))
# NOT deduplicated -- contains duplicate CA_MPA_Name_Short entries
```

This table is then merged with response ratio data using `all = TRUE` (full outer join) on `CA_MPA_Name_Short`:

```r
# 04_pisco_processing.R, lines 854-857:
Swath.join.sub <- merge(Swath.join.sub, sites.short.edit,
                        by.x = c("CA_MPA_Name_Short"),
                        by.y = c("CA_MPA_Name_Short"),
                        all = TRUE)
```

The same pattern occurs in:
- `05_kfm_processing.R` line 632 (KFM.join.ave)
- `05_kfm_processing.R` line 742 (kfm.fish)
- `06_lter_processing.R` line 388 (LTER.join.ave)
- `06_lter_processing.R` line 469 (LTER.lob)
- `06_lter_processing.R` line 571 (Short.lter.macro)
- `06_lter_processing.R` line 671 (LTER.fish)
- `06b_landsat_processing.R` line 157 (Landsat.RR)

If an MPA has N monitoring sites, each response ratio row for that MPA gets replicated N times. The downstream `complete.cases(year)` filter removes rows from `sites.short.edit` that lack matching data, but does NOT deduplicate the inflated rows.

**This is the most likely root cause of the 260 duplicates previously reported in `All.RR.sub`.** The duplicate detection in `07_combine_data.R` (lines 104-123) catches and removes these, but this is a band-aid fix -- the duplication should be prevented at the source.

**Fix:** Add `sites.short.edit <- sites.short.edit[!duplicated(sites.short.edit$CA_MPA_Name_Short), ]` after creating the table, both in `03_data_import.R` and `04_pisco_processing.R`.

---

### C-03: LTER lobster code references `$SIZE` but CSV column is `SIZE_MM` [CRITICAL]

**File:** `code/R/06_lter_processing.R`, line 412
**Impact:** Relies on R's deprecated partial column name matching; may break in future R versions

```r
lter.lob.site$biomass <- lter.lob.site$COUNT * bio_lobster(lter.lob.site$SIZE)
```

The LTER lobster CSV (`Lobster_Abundance_All_Years_20230831.csv`) has column header `SIZE_MM`, not `SIZE`. R's `$` operator performs partial matching, so `$SIZE` matches `SIZE_MM` with a warning. This behavior is deprecated as of R 4.x and may be removed in future versions.

Similarly, line 621 references `lter.fish.sub.site$SIZE`, but this is correct because the fish CSV (`Annual_fish_comb_20240307.csv`) actually has a column named `SIZE`.

**Consequence:** Currently works via partial matching (SIZE_MM values are indeed in mm), so the biomass calculation is numerically correct. But this is fragile and will break when R removes partial matching support.

---

## MODERATE Issues (Could produce wrong results under certain conditions)

### C-04: Landsat data excluded from combined response ratio dataset [MODERATE]

**File:** `code/R/07_combine_data.R`, lines 74-75
**Impact:** Landsat kelp data missing from `All.RR.sub.trans` used in meta-analysis

The `All.RR` rbind in `07_combine_data.R` combines data from LTER, PISCO, KFM, and fish sources:

```r
All.RR <- rbind(LTER.join.ave, Swath.join.sub, KFM.join.ave, LTER.lob,
                Short.lter.macro, LTER.fish, kfm.fish)
```

`Landsat.RR` is conspicuously absent. The Landsat response ratios are only used later in `08_effect_sizes.R` (line 1320-1326) where they are processed separately. However, this means `All.RR.sub.trans` -- the main combined dataset -- does not include Landsat-derived kelp canopy data.

**Consequence:** This may be intentional (Landsat uses different units/methods), but it means the standardized species names, source labels, and sheephead exclusion logic are not applied to Landsat data through the common pathway. The effect size extraction in `08_effect_sizes.R` handles Landsat separately, which is potentially correct but creates an inconsistent pathway that could lead to missed exclusions or naming mismatches.

---

### C-05: `assign_time_from_site_table` produces false-positive warnings [MODERATE]

**File:** `code/R/01_utils.R`, lines 103-109
**Impact:** Spurious warnings clutter output and may mask real problems

```r
# Warn about MPAs that were not found in the site table (defaulting to time=0)
unmatched <- unique(df$CA_MPA_Name_Short[df$time == 0])
```

The warning logic checks for MPAs where `time == 0`, but `time == 0` is the correct value for ALL pre-implementation observations (years <= MPA start). The warning thus fires for every MPA that has any "Before" period data, falsely claiming they were "not found in site table."

**Consequence:** The warning message is misleading. Legitimate unmatched MPAs (truly not in the site table) are buried among the false positives, making it impossible to identify real matching failures.

---

### C-06: PISCO BA assignment differs from other modules [MODERATE]

**File:** `code/R/04_pisco_processing.R`, line 868
**Impact:** Inconsistent Before/After logic between data sources

PISCO uses:
```r
Swath.join.sub$BA <- ifelse(Swath.join.sub$time == 0, "Before", "After")
```

KFM and LTER use:
```r
KFM.join.ave <- assign_ba_from_site_table(KFM.join.ave, Sites2)
```

While the outcome is equivalent (both assign "Before" when year <= MPA start), the PISCO approach goes through an intermediate `time` variable that was already set by `assign_time_from_site_table`. This means PISCO's BA assignment depends on the time assignment being correct. If the time assignment had any bug (e.g., MPA not found, defaulting to time=0), all those rows would be labeled "Before" instead of getting `BA = NA` as the utility function would assign.

The KFM/LTER approach is more robust because `assign_ba_from_site_table` independently looks up MPA start years and assigns NA for unmatched MPAs.

---

### C-07: PISCO `sites.short.edit` overwritten mid-pipeline [MODERATE]

**File:** `code/R/04_pisco_processing.R`, lines 847-852 vs `code/R/03_data_import.R`, lines 212-217
**Impact:** Different definitions of `sites.short.edit` in different scripts

`03_data_import.R` creates `sites.short.edit` from `Site`. Then `04_pisco_processing.R` re-imports Site_List_All.csv (line 727), re-reads MPAfeatures_subset.csv (line 848), re-merges them (line 849), and re-creates `sites.short.edit` (lines 851-852), overwriting the version from `03_data_import.R`.

This means the `sites.short.edit` used by PISCO (and by all downstream scripts 05, 06, 06b that reference this global variable) comes from the re-created version in `04_pisco_processing.R`, not the original from `03_data_import.R`. The two should be identical, but the redundant creation is fragile -- if the data files differ between reads (e.g., edited between import and processing), the results would silently diverge.

---

### C-08: Macrocystis stipe bootstrap uses `rowSums(.)` on data with metadata columns [MODERATE]

**File:** `code/R/05_kfm_processing.R`, lines 415-421
**Impact:** Potential inclusion of non-stipe columns in rowSums

```r
s_bio <- s %>%
  dplyr::mutate(dplyr::across(starts_with("X"), ~ bio_macro(.)))
s_bio <- s_bio %>%
  dplyr::mutate(sum = rowSums(.))
s <- s %>%
  dplyr::mutate(sum = rowSums(.))
```

The `rowSums(.)` on the full dataframe includes ALL columns, not just the bootstrap sample columns. The `s` matrix is created as `data.frame(matrix(NA, nrow = 1000, ncol = n))`, so column names are X1, X2, ..., Xn. The `mutate(across(starts_with("X"), ...))` correctly targets only these columns. But `rowSums(.)` sums ALL columns, including the newly created bio-converted ones AND the original stipe count columns (if the `across` added new columns rather than modifying in place).

Actually, on closer inspection: `across()` with `starts_with("X")` modifies columns in place (they still start with "X"). So after `mutate(across(starts_with("X"), ~ bio_macro(.)))`, the X columns contain biomass values. Then `rowSums(.)` correctly sums all columns (which are all biomass values). For `s` (stipe counts), `rowSums(.)` sums all stipe count columns. This appears correct since the data frame only contains the sampled values. But the `mutate(sum = rowSums(.))` call includes the `sum` column itself in the calculation if it were already present -- which it is not on the first call, so this is actually fine.

**Revised assessment:** After more careful analysis, this is likely correct because the data frame `s` only contains the bootstrap sample columns. However, the code pattern is fragile and would break if any metadata columns were added.

---

### C-09: KFM density filter uses string comparison instead of numeric [MODERATE]

**File:** `code/R/05_kfm_processing.R`, line 560
**Impact:** Potential failure to filter non-numeric "N/A" values

```r
All.den <- subset(All.den, den != "N/A")
```

The variable `den` is created via `safe_divide()` which returns numeric values or NA. The string comparison `den != "N/A"` compares a numeric column to a character string. In R, this coerces the numeric column to character, which may produce unexpected behavior. If `den` contains actual `NA` values (from safe_divide replacing zeros), the comparison `NA != "N/A"` returns `NA`, which `subset()` drops (effectively filtering out NAs).

The same pattern appears at line 580: `All.bio <- subset(All.bio, biomass != "N/A")`.

**Consequence:** Likely works by accident -- the string "N/A" comparison on a numeric column effectively becomes a no-op (always TRUE for non-NA values), while the `NA != "N/A"` evaluating to NA drops true NA rows. But this is confusing and fragile.

---

### C-10: LTER lobster SIZE_MM == -99999 sentinel values not explicitly filtered [MODERATE]

**File:** `code/R/06_lter_processing.R`, lines 412-413
**Impact:** Silent handling of sentinel values via NaN arithmetic

The LTER lobster data uses -99999 as a missing value sentinel for SIZE_MM. The code never explicitly filters these out:

```r
lter.lob.site$biomass <- lter.lob.site$COUNT * bio_lobster(lter.lob.site$SIZE)
lter.lob.site$biomass[is.na(lter.lob.site$biomass)] <- 0
```

`bio_lobster(-99999)` = `0.001352821 * (-99999)^2.913963113` = `NaN` (negative base to non-integer power). Then `COUNT * NaN = NaN`, which `is.na()` catches and converts to 0.

This works because: (1) -99999 rows have COUNT = 0, so they would contribute 0 biomass regardless; (2) even if COUNT > 0, the NaN gets caught. But it relies on implicit NaN arithmetic rather than explicit filtering.

Compare with LTER fish (line 616): `lter.fish.sub.site <- subset(lter.fish.sub.site, COUNT != -99999)`, which explicitly filters the sentinel. LTER kelp (line 508) also filters: `subset(lter.macro.site, FRONDS != -99999)`. The lobster processing is the odd one out.

---

## MINOR Issues (Cosmetic, code quality, or maintainability)

### C-11: Hardcoded survey areas not using constants from `00c_analysis_constants.R` [MINOR]

**Files:** Multiple
**Impact:** Code maintainability

Despite `00c_analysis_constants.R` defining named constants like `PISCO_SWATH_AREA_M2 = 60`, `KFM_QUAD_AREA_M2 = 2`, `LTER_LOBSTER_PLOT_AREA_M2 = 300`, and `LTER_FISH_SURVEY_AREA_M2 = 80`, the processing scripts use hardcoded literals:

- `04_pisco_processing.R` line 185: `count / 60` (should use `PISCO_SWATH_AREA_M2`)
- `04_pisco_processing.R` line 258: `count / 60`
- `04_pisco_processing.R` line 396: `count / 60`
- `04_pisco_processing.R` line 559: `count / 60`
- `04_pisco_processing.R` line 695: `count / 60`
- `04_pisco_processing.R` line 673: `count / 60`
- `05_kfm_processing.R` line 311: `count.ave / 2` (should use `KFM_QUAD_AREA_M2`)
- `06_lter_processing.R` line 479: `COUNT / 300` (should use `LTER_LOBSTER_PLOT_AREA_M2`)
- `06_lter_processing.R` line 480: `biomass / 300`
- `06_lter_processing.R` line 680: `biomass / 80` (should use `LTER_FISH_SURVEY_AREA_M2`)
- `06_lter_processing.R` line 681: `COUNT / 80`

---

### C-12: Fragile column index selection in multiple locations [MINOR]

**File:** `code/R/04_pisco_processing.R`, lines 400, 408, 527, 573, 632
**Impact:** Code fragility

Several places still use numeric column indices for reordering:

```r
VRG.PANINT <- VRG.PANINT[, colnames(VRG.PANINT)[c(1:6, 10, 7:9)]]     # line 400
PANINT.all <- PANINT.all[, colnames(PANINT.all)[c(2, 1, 7, 8, 11, 3:6, 10, 9)]]  # line 408
PISCO.Urchin.site.merge <- PISCO.Urchin.site.merge[, colnames(PISCO.Urchin.site.merge)[c(2, 1, 3, 17, 5:10, 16, 20)]]  # line 527
Swath.PISCO <- Swath.PISCO[, colnames(Swath.PISCO)[c(1:11, 15)]]  # line 573
Fish.sub.site <- Fish.sub.site[, colnames(Fish.sub.site)[c(1:13, 20:22, 25)]]  # line 632
```

If upstream merges change column order (e.g., from different R or dplyr versions), these indices would silently select wrong columns.

---

### C-13: Redundant re-import of `Site_List_All.csv` in `04_pisco_processing.R` [MINOR]

**File:** `code/R/04_pisco_processing.R`, line 727
**Impact:** Unnecessarily re-reads data already loaded in `03_data_import.R`

```r
Site <- read.csv(here::here("data", "Site_List_All.csv"))
```

This overwrites the `Site` object already created and validated in `03_data_import.R` (line 172). The re-import bypasses the schema validation done in `03_data_import.R` and does not use `safe_read_csv()`.

---

### C-14: KFM fish density response uses hardcoded column indices in gather [MINOR]

**File:** `code/R/05_kfm_processing.R`, lines 607, 619, 703
**Impact:** Code fragility

```r
KFM.den.long <- gather(KFM.den, status, value, 5:6)          # line 607
KFM.bio.long <- gather(KFM.bio, status, value, 5:6)          # line 619
KFM.fish.den.long <- gather(KFM.fish.den.sub, status, value, 5:6)  # line 703
```

Using numeric column positions (5:6) instead of column names is fragile and will break if upstream column order changes. The PISCO code (line 815) uses column names: `gather(PISCO.den, status, value, mpa, reference)`.

---

### C-15: LTER bootstrap loop uses O(n^2) `rbind` pattern [MINOR]

**File:** `code/R/06_lter_processing.R`, lines 195-229
**Impact:** Performance only

Unlike the PISCO and KFM bootstrap loops which pre-allocate lists and use `do.call(rbind, results_list)`, the LTER bootstrap loop grows the dataframe row-by-row:

```r
Urchin.site <- rbind(Urchin.site, Urchin.SF)  # line 228
```

This is O(n^2) complexity. For large datasets this could be very slow. The PISCO and KFM scripts were apparently optimized but the LTER script was not.

---

### C-16: LTER urchin bootstrap not wrapped in tryCatch [MINOR]

**File:** `code/R/06_lter_processing.R`, lines 210-214
**Impact:** A single failed bootstrap iteration could crash the entire pipeline

```r
result <- bootstrap_biomass(count = n,
                            size_freq_indices = t2,
                            size_freq_table = SizeFreq.Urch,
                            biomass_fun = bio_fun,
                            n_resamples = 1000)
```

Both PISCO (line 488) and KFM (line 239) wrap their bootstrap calls in `tryCatch`. The LTER call does not, meaning any unexpected error in the bootstrap function will halt the entire pipeline instead of gracefully producing NA for that iteration.

---

### C-17: Missing `set.seed()` in LTER urchin bootstrap [MINOR]

**File:** `code/R/06_lter_processing.R`
**Impact:** Non-reproducible bootstrap results

PISCO uses `set.seed(12346)` (line 471), KFM uses `set.seed(12347)` (line 223), and `00c_analysis_constants.R` defines `SEED_LTER_URCHIN <- 12349`. However, the LTER urchin bootstrap loop (lines 195-229) does not call `set.seed()` before the loop. Results will vary between runs.

---

### C-18: `source(here::here("code", "R", "01_utils.R"))` redundantly called [MINOR]

**Files:** `code/R/06_lter_processing.R` (line 54), `code/R/06b_landsat_processing.R` (line 51)
**Impact:** Minor performance and clarity issue

`01_utils.R` is sourced in the pipeline orchestrator (`run_all.R`) before these scripts execute. Re-sourcing it is harmless but unnecessary and adds confusion about whether these scripts can run standalone.

---

## Summary Table

| ID | Severity | Category | File(s) | Summary |
|----|----------|----------|---------|---------|
| C-01 | CRITICAL | Data integrity | 04_pisco | SizeFreq.Urch.OG saves filtered data, not original |
| C-02 | CRITICAL | Duplicate rows | 03, 04, 05, 06, 06b | sites.short.edit not deduplicated; root cause of 260 duplicates |
| C-03 | CRITICAL | Data integrity | 06_lter | $SIZE references non-existent column (relies on partial matching) |
| C-04 | MODERATE | Data integrity | 07_combine | Landsat data excluded from combined RR dataset |
| C-05 | MODERATE | Code quality | 01_utils | assign_time_from_site_table false-positive warnings |
| C-06 | MODERATE | Consistency | 04_pisco | PISCO BA assignment differs from KFM/LTER approach |
| C-07 | MODERATE | Code quality | 04_pisco | sites.short.edit overwritten mid-pipeline |
| C-08 | MODERATE | Data integrity | 05_kfm | Macro bootstrap rowSums on potentially mixed columns |
| C-09 | MODERATE | Data integrity | 05_kfm | String vs numeric comparison for N/A filtering |
| C-10 | MODERATE | Data integrity | 06_lter | SIZE_MM sentinel -99999 not explicitly filtered |
| C-11 | MINOR | Maintainability | Multiple | Hardcoded survey areas instead of constants |
| C-12 | MINOR | Code fragility | 04_pisco | Numeric column indices for selection |
| C-13 | MINOR | Code quality | 04_pisco | Redundant re-import of Site_List_All.csv |
| C-14 | MINOR | Code fragility | 05_kfm | Hardcoded column indices in gather() |
| C-15 | MINOR | Performance | 06_lter | O(n^2) rbind in bootstrap loop |
| C-16 | MINOR | Robustness | 06_lter | Missing tryCatch around bootstrap |
| C-17 | MINOR | Reproducibility | 06_lter | Missing set.seed() in urchin bootstrap |
| C-18 | MINOR | Code quality | 06_lter, 06b | Redundant source() of utils |

**Totals:** 3 CRITICAL, 7 MODERATE, 8 MINOR

---

## Specific Answers to Audit Questions

### 1. Data integrity: Are joins correct?

**Partially.** The most problematic joins are the `merge(..., all = TRUE)` calls between response ratio data and `sites.short.edit` (see C-02). These create many-to-many joins because `sites.short.edit` has duplicate `CA_MPA_Name_Short` values. The duplicate detection in `07_combine_data.R` catches this post-hoc but does not prevent it.

All other joins appear structurally correct: site table joins use appropriate keys, and the MBON joins on `site_id` are one-to-one (each site_id appears once in Sites2).

### 2. MPA name matching

The `assign_time_from_site_table` and `assign_ba_from_site_table` functions match MPAs by exact string equality on `CA_MPA_Name_Short`. The false-positive warning in `assign_time_from_site_table` (C-05) makes it difficult to identify genuine mismatches. The different site tables used across modules (PISCO's `sites.short` vs MBON's `Sites2` vs the master `Site` table) have different `CA_MPA_Name_Short` values, but within each module the matching appears correct.

### 3. Duplicate rows: Root cause of 260 duplicates

**Root cause identified: C-02.** The `sites.short.edit` table is not deduplicated before being merged with response ratio data. Since multiple monitoring sites share the same MPA, the merge inflates rows. The `07_combine_data.R` deduplication logic (lines 104-123) catches and removes these, but the root cause persists upstream.

### 4. Species name standardization

Generally consistent. The `standardize_species_names()` function in `01_utils.R` handles the main PISCO codes (STRPURAD, MESFRAAD, MACPYRAD, PANINT, SPUL). KFM data uses full scientific names internally and converts to PISCO codes for bootstrapping, then back. LTER follows the same pattern. One inconsistency: KFM fish uses `taxon_name = "SPUL"` (line 692 of 05_kfm_processing.R), which is converted to "Semicossyphus pulcher" in `07_combine_data.R` (line 180).

### 5. Biomass conversions

Allometric equations appear correct and sourced from published literature:
- **Lobster:** `W = 0.001352821 * CL^2.913963113` (LTER source) -- correct for mm input
- **Red urchin:** `W = 0.00059 * TD^2.917` -- consistent with published values
- **Purple urchin:** `W = 0.00059 * TD^2.870` -- consistent with published values
- **Kelp:** `stipe * average_slope * 1000` where slope is seasonal average from LTER
- **Sheephead (LTER):** `W = 0.0144 * TL^3.04` -- standard length-weight relationship
- **Sheephead (PISCO):** Uses species attribute table from PISCO database

One concern: PISCO lobster size is in cm (converted to mm: `size * 10`), while LTER lobster SIZE_MM is already in mm. The LTER code accesses `$SIZE` (partial match for `SIZE_MM`) and passes directly to `bio_lobster()` without conversion, which is correct since SIZE_MM is already in mm.

### 6. Zero-correction

The `calculate_proportions()` function uses an adaptive pseudocount method (half of minimum non-zero proportion) rather than a fixed +0.01 correction. This is applied consistently across all programs (PISCO, KFM, LTER, Landsat). The adaptive method is statistically sound (Aitchison 1986). No issues identified with the zero-correction approach itself.

### 7. Before/After assignment

MPA implementation year is assigned as "Before" across all modules (year <= MPA_Start yields time = 0 and BA = "Before"). This is consistent. The only inconsistency is the mechanism: PISCO uses `ifelse(time == 0, "Before", "After")` while KFM/LTER use `assign_ba_from_site_table()` (see C-06). The substantive result is the same.
