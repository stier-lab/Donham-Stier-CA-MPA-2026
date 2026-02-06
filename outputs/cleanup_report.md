# Repository Cleanup Report

**Generated:** 2026-02-06
**Repository:** Donham-Stier-CA-MPA-2026

---

## Summary

| Category | Size | Action |
|----------|------|--------|
| plots/archive/ | ~25 MB | Delete (superseded figures) |
| code/archive/ | ~888 KB | Delete (old scripts) |
| logs/ | ~156 KB | Delete (pipeline logs) |
| Root temp files | ~34 KB | Delete |
| **Total recoverable** | **~26 MB** | |

---

## 1. Files Safe to Delete Immediately

### Archive Directories

**plots/archive/** (~25 MB) - Obsolete figure versions:
- `fig_01_kelp_timeseries.pdf/png` - Old format
- `fig_01_map_only.pdf/png` - Superseded by `fig_01_mpa_map.pdf`
- `fig_01_mpa_map_composite.pdf/png` - Superseded by `fig_01_mpa_map.pdf`
- `Fig1_plots.pdf`, `Fig2_*.png`, etc. - Old naming convention
- `Rplots.pdf` - Temporary plot file

**code/archive/** (~888 KB) - Historical scripts:
- `0_libraries.R` - Old naming convention
- `MEE_Supplement_example.R` - External example
- `pBACIPS_PISCO_V10.R` - Old monolithic version (now modular)
- `fig01_*.R` files - Old figure generation
- `kelp_mpa_map_insets.R` - References deleted figures
- `utils_spatial.R` - Unused utilities
- `10b_fig01_map.R` - Alternative figure code

### Temporary Files

- `Rplots.pdf` (root) - Auto-generated R graphics device file
- `firebase-debug.log` - Debug logging artifact

### Log Files

**logs/** (~156 KB) - 32 pipeline execution logs:
- All `pipeline_YYYYMMDD_HHMMSS.log` files
- Already gitignored, safe to delete

---

## 2. Files Requiring Review

### Outputs Directory

These analysis artifacts were created during development:

| File | Size | Purpose |
|------|------|---------|
| `crossed_structure_verification.md` | 7.3 KB | Random effects verification |
| `data_flow_summary.csv` | 434 B | Data flow documentation |
| `filter_audit_effect_sizes.csv` | 49 KB | Effect size filtering audit |
| `filter_audit_meta_analysis.csv` | 26 KB | Meta-analysis filtering audit |
| `filter_summary_by_taxa.csv` | 473 B | Taxa filter summary |
| `random_effects_analysis.md` | 19 KB | Random effects justification |

**Decision:** Keep these as supplementary documentation? They're referenced in CHANGELOG_FOR_EMILY.md.

### Diagnostic Plots

**plots/diagnostics/** (~168 KB):
- `dharma_*.pdf` - DHARMa residual diagnostics
- `nls_diag_*.pdf` - NLS fit diagnostics
- `diagnostic_summary.pdf/png` - Summary diagnostics

**Decision:** Keep for methods documentation or delete after publication?

---

## 3. .gitignore Status

### Already Ignored (Good)
- `.DS_Store`, `._*`
- `.Rhistory`, `.RData`, `.Ruserdata`
- `Rplots.pdf`
- `logs/`
- `data/cache/`, `data_cache/`
- `marmap_coord_*.csv`
- `firebase-debug.log`

### Recommended Additions
```gitignore
# R Project files
.Rproj.user/
*.Rproj

# Editor backup files
*.bak
*~
```

---

## 4. Git Status

### Currently Staged (Ready to Commit)
- 52 files with modifications
- Figure consolidation complete
- Documentation updates
- Code improvements

### Untracked Files
- 6 files in outputs/ (analysis artifacts)

### Archive Files in Git
- 3 files in code/archive/ (tracked)
- 7 files in plots/archive/ (tracked)
- **Action:** Use `git rm -r` to remove from tracking

---

## 5. Recommended Cleanup Commands

```bash
# 1. Delete untracked temp files
rm -f Rplots.pdf firebase-debug.log

# 2. Clear logs directory
rm -rf logs/*

# 3. Remove archive directories from git tracking
git rm -r plots/archive/
git rm -r code/archive/

# 4. Delete untracked marmap cache files
rm -f marmap_coord_*.csv

# 5. Commit the cleanup
git add .
git commit -m "Clean up archive directories and temp files"
```

---

## 6. Directory Structure After Cleanup

```
Donham-Stier-CA-MPA-2026/
├── code/
│   └── R/                    # Active pipeline scripts only
├── data/
│   ├── cache/                # Bootstrap results (gitignored)
│   └── *.csv                 # Analysis outputs
├── docs/
│   ├── CHANGELOG_FOR_EMILY.md
│   ├── RESULTS_SUMMARY.md
│   ├── STATISTICAL_REVIEW.md
│   └── ...
├── outputs/                  # Analysis artifacts
│   ├── filter_audit_*.csv
│   ├── random_effects_analysis.md
│   └── ...
├── plots/
│   ├── fig_*.pdf/png        # Publication figures only
│   └── diagnostics/         # Optional: keep for validation
├── CLAUDE.md
└── README.md
```

---

## 7. Post-Cleanup Tasks

1. [ ] Update CLAUDE.md to remove archive file references
2. [ ] Verify all outputs/ files are referenced in documentation
3. [ ] Consider organizing outputs/ into subdirectories:
   - `outputs/audit/` - Filtering audit files
   - `outputs/supplementary/` - Supplementary analysis docs
4. [ ] Create git tag before cleanup: `git tag v1-pre-cleanup`
