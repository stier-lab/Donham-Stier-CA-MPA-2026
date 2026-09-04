# Kumagai et al. 2024 Comparator Inputs

This directory vendors the small public comparator inputs needed for the
resilience methods and environmental-moderator checks.

Source:

- Kumagai, J. A. et al. 2024. "Marine Protected Areas That Preserve Trophic
  Cascades Promote Resilience of Kelp Forests to Marine Heatwaves." *Global
  Change Biology* 30:e17620. https://doi.org/10.1111/gcb.17620
- Code/data archive: Zenodo record 14188853,
  https://doi.org/10.5281/zenodo.14188853, CC-BY-4.0.
- GitHub source used for the local mirror:
  https://github.com/jkumagai96/Kelp_Forests_and_MPAs
  at commit `dacfd7e0b6732cc0ea270a9fb13267fc1a946461`
  on branch `Updated-Master`.

Files:

- `MLPA_data_summarized_wo_siteblocks.csv`: processed South Coast/California
  MPA subtidal data used by `code/R/15_methods_comparison.R`.
- `CS_cummulative_intensity_1km.rds`: cold-spell cumulative-intensity grid used
  by `code/R/16_environmental_moderators.R`,
  `code/R/18_mpa_effectiveness_predictors.R`, and the concordance gate.
- `human_gravity_for_kelp_patches.csv`: population-pressure grid used by
  `code/R/extract_kelp_env_covariates.R` when regenerating
  `data/per_mpa_kelp_env.csv`.

Only these small derived/comparator files are tracked here. The full Kumagai
repository and Zenodo snapshot are not vendored into this analysis repo.
