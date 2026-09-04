# Repository Readiness Audit

Date: 2026-09-04

Scope: public analysis repository for the Donham-Stier California MPA kelp
analysis.

## Status

- Public GitHub repo verified: `stier-lab/sbc-2026-donham-kelp-mpa-cascade`.
- GitHub Pages verified absent; no `docs/` site is being served.
- Main analysis entry point remains `source(here::here("code", "R", "run_all.R"))`.
- Fast figure regeneration remains `source(here::here("code", "R", "run_figures_only.R"))`.
- Resilience concordance gate remains `code/R/24_resilience_pipeline_check.R`.
- Reviewer-readiness rerun on 2026-09-04: the local full analysis pipeline
  completed, the resilience concordance gate passed, and a fresh public GitHub
  clone completed `run_all.R`, reran `24_resilience_pipeline_check.R`, and
  verified all five main figure PDFs.
- Figure 1 no longer requires the ignored raw California MPA shapefile in a
  clean clone. The figure script uses the raw shapefile when present, otherwise
  falls back to tracked cached MPA-boundary objects.
- Small public Kumagai comparator inputs are vendored under
  `data/external/kumagai2024/`, so the main pipeline no longer depends on a
  user-specific mirror for the method-vs-data and environmental-moderator gates.

## Cleanup Decisions

- Removed manuscript drafts and generated manuscript/SI HTML from the public
  tracked tree. Current manuscript construction lives in the private manuscript
  repo.
- Removed `literature/pdfs/` from the public tracked tree. Source-copy PDFs are
  not required for code review and are maintained in the private manuscript
  repo.
- Moved the unlinked Figure 3 icon experiment into local ignored archive:
  `local_archive/2026-09-03_public_repo_cleanup/`.
- Moved ignored root-level scratch files (`Rplots.pdf`, Firebase/excalidraw
  logs, and `.DS_Store`) into the same local archive.
- Removed tracked legacy code archives from the public head and preserved local
  copies under `local_archive/2026-09-03_public_repo_harsh_cleanup/`.
- Removed the old heatwave manuscript-text draft from the public head; the
  current manuscript prose lives in the manuscript repo.
- Replaced the placeholder Dryad/raw-data downloader with an explicit
  harmonized-data availability check.
- Removed home-directory fallbacks from the in-pipeline Kumagai comparison,
  environmental-moderator, effectiveness-predictor, and resilience check code.
  Optional external inputs must now come from repo-controlled files or explicit
  environment variables.
- Corrected the cross-taxa meta-regression implementation so all six Table 3
  `metafor::rma` models use Knapp-Hartung tests (`test = "knha"`), matching the
  manuscript Methods and refreshed Table 3 output.
- Added `docs/CONCORDANCE_AUDIT.md` as a reviewer-facing claim ledger.

## Residual Risk

The current commit removes private/manuscript/PDF artifacts from the public repo
head. It does not rewrite Git history. If the public repository must be treated
as never having distributed those files, perform a deliberate history rewrite
with `git-filter-repo`, force-push, and re-clone downstream working copies.

No open-source reuse license has been assigned. Keep that explicit until the
author team chooses a license/DOI/archive policy for the public release.
