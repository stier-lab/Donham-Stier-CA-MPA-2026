# Repository Readiness Audit

Date: 2026-09-03

Scope: public analysis repository for the Donham-Stier California MPA kelp
analysis.

## Status

- Public GitHub repo verified: `stier-lab/sbc-2026-donham-kelp-mpa-cascade`.
- GitHub Pages verified absent; no `docs/` site is being served.
- Main analysis entry point remains `source(here::here("code", "R", "run_all.R"))`.
- Fast figure regeneration remains `source(here::here("code", "R", "run_figures_only.R"))`.
- Resilience concordance gate remains `code/R/24_resilience_pipeline_check.R`.

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

## Residual Risk

The current commit removes private/manuscript/PDF artifacts from the public repo
head. It does not rewrite Git history. If the public repository must be treated
as never having distributed those files, perform a deliberate history rewrite
with `git-filter-repo`, force-push, and re-clone downstream working copies.
