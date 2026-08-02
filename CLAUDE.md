# CLAUDE.md

Fork of the `synthdid` R package extended with **weighted SDID**
(researcher-specified treated-unit weights), plus the paper and
supplement that develop the method and apply it to the 2014 ACA Medicaid
expansion (county mortality, building on Borgschulte & Vogler 2020).
Target: Elsevier journal submission.

> On the Yale HPC, cluster setup (R module load, Slurm, paths) lives in your user-level `~/.claude/CLAUDE.md`.

## Architecture: results factory vs. PDF press

Heavy computation and document rendering are fully decoupled, sharing
`results/` via git:

- `scripts/analyze_application.R` — ALL heavy computation (estimates,
  bootstrap/jackknife SEs, event studies, placebos, SC diagnostics, LOO,
  DLPS reproduction). Parallel (PSOCK, `SLURM_CPUS_PER_TASK`-aware),
  checkpointed, and **idempotent**: block-level guards skip a group when
  its CSVs exist; set `ANALYZE_FORCE=1` to recompute. Run on the cluster
  via `scripts/slurm_run_analysis.sh`; writes reduced CSVs to `results/`
  (committed, with `results/_manifest.csv` recording seed/SHA/reps).
- `paper/weighted_sdid_paper.Rmd` and `paper/weighted_sdid_supplement.Rmd`
  are **thin**: they read `results/*.csv` (as `../results/...`) and
  `paper/data/*.csv`, then format/plot. They render in ~0.4 min and never
  call synthdid functions in executed code.
- `paper/_build_all.R` — render stage only (`Rscript paper/_build_all.R
  paper|supplement` or no arg for both). Cross-OS pandoc detection +
  TinyTeX bootstrap.
- `scripts/00_run_all.R` — stage-gated orchestrator (data → analysis →
  mc → render); skips stages whose outputs exist.
- Monte Carlo sweeps run outside knitr: `scripts/run_mc_simulations.R` →
  `paper/data/mc_results.csv`; `scripts/run_solutions_mc.R` →
  `paper/data/mc_solutions_results.csv`. Seeding reproduces the
  sequential loop exactly regardless of worker count.

Data pipeline (raw pulls → `paper/data/analysis_data.csv`) is documented
in `scripts/README.md`, including run order and provenance. All inputs
are public; per the CDC WONDER DUA, death **counts** are never committed
(rates + suppression flags only).

## Gotchas

- **knitr caches do not watch data files.** After changing anything in
  `results/` or `paper/data/`, delete `paper/*_cache/` before knitting,
  or stale numbers will silently persist.
- **renv pins R 4.4.2** (the cluster version). The local machine runs
  4.5.2; bypass renv locally with
  `RENV_CONFIG_AUTOLOADER_ENABLED=FALSE` (uses the system library).
- The user frequently edits files externally — **always re-read a file
  before editing it**.
- `tau_sdid` in the Rmds is the *standard* `synthdid_estimate` (the
  robustness "equal" row); the main-results-table "equally weighted"
  cell (`sdid_eq`; Table 2 since the two-unit example became Table 1)
  is the weighted-uniform variant — they differ by ~2e-5
  (Proposition 0), so don't "fix" the discrepancy.
- Bootstrap SEs re-run on a different machine shift ~3% (Monte Carlo
  noise; RNG differs from `vcov()`'s although the resampling algorithm
  is identical). Point estimates are deterministic.
- Hardcoded numbers in prose have gone stale before — verify against
  `results/app_scalars.csv` rather than trusting the text (see
  `docs/econometric_review_2026-06-10.md`, item m6).
- kableExtra captions don't resolve `@citations`; use plain text.
- Use `ntile()` (not `cut(quantile(...))`) for binscatter binning:
  suppression-recovered counties create mass points that make raw
  quantile breaks repeat.
- PowerShell 5.1 `Set-Content -Encoding utf8` writes a BOM that crashes
  Rscript; write files with no-BOM UTF-8.

## Package development

Standard R package layout (`R/`, `man/`, `tests/`, `vignettes/`,
`DESCRIPTION`). Weighted entry points: `synthdid_estimate_weighted()`
(stores `cluster` for automatic use by `vcov()`), weighted
`vcov.synthdid_estimate` methods (bootstrap/jackknife/placebo with
weight renormalization; `cluster = NULL` forces unit-level resampling).
`tmpclaude-*` files at the root are gitignored scratch.

## Docs to read before working

- `docs/TODO.md` — single live todo list.
- `docs/submission_review_checklist.md` — status of every April
  coarse-review item.
- `docs/econometric_review_2026-06-10.md` — referee-style review of the
  math/proofs (major items M1–M5 are pre-submission blockers).
- `docs/coarse_review_2026-04-17.md` — the April external review.
