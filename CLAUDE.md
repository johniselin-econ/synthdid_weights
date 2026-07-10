# CLAUDE.md

Fork of the `synthdid` R package extended with **weighted SDID**
(researcher-specified treated-unit weights), plus the paper and
supplement that develop the method and apply it to the 2014 ACA Medicaid
expansion (county mortality, building on Borgschulte & Vogler 2020).
Target: Elsevier journal submission.

## Repository state (consolidated 2026-07-05)

Everything lives on a **single branch, `master`** (@ `d3e7dfc`,
"single-vintage cluster rebuild"). All feature branches were fast-
forwarded/merged into master and deleted:

- `extend-preperiod-2005` — the 2005 pre-period extension + trend-robust
  variants + DLPS placebo/ladder + MC trend cell; was the furthest-ahead
  branch and became master via a clean fast-forward.
- `pipeline-restructure`, `portable-build-renv` — already contained in
  master's history; deleted as redundant.
- `erica_test` — the **Jan-2026 ancestral HTML working draft**
  (`vignettes/sdid_weights_paper_html.Rmd`). Reviewed 2026-07-05 and
  found to be a **strict subset** of the current `paper/` + supplement:
  its three-solutions framing is already Appendix "Alternative Approaches
  to Weighted Estimation" (with a Sol 1/2/3 Monte Carlo horserace), its
  DLPS/helper code is superseded by the package + `analyze_application.R`,
  and its §2.7 inference claim (rate √N_co, variance driven by N_co) is
  **inconsistent with** the adopted theory (CLT at √N₁ᵉᶠᶠ, driven by the
  effective number of *treated* units) — do NOT reintroduce it. Nothing
  to port; branch deleted (recoverable via reflog). Erica remains a
  credited coauthor.

Only `origin/master` remains on the fork; `upstream/*` (original
synthdid) is untouched.

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

Staggered adoption: `synthdid_estimate_staggered()` in `R/staggered.R`
(structural sibling of `R/stratified.R`) — loops over adoption cohorts,
runs the *block* `synthdid_estimate_weighted()` within each against
not-yet-treated (`control="notyet"`) or never-treated-in-window
(`control="never"`, the Stata `sdid` / golden-test pool) controls, and
aggregates by treated-weight share (`cohort.weight="treated.share"`, the
estimand τ^{ω̃,stag}) or treated-periods (reproduces Clarke et al. 2023 at
uniform weights). `vcov` is cluster/unit bootstrap with per-draw full
refit. Reduces to the block estimator for a single cohort (tested). See
`docs/reorient_plan_2026-07-10.md` Part 2. **Stata golden test DONE
(`docs/gold_test_staggered.md`):** treated-periods aggregation reproduces
Stata `sdid` on the Jones Table-5 panel to 6.4e-08 (1.799401), single-cohort
== block is exact. Jones has ONE never-treated unit -> `control="never"`
gives N0=1 (degenerate SC, no donor pool); the estimator refuses it cleanly
via `min.controls=2` and we did NOT patch core `synthdid.R` for N0=1.
STILL TODO: supplement proofs (per-cohort A5 / effective-N / not-yet-treated
validity / reductions / bootstrap), and a `gamma_trend` generating block in
`scripts/run_mc_simulations.R`.

## Environment gotcha (this Windows box)

`devtools::load_all()` / `roxygen2::roxygenise()` **segfault** here
(pkgload issue). To test package code, `source()` the R files directly
(`utils.R`, `solver.R`, `synthdid.R`, then the target file) in a plain
`Rscript`; NAMESPACE + man pages must be edited by hand. Plain `Rscript`
with `RENV_CONFIG_AUTOLOADER_ENABLED=FALSE` works fine.

## Docs to read before working

- `docs/TODO.md` — single live todo list.
- `docs/submission_review_checklist.md` — status of every April
  coarse-review item.
- `docs/econometric_review_2026-06-10.md` — referee-style review of the
  math/proofs (major items M1–M5 are pre-submission blockers).
- `docs/coarse_review_2026-04-17.md` — the April external review.

## Running R on this cluster

R is **not on `PATH` by default** (Lmod environment modules). Load the module in
the **same shell command** as `Rscript` — module state does not persist across
separate shell invocations (and each Claude Bash tool call is a fresh shell):

```bash
module load R/4.4.2-gfbf-2024a
Rscript your_script.R
```

This gives **R 4.4.2** (matches the Tariff-Rate-Tracker / Tariff-Model manifests)
with `tidyverse`, `scales`, `lubridate`, `here`, `patchwork`, `openxlsx`,
`readxl`, `arrow`, and the rest of the CRAN/Bioconductor bundle already attached.

### Gotchas

- Run the `module load` as the **sole/first** module op in a clean shell. It
  works fine that way (it pulls in `R-bare` + `R-bundle-CRAN` + Bioconductor via
  Lmod `depends_on`). It only appears to "fail" if you `module purge` first or
  stack another R module in the same accumulated shell.
- `~/r_libs_4.4` sits **first** on `.libPaths()`. Its hand-built **Arrow 24 lacks
  zstd** and breaks parquet builds; the module's Arrow 17 has zstd. For
  Arrow/parquet work, don't let the personal lib shadow the module. (CSV-only
  builds are unaffected.)
- Fallback if the meta-module ever won't load: `module load
  R/4.4.2-gfbf-2024a-bare` (Rscript + toolchain only), then `export
  R_LIBS_SITE=/apps/software/2024a/software/R-bundle-CRAN/2024.11-foss-2024a` for
  the CRAN packages. Do **not** `module load` the CRAN bundle on top of `-bare`
  (foss/gfbf toolchain conflict drops Rscript off `PATH`).

### Heavy jobs

Light builds run in seconds interactively. For memory-heavy rebuilds, wrap the
same `module load` line in an `sbatch` script and submit via Slurm.
