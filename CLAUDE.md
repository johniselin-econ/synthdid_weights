# CLAUDE.md

The **`wsynthdid`** R package: a fork of
[`synthdid`](https://github.com/synth-inference/synthdid) (Arkhangelsky et
al. 2021) extended with **weighted SDID** — researcher-specified weights on
treated units — plus staggered adoption, weighted variance estimators, and
event-study decomposition.

The paper, supplement, and full replication package that develop and apply the
method (ACA Medicaid expansion, county mortality) live in a **separate repo**,
[`johniselin-econ/weighted-sdid`](https://github.com/johniselin-econ/weighted-sdid),
which depends on this package (pinned by commit via `renv`). **This repo is the
estimator only.** The GitHub repo was renamed from `synthdid_weights` to
`wsynthdid` on 2026-07-11 (the old URL redirects); `origin` is SSH under
`johniselin-econ`. `upstream` points at the original `synth-inference/synthdid`.

## Layout

Standard R package: `R/`, `man/`, `tests/`, `vignettes/`, `data/`,
`DESCRIPTION`, `NAMESPACE`. Package name is `wsynthdid`; the exported functions
keep their `synthdid_*` names (the original API is fully preserved — weighted
functions are additive).

## Weighted entry points

- `synthdid_estimate_weighted()` — weighted SDID; stores `cluster` for automatic
  use by `vcov()`. `sc_estimate_weighted()`, `did_estimate_weighted()` are the
  SC/DiD boundary cases. Uniform weights reproduce standard SDID (Proposition 0).
- Weighted `vcov.synthdid_estimate_weighted` — bootstrap/jackknife/placebo with
  weight renormalization; `cluster = NULL` forces unit-level resampling.
- `synthdid_event_study()` — event-study decomposition, cluster-robust bands
  when the estimate carries a stored cluster.
- **Staggered:** `synthdid_estimate_staggered()` in `R/staggered.R` (structural
  sibling of `R/stratified.R`) — loops over adoption cohorts, runs the *block*
  `synthdid_estimate_weighted()` within each against not-yet-treated
  (`control="notyet"`) or never-treated-in-window (`control="never"`) controls,
  and aggregates by treated-weight share (`cohort.weight="treated.share"`, the
  estimand τ^{ω̃,stag}) or treated-periods (reproduces Clarke et al. 2023 at
  uniform weights). Cluster/unit bootstrap `vcov` with per-draw full refit.
  Reduces to the block estimator for a single cohort. Requires `min.controls`
  ≥ 2 (a single never-treated unit gives a degenerate SC).
- **Stratified:** `synthdid_estimate_stratified()` — weighted SDID with the
  donor pool restricted to same-stratum controls, aggregated by stratum share.

## Tests

`tests/testthat/`: `test_weighted.R`, `test_staggered.R` (incl. a **golden test**
reproducing Stata `sdid`'s staggered aggregation to 6.4e-8), `test_remedies.R`,
`test_against_reference.R`, `test_synthdid.R`, `test_utils.R`.

## Gotchas (this Windows box)

- `devtools::load_all()` / `roxygen2::roxygenise()` **segfault** here (pkgload).
  To test package code, `source()` the R files directly (`utils.R`, `solver.R`,
  `synthdid.R`, then the target file) in a plain `Rscript`; edit `NAMESPACE` and
  `man/*.Rd` **by hand**. Plain `Rscript` with
  `RENV_CONFIG_AUTOLOADER_ENABLED=FALSE` usually works but R also segfaults
  **intermittently** on plain calls — retry, or verify via a background run's log.
- `R CMD check` and any `library(wsynthdid)` load should be verified on a healthy
  environment (e.g. the cluster), not this box.
- Never `git add -A` blind: a downloaded reference text (`w33719.txt`) was once
  swept in that way; it's now gitignored.
