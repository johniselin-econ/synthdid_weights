# Server build list — pieces to build before the single-vintage run

2026-07-03. Companion to `draft_sketch_2026-07.md` (which marks where each
lands as [SRV-n]). Reference implementations for SRV-2/SRV-3 are committed
as `scripts/dev_*.R` — they reproduce the numbers in
`extend_preperiod_2005.md` exactly and only need porting into guarded
blocks. Delete the dev scripts after porting.

The paper/supplement chunks will be written against the exact CSV schemas
below — please keep column names as specified.

## SRV-1: Monte Carlo trend cell

New script (suggest `scripts/run_mc_trend.R`, same seeding/checkpoint
conventions as `run_mc_simulations.R`) writing
`paper/data/mc_trend_results.csv`.

DGP: as Section 4.1, plus a size-correlated trend in the UNTREATED
outcome of treated units only:
  L[N0+j, t] += gamma_trend * (s_j / s_bar - 1) * (t / T)
(controls unchanged, so no convex combination matches the deviation —
an A5(omega-tilde) violation that leaves A5(1/N1) approximately intact).
Grid: N1 = 20, gamma = 0.5, high HHI (3/N1), gamma_trend in {0, 0.5, 1},
200 sims. Per sim estimate SIX things on the same panel plus one placebo:

Columns: `gamma_trend, sim, true_att_equal, true_att_weighted,
est_eq, est_wt, est_did_wt, est_wt_detrended, est_wt_binned, plac_wt`

- est_wt_detrended: `synthdid_estimate_weighted(..., detrend = TRUE)`
- est_wt_binned: `synthdid_estimate_stratified()` on size quartiles
  (N1 = 20 is small — 4 strata by s_j quartile, donors by matching
  control "size" — assign control sizes from the same lognormal)
- plac_wt: weighted SDID with fake onset at T0/2 on the pre-period only
  (detects the violation; report alongside)

Paper 4.2 will report: bias of each estimator by gamma_trend + the share
of sims where |plac_wt| > 2 * (placebo SE proxy) — keep raw columns, the
Rmd computes the summaries.

## SRV-2: `dlps_placebo` block in analyze_application.R

Guard: `if (!have_all("dlps_placebo"))`, placed inside/after the DLPS
section (reuses its objects for the "fixed" variant). Reference code:
`scripts/dev_dlps_placebo_ladder.R` (fixed-weights variant) and
`scripts/dev_dlps_reestimated.R` (re-estimated variant).

Output `results/dlps_placebo.csv`, columns:
`window` ("2009-2013" | "2005-2013"), `onset` (int),
`variant` ("fixed" | "reestimated"), `spec` ("base" | "controls"),
`estimate`, `se`, `n`.

Rows: fixed x {2009-13: 2011,2012; 2005-13: 2009..2012} x both specs;
reestimated x same grid (2005-13 window uses ALL pre-onset mortality
years as PS predictors — the deliberately generous variant; keep the
NOTE about the fixed-vintage ACS/SAHIE covariates as a code comment).
Needs `paper/data/analysis_data_2005.csv` (committed). SEs are the
feols state-clustered ones (already in the reference code) — no
bootstrap needed.

## SRV-3: `dlps_ladder` block in analyze_application.R

Guard: `if (!have_all("dlps_ladder"))`. Reference:
`scripts/dev_dlps_placebo_ladder.R` part (2).

Output `results/dlps_ladder.csv`, columns:
`rung` ("L1".."L6"), `label`, `n_counties`, `estimate`, `se`.
L1 pop TWFE all; L2 complete-case; L3 PS-trimmed, pop weights;
L4 pop x IPW untrimmed; L5 pop x IPW trimmed (= BV base);
L6 + panel controls (= BV preferred).

## SRV-4: results_2005 run (no new code)

```sh
ANALYZE_PANEL=paper/data/analysis_data_2005.csv \
ANALYZE_RESULTS=results_2005 \
ANALYZE_INTIME_YEARS=2009,2010,2011,2012 \
sbatch scripts/slurm_run_analysis.sh
```

Commit `results_2005/` (new dir; safe). Supplement E.6 consumes:
placebo_intime (grid w/ CIs), event_studies (2005–2017 long event
study), detrended_results, binned_results, app_estimates/app_scalars.
Check the event-study block for hardcoded-year assumptions before
trusting the long figure.

## Then: the single-vintage sequence

1. Build SRV-1..3, push.
2. `ANALYZE_FORCE=1 sbatch scripts/slurm_run_analysis.sh` (fresh
   canonical `results/`, now including dlps_placebo/dlps_ladder/
   omega_size_bins) — one job.
3. SRV-4 job. 4. MC jobs (`run_mc_simulations.R`, `run_solutions_mc.R`,
   `run_mc_trend.R`).
5. Commit all results + manifests (each dir's `_manifest.csv` carries
   seed/SHA/reps = the vintage stamp), push.
6. Locally: pull, delete `paper/*_cache/`, re-knit both docs, run the
   audit-reproducibility pass, commit.

Expected diffs vs current committed results: bootstrap SE digits ~3%
(machine RNG); point estimates identical. NOT server-buildable: the raw
data stage (paper/data/raw is gitignored; committed panels are the data
vintage) and the renders (no pandoc/LaTeX on the cluster — press stays
local).

## NOT on the server list (press-side, mine)

- RR sensitivity chunk in the supplement (reads the committed
  `results/event_study_draws.csv`; reference `scripts/dev_rr_sensitivity.R`).
- Supplement B.6 + E.6 text; Section 2.5/3.6/3.7 prose per the sketch.
