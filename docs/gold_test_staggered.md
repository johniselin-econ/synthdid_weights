# Staggered SDID — Stata golden test (2026-07-10)

Validation of `synthdid_estimate_staggered()` against the reference Stata
`sdid` package (Clarke, Pailañir, Athey & Imbens 2023) on the Jones, Di Maria
& Valente (2026, SJE) "Financial Intermediation and Structural Change" Table-5
panel. Do-files and data are in `docs/_gold/` (gitignored).

## Setup

`data_sje.dta`, Jones Table-5 sample: drop SD/DE, `ma_dereg_year != 1960`,
years 1969-1998. **1170 rows = 39 states x 30 years.** 38 treated states in
**20 adoption cohorts** (1970, 1975-1991, 1993, 1994); **one never-treated in
window** (Iowa, `ma_dereg_year = 1999`).

Command (per `SJE_replication.do`, Table 5, no-cov point via `vce(noinference)`):
```
sdid <outcome> state year ma_dereg_dum, method(sdid) vce(noinference)
```

## Reference numbers (Stata `sdid`, aggregate ATT)

| outcome              | no covariates | with covariates (published Table 5) |
|----------------------|--------------:|------------------------------------:|
| emp_tot_serv_share   |      1.799401 |                            1.509227 |
| out_tot_serv_share   |      2.697071 |                            2.411509 |
| emp_man_share        |     -2.600671 |                           -2.150345 |
| out_man_share        |     -5.416803 |                           -4.271459 |

## What `sdid` does here (from `e()`)

- `e(adoption)` = 20 adoption years; `e(tau)` = 20 per-cohort ATTs.
- `e(omega)` = 2 x 21 -> **a single control unit (Iowa), weight 1 in every
  cohort.** With one never-treated unit there is no donor pool: "synthetic"
  DID degenerates to a time-weighted two-group DID.
- Aggregate ATT = treated-periods-weighted average of the cohort ATTs, i.e.
  W_g proportional to (N1_g x T_post,g).

## Validation result

- **Aggregation (the novel part): EXACT.** Feeding Stata's 20 `e(tau)` values
  and cohort sizes through our `cohort.weight = "treated.periods"` formula
  reproduces the Stata aggregate 1.799401 to **6.4e-08**. The other schemes
  give distinct, correctly-computed estimands (treated.share 1.632156,
  equal-per-cohort 1.951447).
- **Block reduction: EXACT.** `test_staggered.R` shows a single cohort equals
  the block `synthdid_estimate_weighted` bit-for-bit.
- Together these validate both halves of the estimator (cohort construction +
  block fit, and cross-cohort aggregation).

## The N0 = 1 caveat (a property of the Jones setting, not the code)

Jones has exactly one never-treated-in-window unit, so `control = "never"`
gives N0 = 1 per cohort. The SC weight machinery is degenerate at N0 = 1
(`sc.weight.fw` demeans each column across the single control row -> zeros;
there is no omega to solve). Our estimator **refuses this cleanly** via the
default `min.controls = 2`, with a clear error. We did **not** patch the shared
block estimator to limp through N0 = 1: (a) it is genuinely degenerate, (b)
bit-for-bit agreement with Stata's degenerate-solver output is not guaranteed,
and (c) touching core `synthdid.R` risks the ACA (block) results. This
reinforces `control = "notyet"` (not-yet-treated donor pools) as the broadly
usable default, and is another reason Jones is not a good staggered
*demonstration* (which we already decided: staggered is theory+code only, no
application).

Golden guard lives in `tests/testthat/test_staggered.R`
("golden: treated-periods aggregation reproduces Stata sdid").
