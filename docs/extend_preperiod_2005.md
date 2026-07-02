# Extended pre-period panel (2005–2017, T0 = 9)

Branch `extend-preperiod-2005`, 2026-07-01. Response to the in-time placebo
caveat in main-paper Section 3.6: with T0 = 5 the population-weighted
placebos (−7.7 / −7.3, CIs excluding zero) cannot distinguish a secular
differential trend in large urban counties from genuine pre-2014 coverage
effects (CA LIHP 2011–13). The raw WONDER exports already cover 2005+, so
the pre-period extends to 2005 with no new data pulls.

## What exists

- `paper/data/analysis_data_2005.csv` — balanced 2005–2017 panel,
  **identical county set to the canonical panel** (2,824 counties; 1,181
  treated / 1,643 control) and **bit-identical rates/populations on the
  shared 2009–2017 rows** (verified max |diff| = 0). Suppression-recovery
  feasibility holds for all state-years; imputed share is stable at ~8–9%
  in every year 2005–2017. Rates only (DUA-safe).
- `paper/data/seer_laus_covariates.csv` — regenerated with 2005–2008 rows
  added (insertions only; 2009–2017 rows unchanged). 2005–2008 coverage:
  1 county lacks SEER covariates, 11 lack unemp — same pattern as 2009+.
- Parameterization (all defaults reproduce the canonical run exactly):
  - `rebuild_mortality_panel.R`: `PANEL_START_YEAR`, `PANEL_END_YEAR`,
    `PANEL_OUT`.
  - `build_seer_laus_covariates.R`: `COV_START_YEAR`, `COV_END_YEAR`.
  - `analyze_application.R`: `ANALYZE_PANEL` (panel path),
    `ANALYZE_RESULTS` (results dir — keeps the committed `results/`
    vintage untouched), `ANALYZE_INTIME_YEARS` (comma-separated placebo
    onsets). The DLPS block always uses the canonical panel (it mirrors
    BV's design and should not move with the SDID panel).

## Cluster run

```sh
ANALYZE_PANEL=paper/data/analysis_data_2005.csv \
ANALYZE_RESULTS=results_2005 \
ANALYZE_INTIME_YEARS=2009,2010,2011,2012 \
sbatch scripts/slurm_run_analysis.sh
```

(`sbatch` exports the environment by default. `results_2005/` starts empty
so no `ANALYZE_FORCE` needed. T0 = 9 is derived from the panel — nothing
else to set. Skim the event-study block for hardcoded-year assumptions
before trusting its output; the placebo and estimate blocks are clean.)

## Local smoke test (point estimates, no SEs — 2026-07-01)

| Quantity | T0 = 5 (canonical) | T0 = 9 |
|---|---|---|
| SDID equal | −0.71 | +0.51 |
| SDID pop-weighted | −17.45 | **−22.45** |
| DID pop-weighted | −18.39 | −25.96 |
| In-time placebo, popw, onset 2009 (train 2005–08) | — | **−16.7** |
| onset 2010 (train 2005–09) | — | −16.0 |
| onset 2011 (train 2005–10) | −7.68 | −15.1 |
| onset 2012 (train 2005–11) | −7.30 | −14.7 |

**Preliminary read (needs the cluster SEs before concluding):** the longer
pre-period *strengthens* the trend-bias interpretation rather than
rescuing the estimate. Population-weighted placebos are ≈ −15 to −17 at
every onset, including 2009–2010 onsets trained entirely on 2005–2008 —
before ACA passage, so early-coverage (LIHP) cannot explain them. The
T0 = 9 headline *grows* (−22.5), consistent with unremoved trend
accumulating over a longer window rather than being differenced out.
A raw-gap diagnostic points the same way: the population-weighted
treated−control gap widened −5.4 per 100k over 2005–2010 (pre-ACA), and
that widening disappears when California's treated counties are excluded.
Caveats: no SEs yet; the equally-weighted placebos also drift positive
(+2 to +4) on the long panel; short placebo training windows may overfit;
and the placebo post-windows span the Great Recession.

## Drop-all-California robustness (local, point estimates, 2026-07-01)

All CA counties are treated, so dropping CA removes 58 treated counties
(~23% of treated mass; no-CA weight stats: N1 = 1,123, N1_eff = 133,
top-1 share 4.2% — Cook County replaces LA) and leaves the control pool
untouched.

| SDID pop-weighted | full | no-CA |
|---|---|---|
| Headline, T0 = 5 | −17.45 | −13.66 |
| Headline, T0 = 9 | −22.45 | −17.64 |
| Placebo onset 2011, T0 = 5 panel | −7.68 | −6.50 |
| Placebo onset 2012, T0 = 5 panel | −7.30 | −6.35 |
| Placebo onset 2009, T0 = 9 panel | −16.73 | −13.34 |
| Placebo onset 2010, T0 = 9 panel | −15.99 | −12.56 |
| Placebo onset 2011, T0 = 9 panel | −15.08 | −12.58 |
| Placebo onset 2012, T0 = 9 panel | −14.69 | −11.89 |

Read: **the placebo failure is NOT a California artifact.** Dropping CA
shaves ~20–25% off both the headline and the placebos, leaving the
placebo-to-headline ratio roughly unchanged (~45–50% at T0 = 5, ~70% at
T0 = 9). The differential pre-trend generalizes to non-CA large urban
treated counties in the SDID contrast, even though the *raw*
population-weighted gap widening 2005–2010 looked all-CA (the raw
diagnostic compares different objects: pop-weighted control mean vs.
SDID's ω̂-matched combo). Also note the equally-weighted placebos on the
T0 = 9 panel drift positive (+2 to +5, no-CA) — the small-county side has
its own (opposite-signed, opioid-era) differential trend that equal
weighting averages against.

## Follow-ups

1. Cluster run above → `results_2005/` (do NOT commit over `results/`).
2. Long-panel event study (year effects 2005–2017) — dates the divergence
   year by year; the key figure for the Section 3.6 discussion.
3. Drop-all-California: point estimates done locally (table above; B=200
   local SEs in progress). For the paper, add as a Table-3 robustness row
   with the full B=500 state-clustered SE in the cluster run.
4. If the cluster run confirms the smoke test, Section 3.6 and the
   Discussion likely need to shift further toward "upper bound / trend
   bias" framing, and the trend-adjusted ≈ −10 reading may itself be
   generous.
