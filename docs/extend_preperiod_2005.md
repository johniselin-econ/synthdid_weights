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

## Unit-specific linear detrending (local, point estimates, 2026-07-01)

Two-step: per county, fit y ~ time by OLS on pre-onset years only, subtract
the fitted unit trend from all years (extrapolating through the post
period), run weighted SDID on residuals. Slopes are never fit on post-onset
data, so dynamic treatment effects are preserved (avoids the Wolfers-2006
jointly-estimated-trends problem).

| SDID | raw | detrended |
|---|---|---|
| Headline popw, T0 = 5 | −17.45 | **−2.29** |
| Headline popw, T0 = 9 | −22.45 | **−2.78** |
| Headline eq, T0 = 5 / T0 = 9 | −0.71 / +0.51 | +3.20 / −3.76 |
| Placebo popw 2012 (T0 = 5 panel, 3 slope pts) | −7.30 | −0.88 |
| Placebo popw 2009/10/11/12 (T0 = 9, 4–7 slope pts) | −16.7/−16.0/−15.1/−14.7 | −8.8/−2.5/−0.9/−0.3 |
| Placebo eq 2009–2012 (T0 = 9) | +2 to +4 | −4 to −9 |

Read: under the assumption that the differential trend continues linearly,
the population-weighted effect is ≈ **−2.5** (consistent across both
panels) and its placebos pass once slopes have ≥5 fitting points. The raw
popw headline is, on this evidence, almost entirely differential linear
trend. Detrending is NOT uniformly better: the equally-weighted side gets
noisier and worse (small-county slopes are imprecise and the opioid-era
trend is nonlinear; extrapolation injects error) — it is the right surgery
specifically for the big-county linear trend.

Caveats: point estimates only (a proper SE must re-estimate slopes inside
every bootstrap draw); linear continuation is an assumption (Rambachan-Roth
M = 0); if part of the pre-2014 decline is early coverage (LIHP), detrending
over-adjusts — though the pre-ACA-onset placebos argue the trend is mostly
secular.

Precedent (has this been tried?): yes in the SC family, apparently not for
SDID specifically. Powell (2018 wp; 2022 JBES minimum-wage paper; 2026 JAE
"Imperfect Synthetic Controls") proposes exactly this two-step —
unit-specific trends first, synthetic control on the residual/fitted
structure. Ferman & Pinto (2021 QE) analyze SC under imperfect fit and
recommend the demeaned SC (SDID's intercept is that) and assessing fit
after discarding diverging trends. Rambachan & Roth (2023 REStud) give the
inference framework for exactly this counterfactual (linear extrapolation
of pre-trends, their M = 0 case; M > 0 relaxes linearity). Xu (2017)
gsynth/IFE is the estimator-level alternative (latent factors absorb
size-correlated trends). A detrended *weighted SDID* row appears to be
novel and citable to Powell + Ferman-Pinto + Rambachan-Roth.

## Size-binned weighted SDID (local, point estimates, 2026-07-01)

Stratify all counties by pooled 2013 population (quartiles + top-decile
split), run weighted SDID within each bin (donor pool = same-bin controls
only), aggregate with treated-pop shares. Note the aggregate targets the
SAME estimand tau^(omega-tilde) — only the comparisons are restricted.
Extended panel.

| bin (2013 pop) | N1 / N0 | tau_popw | placebo 2010 | placebo 2012 |
|---|---|---|---|---|
| 1: <6.1k | 233 / 473 | +0.5 | +1.9 | −0.8 |
| 2: 6.1–14.5k | 296 / 410 | +3.0 | +8.6 | +1.6 |
| 3: 14.5–38.5k | 292 / 414 | +5.0 | +5.3 | +6.3 |
| 4: 38.5–119k | 198 / 225 | +2.7 | +2.1 | +5.4 |
| 5: >119k (75.9% of treated mass) | 162 / 121 | **−6.5** | −5.0 | −2.7 |
| **Aggregate** | | **−4.1** | **−2.9** | **−0.9** |

vs unrestricted popw: −22.45 with placebos −16.0 / −14.7.

Read: **most of the failure was cross-size donor contamination, not an
unmatchable treated group.** Forcing big-treated vs big-control comparisons
collapses the estimate to −4.1 and the aggregate placebos essentially pass.
Mechanism: SDID's ridge-regularized level-matching spreads omega-hat
diffusely across the (mostly small/mid) donor pool and matching pre-period
LEVELS does not select trend-comparable donors; size-stratification imposes
the comparability the objective doesn't find. Residual top-bin placebo
(−5.0 at onset 2010) suggests some remaining within-class treated-vs-control
divergence (or early-coverage contamination), but small. Middle bins show
their own POSITIVE placebos ≈ their "effects" (opioid-era trend) — nothing
causal there either.

Triangulation across surgeries (all targeting the popw estimand):
raw −22.5/−17.5 → binned −4.1 → detrended −2.5/−2.8 → our DLPS repro −6.2
→ BV published −11.4. The comparability-enforced estimates cluster at
−2.5 to −6.5. Precedent for binning: Abadie-L'Hour (2021 JASA) penalized
SC for disaggregated data (penalizes dissimilar donor matches — the
estimator-level version of stratification); Abadie (2021 JEL) donor-pool
comparability advice. Caveats: no SEs; top bin still weight-concentrated
(LA); bin boundaries arbitrary; N0 = 121 donors for the top bin.

## Status update (2026-07-02)

Both remedies are now package code (committed): `detrend = TRUE` on
`synthdid_estimate_weighted()` and `synthdid_estimate_stratified()`, with
tests (`tests/testthat/test_remedies.R`). `analyze_application.R` gained
three guarded blocks — `detrended_results`, `binned_results` (pooled-2013
quantile bins, 25/50/75/90 breaks), `event_study_draws` (RR sensitivity) —
all with state-clustered bootstrap SEs; the detrended bootstrap re-fits
unit slopes per draw and the binned bootstrap rebuilds strata per draw.

## DLPS placebo, bridge ladder, RR sensitivity (local, 2026-07-03)

**1. BV's DLPS design FAILS the same in-time falsification.** Holding BV's
propensity-score weights fixed (they are functions of pre-2014 data only)
and running the TWFE on pre-2014 windows with simulated onsets
(est (state-clustered SE)):

| window | onset | base | with controls |
|---|---|---|---|
| 2009–13 | 2011 | −9.83 (3.50) | −7.68 (3.08) |
| 2009–13 | 2012 | −7.22 (2.58) | −5.10 (2.05) |
| 2005–13 | 2009 | −11.68 (4.33) | −7.26 (3.40) |
| 2005–13 | 2010 | −13.56 (4.41) | −8.86 (2.88) |
| 2005–13 | 2011 | −13.64 (4.50) | −8.67 (2.76) |
| 2005–13 | 2012 | −12.19 (4.07) | −7.11 (2.42) |

All significant; extended-window placebo magnitudes rival the headline
(−9.6 base / −6.6 controls; BV published −11.36). Lagged-outcome balancing
does NOT deliver trend balance here — the published estimate inherits the
same differential trend. (If anything the test is biased TOWARD passing:
the PS's 2005–09 mortality predictors overlap the placebo windows.)

**2. Bridge ladder (2009–2017, est (SE)) — the asymmetric-weighting
insight.** L1 pop-weighted TWFE, all counties: **−5.24 (5.19)**; L2
complete-case −5.22; L3 + PS trim −4.67; L4 pop x IPW untrimmed −9.24;
L5 + trim (= BV base) −9.56 (5.05); L6 + controls (= BV preferred) −6.62
(4.14). Two lessons. (i) The package's weighted DID (−18.39) weights only
the TREATED side by population against an unweighted control mean;
symmetric population weighting (both sides, as any pop-weighted regression
does) already compares big-to-big and lands at −5.2 — right where the
binned SDID sits (−4.8). Much of the raw −17/−22 reflects the asymmetric
big-treated-vs-diffuse-control comparison. (ii) Within BV's design the
trim is innocuous (−5.2 -> −4.7); it is the IPW reweighting that adds
−5 of magnitude — and the placebo above says that reweighting does not
buy trend balance.

**3. RR-style (M = 0) sensitivity on the popw SDID event study** (500
state-clustered draws): raw post avg −17.45 (CI [−24.8, −7.7]);
pre-period differential slope −0.99/yr (SE 0.45); trend-adjusted effect
**−13.22 (CI [−19.9, −4.5])**. IMPORTANT CAVEAT: this adjusts far less
than unit-level detrending (−2.3) because the event-study pre-coefficients
are computed AFTER omega-hat/lambda-hat were optimized to fit the weighted
pre-trajectory — the estimator absorbs most of the trend into the fit, so
the visible pre-slope understates the violation. RR on a fit-optimized
event study is a LOWER BOUND on the needed adjustment; the in-time
placebos (fit on truncated windows) and unit-level detrending are the
sharper diagnostics. Worth a remark in the paper — this interaction
(pre-fit optimization deflates measured pre-trends) applies to any
SC-family event study fed into HonestDiD.

Updated triangulation (all popw): symmetric-size comparisons −4.7 to −5.2;
binned SDID −4.77 (3.93); detrended SDID −2.29 (4.35); DLPS after its own
placebo evidence is discounted, residual plausibly −2 to −4. Every
level-matched design fails the in-time test, INCLUDING the published one;
comparability-enforced/trend-adjusted estimates cluster at −2 to −5.

## Follow-ups

1. Cluster run above → `results_2005/` (do NOT commit over `results/`).
   A second default-panel run adds only the new remedy CSVs to `results/`
   (all other guards are satisfied), giving the T0 = 5 remedy SEs.
2. Long-panel event study (year effects 2005–2017) — dates the divergence
   year by year; the key figure for the Section 3.6 discussion.
3. Drop-all-California: point estimates done locally (table above; B=200
   local SEs in progress). For the paper, add as a Table-3 robustness row
   with the full B=500 state-clustered SE in the cluster run.
4. If the cluster run confirms the smoke test, Section 3.6 and the
   Discussion likely need to shift further toward "upper bound / trend
   bias" framing, and the trend-adjusted ≈ −10 reading may itself be
   generous.
