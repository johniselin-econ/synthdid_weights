# Rework plan: fold the trend findings into the weights contribution

2026-07-01, branch `extend-preperiod-2005`. Companion to
`docs/extend_preperiod_2005.md` (all numbers cited here are from that file's
local point-estimate runs; cluster SEs pending).

## The framing that prevents the tangent

One sentence, stated early and repeated at each turn of the paper:

> Researcher-chosen treated-unit weights select an estimand — and with it an
> identification burden. A5 is a condition on the omega-tilde-WEIGHTED
> treated factor structure: it can hold at 1/N1 and fail at omega-tilde in
> the same panel with the same estimator.

Under this framing nothing new is an ACA detour. The diagnostics
(estimand-specific placebos), the remedies (detrended and stratified
weighted SDID), and the ACA findings are all consequences of taking the
weights seriously: if you change what you estimate, you must re-verify what
must be true — at the estimand you chose. The application is then the
workflow demonstrated end-to-end, not a policy paper bolted on:

1. Divergence test: estimands differ (0 vs −17.5).
2. Estimand-specific placebos: the equal-weight target passes, the
   population-weight target fails — same panel, same estimator.
3. Diagnosis: a size-correlated secular trend; equal weighting doesn't
   solve it, it CONCEALS it (small-county positive trend nets against
   big-county negative trend).
4. Remedies within the weighted framework converge: detrended −2.5/−2.8,
   stratified −4.1, both placebo-clean; same estimand throughout.
5. Reconciliation: comparability-enforced estimates cluster at −2.5 to
   −6.5; our DLPS repro (−6.2, SE 4.3) covers the cluster and BV's −11.4
   sits at its edge.

The application's headline lesson changes from "Medicaid saved lives in
big counties" to: **the weighting choice changed the answer twice — once
through the estimand, once through the identification burden — and the
weighted toolkit is what caught it.** That is a WEIGHTS lesson.

## Decisions to settle first

- **Primary panel stays 2009–2017** (BV comparability; keeps existing
  results/ vintage valid). The 2005–2017 panel enters as the *diagnostic*
  panel, introduced inside the diagnostics section and detailed in the
  supplement. Do not re-base Tables 2–4 on it.
- Main text gets at most ONE new exhibit (the diagnostics-and-remedies
  table). Everything else — extended-panel construction, drop-CA, per-bin
  results, DLPS placebo, bridge ladder, sensitivity — goes to the
  supplement.
- Tone for the popw headline everywhere: "the estimand under SDID's
  assumptions," with a forward pointer to 3.6. Never "the effect."

## Phase 1 — Theory (Section 2; ~1 page of additions)

- **2.2 two-unit example**: add a short paragraph (or a second panel row
  under Table 1): give the large unit a differential trend delta. Bias of
  the weighted estimand = cov(omega-tilde, delta_j) x (post−pre time gap).
  Trend heterogeneity correlated with weights biases the weighted target
  exactly where effect heterogeneity moves it. Sets up 3.6 from page 6.
- **2.3, after the A5 gloss**: named remark *"Identification is
  estimand-specific."* A5 = A5(omega-tilde); the parallel of the existing
  "when the asymptotics lose their bite" paragraph, but for bias rather
  than inference. Two or three sentences + forward refs to 2.5 and 3.6.
- **2.4 diagnostics paragraph** becomes two-step workflow: (i) paired
  bootstrap divergence test (as now); (ii) if estimands diverge, rerun the
  falsification battery UNDER omega-tilde — placebos at the default
  weights certify the default estimand only.
- **NEW 2.5 "Trend-robust variants"** (half page, one display each):
  - *Detrended weighted SDID*: per-unit OLS trend on pre-onset years only,
    subtract from all years, run weighted SDID on residuals. Pre-only
    fitting preserves dynamic effects (contrast Wolfers 2006). Cite
    Powell (2018 wp / 2022 JBES / 2026 JAE), Ferman & Pinto (2021 QE)
    (demeaning = SDID's intercept; detrending is the next order),
    Rambachan & Roth (2023) (this is their M = 0 counterfactual).
  - *Stratified weighted SDID*: partition units on an observable that
    predicts trends; weighted SDID within strata (donor pool = same
    stratum); aggregate with omega-tilde stratum shares. SAME estimand —
    only the comparisons are restricted. Cite Abadie & L'Hour (2021 JASA),
    Abadie (2021 JEL). Note why it can succeed where omega-hat fails:
    ridge-regularized LEVEL matching does not select trend-comparable
    donors; stratification imposes comparability.
  - Formal statements (props inherit within stratum / under the detrended
    factor model) -> supplement B.6. Keep 2.5 practitioner-level.
- Bib adds: powell2022 (JBES), powell2026 (JAE), ferman2021 (QE),
  rambachan2023 (REStud), abadie_lhour2021 (JASA), abadie2021 (JEL),
  wolfers2006 (AER).

## Phase 2 — Application (Section 3)

- **3.2/intro framing sentence** after the Table 2 discussion: "Both
  numbers are estimates of different estimands under the same assumptions;
  Section 3.6 asks whether those assumptions hold for each estimand."
- **3.6 rewritten** as *"Estimand-Specific Diagnostics and Trend-Robust
  Estimates"* (replaces the 2026-07-01 uncommitted rewrite, whose
  "trend-adjusted ~ −10" reading is now stale). Structure:
  1. In-time placebos: eq passes, popw fails (as computed). One sentence
     on the 2x2: weighted DID fails identically, SDID mitigates — failure
     tracks the WEIGHTS, not the estimator.
  2. Extended pre-period (one paragraph + supplement pointer): placebos at
     pre-ACA onsets are −15 to −17 (kills the early-coverage-only story);
     the T0 = 9 headline GROWS to −22.5 (unremoved trend accumulates —
     the simplex time weights cannot extrapolate); not CA-specific
     (drop-CA leaves the placebo/headline ratio unchanged).
  3. Remedies table (the one new main-text exhibit). Rows: raw popw,
     detrended popw, stratified popw (+ eq comparisons); columns:
     estimate, placebo(s), [SEs when cluster run lands]. Point estimates:
     −17.5/−22.5 raw -> −2.3/−2.8 detrended -> −4.1 stratified, remedies
     placebo-clean. One caveat sentence each (linear-continuation
     assumption; bin-boundary choice + thin top-bin donors).
  4. Reconciliation with BV (one paragraph): triangulation −2.5 to −6.5;
     our public-data DLPS repro −6.2 (SE 4.3) covers the cluster and
     nearly covers BV; frame the residual as a sensitivity statement
     (Rambachan-Roth: how much trend nonlinearity must one allow before
     −11.4 or 0 is admissible) rather than a point-estimate contest.
     Full DLPS placebo + DID->DLPS bridge ladder -> supplement.
  5. Closing lesson: equal weighting passed every test here not because
     the panel is clean but because opposing trends net out — the
     diagnostics must be run at the estimand you report.
- In-space placebo paragraph: keep, with its limitations text (already
  drafted).
- Table numbering ripples again (remedies table lands after Table 4);
  MC tables are \ref'd and safe — sweep hardcoded refs in both docs.

## Phase 3 — Monte Carlo (Section 4)

- Add a **trend cell**: gamma_trend x (s_j/s-bar − 1) x t added to L,
  violating A5 for the popw target only. Report: weighted SDID and
  weighted DID inherit the same bias; equal-weight variants unbiased;
  the in-time placebo detects it; detrended weighted SDID recovers the
  target; stratified recovers most. This is the structural proof that the
  failure is the DGP, not the weighted machinery. ~3 sentences + one
  table block; script change in run_mc_simulations.R (cluster, cheap
  relative to the existing grid).

## Phase 4 — Discussion (Section 5)

- Rebuild paragraph 2 around the workflow (choose omega-tilde -> divergence
  test -> estimand-specific placebos -> trend-robust variants). Replace the
  stale "upper bound / trend-adjusted −10" caveat with the remedies result.
- BV sentence: comparability-enforced estimates sit below the published
  point; the design cannot trend-robustly distinguish 0 from −11 — cast as
  what the weighted toolkit revealed, not as a refutation.
- Caveats list: add "stratification/detrending as estimand-preserving
  remedies chosen ex post; pre-registering the diagnostic workflow is the
  discipline" (referee-proofing).

## Phase 5 — Supplement

- **B.6** (new): formal notes for the two variants. Detrended: A1 with
  unit-specific linear trends; slopes fit pre-only; Props 1–3 carry with
  an added slope-estimation term (vanishes as T0 grows; bootstrap must
  re-fit slopes per draw). Stratified: A5-within-stratum implies the
  aggregate result; aggregation = cluster structure of Prop 3.
- **New empirical section** (E.6 or F): extended-panel construction +
  validation (identical county set, bit-identical shared rows, stable
  imputation share); full placebo grids (onsets 2009–2012, both panels);
  drop-CA table; per-bin table; DLPS in-time placebo; DID->DLPS bridge
  ladder; Rambachan-Roth sensitivity. Content largely exists in
  docs/extend_preperiod_2005.md — this is a port once SEs land.
- Update E.1 cross-references (in-time contrast sentence now cites 3.6).

## Phase 6 — Computation before re-knit

Local (can run now): DLPS in-time placebo; DID->DLPS bridge ladder.
Cluster (one results_2005 vintage): full extended-panel run (existing
blocks); NEW analyze_application.R blocks writing `detrended_results.csv`
and `binned_results.csv` with state-clustered bootstrap SEs — detrended
bootstrap MUST re-estimate unit slopes inside every draw; binned bootstrap
resamples states within stratum structure; MC trend cell; event-study
bootstrap draws saved for the RR sensitivity. Then delete paper caches and
re-knit both documents.

## Phase 7 — Consistency sweep

- Abstract: keep contributions 1–3; replace the application sentence pair:
  divergence stays the headline; add one clause that estimand-specific
  diagnostics reveal the population-weighted estimate rests on a failing
  assumption and trend-robust variants of the weighted estimator give
  −2.5 to −6.5. (User to finalize language.)
- Intro contribution paragraph: consider promoting the
  diagnostics-and-remedies workflow to a named contribution (3a) rather
  than a subsection surprise.
- Sweep stale claims: "significant at the 1% level" (qualify), yesterday's
  "upper bound ~ −10" texts in intro/3.6/discussion, E.1 contrast
  sentence, any "robust to placebo tests" residue.
- audit-reproducibility pass after re-knit.

## Addendum 2026-07-03: reconciliation branch resolved

The Phase-2 reconciliation paragraph was written with two branches
("if DLPS placebos pass / fail"). They FAIL — significantly, at every
onset, both windows (see extend_preperiod_2005.md). Consequences for the
plan:

- 3.6's reconciliation paragraph takes the stronger form: every
  level-matched design applied to this panel — weighted SDID, weighted
  DID, and BV's DLPS itself — fails the same in-time falsification; the
  trend-robust cluster is −2 to −5. Frame respectfully (BV lacked the
  2005-08 window's diagnostic power on their sample vintage; our repro of
  their design is within 1 SE of theirs) and as a demonstration of the
  workflow's value, not a refutation exercise.
- New remark needed (2.4 or 3.6): the asymmetric-weighting point. The
  weighted DID/SDID compare a pop-weighted treated average against a
  (diffuse) control combination; symmetric pop-weighted TWFE already
  compares big-to-big and lands at −5.2, near the stratified estimate.
  Practitioners reading "weighted DID −18" vs "pop-weighted regression
  −5" need this distinction spelled out — it is an estimand/comparison
  statement, not a contradiction.
- RR/HonestDiD caveat (one paragraph, likely supplement): event-study
  pre-coefficients computed after omega/lambda fit the pre-period
  UNDERSTATE trend violations (fit absorbs the trend), so M = 0 on the
  SDID event study (−13.2) is a lower-bound adjustment; placebos on
  truncated windows and unit detrending are the sharper tools. This is a
  publishable observation about SC-family + HonestDiD interactions.

## Risks / referee anticipation

- Wide SEs on remedies -> phrase as "cannot trend-robustly distinguish
  from zero," never as a precise null.
- "Too many estimators" -> main text carries exactly two remedies; the
  rest stays in the supplement.
- "Remedies chosen after seeing the data" -> own it in Discussion; the
  workflow is the prescription for doing this ex ante.
- Top-bin donor thinness (N0 = 121) and within-bin weight concentration
  (LA) -> report bin-boundary sensitivity in supplement.
