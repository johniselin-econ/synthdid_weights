# Draft sketch — next revision

2026-07-03, branch `extend-preperiod-2005`. Converts `rework_plan_2026-07.md`
(+ addendum) into a section-level sketch. Numbers marked [T0=5] exist in
`results/` today; [SRV] marks placeholders awaiting the server build list
(`server_build_list_2026-07.md`). Voice-drafted passages are quoted; the
rest is paragraph-level bullets.

## Thesis sentence (recurs at every turn)

> Researcher-chosen treated-unit weights select an estimand — and with it an
> identification burden. Assumption A5 is a statement about the
> omega-tilde-weighted treated factor structure: it can hold at 1/N1 and
> fail at omega-tilde in the same panel with the same estimator.

## Frozen exhibit list

Main paper (9 exhibits, +2 vs current):
- T1 two-unit example [exists] + trend row/paragraph extension
- T2 main results, T3 weight-spec robustness, T4 BV comparison [exist]
- T5 omega size bins [exists, added 2026-07-03]
- T6 (NEW) diagnostics & remedies [T0=5: app_estimates + detrended_results
  + binned_results — all present]
- F1 event studies, F2 size binscatter, F3 urban binscatter, F4 placebos
  [exist]
- MC tables (\ref'd) + trend-cell block [SRV-1]

Supplement adds: B.6 (variant theory notes), E.6 (extended-panel and
reconciliation section: placebo grid [SRV-4], drop-CA [docs], per-bin
detail [exists], DLPS placebo [SRV-2], ladder [SRV-3], RR sensitivity
[press-side chunk], long event study [SRV-4]).

## Title / abstract

Title: keep "Weighted Synthetic Difference-in-Differences". (Option if a
subtitle is wanted: "…: Estimands, Diagnostics, and Trend-Robust
Variants" — decide late.)

Abstract draft (replaces sentence 5 onward of current):

> We extend synthetic difference-in-differences (SDID) to researcher-
> specified treated-unit weights, implemented via a single modification of
> the collapsed form of Arkhangelsky et al. (2021), with adapted bootstrap,
> jackknife, and placebo variance estimators (including cluster-robust
> versions) and asymptotic theory under an effective-treated-sample
> condition. Because the weights define both the estimand and the
> identification burden — the oracle-weights assumption must hold for the
> weighted treated composite — we pair the estimator with estimand-specific
> diagnostics and two trend-robust variants: a unit-detrended and a
> size-stratified weighted SDID. Applying the toolkit to the county-level
> mortality effects of the 2014 ACA Medicaid expansion, the equally
> weighted estimate is a precise null while the population-weighted
> estimate is a reduction of roughly 17 deaths per 100,000 — but the
> population-weighted target fails its in-time falsification: the estimate
> loads on a secular mortality decline in large urban counties that
> convex, level-matched comparisons cannot remove, and that also
> contaminates the published propensity-score benchmark. Trend-robust
> estimates cluster at 2–5 fewer deaths per 100,000 and cannot be
> distinguished from zero. The weighting choice thus moves the answer
> twice — once through the estimand, once through the assumptions — and
> the weighted toolkit is what makes both moves visible.

## Section 1 — Introduction

- P1–P2 (SDID background, aggregation literature): as-is.
- P3 contributions: renumber to four. (1) estimator; (2) asymptotics
  (W1–W3, cluster version); (3) inference + the paired-bootstrap
  divergence diagnostic; (4) NEW: estimand-specific diagnostics and two
  trend-robust variants (detrended; stratified), motivated by the
  observation that A5 is indexed by the weights. Keep each to 2–3
  sentences; (4) cites Powell 2022/2026, Ferman–Pinto 2021,
  Rambachan–Roth 2023, Abadie–L'Hour 2021.
- P4 application: rewrite. Sequence: divergence (0 vs −17.5, both SEs);
  the divergence is real (paired bootstrap) — but the population-weighted
  target fails estimand-specific placebos; diagnostics localize the
  failure (size-blind level matching, Table 5); trend-robust variants
  give −2 to −5, indistinguishable from zero; the same falsification
  fails for the published DLPS benchmark, whose −11.36 sits at the edge
  of the trend-robust cluster. One sentence of the lesson (thesis).
- P5 roadmap: add Section 3.7.
- DELETE from current draft: the 2026-07-01 "upper bound ≈ −10" sentence
  (intro) — superseded.

## Section 2 — The Weighted SDID Estimator

- 2.1 Setup: as-is.
- 2.2 Weighted modification + two-unit example: as-is, then EXTEND the
  example paragraph/Table 1 with the trend case: give unit 1 a
  differential trend delta per period; bias of the weighted estimand =
  cov(omega-tilde_j, delta_j) x (post−pre time gap), by the same
  covariance algebra. One display or one table row. Purpose: the reader
  meets "trend heterogeneity correlated with weights" on page 6, so 3.6
  is foreshadowed, not a surprise.
- 2.3 Asymptotics: unchanged through Prop 3; after the A5 gloss insert
  named remark, draft:

  > **Remark (Identification is estimand-specific).** A5 is a statement
  > about the omega-tilde-weighted treated composite. Changing
  > omega-tilde changes which factor loadings must be matched by a convex
  > combination of controls, so A5 can hold for one weight vector and fail
  > for another in the same panel: the identification burden travels with
  > the estimand. Section 2.4's diagnostics should therefore be run at the
  > weights the researcher reports, and Section 2.5 provides estimators
  > for the leading failure mode, a trend factor whose loadings are
  > correlated with the weights.

- 2.4 Inference: bootstrap/jackknife/placebo paragraphs as-is. Expand
  "Diagnosing estimand divergence" into "Diagnosing estimand divergence
  and estimand-specific validity": step (i) paired bootstrap (as now);
  step (ii) rerun in-time placebos UNDER omega-tilde. One sentence on why
  equal-weight placebos certify only the equal-weight estimand.
- 2.5 (NEW) Trend-robust variants (~0.75 pp). Opening draft:

  > When the estimand-specific placebo fails, the leading interpretation
  > is a trend factor whose loadings covary with the weights. Two
  > variants of the weighted estimator target it directly, each keeping
  > the estimand tau^(omega-tilde) unchanged.

  Then one display + 3–4 sentences each:
  - Detrended: y-tilde_it = y_it − a-hat_i − b-hat_i t with (a,b) fit by
    unit OLS on PRE-onset years only; weighted SDID on y-tilde. Pre-only
    fitting preserves dynamic effects (contrast Wolfers 2006); this is
    the Rambachan–Roth M = 0 counterfactual and the SDID analog of
    Powell's two-step SC; Ferman–Pinto's demeaning is the intercept SDID
    already has. Bootstrap re-fits slopes per draw. API:
    `synthdid_estimate_weighted(..., detrend = TRUE)`.
  - Stratified: partition units on an observable that predicts trend
    loadings; weighted SDID within strata (donors restricted to the
    stratum); aggregate with omega-tilde stratum shares — algebraically
    the same estimand. Why it can succeed where omega-hat fails:
    ridge-regularized level matching is size-blind (forward-ref Table 5);
    stratification imposes the comparability the objective does not seek.
    Abadie–L'Hour penalization is the estimator-level analog. API:
    `synthdid_estimate_stratified()`.
  - Closing sentence: both are diagnostics-driven; formal statements and
    the slope-noise/bin-boundary caveats in Supplement B.6.

## Section 3 — Application

- 3.1 Setting/data: as-is. 3.2 Main results: as-is + closing framing
  sentence: "Both columns estimate well-defined estimands under the same
  assumptions; Sections 3.6–3.7 test whether those assumptions hold for
  each estimand." (Delete nothing else.)
- 3.3 heterogeneity, 3.4 robustness, 3.5 BV comparison: as-is; in 3.5 add
  a forward reference: "Section 3.7 subjects this benchmark to the same
  falsification as our own estimates."
- 3.6 RENAMED "Estimand-Specific Diagnostics". Paragraph flow (P1–P3
  exist in working tree, P4 replaces the stale '≈ −10' paragraph):
  1. [keep] In-time placebo setup; eq passes; popw fails (−7.7/−7.3, CIs
     exclude zero).
  2. [keep] Mechanism: Figure 1 pre-gap; convex hull + simplex
     no-extrapolation; the ≥2.5-year time gap; A5(omega-tilde) violation
     specific to the popw target.
  3. [keep] Table 5 omega size bins: emergent asymmetry (added
     2026-07-03).
  4. [NEW, replaces '≈ −10' para] Corroboration triad, one compact
     paragraph + supplement pointer: (i) weighted DID fails identically
     and SDID mitigates — failure tracks weights, not estimator; (ii)
     extended 2005-pre-period placebos are −15 to −17 at every onset
     including pre-ACA ones [SRV-4 for CIs; points in docs] — the
     early-coverage (LIHP) story cannot carry them, and the T0 = 9
     headline GROWS (trend accumulates); (iii) not CA-specific (drop-CA
     ratio unchanged). Keep @frean2017 for the LIHP sentence.
  5. [keep] In-space placebo + its two limitations.
- 3.7 (NEW) "Trend-Robust Estimates and Reconciliation" (~1 page):
  1. Table 6: raw popw −17.45 (4.66) / detrended −2.29 (4.35) /
     stratified −4.77 (3.93), with their placebos [T0=5]. Topic sentence
     draft: "Applying the Section 2.5 variants, the population-weighted
     effect that survives trend adjustment is a fifth to a quarter of the
     raw estimate and statistically indistinguishable from zero."
  2. Per-bin sentence: effect concentrated in the top decile (−6.9
     (4.09), marginal); middle bins' positive estimates match their own
     positive placebos.
  3. Reconciliation: BV's DLPS fails the same falsification — fixed
     weights AND full re-estimation at each onset [SRV-2 for committed
     CSV; numbers in docs]; the ladder shows symmetric population
     weighting alone lands at −5.2 and IPW adds −5 without buying trend
     balance [SRV-3]; our public-data DLPS repro (−6.62 (4.14)) covers
     the trend-robust cluster; BV's published −11.36 sits at its edge.
     Frame respectfully; one sentence that this demonstrates the
     workflow, not a refutation.
  4. RR sensitivity + lower-bound remark: M = 0 on the SDID event study
     gives −13.2 [press-side chunk], but pre-fit optimization absorbs
     trend into omega-hat/lambda-hat, so event-study-based sensitivity
     is a lower bound on the needed adjustment — placebos on truncated
     windows and unit detrending are the sharper tools. (Candidate for
     supplement if 3.7 runs long; keep one sentence in main text either
     way.)
  5. Closing lesson: "Equal weighting passed every test here not because
     the panel is clean but because two opposing size-correlated trends
     net out. The diagnostics must be run at the estimand one reports."

## Section 4 — Monte Carlo

- 4.1 design: add gamma_trend to the DGP paragraph [SRV-1].
- 4.2 results: existing patterns + one new block/table rows: under
  gamma_trend > 0, raw weighted SDID and weighted DID inherit the same
  bias, equal-weight variants do not, the in-time placebo detects it,
  detrended recovers the target, stratified recovers most [SRV-1
  placeholders]. 3–4 sentences.

## Section 5 — Discussion

- P1 (asymptotics lose bite): as-is.
- P2: rewrite around the workflow: choose omega-tilde -> paired-bootstrap
  divergence -> estimand-specific placebos -> trend-robust variants.
  Replace the 2026-07-01 "upper bound" insert with the remedies result and
  the twice-moved-answer lesson.
- P3 caveats: keep three existing; add fourth: remedies were chosen after
  seeing the diagnostics; the workflow is the ex-ante prescription. Add
  BV sentence: the design cannot trend-robustly distinguish 0 from −11 —
  what the toolkit revealed, not a refutation.

## Supplement

- B.6 (NEW): formal notes. Detrended: A1 with unit-specific linear
  trends; slope estimation adds an O_p(extrapolation-window / T0^{3/2})
  noise term absorbed by the per-draw-slope bootstrap; Props 1–3 carry.
  Stratified: A5-within-stratum implies the aggregate result; aggregation
  is the Prop 3 cluster algebra with strata as clusters.
- E.6 (NEW) empirical section, in order: extended-panel construction +
  validation; full placebo grid w/ CIs [SRV-4]; drop-CA table; per-bin
  table (binned_results.csv); DLPS placebo (fixed + re-estimated)
  [SRV-2]; ladder [SRV-3]; RR chunk from event_study_draws.csv with the
  lower-bound remark; long event study figure [SRV-4].
- E.1: cross-ref sentence already updated 2026-07-01; renumber "Section
  3.6" refs to 3.6/3.7 as applicable.
- A/B: no changes beyond B.6.

## Consistency sweep checklist (Phase 7, unchanged from plan)

Stale "≈ −10" texts (intro, 3.6 P4, Discussion) — all replaced above;
"significant at the 1% level" qualified; Table numbering re-swept after
T6 lands; audit-reproducibility after the single-vintage re-knit.
