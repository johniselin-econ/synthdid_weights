# Reorientation plan (2026-07-10)

Three decisions, taken together (**settled 2026-07-10**):

1. **Paper**: pull back from V2's trend-robustness immersion to a
   *method-forward* draft organized around four beats —
   (A) re-weighting matters, (B) new SDID capacity, (C) ACA context/
   findings, (D) identification here is shaky (in-time placebo), but that
   is a property of the **test case**, not the weighting scheme.
   **ACA-only; method-first.**
2. **Code**: add **staggered** weighted SDID as **theory + code only** —
   package capability, supplement proofs, MC, and a validation golden test.
   **No Jones paper, no empirical application** (dropped).
3. **Repos**: split the R package (**`wsynthdid`**) out from the ACA paper.
   Two repos, not three.

Locked-in answers to the earlier open decisions:
- Jones application paper — **dropped**; staggered is capability only.
- §4 MC trend cell — **add, diagnosis-only**.
- Package name — **`wsynthdid`**.
- Staggered in paper #1 — appears as a *brief capability note in §2 + full
  treatment in the supplement* (it has nowhere else to live and it
  strengthens contribution B); **no main-text application**.

---

## Part 1 — Paper: the "part-way" draft

### The spine (A/B/C/D), one sentence each

- **A. Re-weighting matters.** Researcher-chosen treated-unit weights
  select an *estimand*; with size-correlated heterogeneity the equal-weight
  estimand and the population-weight estimand diverge sharply (ACA: 0 vs
  −17.5).
- **B. New SDID capacity.** A parsimonious weighted SDID (one change to the
  collapsed form) with full inference (bootstrap / jackknife / placebo,
  cluster-robust) and asymptotics under an effective-treated-sample
  condition. [Companion capability, paper #2: staggered adoption.]
- **C. ACA context / findings.** Applied to county mortality under the 2014
  Medicaid expansion (BV design); the divergence is real, driven by the
  weight not the estimator, concentrated in the largest urban counties.
- **D. Shaky identification is the test case, not the method.** The
  population-weighted target fails an in-time placebo: a treated-specific
  differential mortality trend in large counties that *no* convex,
  level-matched estimator removes — weighted SDID, weighted DID, **and BV's
  own DLPS**. This is a feature of the ACA county design, not of weighting.
  Recent triple-difference / administrative-data designs (@miller2021,
  @wyse2025) are built to overcome exactly this; @wyse2025's population-
  weighted, state-identified estimand faces the *same* concentration we
  study and *survives* its trend check via a triple difference — the
  contrast that makes the point.

### What this is NOT (the de-scope vs V2)

Drop from scope (leave in `paper/..._v2.Rmd` and the supplement, do not
promote):

- §2.5 "Trend-Robust Variants" (detrended + stratified as main-text
  estimators).
- §3.7 "Trend-Robust Estimates and Reconciliation" — the DLPS bridge
  ladder, per-bin remedy table, Rambachan–Roth sensitivity.
- The intro/3.6/discussion "trend-adjusted ≈ −10 / upper bound" rescue of a
  point estimate.

Rationale: those turn a *methods* paper into an ACA-identification paper.
The part-way draft **keeps the diagnosis** (placebo fails, Table 5
mechanism, triple-diff literature) but **refuses to manufacture a rescued
number or a remedy menu** — that honesty is cleaner and shorter, and it is
what "don't get lost in ACA" means. Detrended/stratified remain in the
*package* and the *supplement* as available tools; they are not the story.

Practical bonus: this removes `detrended_results.csv`, `binned_results.csv`
(main-text use), the DLPS ladder, and the RR chunk from the critical-path
build — shrinks `server_build_list_2026-07.md`.

### Section-by-section edits, starting from `weighted_sdid_paper.Rmd` (working tree)

Base = current working tree (already has the Table 5 mechanism, the placebo
section, and the uncommitted Miller/Wyse paragraph at intro L73). Edits:

- **Abstract (YAML).** Keep the three method contributions. Replace the
  application sentences: (i) divergence 0 vs −17.5; (ii) one D sentence:
  *"The population-weighted estimand fails an in-time placebo — a
  differential-trend problem intrinsic to the ACA county design that
  equally afflicts level-matched benchmarks (Borgschulte & Vogler 2020),
  not the weighting method — which triple-difference and administrative-data
  designs (Miller 2021; Wyse & Meyer 2025) are built to resolve."* No
  "−10", no remedies. (User finalizes wording.)
- **Intro P1–P3** (L65–69): keep. Optional one clause in P3 that the
  estimator "extends naturally to staggered adoption (developed in
  companion work)" — only if we want to advertise paper #2 here.
  Recommend: leave out for now.
- **Intro P4** (L71, application paragraph): keep through the in-space
  sentence; **cut the tail** ("An in-time placebo, however … upper bound …
  trend-adjusted effect nearer −10 …"). Replace with a 1–2 sentence D lead:
  the population-weighted target fails an in-time placebo; the failure is a
  treated-specific differential trend that level-matched convex estimators
  (incl. weighted DID and BV's DLPS) cannot remove, i.e. a property of this
  application; forward-ref §3.6 and the next paragraph.
- **Intro P5** (L73, the Miller/Wyse paragraph, currently uncommitted):
  **keep and sharpen** — this is the D payload. Make explicit: triple
  difference nets out the common trend and survives; our single-difference
  county design cannot; therefore *estimand choice and a formal placebo*,
  not eyeballed pre-trends, decide whether population weighting also moves
  the identification burden. Commit it.
- **Intro roadmap** (L75): unchanged (no §3.7).
- **§2 (estimator).** Setup / Weighted Modification / Asymptotics /
  Inference — unchanged. **Do not add the V2 §2.5 "Trend-Robust Variants".**
  *Do* add a short **§2.5 "Staggered adoption" capability note** (≈1 short
  paragraph): the weighted collapsed form applies cohort-by-cohort to
  not-yet-treated controls and aggregates to `τ^{ω̃,stag}`; formal statement
  and reductions in the supplement; no application in this paper. This is
  the only new §2 subsection, and it is a pointer, not a development.
- **§3.1–3.5.** Unchanged (Setting/Data, Main Results, Heterogeneity,
  Robustness, BV comparison). In §3.5 add one honest sentence: BV's DLPS
  identification rests on the same parallel-trend assumption the in-time
  placebo stresses, so their −11.36 and our −17.5 share the exposure — the
  comparison is level-to-level, not a validation.
- **§3.6 Placebo Tests.** Keep in-time (popw fails), mechanism paragraph,
  and **Table 5** (omega size bins — strong diagnostic, keep). **Rewrite the
  closing paragraph** (L607): delete "how much bias / trend-adjusted −10 /
  LIHP could rescue it". Replace with the D close: the failure localizes to
  a treated-specific differential trend that convex level-matching cannot
  span (SDID's λ lives on the simplex, cannot extrapolate); the same
  falsification would fail for weighted DID and for BV's DLPS; we therefore
  *do not* report a trend-corrected point estimate for this design and read
  the population-weighted number as the estimand SDID identifies *only if*
  parallel trends hold at that weighting — which here they do not. Point to
  triple-diff/admin designs. Keep in-space placebo + its two limitations.
- **§4 Monte Carlo.** Keep, **and add the trend cell — diagnosis rows
  only** (settled): show weighted SDID (and weighted DID) inherit bias when
  a weight-correlated trend exists and the in-time placebo detects it;
  equal-weight variants are unbiased. **Omit** the detrended/stratified
  *recovery* rows (that is the remedy story). This gives D a structural
  backstop ("it's the DGP, not the machinery") without reopening remedies.
  Script: add `gamma_trend` to the DGP in `run_mc_simulations.R` (one
  cluster, cheap vs the existing grid) → `mc_trend_results.csv` (already
  present from the V2 branch — reuse if the config matches).
- **§5 Discussion.** Rebuild P2 around A→B→D→C: weights select an estimand
  (A); we supply the estimator + inference (B); but the estimand can carry
  its own identification burden, checkable with an estimand-specific placebo
  (D); in ACA the population-weighted target fails that check as the
  level-matched literature does, so a credible population-weighted mortality
  effect needs a triple-difference/admin design (C→forward). Delete the
  "upper bound" and remedy caveats.

Net: keep exactly the current exhibits (T1–T5, F1–F4, MC) — **no new
main-text exhibit**. The change is framing + one rewritten paragraph +
cut tails, not new computation.

---

## Part 2 — Staggered weighted SDID (theory + code only)

**Scope:** ships as a package capability with supplement proofs, MC
coverage, and a validation golden test. **No empirical application / no
Jones paper.** The Jones `data_sje.dta` panel is used *only* as the
staggered golden-test fixture (real staggered data), not as a paper result.

**Finding (from exploration):** the fork is **block-only**. Every path
takes a scalar `T0` / single adoption (`synthdid.R:253`,
`analyze_application.R:173` `post = time>=2014`). Staggered has never run.

**Good news:** it is a *wrapper over the existing block estimator*, and
`R/stratified.R` is a near-exact structural template.

### Design

`synthdid_estimate_staggered(Y, adoption_time, treated.weights = NULL,
cluster = NULL, control = c("notyet","never"), ...)`

- `adoption_time`: length-N, finite for treated (adoption period), `Inf`
  for never-treated. (Alt long-format entry mirroring Stata `sdid`.)
- For each distinct finite cohort g:
  - controls = never-treated ∪ not-yet-treated through g's post-window
    (`control` toggle);
  - treated = units with `adoption_time == g`;
  - window = pre `< g`, post = g's event window;
  - build `Y_g`, call **block** `synthdid_estimate_weighted(Y_g, N0_g, T0_g,
    treated.weights = normalize(π restricted to g), cluster = ..., ...)`.
- Aggregate: `τ = Σ_g W_g τ_g / Σ_g W_g`, with `W_g = Σ_{i∈g} π_i`
  (importance-weight share) optionally × post-length (Clarke's
  treated-period weighting when `π = 1/N1`).
- `vcov`: cluster bootstrap resampling states → full refit (copy
  `stratified.R`'s `vcov` + `*_boot_rep`, swapping cohort for stratum).

Mapping to `stratified.R`: grouping var stratum→cohort; donor pool
same-stratum→not-yet-treated; **new piece** = per-cohort *time window* (the
one thing stratified does not do). ~200–300 lines + tests.

### Proofs / theory needed (supplement)

1. **Estimand + consistency.** Define `τ^{ω̃,stag} = Σ_g W_g τ_g^{ω̃}` as a
   convex combination; each cohort block is consistent under A1–A5 restricted
   to its window/donor pool, with **A5 estimand-specific per cohort**
   (`A5(ω̃_g)`) — the existing "identification is estimand-specific" remark,
   applied cohort-wise.
2. **Asymptotic normality / effective-N.** Aggregate CLT via Prop 3's
   cluster algebra with cohorts as blocks. The `√N₁ᵉᶠᶠ` condition now binds
   **per cohort**: a cohort with 1–2 treated units (Jones has several) or a
   single dominant state collapses the aggregate's effective N. State this
   explicitly — it is the honest limit for Jones.
3. **Not-yet-treated controls.** No-anticipation + parallel-trends for using
   not-yet-treated as controls (Callaway–Sant'Anna flavor) in the SDID
   factor model; a control used for cohort g must be untreated through g's
   post-window. **Required for Jones** (no never-treated states).
4. **Reductions.** (i) one cohort ⇒ block weighted SDID (exact);
   (ii) `ω̃ = 1/N₁` + treated-period aggregation ⇒ Clarke et al. staggered
   `sdid` — this is the golden-test bridge to Jones Table 5.
5. **Bootstrap validity** for the aggregate (states = clusters; a state is
   treated in one cohort, a control in others' windows — cluster resampling
   handles the dependence). `stratified.R`'s per-draw full refit already
   embodies the pattern.

### Tests

- **Golden test (validation only)**: reproduce Jones et al. Table 5
  (unweighted) on `data_sje.dta` via Stata `sdid` (have the do-file + data +
  Stata MCP); match our staggered estimator at `ω̃ = 1/N₁` to tolerance.
  Lives in `tests/`, not in the paper — the credibility anchor for the
  staggered code, mirroring the ETR golden test.
- Unit: single-cohort ≡ block (exact); uniform+Clarke-agg ≡ Stata sdid
  (tol); notyet/never toggle; weight renormalization sums to 1;
  degenerate-cohort skip (like `drop.infeasible`); bootstrap SE sane.
- MC: staggered DGP with weight-correlated trend → aggregate inherits bias,
  placebo detects, coverage of cluster bootstrap CI. Extend `test_weighted.R`
  / `test_remedies.R`.

**Effort:** wrapper + vcov ≈ bounded (stratified.R template); the real work
is per-cohort window/not-yet-treated construction, the aggregation-estimand
definition, and the Stata golden harness. Proofs are *extensions* of Props
1–3, not new theory. Scope ≈ a methods section for paper #2.

---

## Part 3 — Two-repo split

**Recommend: yes, split** — the package now carries weighted **and**
staggered capability that outlives this one paper, and a clean package is
worth publishing on its own. The "results factory / PDF press" decoupling
(CLAUDE.md) already means the Rmds read `results/*.csv` and **call no
package functions at render time**, so the paper repo needs the package
only at the *analysis* stage. Clean seam.

### Layout (two repos)

- **Repo A — `wsynthdid`** (new repo, clean carve): `R/ man/ tests/
  vignettes/ DESCRIPTION NAMESPACE`. The fork of `synthdid` + weighted +
  staggered. Distinct name so it installs side-by-side with upstream
  `synthdid` and via a clean `remotes::install_github`. Tagged releases,
  `R CMD check` clean, the staggered golden test lives here.
- **Repo B — the ACA paper** = the *current* repo minus package internals:
  `paper/ scripts/ results/ docs/ paper/data/`. Keeps results/ history,
  issues, and CLAUDE.md. Depends on `wsynthdid` pinned in `renv.lock` by git
  SHA.

### Mechanics

- Carve with `git subtree split` (or `git filter-repo`) on `R/ man/ tests/
  DESCRIPTION NAMESPACE` (+ pkg vignettes) → Repo A, history preserved.
  Remove them from the current repo; add Repo A as a dependency.
- Reproducibility: each paper's `renv.lock` records the exact package SHA
  used for its `results/` vintage (you already pin renv to R 4.4.2). This is
  strictly better than today — the estimator version behind each paper
  becomes explicit.
- Dev loop friction (edit package → reinstall to test in paper) mitigated
  with `renv` + `remotes::install_local()` / a dev install during active
  work.
- Split `CLAUDE.md`, `scripts/README.md`, renv between package and paper.

### Caveats

- The `upstream/*` fork remotes go with **`wsynthdid`** (it is the fork).
  The Borgschulte–Vogler reproduction data, `tmpclaude-*` gitignore, and
  `results/` stay with the **paper** repo.
- Dev-loop friction (edit package → reinstall to test) handled with
  `remotes::install_local()` during active work; `renv.lock` pins the SHA
  for the committed `results/` vintage.

---

## Execution order (all decisions settled)

1. **Paper edits** on `weighted_sdid_paper.Rmd` (Part 1): abstract, intro
   P4 tail cut + P5 commit/sharpen, §2.5 staggered capability note, §3.5
   BV sentence, §3.6 closing rewrite, §4 MC trend-cell diagnosis rows, §5
   discussion. No new main-text exhibit; reuses existing `results/`.
2. **`wsynthdid` staggered code** (Part 2): `synthdid_estimate_staggered()`
   + `vcov`, off the `stratified.R` template; supplement proofs; Stata
   golden test on `data_sje.dta`; MC trend cell.
3. **Repo split** (Part 3): carve `wsynthdid`, repoint the paper repo via
   renv. Do this *after* the paper edits land, so the ACA `results/` vintage
   is stable before the SHA pin.
