# Pre-submission review checklist (regenerated 2026-06-06)

Supersedes the earlier version. Includes status of every item from the April
coarse review (docs/coarse_review_2026-04-17.md). Legend: [x] done,
[~] partial / verify, [ ] open.

## Fixed today (commit pending)

- [x] "Appendix F" -> "Supplementary Appendix E.4" (4 spots in main paper)
- [x] Economics Letters ghost reference + "for this draft" language removed
      from Section D; MC machinery is now REAL: scripts/run_mc_simulations.R
      (parallelized, seeding identical to sequential) writes
      paper/data/mc_results.csv; Tables D.1-D.2 render when it exists
- [x] highlights.txt + README updated (no sign-flip claim; 17 not 13.5)
- [x] secnumdepth=3 (kills "D.2.0.1"-style heading artifacts)
- [x] Sign-change/flip language sweep -- clean in both Rmds
- [x] deaths counts dropped from public CSVs (WONDER DUA); rates + death_supp
      flag remain; no downstream code used the counts
- [x] AI-use Disclosure section added to paper + README

## Open blockers

- [ ] **MC results**: the parallel run must finish and the supplement re-knit.
      First smoke test took 2.16h for 18 sims (vs ~105s/sim when profiled
      solo) -- timing under investigation; full run feasibility TBD. After
      results land: **verify the three claims in the Section D bullet list
      against the actual tables** (unbiasedness; boot/placebo coverage with
      placebo undercoverage at high HHI; conservative jackknife).
- [ ] **Abstract ~370 words** -> cut to 100-150, finding first (your edit).
- [ ] **Paper PDF locked by your viewer** -- close it; rebuild is queued.
- [ ] Small-county framing tension: pop bins "oscillate around zero" vs
      most-rural urban quintile "+21" vs "near zero or positive" (Discussion).
      Unify (your call on framing).

## April coarse review: major concerns

- [x] **Formal identification/asymptotics absent** -> Props 1-3, weight
      regularity (W2), proofs in Supp B; "when the asymptotics fail" in main
      text. *Re-read B.5 with fresh eyes.*
- [~] **Monte Carlo evidence absent** -> infrastructure now real (above), but
      the reviewer asked for a main-text summary (small coverage table or a
      results paragraph in Sec 2/Discussion). Decide after results land.
- [~] **Placebo variance inconsistency** (uniform pseudo-weights vs weighted
      estimand) -> MC reports placebo coverage; paper recommends cluster
      bootstrap as primary. Verify Supp C.3 states a clear recommendation
      and that MC results back the "undercoverage at high HHI" claim.
- [x] **Effective-N collapse / "modestly larger" SEs** -> N1eff (Kish) in
      Table 2, dominant-weight discussion, "modestly" wording gone, cluster
      bootstrap rationale explicit.
- [x] **Single application, no weight robustness** -> Table 2 (2013/2010/
      sqrt/log/equal) + leave-one-out (E.3).
- [ ] **SUTVA / cross-border spillovers** -> still no discussion. Add a
      paragraph (cite Medicaid spillover lit); optional border-county check.
- [~] **Staggered adoption / late expanders** -> AK/IN/LA/MT/NH/PA are
      excluded (now correctly after the state-fips fix). Verify the data
      section names them; consider a sensitivity footnote (e.g., late
      expanders as never-treated controls through 2015).
- [ ] **Worked two-unit parametric example** -> not in draft. Cheap to add
      (half-page algebra box in Sec 2 or supplement).
- [~] **Practitioner decision rule** -> the paired-bootstrap divergence test
      (p < 0.001) IS a Hausman-style diagnostic -- frame it explicitly as
      the recommended ex-ante check ("estimate both, test the gap").
- [x] **Relation to heterogeneity-robust DiD** (Callaway-Sant'Anna etc.) ->
      cited and positioned in intro + supplement. Verify depth satisfies you.
- [x] **Decompose the divergence** -> corrected Figs 3-4 binscatters, LOO,
      weight-compression gradient; pop-weighted mean of unit effects
      recovers tau-w (internal-consistency check in text).

## April coarse review: detailed comments 1-11

1. [~] Asymmetric SE methods in Table 1 (jackknife for equal, bootstrap for
   weighted) -- still asymmetric; add one rationale sentence or harmonize.
2. [x] Cluster-robust inference now core (Prop 3), not "future extension."
3. [~] DiD pre-trend failure vs robustness claim -- now framed as motivation
   for SDID's time weights (Fig 2 text); confirm you're happy with framing.
4. [ ] Hat/tilde notation consistency intro vs formal defs -- needs your read.
5. [ ] tilde-omega used for researcher weights AND control weights -- check
   Sec 2 notation table.
6. [x] Sign-change claim -- moot; no sign-change claim remains (-0.71 vs
   -17.45, framed as divergence).
7. [x] "Modestly larger" SE wording removed (ratio is now ~2x: 2.32->4.58).
8. [x] "Order of magnitude" overstatement removed from SE/estimate claims
   (one remaining literal use re county populations, 185k vs 3.8k, is
   factually correct).
9. [x] Sample construction replicability -- public-data pipeline, E.4, data
   statement.
10. [~] Clusters now 45 -- check the excluded-state list appears in the data
    section.
11. [x] Bootstrap replications 200 -> 500.

## Hygiene (unchanged from previous list)

- [ ] Setup chunks `source("../R/*.R")` -- document in replication README or
      switch to library(synthdid)
- [ ] Erica's Amazon affiliation + any employer disclaimer
- [ ] lit-counts figure caption length (move method detail to a footnote?)
- [ ] Final knit of BOTH documents after MC results + any text edits
      (remember: caches don't watch data; clear or use dependson)
