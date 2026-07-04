# TODO — remaining items (consolidated 2026-06-07)

Single live todo list. Consolidates and supersedes `docs/todo.txt`,
`docs/bv_replication_todo.md`, and `paper/DLPS_REPLICATION_PLAN.md` (all
deleted; full history in git). Done/partial status of every April
coarse-review item is tracked in `docs/submission_review_checklist.md`; only
genuinely open work is listed here.

The BV2020/DLPS replication is **complete** — the record lives in supplement
Appendix E.5 and `scripts/README.md` documents the pipeline. No open BV items.

## Monte Carlo — COMPLETE (2026-06-08)

- [x] Main sweep finished (1800 sims) + `run_solutions_mc.R` (3600 sims);
      both CSVs committed (b0eebe1)
- [x] Three hardcoded MC claims verified against final `mc_results.csv`;
      jackknife coverage claim corrected (b0eebe1)
- [x] Supplement F.3 Solutions table renders; PDFs re-rendered from server
      results (f45a416)
- [x] Final clean knit of BOTH documents (2026-06-10: caches cleared,
      paper + supplement rendered clean). Re-knit again if further text
      edits land (clear `paper/*_cache/` first — caches do not watch data
      files)

## Math/proof fixes (econometric review 2026-06-10) — DONE 2026-06-10

See `docs/econometric_review_2026-06-10.md` for full statements (NOTE
2026-06-10: this file is not in the repo — commit it from the machine it
was written on, or drop this pointer); all items fixed in both Rmds and
re-rendered:

- [x] **M1**: A4 replaced with internally consistent log N0 = o(T0),
      log T0 = o(N0) + pointer to AHIW's primitive rates
- [x] **M2**: A5 restated as oracle-weights assumption at the
      factor-loading level (convexity + intercepts + λ* oracle +
      diffuseness ‖ω*‖₂ → 0)
- [x] **M3**: A2 now states cross-unit independence; Props 1–2 restated
      on post-minus-pre contrasts ξ_j (serial covariance absorbed in σ²_ξ)
- [x] **M4**: Prop 2 unified at √N1eff rate, limit N(0, σ²_ξ) +
      negligibility condition (R) with two remarks (primitive conditions;
      when (R) is plausible — T0=5 makes it doubtful, motivating bootstrap)
- [x] **M5**: Prop 3 restated with cluster-aggregated contrasts ξ^(k),
      target σ²_cl ≠ σ²_ξ; bounded-cluster-sizes dropped (weight-share
      condition instead)
- [x] **m6**: B.5 + Prop 2 remark numbers now inline from
      `load-scalars-early` chunk (app_scalars.csv + heterogeneity.csv)
- [x] m7 ("W2 fails" framing unified to "loses practical bite"; CA-23%
      threshold sentence rewritten), m8 (ω/λ programs displayed in §2.1),
      m9 (assumption gloss aligned, no longer attributed to AHIW labels),
      m10 (cross term bounded via rank-r expansion, fixed-effect parts
      vanish exactly), n11 (B−1), n12 (ℓ index), n13 (estimand-drift
      remark in C.1), n14 (W3 labeled in main text), n15 (period-weights
      T1→∞ caveat in G)

## Paper revisions (April coarse review) — DONE 2026-06-10

- [x] SUTVA / spillovers paragraph in §3.1 (cites Frean–Gruber–Sommers
      2017 + Sommers 2016; attenuation + within-treated-group arguments).
      Border-county robustness check DISCARDED (paragraph suffices)
- [x] Two-unit worked example in §2.2 (covariance identity + 0.9/0.1
      numeric example)
- [x] Paired-bootstrap divergence test framed as recommended ex-ante
      diagnostic (new ¶ in §2.4; §3.2 references it)
- [x] Notation pass: k-index clash in Prop 3/C.1 fixed (ℓ over treated
      units); intro hat/tilde already consistent
- [x] Late expanders named in §3.1 (AK/IN/LA/MT/NH/PA); 45 clusters =
      44 states + DC explained; ALSO FIXED: the old "excluding states with
      substantial pre-2014 expansions" sentence misdescribed the sample
      (early expanders MA/NY/VT/DE/DC are *treated*)
- [x] Abstract cut to ~140 words, finding first
- [x] Small-county framing unified ("oscillate around zero"; rural-quintile
      positive flagged as noisy; Discussion phrasing aligned)
- [x] Data availability section now holds the data text; duplicate AI
      disclosure (Disclosure section) merged into the Elsevier AI
      declaration

Still open from the April list (deliberate): permute-placebo MC comparison
and cluster-robust placebo SE (discarded / deferred as documented in the
paper). The ACA-calibrated MC cell is now implemented (see below) and
awaits a server run.

## ACA-calibrated MC cell + clustered event-study SEs (2026-06-10)

- [x] Config 19 added to `run_mc_simulations.R`: ACA-calibrated cell with
      the application's panel dimensions (N0=1643, N1=1181, T0=5, T1=4)
      and treated-unit sizes set to the application's ACTUAL 2013 county
      populations (`results/heterogeneity.csv`), so every sim has exactly
      HHI=0.009829 / N1eff~102. (The old TODO numbers "N1=118, HHI=0.0085"
      were the pre-rebuild N1eff and its reciprocal — stale.) gamma=1.
      Resume-safe: configs 1–18 keep their numbers; a completed grid run
      computes only the new cell. Verified: one full-scale sim runs clean
      (~1.4 min at B=4, est_wt 11.65 vs target 11.61; placebo SE 0.005 vs
      bootstrap 2.8 — the uniform-placebo understatement is dramatic at
      this concentration).
- [x] **Server run DONE** (slurm job 14602357, completed 2026-06-10
      15:16, 12 min — resume logic skipped configs 1–18):
      `paper/data/mc_results.csv` now has all 1900 rows (19 × 100) and is
      committed. Re-knit on the Windows machine DONE (2026-06-10:
      caches cleared, both docs rendered clean): the `mc_has_aca`
      conditional now fires — §4 bias/coverage tables carry the
      N1=1,181 row (pop-wt bias 0.005, RMSE 0.053; coverage boot 1.00,
      jk 1.00, placebo 0.74) and the ACA-calibrated paragraph renders
      (target ATT 11.6; the 0.74 placebo coverage is the uniform-placebo
      undercoverage that motivates the cluster bootstrap).
- [x] Cluster support in `synthdid_event_study()` /
      `.event_study_replications()`: new `cluster` argument defaulting to
      the estimate's stored cluster (vcov-consistent); bootstrap resamples
      whole clusters (mirrors `cluster_bootstrap_se_weighted`); placebo
      warns and falls back to unit-level. Tests + man page updated.
- [x] `analyze_application.R`: event studies (main + SC) moved to their own
      idempotency guards (like placebo_intime) and now use the
      state-clustered bootstrap; figure captions updated to say
      "state-clustered". Curves are unchanged (deterministic); only the
      CI bands change.

## 2005 pre-period extension (branch `extend-preperiod-2005`, 2026-07-04)

Full-quality run committed as `results_2005/` (seed 20240101, B=500,
fast_mode=FALSE, SHA ea9b090; slurm job 17073601). Same 2,824-county
sample; window 2005–2017 (T0=9 vs 5).

Event-study comparison vs baseline `results/` (2026-07-04):

- DID curves shift by a **constant** across all overlapping event times
  (eq +1.80, pop-wt −7.57) — purely the longer pre-period baseline mean
  (matches the obs_diff change in app_scalars). SDID shifts are not
  constant: ω/λ genuinely re-estimated.
- New 2005–2008 leads reveal a large positive pre-trend in the
  pop-weighted comparison (DID leads +10 to +12 at −9..−7). Weighted
  SDID absorbs much of it (+6.4 to +7.2 far out, ~0 by −5/−4), but
  near-treatment leads deteriorate: −1 lead goes −2.7 → −6.3 (SE 1.6),
  −3/−2 flip from ~0 to −1.7/−1.9. Eq-SDID leads stay small and noisy
  all the way back (−3.9 … +0.4).
- Post effects deepen ~5 pts uniformly (sdid_w headline −17.45 → −22.45;
  post path −16/−16/−14/−23 → −22/−21/−19/−29), SEs up ~15–20%. Eq specs
  remain near zero / insignificant in both vintages, so the eq-vs-weighted
  contrast *widens* with the longer window.

Open items:

- [ ] Decide whether the 2005 window goes in as a robustness exhibit
      (event-study figure + estimates row) in paper or supplement
- [ ] If it goes in: the weighted-SDID −1 lead of −6.3 (~4 SEs, ~29% of
      the on-impact effect vs 17% in the 2009 window) is the obvious
      referee target — add a sentence connecting the 2013 dip to
      anticipation/woodwork and the existing placebo-bias discussion

## Submission logistics (Elsevier)

- [ ] `devtools::check()` clean pass
- [ ] Competing-interests declaration (Word doc); Erica's Amazon affiliation
      / employer disclaimer
- [ ] Deposit data in a repository per Elsevier Option C
- [x] Document the package-loading setup in a replication README (done
      2026-06-10: `scripts/README.md` now documents the analysis stage
      (`analyze_application.R` sources `R/` directly), the orchestrator,
      and which scripts need the installed package; the Rmds themselves
      no longer source `R/` and call no synthdid functions)
- [ ] Submission fee + Editorial Manager upload

## Optional / deferred

- [ ] Cluster-robust placebo SE (currently warns and falls back to
      unit-level; applies to vcov and to event-study placebo SEs)
- [x] Cluster support in `.event_study_replications()` for clustered
      event-study SEs (done 2026-06-10, see above)
- [ ] "permute" placebo pseudo-weight variant: `vcov.R` implements it but
      the MC only simulates "uniform" — decide whether to add the comparison
      or defer
- [x] ACA-calibrated MC cell (done 2026-06-10, see above; server run
      pending)
- [x] CLAUDE.md for the repo (added 2026-06-10)
