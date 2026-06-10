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

See `docs/econometric_review_2026-06-10.md` for full statements; all items
fixed in both Rmds and re-rendered:

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

Still open from the April list (deliberate): ACA-calibrated MC cell
(recommended tackle — one config addition to run_mc_simulations.R);
permute-placebo MC comparison and cluster-robust placebo SE (discarded /
deferred as documented in the paper)

## Submission logistics (Elsevier)

- [ ] `devtools::check()` clean pass
- [ ] Competing-interests declaration (Word doc); Erica's Amazon affiliation
      / employer disclaimer
- [ ] Deposit data in a repository per Elsevier Option C
- [ ] Document the `source("../R/*.R")` setup chunks in a replication README
      (or switch the Rmds to `library(synthdid)`)
- [ ] Submission fee + Editorial Manager upload

## Optional / deferred

- [ ] Cluster-robust placebo SE (currently warns and falls back to
      unit-level)
- [ ] Cluster support in `.event_study_replications()` for clustered
      event-study SEs
- [ ] "permute" placebo pseudo-weight variant: `vcov.R` implements it but
      the MC only simulates "uniform" — decide whether to add the comparison
      or defer
- [ ] ACA-calibrated MC cell (N1 = 118, HHI ≈ 0.0085) so Table 1's SEs have
      direct simulation backing — decide whether to add
- [x] CLAUDE.md for the repo (added 2026-06-10)
