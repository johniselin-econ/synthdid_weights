# TODO — remaining items (consolidated 2026-06-07)

Single live todo list. Consolidates and supersedes `docs/todo.txt`,
`docs/bv_replication_todo.md`, and `paper/DLPS_REPLICATION_PLAN.md` (all
deleted; full history in git). Done/partial status of every April
coarse-review item is tracked in `docs/submission_review_checklist.md`; only
genuinely open work is listed here.

The BV2020/DLPS replication is **complete** — the record lives in supplement
Appendix E.5 and `scripts/README.md` documents the pipeline. No open BV items.

## In flight: Monte Carlo (Slurm job 13958206, submitted 2026-06-07)

- [ ] Main sweep finishes (resumed at 270/1800 tasks, 32 workers), then
      `run_solutions_mc.R` runs in the same job →
      `paper/data/mc_solutions_results.csv`
- [ ] **Verify the three hardcoded claims** in the main-text MC Results
      bullets against the final `mc_results.csv`: (i) approximate
      unbiasedness for each target; (ii) bootstrap/placebo near-nominal
      coverage with placebo undercoverage at high HHI; (iii) jackknife
      conservative at small N1. These were written against the partial
      270-sim checkpoint.
- [ ] Confirm the supplement F.3 Solutions 1–3 RMSE table renders (chunk is
      conditional on the CSV existing) and the generated comparison text
      reads sensibly
- [ ] Clear all knitr caches and do a clean knit of BOTH documents
      (caches do not watch data files)

## Paper revisions still open (April coarse review)

- [ ] SUTVA / cross-border spillovers: add a paragraph (cite Medicaid
      spillover literature); optional border-county robustness check
- [ ] Worked two-treated-unit algebra example showing when weighting moves
      the estimate (half-page box, Section 2 or supplement)
- [ ] Frame the paired-bootstrap divergence test explicitly as the
      recommended ex-ante practitioner diagnostic ("estimate both, test the
      gap")
- [ ] Notation pass: hat/tilde consistency intro vs. Section 2;
      $\tilde\omega$ currently used for both researcher weights and control
      weights
- [ ] Verify the data section names the excluded late expanders
      (AK/IN/LA/MT/NH/PA) and the excluded clusters behind the state count
- [ ] Abstract: cut to 100–150 words, finding first
- [ ] Small-county framing tension: pop bins "oscillate around zero" vs.
      most-rural urban quintile "+21" vs. "near zero or positive"
      (Discussion) — unify
- [ ] Empty `\section*{Data availability}` in the paper: the data text
      currently sits under `\section*{Disclosure}` — split or merge

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
- [ ] CLAUDE.md for the repo (init skill)
