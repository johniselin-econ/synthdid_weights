# BV2020 thorough-replication TODO

Status as of 2026-06-04. Goal: replicate the all-cause-mortality parts of
Borgschulte & Vogler (2020) Tables 1 and 2 (paper: `docs/Borgschulte and
Vogler (2020).pdf`; extracted text: `docs/bv2020_text.txt`).

Pipeline: `scripts/county_dataset_creation.R` (raw data, Mac) →
`paper/data/analysis_data.csv` → `scripts/pull_bv_covariates.R` +
`scripts/pull_county_votes.R` → `scripts/bv_replication_tables.R` →
`paper/data/bv_table{1,2}_comparison.csv`.

---

## 1. Data to download (you)

- [ ] **CDC WONDER county re-export** (web UI; API blocks county geography).
      All-cause deaths, ages 20–64, by county × year, **2005–2017** (extends
      the current 2009–2017 pull so the PS model can use 2005–09 mortality
      year-by-year and hold out 2010–13, exactly as BV do). Make sure
      suppressed cells appear as "Suppressed" rows. Save over
      `<dropbox>/sdid_weights/data/cdc_wonder_data.csv`.
- [ ] **CDC WONDER state-level export**, same query at state × year →
      `<dropbox>/sdid_weights/data/cdc_wonder_states.csv`. Needed for the
      "residual" suppression imputation (state totals are exact).
- [ ] *(Phase 2, optional)* WONDER exports **by sex** and **by the five age
      groups** to replicate Table 2 columns (3)–(9).
- [ ] **MIT Election Lab** `countypres_2000-2024.csv` — guestbook-protected,
      manual download from
      https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/VOQCHQ
      → save as `paper/data/raw/countypres_2000-2024.csv`.
      (Interim: `pull_county_votes.R` auto-falls-back to the tonmcg GitHub
      mirror — fine for now, switch to MIT before publication.)
- [ ] **Census API key on this (Windows) machine**: free at
      https://api.census.gov/data/key_signup.html → add
      `CENSUS_API_KEY=...` to `~/.Renviron`. SAHIE now requires it, and the
      uninsured variable must be re-pulled (see bug 3 below).
- [ ] **Census of Governments / Annual Survey of State Government Finances**:
      state health + welfare expenditures 2005–2013 (BV use the logged
      average). Not yet scripted — last missing PS candidate.
- [ ] Get the raw Dropbox folder (`sdid_weights/data/`) reachable from this
      machine, or plan to run `county_dataset_creation.R` on the Mac.

## 2. Bugs found and fixed in code — review these (data not yet regenerated)

1. **`log_20_64` was log(pop 20–24)** — `across(pop_20_24:pop_20_24)` is a
   one-column range (`county_dataset_creation.R`, fixed). It's one of the six
   DID controls, so Appendix F numbers in the current supplement were
   estimated with a wrong control.
2. **`pct_55_64` used total population as denominator** — BV's five age
   shares sum to 1 within the 20–64 population (their Table 1), so the
   denominator is working-age pop (fixed in the same script).
3. **SAHIE `AGECAT=3` is ages 50–64, not 18–64** (verified against the API
   metadata). The `uninsured_rate` in `bv_covariates.csv` is currently the
   50–64 uninsured rate, understating BV's non-elderly rate by ~7pp
   (ours 0.12/0.15 vs BV 0.19/0.24). `pull_bv_covariates.R` now uses
   `AGECAT=1` — **needs a re-run with a Census key** to take effect.
4. (Data source, not ours) tonmcg vote file: Laclede County MO `total_2008`
   is a transcription error; script computes denominators from dem+gop+oth.

**Until `analysis_data.csv` and `bv_covariates.csv` are regenerated, every
DLPS number (supplement Appendix F and `bv_replication_tables.R` output) is
provisional.** The main-paper SDID results are NOT affected: they use only
`crude_rate` and `population`, which are untouched by bugs 1–3. The
suppression rework (below) WILL change the SDID sample, though.

## 3. Suppression imputation — review the design

- WONDER distinguishes "Suppressed" (1–9 deaths) from true zeros, so yes, we
  can tell them apart. The old pipeline dropped any county with a suppressed
  year (~380 counties, small/rural/high-mortality, mostly non-expansion) —
  the likely driver of our control-mortality gap vs BV (323 vs 359).
- New code in `county_dataset_creation.R` (untested until raw data arrives):
  `impute_method = "residual"` allocates the exact state-year residual
  (state total − unsuppressed county sum) across suppressed cells by
  population, clamped to [1,9]; falls back to truncated-Poisson; `"lower"` /
  `"upper"` give deaths=1/deaths=9 bound runs for a robustness footnote.
- [ ] Decide: residual allocation as headline + bounds in appendix?
- [ ] After rebuild: check county count vs BV's 2,823 pre-trim sample.

## 4. Empirical-approach alignment (your Q3) — remaining items

- [x] Logit PS on union of double-lasso selections; trim [0.038, 0.971] — matches.
- [x] DID: county+year FE, state clusters, weights = pop 20–64 × IPW — matches.
- [x] Square-root lasso (BV fn. 8) — implemented via `RPtests::sqrt_lasso`
      in `bv_replication_tables.R` (falls back to hdm if absent).
- [ ] Outcome-lasso DV: BV use 2005–09 mortality (2010–13 held out); we use
      2009–13 average until the WONDER re-export lands. **Holdout violated
      until then.**
- [ ] PS candidate set still missing: state H&W expenditures (item 1),
      log total population (have components, trivial to add).
- [ ] BV's panel-lasso (Belloni et al. 2016) re-selection of DID controls —
      we adopt their published 6 controls instead; fine, but note it.

## 5. Open questions from the first comparison run (investigate after rebuild)

Current output (`Rscript scripts/bv_replication_tables.R`):

| | ours | BV |
|---|---|---|
| Table 2 col (2) | **+4.05 (4.98)** | −11.36 (3.59) |
| Final sample | 1,730 counties | 2,260 |
| Trim | 696 (28 exp / 668 ctrl) | 563 (78 / 485) |

- [ ] **Sign flip + trim blowout**: our PS separates far more than BV's —
      Harris, Miami-Dade, Dallas, Tarrant, Broward all get p̂ ≈ 0 and are
      trimmed (36.4% of control population!). Expected to improve with the
      corrected uninsured rate, MIT vote shares, expenditures, and the
      suppression-expanded sample. If mega-counties still trim out after the
      rebuild, revisit sqrt-lasso penalty level and covariate scaling.
- [ ] **Obama-share levels**: our pop-weighted means (0.58 exp / 0.53 ctrl)
      sit well above BV's Table 1 (0.46 / 0.38) for both groups; same
      direction, different level. Check BV Appendix Table A.2 / variable
      definition — possibly share of voting-eligible population or a
      different weighting. Worth resolving before trusting the PS.
- [ ] **% Hispanic flipped** (ours exp 0.12 / ctrl 0.22 vs BV 0.19 / 0.16) —
      mostly a sample-composition artifact of the trim; recheck after rebuild.
- [ ] Compare selected covariates against BV's described selection (age
      45–54/55–64 shares, male, Hispanic, black, income, unemployment,
      uninsured, political vars; Appendix Table A.2 has the full list).

## 6. Run order once data is in place

1. `county_dataset_creation.R` (Mac or wherever raw data lives) — rebuilds
   `analysis_data.csv` with imputation + covariate fixes
2. `pull_county_votes.R` (uses MIT file if present)
3. `pull_bv_covariates.R` (needs CENSUS_API_KEY; picks up SAHIE fix +
   2005–09 mortality once the panel extends back)
4. `bv_replication_tables.R` → compare against published values
5. Sensitivity: rerun 1+4 with `impute_method = "lower"` / `"upper"`
6. Update supplement Appendix F to use the corrected controls and (if the
   match is close) the new Table 1/Table 2 comparison tables
