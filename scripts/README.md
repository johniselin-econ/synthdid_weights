# scripts/ — data pipeline and simulation runners

All scripts are run from the **repository root** (`Rscript scripts/<name>.R`).
They fall into two groups: the ACA application data pipeline (steps 1–6) and
the Monte Carlo runners (step 7).

## Data pipeline run order

| # | Script | Inputs | Outputs |
|---|--------|--------|---------|
| 1 | `build_seer_laus_covariates.R` | `paper/data/raw/seer_us.1990_2024.20ages.adjusted.txt.gz` (SEER popdata), `paper/data/raw/la.data.64.County` (BLS LAUS) | `paper/data/seer_laus_covariates.csv` |
| 2 | `rebuild_mortality_panel.R` | `paper/data/raw/cdc_wonder_data_2005-2011.csv`, `…_2012-2017.csv`, `…_states.csv` (CDC WONDER exports), output of step 1 | `paper/data/analysis_data.csv`, `paper/data/county_mortality_pre.csv` |
| 3 | `pull_county_votes.R` | `paper/data/raw/countypres_2000-2024.csv` (MIT Election Lab, manual download) or tonmcg GitHub fallback | `paper/data/county_votes.csv` |
| 4 | `pull_bv_covariates.R` | ACS/tigris/SAHIE via Census API (**needs `CENSUS_API_KEY` in `~/.Renviron`**), outputs of steps 2–3 | `paper/data/bv_covariates.csv` |
| 5 | `pull_urban_share.R` | 2010 Decennial SF1 via Census API | `paper/data/urban_share.csv` |
| 6 | `bv_replication_tables.R` | outputs of steps 2 and 4 | `paper/data/bv_table1_comparison.csv`, `paper/data/bv_table2_comparison.csv` |

Steps 3–5 are independent of one another; step 4 requires steps 2 and 3 first.
The paper and supplement (`paper/*.Rmd`) read only the committed CSVs in
`paper/data/`, so a knit never needs the raw inputs.

## Monte Carlo runners

| Script | Purpose | Output |
|--------|---------|--------|
| `run_mc_simulations.R` | Main MC sweep (18 configs × 100 sims, B = 50): bias/RMSE and coverage for the weighted estimator (main-paper Section 4). Parallelized (`SLURM_CPUS_PER_TASK`-aware); checkpoints in batches and **resumes** from the output file. | `paper/data/mc_results.csv` |
| `run_solutions_mc.R` | RMSE comparison of the three weighting strategies, Solutions 1–3 (supplement Appendix F). Same DGP and seeding as the main sweep. | `paper/data/mc_solutions_results.csv` |
| `slurm_run_mc.sh` | Slurm batch script running both sweeps in sequence (`sbatch scripts/slurm_run_mc.sh` from the repo root). | log in `scripts/slurm_mc_<jobid>.log` |

Seeding in both runners reproduces the equivalent sequential loop exactly
(`set.seed(20240101 + sim)` per task), so results are independent of worker
count.

## Superseded

`county_dataset_creation.R` is the original from-raw-files builder (Mac,
Dropbox paths); it is retained for provenance but superseded by steps 1–2,
which use only public re-downloadable inputs. See its header note.

## Data provenance

All inputs are public; no restricted microdata is used or committed.

- **CDC WONDER** (mortality): public county-aggregated exports, all-cause
  deaths ages 20–64. Cells with 1–9 deaths are suppressed in the public file
  and imputed by allocating the exact state-year residual (step 2). Per the
  WONDER data-use agreement, **death counts are dropped from the committed
  CSVs** — only rates and suppression flags are retained.
- **SEER** (county population by age/sex), **BLS LAUS** (county
  unemployment), **Census ACS / SAHIE / 2010 SF1**, and **MIT Election Lab**
  county presidential returns are all public; the large raw files are
  gitignored and re-downloadable (URLs in the script headers; the MIT file is
  guestbook-protected and must be fetched manually).
