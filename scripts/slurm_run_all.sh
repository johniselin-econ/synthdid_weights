#!/bin/bash
#SBATCH --job-name=synthdid_all
#SBATCH --partition=day
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=23:00:00
#SBATCH --output=scripts/slurm_all_%j.log

# One-shot cluster rebuild of the ENTIRE results factory, in dependency order:
#   1. canonical application analysis          -> results/*.csv
#   2. extended 2005 pre-period run            -> results_2005/*.csv
#   3. Monte Carlo sweeps (main+solutions+trend)-> paper/data/mc_*.csv
#   4. commit the results + manifests and (best-effort) push
#
# The Rmd render is deliberately NOT here: there is no pandoc/LaTeX on the
# cluster and rendering is a seconds-long LOCAL step. After this job lands:
# pull locally, delete paper/*_cache/, then Rscript paper/_build_all.R.
#
# Submit from the repository root:  sbatch scripts/slurm_run_all.sh
#
# Env knobs (all optional):
#   RUN_ALL_FORCE=0   skip the ANALYZE_FORCE recompute and only fill missing
#                     CSVs. Default 1: recompute BOTH results dirs into ONE
#                     coherent vintage (the whole point of a single-shot run).
#                     Affects the two analysis stages only -- the MC sweeps use
#                     their own resume logic and are never force-rerun here.
#   RUN_ALL_COMMIT=0  do not git add/commit/push the outputs. Default 1.
#
# Prerequisite (once): renv restored on the cluster so the PSOCK workers find
# the CRAN packages:  module load R/4.4.2-gfbf-2024a && Rscript scripts/setup.R
# (synthdid itself is sourced from R/, so it is never installed.)

set -euo pipefail
module load R/4.4.2-gfbf-2024a

FORCE="${RUN_ALL_FORCE:-1}"
if [ "$FORCE" = "1" ]; then export ANALYZE_FORCE=1; fi
echo "[$(date)] run_all starting on ${SLURM_CPUS_PER_TASK:-?} cpus (force=$FORCE)"

# --- 1. canonical application analysis (results/) --------------------------
echo "[$(date)] (1/4) application analysis -> results/"
Rscript scripts/analyze_application.R

# --- 2. extended 2005 pre-period run (results_2005/) -----------------------
echo "[$(date)] (2/4) extended 2005 analysis -> results_2005/"
ANALYZE_PANEL=paper/data/analysis_data_2005.csv \
ANALYZE_RESULTS=results_2005 \
ANALYZE_INTIME_YEARS=2009,2010,2011,2012 \
Rscript scripts/analyze_application.R

# --- 3. Monte Carlo sweeps (resume-skip if already committed) --------------
echo "[$(date)] (3/4) Monte Carlo sweeps -> paper/data/mc_*.csv"
Rscript scripts/run_mc_simulations.R
Rscript scripts/run_solutions_mc.R
Rscript scripts/run_mc_trend.R

# --- 4. commit + push (best-effort) ----------------------------------------
if [ "${RUN_ALL_COMMIT:-1}" = "1" ]; then
  echo "[$(date)] (4/4) committing results"
  git add results results_2005 \
          paper/data/mc_results.csv paper/data/mc_solutions_results.csv \
          paper/data/mc_trend_results.csv
  if git diff --cached --quiet; then
    echo "[info] no result changes to commit"
  else
    git commit -m "results: single-vintage cluster rebuild (analysis + 2005 + MC)"
    git push || echo "[warn] push failed (compute node egress?); push from a login node"
  fi
else
  echo "[$(date)] (4/4) commit skipped (RUN_ALL_COMMIT=0)"
fi

echo "[$(date)] run_all complete"
