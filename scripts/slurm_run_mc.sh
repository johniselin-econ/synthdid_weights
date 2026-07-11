#!/bin/bash
#SBATCH --job-name=synthdid_mc
#SBATCH --partition=day
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=23:00:00
#SBATCH --output=scripts/slurm_mc_%j.log

# Complete the main Monte Carlo sweep (resumes from paper/data/mc_results.csv),
# then the Solutions 1-3 RMSE comparison for supplement Section F, then the
# trend cell (rework plan Phase 3). Completed sweeps resume-skip instantly.
# Submit from the repository root: sbatch scripts/slurm_run_mc.sh

set -euo pipefail
module load R/4.4.2-gfbf-2024a

echo "[$(date)] main MC sweep starting"
Rscript scripts/run_mc_simulations.R
echo "[$(date)] main MC sweep done; solutions MC starting"
Rscript scripts/run_solutions_mc.R
echo "[$(date)] solutions MC done; trend MC starting"
Rscript scripts/run_mc_trend.R
echo "[$(date)] all sweeps complete"
