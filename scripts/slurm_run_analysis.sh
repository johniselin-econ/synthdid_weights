#!/bin/bash
#SBATCH --job-name=synthdid_analysis
#SBATCH --partition=day
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=scripts/slurm_analysis_%j.log

# Run the ACA application analysis (estimates, bootstrap/cluster SEs, event
# studies, robustness, placebos) and write the reduced result tables to
# results/. The covariate-adjusted bootstrap SEs dominate runtime (~130
# CPU-hours); with 32 cpus this is roughly 4-6h. analyze_application.R reads
# SLURM_CPUS_PER_TASK for its worker count and checkpoints each (spec, method)
# SE to results/_ses_partial.csv, so a re-submit resumes where it left off.
#
# Submit from the repository root:  sbatch scripts/slurm_run_analysis.sh
#
# Prerequisite (once): the renv environment must be restored on the cluster so
# the worker R sessions find dplyr/tidyr/readr/tibble/purrr. From the repo root:
#   module load R/4.4.2-gfbf-2024a && Rscript scripts/setup.R
# (synthdid itself is sourced from R/, so it does not need to be installed.)

set -euo pipefail
module load R/4.4.2-gfbf-2024a

echo "[$(date)] application analysis starting on ${SLURM_CPUS_PER_TASK:-?} cpus"
Rscript scripts/analyze_application.R
echo "[$(date)] application analysis complete; results/ written"
