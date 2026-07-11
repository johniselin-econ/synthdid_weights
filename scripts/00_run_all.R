#!/usr/bin/env Rscript
# =============================================================================
# 00_run_all.R — top-level pipeline orchestrator
# =============================================================================
# Runs the full project in dependency order and, by default, SKIPS any stage
# whose outputs already exist. Because the data and analysis artifacts are
# committed, a fresh clone goes straight to `render` and produces the PDFs in
# seconds, without re-running the (multi-hour) analysis or Monte Carlo.
#
# Stages (in order):
#   data      build cleaned inputs in paper/data/ from raw + public APIs
#             (build_seer_laus -> rebuild_mortality -> pull_votes ->
#              pull_bv_covariates -> pull_urban_share -> bv_replication_tables;
#              note: pull_bv_covariates depends on the cleaned mortality panel,
#              so "download" and "clean" cannot be cleanly separated)
#   analysis  run all expensive ACA estimators once -> results/*.csv
#   mc        run the Monte Carlo sweeps          -> paper/data/mc_*.csv
#   render    knit the paper and supplement       -> paper/*.pdf
#
# Usage (from anywhere; the script cd's to the repo root):
#   Rscript scripts/00_run_all.R                  # run missing stages, then render
#   Rscript scripts/00_run_all.R render           # one stage only
#   Rscript scripts/00_run_all.R data analysis    # a subset of stages
#   Rscript scripts/00_run_all.R --from analysis  # this stage and everything after
#   Rscript scripts/00_run_all.R --force          # rerun every stage
#
# Notes:
#   * `data` and `mc` need their raw/API inputs to (re)run; with the committed
#     CSVs present they are skipped, so no API keys or downloads are needed to
#     render. Forcing `data` requires CENSUS_API_KEY (see scripts/README.md).
#   * `render` always runs when it is reached (it is the deliverable), unlike
#     the upstream stages which are presence-gated.
# =============================================================================

# ---- locate the repo root (parent of this script's directory) --------------
args0    <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args0, value = TRUE))
if (length(file_arg) == 1) {
  repo_root <- normalizePath(file.path(dirname(file_arg), ".."))
} else {
  repo_root <- normalizePath(getwd())  # assume invoked from repo root
}
setwd(repo_root)
stopifnot(dir.exists("scripts"), dir.exists("paper"))

# ---- stage definitions ------------------------------------------------------
# `scripts`: run in order, each as its own Rscript subprocess (full isolation,
#            exactly as if invoked by hand).
# `outputs`: files whose presence means the stage can be skipped.
STAGES <- list(
  data = list(
    scripts = file.path("scripts", c(
      "build_seer_laus_covariates.R",
      "rebuild_mortality_panel.R",
      "pull_county_votes.R",
      "pull_bv_covariates.R",
      "pull_urban_share.R",
      "bv_replication_tables.R"
    )),
    outputs = file.path("paper", "data", c(
      "seer_laus_covariates.csv", "analysis_data.csv", "county_mortality_pre.csv",
      "county_votes.csv", "bv_covariates.csv", "urban_share.csv",
      "bv_table1_comparison.csv", "bv_table2_comparison.csv"
    ))
  ),
  analysis = list(
    scripts = file.path("scripts", "analyze_application.R"),
    # Terminal artifacts of the application analysis (see analyze_application.R):
    # paper (app_*) and supplement (sc_*, loo_results, dlps_results).
    outputs = file.path("results", c(
      "app_estimates.csv", "app_scalars.csv", "event_studies.csv",
      "heterogeneity.csv", "robustness.csv", "headline_comparison.csv",
      "placebo_intime.csv", "placebo_distribution.csv",
      "detrended_results.csv", "binned_results.csv", "event_study_draws.csv",
      "omega_size_bins.csv",
      "sc_results.csv", "sc_diagnostics.csv", "sc_loo.csv", "sc_intime.csv",
      "sc_event_studies.csv", "omega_weights.csv", "loo_results.csv",
      "dlps_results.csv", "dlps_placebo.csv", "dlps_ladder.csv", "sc_scalars.csv"
    ))
  ),
  mc = list(
    scripts = file.path("scripts", c("run_mc_simulations.R", "run_solutions_mc.R",
                                     "run_mc_trend.R")),
    outputs = file.path("paper", "data", c("mc_results.csv", "mc_solutions_results.csv",
                                           "mc_trend_results.csv"))
  ),
  render = list(
    scripts = file.path("paper", "_build_all.R"),
    outputs = character(0)  # always run when reached
  )
)
STAGE_ORDER <- names(STAGES)

# ---- argument parsing -------------------------------------------------------
cli      <- commandArgs(trailingOnly = TRUE)
force    <- "--force" %in% cli
from_ix  <- grep("^--from$", cli)
from_stage <- if (length(from_ix)) cli[from_ix + 1L] else NA_character_
named_stages <- setdiff(cli, c("--force", "--from", from_stage))

if (!is.na(from_stage) && !from_stage %in% STAGE_ORDER)
  stop("--from must be one of: ", paste(STAGE_ORDER, collapse = ", "), call. = FALSE)
if (length(named_stages) && !all(named_stages %in% STAGE_ORDER))
  stop("unknown stage(s): ", paste(setdiff(named_stages, STAGE_ORDER), collapse = ", "),
       "\nvalid stages: ", paste(STAGE_ORDER, collapse = ", "), call. = FALSE)

# Which stages are in scope this run?
selected <- if (length(named_stages)) {
  STAGE_ORDER[STAGE_ORDER %in% named_stages]
} else if (!is.na(from_stage)) {
  STAGE_ORDER[seq(match(from_stage, STAGE_ORDER), length(STAGE_ORDER))]
} else {
  STAGE_ORDER
}
# --from and --force imply "run regardless of existing outputs".
force_run <- force || !is.na(from_stage) || length(named_stages) > 0

# ---- runner -----------------------------------------------------------------
run_script <- function(path) {
  if (!file.exists(path)) stop("missing script: ", path, call. = FALSE)
  cat(sprintf("  -> Rscript %s\n", path))
  status <- system2("Rscript", shQuote(path))
  if (status != 0L) stop(sprintf("'%s' exited with status %d", path, status), call. = FALSE)
  invisible(TRUE)
}

stage_complete <- function(st) {
  out <- STAGES[[st]]$outputs
  length(out) > 0 && all(file.exists(out))
}

t_start <- Sys.time()
for (st in STAGE_ORDER) {
  if (!st %in% selected) next
  is_render <- st == "render"
  if (!force_run && !is_render && stage_complete(st)) {
    cat(sprintf("[skip] %-9s outputs present\n", st))
    next
  }
  cat(sprintf("[run ] %-9s\n", st))
  for (s in STAGES[[st]]$scripts) run_script(s)
}
cat(sprintf("\n==== run_all complete in %.1f min ====\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
