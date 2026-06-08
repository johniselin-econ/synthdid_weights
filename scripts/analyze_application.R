#!/usr/bin/env Rscript
# =============================================================================
# analyze_application.R — run the expensive ACA application analysis ONCE and
# write reduced result tables to results/, so the paper/supplement Rmds can
# render without re-computing anything heavy.
#
# This is the canonical home for the estimation, bootstrap/jackknife/placebo
# standard errors, event studies, robustness sweeps, and placebo distributions
# that previously ran inside the paper's knitr chunks. The Rmds now read the
# CSVs this script produces and only format/plot them.
#
# Run from the repository root (the orchestrator does this for you):
#   Rscript scripts/analyze_application.R
#
# Fast smoke mode (tiny replication counts; point estimates are unaffected and
# still exact, only the SEs/placebo p-values are rough) for validating that the
# pipeline runs end to end:
#   ANALYZE_FAST=1 Rscript scripts/analyze_application.R
#
# Inputs : paper/data/analysis_data.csv, paper/data/urban_share.csv
# Outputs: results/{app_estimates,app_scalars,event_studies,robustness,
#                   headline_comparison,heterogeneity,placebo_intime,
#                   placebo_distribution}.csv  (+ _manifest.csv)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(purrr)
  library(parallel)
})

# ---- load the weighted-SDID functions from source (no installed package) ----
invisible(lapply(list.files("R", pattern = "[.][Rr]$", full.names = TRUE), source))

set.seed(20240101)  # matches the Rmd setup chunk

# Replication counts: full for the committed run, tiny under ANALYZE_FAST=1.
FAST   <- nzchar(Sys.getenv("ANALYZE_FAST"))
B_SE   <- if (FAST)   5L else 500L   # bootstrap reps for SEs / event studies
B_PB   <- if (FAST)  10L else 500L   # paired-bootstrap reps
B_IT   <- if (FAST)   5L else 200L   # in-time placebo bootstrap reps
N_PLAC <- if (FAST)  20L else 500L   # in-space placebo draws
if (FAST) message("** ANALYZE_FAST: reduced replications; SEs/p-values are rough **")

dir.create("results", showWarnings = FALSE)
out <- function(x, name) write_csv(x, file.path("results", paste0(name, ".csv")))

# =============================================================================
# Parallel bootstrap machinery
# =============================================================================
# The covariate-adjusted estimator costs ~112s per fit, so bootstrap SEs for
# the X-variant Table-1 cells are ~100 CPU-hours. We parallelize at the
# replication level across a PSOCK cluster (Windows-safe; honours a Slurm
# allocation). `boot_rep` copies the exact resample/renormalize/refit logic of
# vcov.R's unit and cluster bootstrap theta() closures, so the *algorithm* is
# identical to vcov(); only the RNG stream differs (each rep is seeded
# independently), which shifts the SEs by Monte-Carlo noise (point estimates
# are unaffected). SE = sqrt((R-1)/R) * sd(estimates), matching bootstrap_se().
boot_rep <- function(estimate, method, cluster_vec, seed) {
  set.seed(seed)
  setup <- attr(estimate, 'setup'); opts <- attr(estimate, 'opts')
  weights <- attr(estimate, 'weights')
  treated.weights <- attr(estimate, 'treated.weights')
  period.weights <- attr(estimate, 'period.weights')
  N0 <- setup$N0; N <- nrow(setup$Y)
  if (method == "unit") {
    ind <- sample(1:N, replace = TRUE)
    control.ind <- sort(ind[ind <= N0]); treated.ind <- sort(ind[ind > N0])
    treated.local <- treated.ind - N0
    if (length(control.ind) == 0 || length(treated.ind) == 0) return(NA_real_)
    weights.boot <- weights; weights.boot$omega <- sum_normalize(weights$omega[control.ind])
    tw.boot <- sum_normalize(treated.weights[treated.local])
    all.ind <- c(control.ind, treated.ind); N0.boot <- length(control.ind)
  } else { # cluster
    cc <- cluster_vec[1:N0]; ct <- cluster_vec[(N0 + 1):N]
    drawn <- sample(unique(cluster_vec), replace = TRUE)
    control.ind  <- unlist(lapply(drawn, function(g) which(cc == g)))
    treated.local <- unlist(lapply(drawn, function(g) which(ct == g)))
    if (length(control.ind) == 0 || length(treated.local) == 0) return(NA_real_)
    weights.boot <- weights; weights.boot$omega <- sum_normalize(weights$omega[control.ind])
    tw.boot <- sum_normalize(treated.weights[treated.local])
    all.ind <- c(control.ind, N0 + treated.local); N0.boot <- length(control.ind)
  }
  Yb <- setup$Y[all.ind, , drop = FALSE]; Xb <- setup$X[all.ind, , , drop = FALSE]
  tryCatch(as.numeric(c(do.call(synthdid_estimate_weighted,
    c(list(Y = Yb, N0 = N0.boot, T0 = setup$T0, X = Xb,
           treated.weights = tw.boot, period.weights = period.weights,
           weights = weights.boot), opts)))),
    error = function(e) NA_real_)
}

n_cores <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", detectCores() - 2L)))
message("Parallel workers: ", n_cores)
CL <- makeCluster(n_cores)
on.exit(stopCluster(CL), add = TRUE)
# PSOCK workers do not reliably inherit the master's working directory, so pass
# the repo root explicitly and source R/ from there (else workers lack the
# synthdid functions and every bootstrap rep silently fails).
ROOT_ <- normalizePath(getwd())
clusterExport(CL, "ROOT_")
worker_ok <- clusterEvalQ(CL, {
  setwd(ROOT_)
  suppressPackageStartupMessages({library(dplyr)})
  invisible(lapply(list.files("R", pattern = "[.][Rr]$", full.names = TRUE), source))
  exists("synthdid_estimate_weighted") && exists("sum_normalize")
})
if (!all(unlist(worker_ok))) stop("workers failed to load R/ functions", call. = FALSE)
clusterExport(CL, "boot_rep")

# Distribute `replications` reps across CL (load-balanced) and pool the valid
# estimates. Rejected draws (all-one-side / empty side) are vanishingly rare on
# this panel, so we run exactly `replications` first and only top up the few
# NAs, rather than over-generating up front.
par_boot_se <- function(estimate, method, cluster_vec, replications, seed_base) {
  EST_ <- estimate; METH_ <- method; CLU_ <- cluster_vec
  clusterExport(CL, c("EST_", "METH_", "CLU_"), envir = environment())
  run_seeds <- function(seeds)
    unlist(parLapplyLB(CL, seeds, function(sd) boot_rep(EST_, METH_, CLU_, sd)))
  s <- seed_base
  seeds <- s + seq_len(replications); s <- s + replications
  collected <- run_seeds(seeds); collected <- collected[!is.na(collected)]
  while (length(collected) < replications) {        # top up only the shortfall
    need <- replications - length(collected)
    seeds <- s + seq_len(need); s <- s + need
    more <- run_seeds(seeds); collected <- c(collected, more[!is.na(more)])
  }
  collected <- collected[seq_len(replications)]
  sqrt((replications - 1) / replications) * sd(collected)
}

# =============================================================================
# 1. Load + prepare data  (paper chunks: load-data, prepare-sdid-panel)
# =============================================================================
panel <- read_csv(file.path("paper", "data", "analysis_data.csv"),
                  show_col_types = FALSE) %>%
  select(unit = fips, time = year, y = crude_rate,
         treated_unit = expansion, pop = population, state_fips,
         pct_white, pct_55_64, log_20_64, log_35_44, log_f_20_64, unemp) %>%
  mutate(post = as.integer(time >= 2014))

n_counties <- n_distinct(panel$unit)
n_treated  <- n_distinct(panel$unit[panel$treated_unit == 1])
n_control  <- n_distinct(panel$unit[panel$treated_unit == 0])
year_range <- range(panel$time)
n_years    <- n_distinct(panel$time)
T0_desc    <- sum(sort(unique(panel$time)) < 2014)
T1_desc    <- n_years - T0_desc

treated_pop <- panel %>%
  filter(treated_unit == 1, time == 2013) %>%
  summarise(total = sum(pop, na.rm = TRUE), median = median(pop, na.rm = TRUE),
            p10 = quantile(pop, 0.10, na.rm = TRUE),
            p90 = quantile(pop, 0.90, na.rm = TRUE))

panel_sdid <- panel %>%
  arrange(unit, time) %>%
  select(unit, time, y, treated_unit, pop, state_fips) %>%
  mutate(.unit = as.factor(unit), .time = time,
         .W = as.integer(treated_unit == 1 & time >= 2014))

dfp   <- panel_sdid %>% select(.unit, .time, y, .W) %>% as.data.frame()
setup <- panel.matrices(dfp)

pop_2013      <- panel %>% filter(treated_unit == 1, time == 2013) %>% select(unit, pop)
treated_names <- rownames(setup$Y)[(setup$N0 + 1):nrow(setup$Y)]
pop_ordered   <- pop_2013 %>% mutate(unit_char = as.character(unit))
weight_order  <- match(treated_names, pop_ordered$unit_char)
treated_weights <- pop_ordered$pop[weight_order]
treated_weights <- treated_weights / sum(treated_weights, na.rm = TRUE)

unit_state  <- panel_sdid %>% distinct(unit, state_fips) %>% mutate(unit_char = as.character(unit))
all_names   <- rownames(setup$Y)
cluster_vec <- unit_state$state_fips[match(all_names, unit_state$unit_char)]
n_states    <- n_distinct(cluster_vec)

# Alternative weight specifications for robustness
pop_2010              <- panel %>% filter(treated_unit == 1, time == 2010) %>% select(unit, pop)
pop_ordered_2010      <- pop_2010 %>% mutate(unit_char = as.character(unit))
weight_order_2010     <- match(treated_names, pop_ordered_2010$unit_char)
treated_weights_2010  <- pop_ordered_2010$pop[weight_order_2010]
treated_weights_2010  <- treated_weights_2010 / sum(treated_weights_2010, na.rm = TRUE)
treated_weights_log   <- log(pop_ordered$pop[weight_order])
treated_weights_log   <- treated_weights_log / sum(treated_weights_log, na.rm = TRUE)
treated_weights_sqrt  <- sqrt(pop_ordered$pop[weight_order])
treated_weights_sqrt  <- treated_weights_sqrt / sum(treated_weights_sqrt, na.rm = TRUE)

# Time-varying covariates (BV2020's six controls) as an N x T x C array
covar_names <- c("pct_white", "pct_55_64", "log_20_64", "log_35_44",
                 "log_f_20_64", "unemp")
X_arr <- array(NA_real_, dim = c(nrow(setup$Y), ncol(setup$Y), length(covar_names)),
               dimnames = list(rownames(setup$Y), colnames(setup$Y), covar_names))
for (k in seq_along(covar_names)) {
  Xk_wide <- panel %>% select(unit, time, val = dplyr::all_of(covar_names[k])) %>%
    tidyr::pivot_wider(names_from = time, values_from = val)
  Xk_mat <- as.matrix(Xk_wide[, -1])
  rownames(Xk_mat) <- as.character(Xk_wide$unit)
  X_arr[, , k] <- Xk_mat[rownames(setup$Y), colnames(setup$Y)]
}
if (anyNA(X_arr)) {
  state_of <- unit_state$state_fips[match(rownames(setup$Y), unit_state$unit_char)]
  for (k in seq_along(covar_names)) {
    Xk <- X_arr[, , k]
    na_idx <- which(is.na(Xk), arr.ind = TRUE)
    for (r in seq_len(nrow(na_idx))) {
      i <- na_idx[r, 1]; t <- na_idx[r, 2]
      Xk[i, t] <- mean(X_arr[state_of == state_of[i], t, k], na.rm = TRUE)
    }
    X_arr[, , k] <- Xk
  }
}
stopifnot(!anyNA(X_arr))

# weight-stats
hhi         <- sum(treated_weights^2)
n1_eff      <- 1 / hhi
top10_share <- sum(sort(treated_weights, decreasing = TRUE)[1:10])
top1_share  <- max(treated_weights)
N1          <- nrow(setup$Y) - setup$N0

# =============================================================================
# 2. Estimates  (paper chunk: run-estimates)
# =============================================================================
uniform_weights <- rep(1 / N1, N1)
tau_did_eq    <- did_estimate_weighted(setup$Y, setup$N0, setup$T0, treated.weights = uniform_weights, cluster = cluster_vec)
tau_did_eq_x  <- did_estimate_weighted(setup$Y, setup$N0, setup$T0, X = X_arr, treated.weights = uniform_weights, cluster = cluster_vec)
tau_sdid_eq   <- synthdid_estimate_weighted(setup$Y, setup$N0, setup$T0, treated.weights = uniform_weights, cluster = cluster_vec)
tau_sdid_eq_x <- synthdid_estimate_weighted(setup$Y, setup$N0, setup$T0, X = X_arr, treated.weights = uniform_weights, cluster = cluster_vec)
tau_did_w     <- did_estimate_weighted(setup$Y, setup$N0, setup$T0, treated.weights = treated_weights, cluster = cluster_vec)
tau_did_w_x   <- did_estimate_weighted(setup$Y, setup$N0, setup$T0, X = X_arr, treated.weights = treated_weights, cluster = cluster_vec)
tau_sdid_w    <- synthdid_estimate_weighted(setup$Y, setup$N0, setup$T0, treated.weights = treated_weights, cluster = cluster_vec)
tau_sdid_w_x  <- synthdid_estimate_weighted(setup$Y, setup$N0, setup$T0, X = X_arr, treated.weights = treated_weights, cluster = cluster_vec)
tau_sdid <- synthdid_estimate(setup$Y, setup$N0, setup$T0)
tau_did  <- did_estimate(setup$Y, setup$N0, setup$T0)

# =============================================================================
# 3. Standard errors  (paper chunk: compute-ses)
# =============================================================================
estimates_t1 <- list(
  did_eq  = tau_did_eq,  did_eq_x  = tau_did_eq_x,
  sdid_eq = tau_sdid_eq, sdid_eq_x = tau_sdid_eq_x,
  did_w   = tau_did_w,   did_w_x   = tau_did_w_x,
  sdid_w  = tau_sdid_w,  sdid_w_x  = tau_sdid_w_x)
specs <- names(estimates_t1)

# Checkpointed parallel SEs. Each (spec, method) is appended to a partial file
# as it finishes, so a crash only costs the in-flight cell (the X-variant cells
# are the ~hour-each ones). se_unit uses unit resampling; se_cluster resamples
# states (the estimates carry a stored state cluster, as in the paper).
ses_ckpt <- file.path("results", "_ses_partial.csv")
done_ses <- if (file.exists(ses_ckpt)) read_csv(ses_ckpt, show_col_types = FALSE) else
  tibble(spec = character(), method = character(), se = numeric())
se_grid <- expand.grid(method = c("unit", "cluster"), spec = specs, stringsAsFactors = FALSE)
for (i in seq_len(nrow(se_grid))) {
  sp <- se_grid$spec[i]; mth <- se_grid$method[i]
  if (nrow(done_ses) && any(done_ses$spec == sp & done_ses$method == mth)) next
  cl_arg <- if (mth == "cluster") cluster_vec else NULL
  t0 <- Sys.time()
  se_val <- tryCatch(par_boot_se(estimates_t1[[sp]], mth, cl_arg, B_SE,
                                 seed_base = 20240101L + 1000L * i),
                     error = function(e) NA_real_)
  row <- tibble(spec = sp, method = mth, se = se_val)
  write_csv(row, ses_ckpt, append = file.exists(ses_ckpt))
  done_ses <- bind_rows(done_ses, row)
  message(sprintf("  SE %-9s %-7s = %8.3f  (%.1f min)", sp, mth, se_val,
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
get_se <- function(sp, mth) {
  v <- done_ses$se[done_ses$spec == sp & done_ses$method == mth]
  if (length(v)) v[1] else NA_real_
}
se_unit    <- sapply(specs, get_se, mth = "unit")
se_cluster <- sapply(specs, get_se, mth = "cluster")

app_estimates <- tibble(
  spec       = specs,
  estimate   = sapply(estimates_t1, as.numeric),
  se_unit    = as.numeric(se_unit),
  se_cluster = as.numeric(se_cluster)
)
out(app_estimates, "app_estimates")

se_sdid   <- se_cluster[["sdid_eq"]]; se_did   <- se_cluster[["did_eq"]]
se_sdid_w <- se_cluster[["sdid_w"]];  se_did_w <- se_cluster[["did_w"]]

# =============================================================================
# 4. Paired bootstrap divergence test  (paper chunk: paired-bootstrap-divergence)
# =============================================================================
paired_bootstrap <- function(Y, N0, T0, treated.weights, replications = B_PB, seed = 42) {
  set.seed(seed)
  N <- nrow(Y); diffs <- rep(NA, replications); count <- 0
  while (count < replications) {
    ind <- sample(1:N, replace = TRUE)
    if (all(ind <= N0) || all(ind > N0)) next
    ind_s <- sort(ind); N0_b <- sum(ind <= N0); Y_b <- Y[ind_s, ]
    tw_b  <- sum_normalize(treated.weights[sort(ind[ind > N0]) - N0])
    est_eq <- tryCatch(as.numeric(synthdid_estimate(Y_b, N0_b, T0)), error = function(e) NA)
    est_wt <- tryCatch(as.numeric(synthdid_estimate_weighted(Y_b, N0_b, T0, treated.weights = tw_b)), error = function(e) NA)
    if (is.na(est_eq) || is.na(est_wt)) next
    diffs[count + 1] <- est_wt - est_eq; count <- count + 1
  }
  diffs
}
boot_diffs <- paired_bootstrap(setup$Y, setup$N0, setup$T0, treated_weights, replications = B_PB)
obs_diff   <- as.numeric(tau_sdid_w) - as.numeric(tau_sdid)
boot_null  <- boot_diffs - mean(boot_diffs)
p_diverge  <- mean(abs(boot_null) >= abs(obs_diff))

# =============================================================================
# 5. Event studies  (paper chunk: event-study-both)
# =============================================================================
es_sdid   <- synthdid_event_study(tau_sdid,   se.method = "bootstrap", replications = B_SE)
es_sdid_w <- synthdid_event_study(tau_sdid_w, se.method = "bootstrap", replications = B_SE)
es_did    <- synthdid_event_study(tau_did,    se.method = "bootstrap", replications = B_SE)
es_did_w  <- synthdid_event_study(tau_did_w,  se.method = "bootstrap", replications = B_SE)
es_df <- bind_rows(
  es_did    %>% mutate(Weighting = "Equally weighted",    Estimator = "DID"),
  es_did_w  %>% mutate(Weighting = "Population weighted", Estimator = "DID"),
  es_sdid   %>% mutate(Weighting = "Equally weighted",    Estimator = "SDID"),
  es_sdid_w %>% mutate(Weighting = "Population weighted", Estimator = "SDID")
) %>% mutate(year = as.numeric(time))
out(es_df, "event_studies")

# =============================================================================
# 6. Unit-level heterogeneity  (paper chunks: heterogeneity-plot, urban-binscatter)
#    Per-unit effects are the heavy part (need tau_sdid's omega/lambda); the
#    Rmd does the cheap binning/plotting from this CSV.
# =============================================================================
omega  <- attr(tau_sdid, 'weights')$omega
lambda <- attr(tau_sdid, 'weights')$lambda
Y      <- attr(tau_sdid, 'setup')$Y
N0_loc <- attr(tau_sdid, 'setup')$N0
T0_loc <- attr(tau_sdid, 'setup')$T0
Y_control <- Y[1:N0_loc, ]
Y_treated <- Y[(N0_loc + 1):nrow(Y), ]
Y_hat     <- matrix(rep(as.vector(omega %*% Y_control), nrow(Y_treated)),
                    nrow = nrow(Y_treated), byrow = TRUE)
gap          <- Y_treated - Y_hat
post_cols    <- (T0_loc + 1):ncol(Y)
unit_effects <- rowMeans(gap[, post_cols]) - as.vector(gap[, 1:T0_loc, drop = FALSE] %*% lambda)
stopifnot(abs(mean(unit_effects) - as.numeric(tau_sdid)) < 1e-8)

heterogeneity <- tibble(
  fips   = as.numeric(pop_ordered$unit_char[weight_order]),
  effect = unit_effects,
  pop    = pop_ordered$pop[weight_order]
)
out(heterogeneity, "heterogeneity")

# =============================================================================
# 7. Robustness to weight specifications  (paper chunk: robustness-table)
# =============================================================================
tau_sdid_log_app  <- synthdid_estimate_weighted(setup$Y, setup$N0, setup$T0, treated.weights = treated_weights_log,  cluster = cluster_vec)
tau_sdid_sqrt_app <- synthdid_estimate_weighted(setup$Y, setup$N0, setup$T0, treated.weights = treated_weights_sqrt, cluster = cluster_vec)
tau_sdid_2010_app <- synthdid_estimate_weighted(setup$Y, setup$N0, setup$T0, treated.weights = treated_weights_2010, cluster = cluster_vec)
# Estimates carry a stored state cluster, so the paper's vcov() defaulted to
# the cluster bootstrap here; preserve that.
se_sdid_log_app  <- tryCatch(par_boot_se(tau_sdid_log_app,  "cluster", cluster_vec, B_SE, 20240101L + 91000L), error = function(e) NA)
se_sdid_sqrt_app <- tryCatch(par_boot_se(tau_sdid_sqrt_app, "cluster", cluster_vec, B_SE, 20240101L + 92000L), error = function(e) NA)
se_sdid_2010_app <- tryCatch(par_boot_se(tau_sdid_2010_app, "cluster", cluster_vec, B_SE, 20240101L + 93000L), error = function(e) NA)

robustness <- tibble(
  spec     = c("pop_2013_baseline", "pop_2010", "sqrt_pop", "log_pop", "equal"),
  label    = c("Population (2013, baseline)", "Population (2010 base year)",
               "Square root of population", "Log population", "Equal ($1/N_1$)"),
  estimate = c(as.numeric(tau_sdid_w), as.numeric(tau_sdid_2010_app),
               as.numeric(tau_sdid_sqrt_app), as.numeric(tau_sdid_log_app), as.numeric(tau_sdid)),
  se       = c(se_sdid_w, se_sdid_2010_app, se_sdid_sqrt_app, se_sdid_log_app, se_sdid),
  n1_eff   = c(round(1/sum(treated_weights^2), 0), round(1/sum(treated_weights_2010^2), 0),
               round(1/sum(treated_weights_sqrt^2), 0), round(1/sum(treated_weights_log^2), 0),
               nrow(setup$Y) - setup$N0)
)
out(robustness, "robustness")

# =============================================================================
# 8. Headline comparison with BV2020  (paper chunk: headline-comparison-table)
# =============================================================================
headline_comparison <- tibble(
  method   = c("Standard SDID (equally weighted)", "Weighted SDID (population)",
               "Weighted DID (population)", "Borgschulte \\& Vogler (2020), Table 2 col.\\ 2"),
  estimate = c(as.numeric(tau_sdid), as.numeric(tau_sdid_w), as.numeric(tau_did_w), -11.36),
  se       = c(se_sdid, se_sdid_w, se_did_w, 3.59)  # BV2020 SE 3.59 (state-clustered, 45 clusters)
)
out(headline_comparison, "headline_comparison")

# =============================================================================
# 9. Placebo tests  (paper chunks: placebo-intime-app, placebo-inspace-app)
# =============================================================================
run_intime_placebo <- function(placebo_year, B = B_IT) {
  panel_pb <- panel %>% filter(time <= 2013) %>% arrange(unit, time) %>%
    mutate(.unit = as.factor(unit), .time = time,
           .W = as.integer(treated_unit == 1 & time >= placebo_year)) %>%
    select(.unit, .time, y, .W) %>% as.data.frame()
  setup_pb   <- panel.matrices(panel_pb)
  cl_pb      <- unit_state$state_fips[match(rownames(setup_pb$Y), unit_state$unit_char)]
  treated_pb <- rownames(setup_pb$Y)[(setup_pb$N0 + 1):nrow(setup_pb$Y)]
  N1_pb      <- nrow(setup_pb$Y) - setup_pb$N0
  tw_pb      <- pop_ordered$pop[match(treated_pb, pop_ordered$unit_char)]; tw_pb <- tw_pb / sum(tw_pb)
  fit_eq <- tryCatch(synthdid_estimate_weighted(setup_pb$Y, setup_pb$N0, setup_pb$T0, treated.weights = rep(1 / N1_pb, N1_pb), cluster = cl_pb), error = function(e) NULL)
  fit_wt <- tryCatch(synthdid_estimate_weighted(setup_pb$Y, setup_pb$N0, setup_pb$T0, treated.weights = tw_pb, cluster = cl_pb), error = function(e) NULL)
  se_of <- function(fit, sb) if (is.null(fit)) NA_real_ else
    tryCatch(par_boot_se(fit, "cluster", cl_pb, B, seed_base = sb), error = function(e) NA_real_)
  yr <- as.integer(placebo_year)
  tibble(placebo_year = as.character(placebo_year),
         estimate_eq = if (is.null(fit_eq)) NA_real_ else as.numeric(fit_eq), se_eq = se_of(fit_eq, 20240101L + 80000L + yr * 10L),
         estimate_wt = if (is.null(fit_wt)) NA_real_ else as.numeric(fit_wt), se_wt = se_of(fit_wt, 20240101L + 80000L + yr * 10L + 5L))
}
placebo_intime <- purrr::map_dfr(c(2011L, 2012L), run_intime_placebo)
out(placebo_intime, "placebo_intime")

set.seed(20240101)
placebo_ests <- rep(NA_real_, N_PLAC)
for (i in seq_len(N_PLAC)) {
  placebo_treated <- sample(seq_len(setup$N0), N1, replace = FALSE)
  placebo_control <- setdiff(seq_len(setup$N0), placebo_treated)
  Y_placebo  <- setup$Y[c(placebo_control, placebo_treated), ]
  N0_placebo <- length(placebo_control)
  placebo_ests[i] <- tryCatch(
    as.numeric(synthdid_estimate_weighted(Y_placebo, N0_placebo, setup$T0, treated.weights = rep(1/N1, N1))),
    error = function(e) NA_real_)
}
p_placebo <- mean(abs(placebo_ests) >= abs(as.numeric(tau_sdid_w)), na.rm = TRUE)
out(tibble(est = placebo_ests), "placebo_distribution")

# =============================================================================
# 10. Scalars that depend on setup/weights or heavy compute (inline numbers)
# =============================================================================
app_scalars <- tibble::tribble(
  ~name,               ~value,
  "n_counties",        n_counties,
  "n_treated",         n_treated,
  "n_control",         n_control,
  "year_min",          year_range[1],
  "year_max",          year_range[2],
  "n_years",           n_years,
  "T0",                T0_desc,
  "T1",                T1_desc,
  "treated_pop_total", treated_pop$total,
  "treated_pop_median",treated_pop$median,
  "treated_pop_p10",   treated_pop$p10,
  "treated_pop_p90",   treated_pop$p90,
  "hhi",               hhi,
  "n1_eff",            n1_eff,
  "N1",                N1,
  "top1_share",        top1_share,
  "top10_share",       top10_share,
  "n_states",          n_states,
  "obs_diff",          obs_diff,
  "p_diverge",         p_diverge,
  "p_placebo",         p_placebo
)
out(app_scalars, "app_scalars")

# =============================================================================
# manifest: seed, git SHA, build date, replication counts
# =============================================================================
git_sha <- tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE),
                    error = function(e) NA_character_)
manifest <- tibble(
  key   = c("seed", "git_sha", "built", "fast_mode", "B_SE", "B_PB", "B_IT", "N_PLAC"),
  value = c("20240101", as.character(git_sha[1]), as.character(Sys.time()),
            as.character(FAST), B_SE, B_PB, B_IT, N_PLAC)
)
out(manifest, "_manifest")

message(sprintf("\nDone. tau_sdid=%.2f  tau_sdid_w=%.2f  p_diverge=%.3f  p_placebo=%.3f",
                as.numeric(tau_sdid), as.numeric(tau_sdid_w), p_diverge, p_placebo))
