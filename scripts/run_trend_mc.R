## Monte Carlo trend cell (rework plan Phase 3, docs/rework_plan_2026-07.md):
## a size-correlated differential trend gamma_trend * (s_i/s_bar - 1) * t/T is
## added to L, violating assumption A5 for the population-weighted target only.
## The cell is the structural proof that the application's placebo failure is
## the DGP, not the weighted machinery: weighted SDID and weighted DID inherit
## the same bias, the equal-weight variants stay unbiased, the in-time placebo
## detects the violation, detrended weighted SDID recovers the target, and
## stratified weighted SDID recovers most of it.
##
## Design notes (where this cell deliberately departs from the main grid):
## - Sizes are drawn for ALL units (controls too, same lognormal), because the
##   application's trend is a property of large counties on both sides of the
##   treatment divide -- stratification can only work if same-size donors share
##   the trend. Deviations d_i = s_i / s_bar_group - 1 are centered within the
##   treated and control groups separately, so the equal-weighted treated
##   loading is exactly zero every draw and the population-weighted loading is
##   exactly N1*HHI - 1.
## - Panel dimensions match the application design (T0 = 5, T1 = 4), not the
##   grid's T0 = 40: the violation only binds when the pre-period is too short
##   for level matching to select trend-comparable donors, which is the regime
##   the application sits in (and T0 = 5 supports detrend = TRUE and an
##   in-time placebo with a 3-period training window, like onset 2012).
## - Point estimates only (bias/RMSE comparison, like run_solutions_mc.R), so
##   the sweep is cheap: 3 gamma_trend values x 200 sims, minutes on a node.
##
## Same seeding scheme as the other MC scripts: each (config, sim) task runs
## under set.seed(20240101 + sim). The weighted-SDID functions are sourced
## from R/ (like analyze_application.R) because the remedies (detrend,
## stratified) may postdate the installed synthdid.
##
## Writes paper/data/mc_trend_results.csv, read by the supplement's
## trend-mc chunk.
##
## Usage: Rscript scripts/run_trend_mc.R          (MC_NSIM=5 for a smoke test)

suppressPackageStartupMessages(library(parallel))
invisible(lapply(list.files("R", pattern = "[.][Rr]$", full.names = TRUE), source))

OUT_CSV <- "paper/data/mc_trend_results.csv"

# ---------------------------------------------------------------------------
# DGP: the main MC's factor + fixed-effects skeleton (run_mc_simulations.R)
# plus sizes for all units and the size-correlated differential trend.
# ---------------------------------------------------------------------------

## Lognormal sizes for all N0 + N1 units, with sdlog calibrated so a treated
## sample of N1 hits the target HHI on average (the same uniroot calibration
## as generate_sizes in run_mc_simulations.R, drawing N0 + N1 values).
generate_sizes_all <- function(N0, N1, target_hhi) {
  if (target_hhi <= 1/N1 + 0.005) return(rep(1, N0 + N1))
  saved_seed <- .Random.seed
  obj <- function(log_sd) {
    set.seed(999)
    s <- rlnorm(N1 * 100, meanlog = 0, sdlog = exp(log_sd))
    hhis <- sapply(split(s, rep(1:100, each = N1)), function(x) sum((x/sum(x))^2))
    mean(hhis) - target_hhi
  }
  opt <- tryCatch(uniroot(obj, c(-3, 3), tol = 0.001), error = function(e) list(root = 0.5))
  .Random.seed <<- saved_seed
  rlnorm(N0 + N1, meanlog = 0, sdlog = exp(opt$root))
}

run_one_sim_trend <- function(N0, N1, T0, T1, gamma, target_hhi, gamma_trend,
                              sigma = 0.5, tau_bar = 1) {
  N <- N0 + N1; T_ <- T0 + T1; rank <- 2; rho <- 0.7
  U <- matrix(rpois(rank * N, sqrt(sample(1:N)) / sqrt(N)), N, rank)
  V <- matrix(rpois(rank * T_, sqrt(1:T_) / sqrt(T_)), T_, rank)
  alpha <- outer(10 * sample(1:N) / N, rep(1, T_))
  beta  <- outer(rep(1, N), 10 * (1:T_) / T_)
  L     <- U %*% t(V) + alpha + beta
  var_mat <- outer(1:T_, 1:T_, FUN = function(x, y) rho^(abs(x - y)))
  error   <- mvtnorm::rmvnorm(N, sigma = var_mat, method = "chol")

  sizes <- generate_sizes_all(N0, N1, target_hhi)   # controls first, like Y
  sizes_ctl <- sizes[1:N0]; sizes_trt <- sizes[(N0+1):N]
  tw <- sizes_trt / sum(sizes_trt); s_bar <- mean(sizes_trt)
  ## Group-centered size deviations: equal-weighted treated loading is exactly
  ## 0, population-weighted loading is exactly N1*HHI - 1, control average 0.
  d <- c(sizes_ctl / mean(sizes_ctl) - 1, sizes_trt / s_bar - 1)
  L <- L + gamma_trend * outer(d, (1:T_) / T_)

  tau_vec <- tau_bar + gamma * (sizes_trt / s_bar - 1)
  true_att_equal    <- mean(tau_vec)
  true_att_weighted <- sum(tw * tau_vec)
  tau_mat <- matrix(0, N, T_)
  for (j in 1:N1) tau_mat[N0 + j, (T0+1):T_] <- tau_vec[j]
  Y <- L + tau_mat + sigma * error
  rownames(Y) <- 1:N; colnames(Y) <- 1:T_

  eqw <- rep(1 / N1, N1)
  est <- function(expr) tryCatch(as.numeric(expr), error = function(e) NA_real_)

  est_sdid_eq <- est(synthdid_estimate_weighted(Y, N0, T0, treated.weights = eqw))
  est_sdid_wt <- est(synthdid_estimate_weighted(Y, N0, T0, treated.weights = tw))
  est_did_eq  <- est(did_estimate_weighted(Y, N0, T0, treated.weights = eqw))
  est_did_wt  <- est(did_estimate_weighted(Y, N0, T0, treated.weights = tw))
  est_det_wt  <- est(synthdid_estimate_weighted(Y, N0, T0, treated.weights = tw,
                                                detrend = TRUE))
  ## Stratified: pooled-size quartile bins (the application uses pooled-2013
  ## population bins); drop.infeasible in case a bin draws < 2 controls.
  bins <- cut(sizes, quantile(sizes, c(0, .25, .5, .75, 1)),
              include.lowest = TRUE, labels = FALSE)
  est_strat_wt <- est(synthdid_estimate_stratified(Y, N0, T0, strata = bins,
                                                   treated.weights = tw,
                                                   drop.infeasible = TRUE))
  ## In-time placebo: pre-period only, onset moved to T0 - 1 (3-period
  ## training window at T0 = 5, like the application's 2012 onset); truth 0.
  est_plac_wt <- est(synthdid_estimate_weighted(Y[, 1:T0, drop = FALSE], N0,
                                                T0 - 2, treated.weights = tw))

  data.frame(est_sdid_eq = est_sdid_eq, est_sdid_wt = est_sdid_wt,
             est_did_eq = est_did_eq,   est_did_wt = est_did_wt,
             est_det_wt = est_det_wt,   est_strat_wt = est_strat_wt,
             est_plac_wt = est_plac_wt,
             true_att_equal = true_att_equal,
             true_att_weighted = true_att_weighted)
}

# ---------------------------------------------------------------------------
# Sweep: gamma_trend in {0, 0.5, 1} at the application-like design point.
# gamma_trend = 0 is the control cell (A5 holds; everything unbiased);
# with the popw loading N1*HHI - 1 = 2 and the post-pre gap of 0.5 in t/T
# units, the popw DID/SDID bias is approximately gamma_trend itself.
# ---------------------------------------------------------------------------

configs <- data.frame(gamma_trend = c(0, 0.5, 1.0))
configs$N1 <- 20; configs$gamma <- 0.5; configs$hhi <- 3 / configs$N1
configs$N0 <- 100; configs$T0 <- 5; configs$T1 <- 4
configs$config <- seq_len(nrow(configs))

N_SIM <- as.integer(Sys.getenv("MC_NSIM", "200"))  # override for smoke tests
if (N_SIM < 200) OUT_CSV <- sub("\\.csv$", "_smoke.csv", OUT_CSV)

tasks <- merge(configs, data.frame(sim = seq_len(N_SIM)))
tasks <- tasks[order(tasks$config, tasks$sim), ]

## Resume support: skip (config, sim) pairs already in the output file
if (file.exists(OUT_CSV)) {
  done <- read.csv(OUT_CSV)[, c("config", "sim")]
  before <- nrow(tasks)
  tasks <- tasks[!paste(tasks$config, tasks$sim) %in%
                   paste(done$config, done$sim), ]
  message(sprintf("Resuming: %d of %d tasks already done", before - nrow(tasks),
                  before))
}

message(sprintf("Trend MC sweep: %d configs x %d sims; %d tasks remaining",
                nrow(configs), N_SIM, nrow(tasks)))

if (nrow(tasks) > 0) {
  n_cores <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK",
                                          detectCores() - 2)))
  message("Workers: ", n_cores)
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)
  invisible(clusterEvalQ(cl, invisible(
    lapply(list.files("R", pattern = "[.][Rr]$", full.names = TRUE), source))))
  clusterExport(cl, c("generate_sizes_all", "run_one_sim_trend"))

  batch_size <- n_cores * 9
  batches <- split(seq_len(nrow(tasks)), ceiling(seq_len(nrow(tasks)) / batch_size))
  t0 <- Sys.time()
  for (b in seq_along(batches)) {
    idx <- batches[[b]]
    results <- parLapply(cl, split(tasks[idx, ], seq_along(idx)), function(tk) {
      set.seed(20240101 + tk$sim)   # matches the main MC's seeding scheme
      out <- run_one_sim_trend(tk$N0, tk$N1, tk$T0, tk$T1, tk$gamma, tk$hhi,
                               gamma_trend = tk$gamma_trend)
      out$config <- tk$config; out$gamma_trend <- tk$gamma_trend
      out$N1 <- tk$N1; out$gamma <- tk$gamma; out$hhi <- tk$hhi; out$sim <- tk$sim
      out
    })
    batch_df <- do.call(rbind, results)
    write.table(batch_df, OUT_CSV, sep = ",", row.names = FALSE,
                col.names = !file.exists(OUT_CSV), append = file.exists(OUT_CSV))
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    eta <- el / b * (length(batches) - b)
    message(sprintf("batch %d/%d done (%.1f min elapsed, ~%.1f min remaining)",
                    b, length(batches), el, eta))
  }
}

trend_df <- read.csv(OUT_CSV)
message(sprintf("%s now has %d rows", OUT_CSV, nrow(trend_df)))
