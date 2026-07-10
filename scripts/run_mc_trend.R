## Monte Carlo "trend cell" for the main paper's Section 4 (diagnosis only).
##
## Adds a size-correlated linear trend to the factor structure of the TREATED
## units, L_it -> L_it + gamma_trend * (s_j / s_bar - 1) * t, holding the
## treatment effect fixed. This violates the oracle-weights assumption for the
## population-weighted target while leaving the equally-weighted target
## essentially unexposed (the treated-average trend is zero because
## mean(s_j/s_bar - 1) = 0, but the population-weighted average is positive).
## The cell shows that weighted SDID and weighted DID inherit the resulting
## bias while the equally-weighted estimators do not, and that the in-time
## placebo detects it. Detrended/stratified REMEDIES are deliberately not
## computed here: the paper's use of this cell is diagnostic only.
##
## Same DGP and seeding scheme as scripts/run_mc_simulations.R: each
## (gamma_trend, sim) task runs under set.seed(20240101 + sim). Because the
## seed depends only on `sim` and gamma_trend enters deterministically AFTER
## all random draws, the three gamma_trend cells share an identical base panel
## within each sim -- the only thing that changes across them is the trend.
##
## Writes paper/data/mc_trend_results.csv, read by the main paper's
## mc-trend-run / mc-trend-table chunks.
##
## Usage: Rscript scripts/run_mc_trend.R

suppressPackageStartupMessages({
  library(parallel)
  library(synthdid)
})

OUT_CSV <- "paper/data/mc_trend_results.csv"

# ---------------------------------------------------------------------------
# DGP (identical draws to scripts/run_mc_simulations.R, plus the trend term)
# ---------------------------------------------------------------------------

generate_sizes <- function(N1, target_hhi) {
  if (target_hhi <= 1/N1 + 0.005) return(rep(1, N1))
  saved_seed <- .Random.seed
  obj <- function(log_sd) {
    set.seed(999)
    s <- rlnorm(N1 * 100, meanlog = 0, sdlog = exp(log_sd))
    hhis <- sapply(split(s, rep(1:100, each = N1)), function(x) sum((x/sum(x))^2))
    mean(hhis) - target_hhi
  }
  opt <- tryCatch(uniroot(obj, c(-3, 3), tol = 0.001), error = function(e) list(root = 0.5))
  .Random.seed <<- saved_seed
  rlnorm(N1, meanlog = 0, sdlog = exp(opt$root))
}

run_one_trend_sim <- function(N0, N1, T0, T1, gamma, target_hhi, gamma_trend,
                              sigma = 0.5, tau_bar = 1) {
  N <- N0 + N1; T_ <- T0 + T1; rank <- 2; rho <- 0.7
  U <- matrix(rpois(rank * N, sqrt(sample(1:N)) / sqrt(N)), N, rank)
  V <- matrix(rpois(rank * T_, sqrt(1:T_) / sqrt(T_)), T_, rank)
  alpha <- outer(10 * sample(1:N) / N, rep(1, T_))
  beta  <- outer(rep(1, N), 10 * (1:T_) / T_)
  L     <- U %*% t(V) + alpha + beta
  var_mat <- outer(1:T_, 1:T_, FUN = function(x, y) rho^(abs(x - y)))
  error   <- mvtnorm::rmvnorm(N, sigma = var_mat, method = "chol")
  sizes <- generate_sizes(N1, target_hhi)
  tw <- sizes / sum(sizes); s_bar <- mean(sizes)
  tau_vec <- tau_bar + gamma * (sizes / s_bar - 1)
  true_att_equal    <- mean(tau_vec)
  true_att_weighted <- sum(tw * tau_vec)
  tau_mat <- matrix(0, N, T_)
  for (j in 1:N1) tau_mat[N0 + j, (T0+1):T_] <- tau_vec[j]

  ## Size-correlated trend on the treated units, linear in normalized time
  ## t/T_ in [0,1] so its total accumulation over the panel is on the same
  ## scale as tau_bar (a trend comparable to the effect, not one that swamps
  ## it). Perturbs the factor structure L, not the treatment effect tau, so
  ## the target ATT is unchanged across gamma_trend cells.
  trend_mat <- matrix(0, N, T_)
  for (j in 1:N1) trend_mat[N0 + j, ] <- gamma_trend * (sizes[j] / s_bar - 1) * ((1:T_) / T_)

  Y <- L + tau_mat + trend_mat + sigma * error
  rownames(Y) <- 1:N; colnames(Y) <- 1:T_

  est_sdid_eq <- tryCatch(as.numeric(synthdid_estimate(Y, N0, T0)), error = function(e) NA_real_)
  est_wt_obj  <- tryCatch(synthdid_estimate_weighted(Y, N0, T0, treated.weights = tw), error = function(e) NULL)
  est_sdid_wt <- if (is.null(est_wt_obj)) NA_real_ else as.numeric(est_wt_obj)
  est_did_eq  <- tryCatch(as.numeric(did_estimate(Y, N0, T0)), error = function(e) NA_real_)
  est_did_wt  <- tryCatch(as.numeric(did_estimate_weighted(Y, N0, T0, treated.weights = tw)), error = function(e) NA_real_)
  est_plac_wt <- if (is.null(est_wt_obj)) NA_real_ else
                 tryCatch(as.numeric(synthdid_placebo_weighted(est_wt_obj)), error = function(e) NA_real_)

  data.frame(est_sdid_eq = est_sdid_eq, est_sdid_wt = est_sdid_wt,
             est_did_eq = est_did_eq, est_did_wt = est_did_wt,
             est_plac_wt = est_plac_wt,
             true_att_equal = true_att_equal, true_att_weighted = true_att_weighted)
}

# ---------------------------------------------------------------------------
# Sweep: three gamma_trend cells at a single (N1, gamma, HHI) configuration
# ---------------------------------------------------------------------------

configs <- data.frame(gamma_trend = c(0, 0.5, 1.0))
configs$N1 <- 20; configs$gamma <- 0.5; configs$hhi <- 0.15
configs$config <- seq_len(nrow(configs))
N0 <- 100; T0 <- 40; T1 <- 10
N_SIM <- as.integer(Sys.getenv("MC_NSIM", "200"))   # override for smoke tests
if (N_SIM < 200) OUT_CSV <- sub("\\.csv$", "_smoke.csv", OUT_CSV)

tasks <- merge(configs, data.frame(sim = seq_len(N_SIM)))
tasks <- tasks[order(tasks$config, tasks$sim), ]

## Resume support: skip (config, sim) pairs already in the output file
if (file.exists(OUT_CSV)) {
  done <- read.csv(OUT_CSV)[, c("config", "sim")]
  before <- nrow(tasks)
  tasks <- tasks[!paste(tasks$config, tasks$sim) %in% paste(done$config, done$sim), ]
  message(sprintf("Resuming: %d of %d tasks already done", before - nrow(tasks), before))
}

message(sprintf("Trend MC sweep: %d gamma_trend cells x %d sims; %d tasks remaining",
                nrow(configs), N_SIM, nrow(tasks)))

if (nrow(tasks) > 0) {
  n_cores <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", detectCores() - 2)))
  message("Workers: ", n_cores)
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)
  invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(synthdid))))
  clusterExport(cl, c("generate_sizes", "run_one_trend_sim", "N0", "T0", "T1"))

  batch_size <- n_cores * 9
  batches <- split(seq_len(nrow(tasks)), ceiling(seq_len(nrow(tasks)) / batch_size))
  for (b in seq_along(batches)) {
    idx <- batches[[b]]
    results <- parLapply(cl, split(tasks[idx, ], seq_along(idx)), function(tk) {
      set.seed(20240101 + tk$sim)   # matches the main MC's seeding scheme
      out <- run_one_trend_sim(N0, tk$N1, T0, T1, tk$gamma, tk$hhi, tk$gamma_trend)
      out$config <- tk$config; out$gamma_trend <- tk$gamma_trend; out$N1 <- tk$N1
      out$gamma <- tk$gamma; out$hhi <- tk$hhi; out$sim <- tk$sim
      out
    })
    batch_df <- do.call(rbind, results)
    write.table(batch_df, OUT_CSV, sep = ",", row.names = FALSE,
                col.names = !file.exists(OUT_CSV), append = file.exists(OUT_CSV))
    message(sprintf("batch %d/%d done", b, length(batches)))
  }
}

trend_df <- read.csv(OUT_CSV)
message(sprintf("%s now has %d rows", OUT_CSV, nrow(trend_df)))
