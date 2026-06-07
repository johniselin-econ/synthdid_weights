## Monte Carlo comparison of the three strategies for incorporating
## treated-unit weights into SDID (supplement Section F):
##   Solution 1: collapse treated units to a single weighted-average row,
##               then run standard SDID with one treated unit;
##   Solution 2: run a separate single-treated-unit SDID for each treated
##               unit and average the estimates with weights tilde-omega;
##   Solution 3: the modified-collapsed-form estimator adopted in the paper
##               (synthdid_estimate_weighted).
##
## Same DGP, configuration grid, and seeding scheme as the main Monte Carlo
## (scripts/run_mc_simulations.R): each (config, sim) task runs under
## set.seed(20240101 + sim). No variance estimation is needed (the
## comparison is RMSE-only), so the sweep is far cheaper than the main MC.
##
## Writes paper/data/mc_solutions_results.csv, read by the supplement's
## solutions-mc-run chunk.
##
## Usage: Rscript scripts/run_solutions_mc.R

suppressPackageStartupMessages({
  library(parallel)
  library(synthdid)
})

OUT_CSV <- "paper/data/mc_solutions_results.csv"

# ---------------------------------------------------------------------------
# DGP (identical to scripts/run_mc_simulations.R)
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

run_one_sim_solutions <- function(N0, N1, T0, T1, gamma, target_hhi,
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
  true_att_weighted <- sum(tw * tau_vec)
  tau_mat <- matrix(0, N, T_)
  for (j in 1:N1) tau_mat[N0 + j, (T0+1):T_] <- tau_vec[j]
  Y <- L + tau_mat + sigma * error
  rownames(Y) <- 1:N; colnames(Y) <- 1:T_

  Y_control <- Y[1:N0, , drop = FALSE]
  Y_treated <- Y[(N0+1):N, , drop = FALSE]

  ## Solution 1: collapse treated rows to one weighted-average unit
  Y_sol1 <- rbind(Y_control, `wavg` = as.vector(tw %*% Y_treated))
  est_sol1 <- tryCatch(as.numeric(synthdid_estimate(Y_sol1, N0, T0)),
                       error = function(e) NA_real_)

  ## Solution 2: unit-by-unit SDID, weighted average of unit estimates
  est_sol2 <- tryCatch({
    taus <- vapply(1:N1, function(j) {
      Yj <- rbind(Y_control, Y_treated[j, , drop = FALSE])
      as.numeric(synthdid_estimate(Yj, N0, T0))
    }, numeric(1))
    sum(tw * taus)
  }, error = function(e) NA_real_)

  ## Solution 3: modified collapsed form (the paper's estimator)
  est_sol3 <- tryCatch(as.numeric(
    synthdid_estimate_weighted(Y, N0, T0, treated.weights = tw)),
    error = function(e) NA_real_)

  data.frame(est_sol1 = est_sol1, est_sol2 = est_sol2, est_sol3 = est_sol3,
             true_att_weighted = true_att_weighted)
}

# ---------------------------------------------------------------------------
# Parallel sweep (same grid as the main MC; 200 sims per configuration)
# ---------------------------------------------------------------------------

configs <- expand.grid(N1 = c(5, 10, 20), gamma = c(0, 0.5, 1.0), hhi_mult = c(1.5, 3))
configs$hhi <- configs$hhi_mult / configs$N1
configs$config <- seq_len(nrow(configs))
N_SIM <- as.integer(Sys.getenv("MC_NSIM", "200"))  # override for smoke tests
N0 <- 100; T0 <- 40; T1 <- 10
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

message(sprintf("Solutions MC sweep: %d configs x %d sims; %d tasks remaining",
                nrow(configs), N_SIM, nrow(tasks)))

if (nrow(tasks) > 0) {
  n_cores <- max(1, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK",
                                          detectCores() - 2)))
  message("Workers: ", n_cores)
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)
  invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(synthdid))))
  clusterExport(cl, c("generate_sizes", "run_one_sim_solutions",
                      "N0", "T0", "T1"))

  batch_size <- n_cores * 9
  batches <- split(seq_len(nrow(tasks)), ceiling(seq_len(nrow(tasks)) / batch_size))
  t0 <- Sys.time()
  for (b in seq_along(batches)) {
    idx <- batches[[b]]
    results <- parLapply(cl, split(tasks[idx, ], seq_along(idx)), function(tk) {
      set.seed(20240101 + tk$sim)   # matches the main MC's seeding scheme
      out <- run_one_sim_solutions(N0, tk$N1, T0, T1, tk$gamma, tk$hhi)
      out$config <- tk$config; out$N1 <- tk$N1
      out$gamma <- tk$gamma; out$hhi <- tk$hhi; out$sim <- tk$sim
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

sol_df <- read.csv(OUT_CSV)
message(sprintf("%s now has %d rows", OUT_CSV, nrow(sol_df)))
