## Run the supplement's Monte Carlo study (Section D) in parallel and write
## paper/data/mc_results.csv, which the supplement's mc-run chunk reads.
##
## This is the canonical copy of the simulation code; the mc-functions chunk
## in paper/weighted_sdid_supplement.Rmd mirrors it for documentation (kept
## eval=FALSE there). Seeding matches the original sequential design exactly:
## each (config, sim) task runs under set.seed(20240101 + sim), so results
## are identical to the single-threaded loop.
##
## Design: 18 configurations (N1 in {5,10,20} x gamma in {0,0.5,1} x
## HHI in {1.5/N1, 3/N1}), 100 sims each, B = 50 bootstrap/placebo reps.
## Single-threaded estimate ~50h; wall time here is roughly that divided by
## the number of workers.
##
## Usage: Rscript scripts/run_mc_simulations.R

suppressPackageStartupMessages({
  library(parallel)
  library(synthdid)
})

OUT_CSV <- "paper/data/mc_results.csv"

# ---------------------------------------------------------------------------
# Simulation functions (mirrored in the supplement's mc-functions chunk)
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

run_one_sim <- function(N0, N1, T0, T1, gamma, target_hhi,
                        sigma = 0.5, tau_bar = 1, B = 200) {
  N <- N0 + N1; T_ <- T0 + T1; rank <- 2; rho <- 0.7
  U <- matrix(rpois(rank * N, sqrt(sample(1:N)) / sqrt(N)), N, rank)
  V <- matrix(rpois(rank * T_, sqrt(1:T_) / sqrt(T_)), T_, rank)
  alpha <- outer(10 * sample(1:N) / N, rep(1, T_))
  beta  <- outer(rep(1, N), 10 * (1:T_) / T_)
  L     <- U %*% t(V) + alpha + beta
  var_mat <- outer(1:T_, 1:T_, FUN = function(x, y) rho^(abs(x - y)))
  error   <- mvtnorm::rmvnorm(N, sigma = var_mat, method = "chol")
  W <- matrix(0, N, T_); W[(N0+1):N, (T0+1):T_] <- 1
  sizes <- generate_sizes(N1, target_hhi)
  tw <- sizes / sum(sizes); s_bar <- mean(sizes)
  tau_vec <- tau_bar + gamma * (sizes / s_bar - 1)
  true_att_equal    <- mean(tau_vec)
  true_att_weighted <- sum(tw * tau_vec)
  tau_mat <- matrix(0, N, T_)
  for (j in 1:N1) tau_mat[N0 + j, (T0+1):T_] <- tau_vec[j]
  Y <- L + tau_mat + sigma * error
  rownames(Y) <- 1:N; colnames(Y) <- 1:T_
  est_unwt <- tryCatch(synthdid_estimate(Y, N0, T0), error = function(e) NA)
  est_wt   <- tryCatch(synthdid_estimate_weighted(Y, N0, T0, treated.weights = tw),
                        error = function(e) NA)
  if (is.na(est_unwt[1]) || is.na(est_wt[1])) {
    return(data.frame(est_unwt = NA, est_wt = NA,
                      se_boot_unwt = NA, se_jk_unwt = NA, se_plac_unwt = NA,
                      se_boot_wt = NA,   se_jk_wt = NA,   se_plac_wt = NA,
                      true_att_equal = true_att_equal,
                      true_att_weighted = true_att_weighted))
  }
  se_boot_unwt <- tryCatch(sqrt(vcov(est_unwt, method = "bootstrap", replications = B)), error = function(e) NA)
  se_jk_unwt   <- tryCatch(sqrt(vcov(est_unwt, method = "jackknife")), error = function(e) NA)
  se_plac_unwt <- tryCatch(sqrt(vcov(est_unwt, method = "placebo", replications = B)), error = function(e) NA)
  se_boot_wt   <- tryCatch(sqrt(vcov(est_wt, method = "bootstrap", replications = B)), error = function(e) NA)
  se_jk_wt     <- tryCatch(sqrt(vcov(est_wt, method = "jackknife")), error = function(e) NA)
  se_plac_wt   <- tryCatch(sqrt(vcov(est_wt, method = "placebo", replications = B)), error = function(e) NA)
  data.frame(
    est_unwt = c(est_unwt), est_wt = c(est_wt),
    se_boot_unwt = c(se_boot_unwt), se_jk_unwt = c(se_jk_unwt), se_plac_unwt = c(se_plac_unwt),
    se_boot_wt   = c(se_boot_wt),   se_jk_wt   = c(se_jk_wt),   se_plac_wt   = c(se_plac_wt),
    true_att_equal = true_att_equal, true_att_weighted = true_att_weighted)
}

# ---------------------------------------------------------------------------
# Parallel sweep
# ---------------------------------------------------------------------------

configs <- expand.grid(N1 = c(5, 10, 20), gamma = c(0, 0.5, 1.0), hhi_mult = c(1.5, 3))
configs$hhi <- configs$hhi_mult / configs$N1
configs$config <- seq_len(nrow(configs))
N_SIM <- as.integer(Sys.getenv("MC_NSIM", "100"))  # override for smoke tests
B_SE  <- 50
N0 <- 100; T0 <- 40; T1 <- 10
if (N_SIM < 100) OUT_CSV <- sub("\\.csv$", "_smoke.csv", OUT_CSV)

tasks <- merge(configs, data.frame(sim = seq_len(N_SIM)))
tasks <- tasks[order(tasks$config, tasks$sim), ]

## Resume support: skip (config, sim) pairs already in the output file, so a
## crash or shutdown only costs the in-flight batch.
if (file.exists(OUT_CSV)) {
  done <- read.csv(OUT_CSV)[, c("config", "sim")]
  before <- nrow(tasks)
  tasks <- tasks[!paste(tasks$config, tasks$sim) %in%
                   paste(done$config, done$sim), ]
  message(sprintf("Resuming: %d of %d tasks already done", before - nrow(tasks),
                  before))
}

message(sprintf("MC sweep: %d configs x %d sims; %d tasks remaining, B = %d",
                nrow(configs), N_SIM, nrow(tasks), B_SE))

if (nrow(tasks) > 0) {
  n_cores <- max(1, detectCores() - 2)
  message("Workers: ", n_cores)
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)
  invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(synthdid))))
  clusterExport(cl, c("generate_sizes", "run_one_sim", "N0", "T0", "T1", "B_SE"))

  progress_file <- "paper/_mc_progress.log"
  cat(sprintf("[%s] parallel MC start: %d tasks on %d workers\n",
              format(Sys.time()), nrow(tasks), n_cores),
      file = progress_file, append = TRUE)

  ## Run in batches: progress + incremental checkpointing
  batch_size <- n_cores * 9
  batches <- split(seq_len(nrow(tasks)), ceiling(seq_len(nrow(tasks)) / batch_size))
  t0 <- Sys.time()
  for (b in seq_along(batches)) {
    idx <- batches[[b]]
    results <- parLapply(cl, split(tasks[idx, ], seq_along(idx)), function(tk) {
      set.seed(20240101 + tk$sim)   # matches the original sequential seeding
      out <- run_one_sim(N0, tk$N1, T0, T1, tk$gamma, tk$hhi, B = B_SE)
      out$config <- tk$config; out$N1 <- tk$N1
      out$gamma <- tk$gamma; out$hhi <- tk$hhi; out$sim <- tk$sim
      out
    })
    batch_df <- do.call(rbind, results)
    write.table(batch_df, OUT_CSV, sep = ",", row.names = FALSE,
                col.names = !file.exists(OUT_CSV), append = file.exists(OUT_CSV))
    el <- as.numeric(difftime(Sys.time(), t0, units = "hours"))
    eta <- el / b * (length(batches) - b)
    cat(sprintf("[%s] batch %d/%d done (%.2f h elapsed, ~%.1f h remaining)\n",
                format(Sys.time()), b, length(batches), el, eta),
        file = progress_file, append = TRUE)
  }
  cat(sprintf("[%s] parallel MC complete\n", format(Sys.time())),
      file = progress_file, append = TRUE)
}

mc_df <- read.csv(OUT_CSV)
message(sprintf("%s now has %d rows (%d complete configs x sims)",
                OUT_CSV, nrow(mc_df), nrow(unique(mc_df[, c("config", "sim")]))))
