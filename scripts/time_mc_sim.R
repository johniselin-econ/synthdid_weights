suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library(mvtnorm)
})

source_sim_fun <- function() {
  generate_sizes <<- function(N1, target_hhi) {
    if (target_hhi <= 1/N1 + 0.005) return(rep(1, N1))
    saved_seed <- .Random.seed
    obj <- function(log_sd) {
      set.seed(999)
      s <- rlnorm(N1 * 100, meanlog = 0, sdlog = exp(log_sd))
      hhis <- sapply(split(s, rep(1:100, each = N1)),
                     function(x) sum((x/sum(x))^2))
      mean(hhis) - target_hhi
    }
    opt <- tryCatch(uniroot(obj, c(-3, 3), tol = 0.001),
                    error = function(e) list(root = 0.5))
    .Random.seed <<- saved_seed
    rlnorm(N1, meanlog = 0, sdlog = exp(opt$root))
  }

  run_one_sim <<- function(N0, N1, T0, T1, gamma, target_hhi,
                           sigma = 0.5, tau_bar = 1, B = 200) {
    N <- N0 + N1; T_ <- T0 + T1; rank <- 2; rho <- 0.7
    U <- matrix(rpois(rank * N, sqrt(sample(1:N)) / sqrt(N)), N, rank)
    V <- matrix(rpois(rank * T_, sqrt(1:T_) / sqrt(T_)), T_, rank)
    alpha <- outer(10 * sample(1:N) / N, rep(1, T_))
    beta  <- outer(rep(1, N), 10 * (1:T_) / T_)
    L <- U %*% t(V) + alpha + beta
    var_mat <- outer(1:T_, 1:T_, FUN = function(x, y) rho^(abs(x - y)))
    error <- mvtnorm::rmvnorm(N, sigma = var_mat, method = "chol")
    W <- matrix(0, N, T_); W[(N0+1):N, (T0+1):T_] <- 1
    sizes <- generate_sizes(N1, target_hhi)
    tw <- sizes / sum(sizes)
    tau_vec <- tau_bar + gamma * (sizes / mean(sizes) - 1)
    tau_mat <- matrix(0, N, T_)
    for (j in 1:N1) tau_mat[N0 + j, (T0+1):T_] <- tau_vec[j]
    Y <- L + tau_mat + sigma * error
    rownames(Y) <- 1:N; colnames(Y) <- 1:T_
    est_wt <- tryCatch(
      synthdid_estimate_weighted(Y, N0, T0, treated.weights = tw),
      error = function(e) NA)
    if (is.na(est_wt[1])) return(NA_real_)
    se_boot <- sqrt(vcov(est_wt, method = "bootstrap", replications = B))
    se_jk   <- sqrt(vcov(est_wt, method = "jackknife"))
    se_plac_u <- sqrt(vcov(est_wt, method = "placebo",
                           replications = B, placebo.weights = "uniform"))
    se_plac_p <- sqrt(vcov(est_wt, method = "placebo",
                           replications = B, placebo.weights = "permute"))
    se_plac_s <- sqrt(vcov(est_wt, method = "placebo",
                           replications = B, placebo.weights = "size_match"))
    c(se_boot, se_jk, se_plac_u, se_plac_p, se_plac_s)
  }
}

source_sim_fun()

N_TIMING <- 5
N0 <- 100; T0 <- 40; T1 <- 10
cat("Timing one sim replicate (N1=10, gamma=0.5, hhi=0.15, B=200) x",
    N_TIMING, "reps...\n")

t0 <- proc.time()
for (s in 1:N_TIMING) {
  set.seed(20240101 + s)
  run_one_sim(N0, 10, T0, T1, 0.5, 0.15, B = 200)
}
elapsed <- (proc.time() - t0)["elapsed"]
per_rep <- elapsed / N_TIMING
cat(sprintf("Per replicate: %.2f sec\n", per_rep))
cat(sprintf("Projected full run (18 cells x 500 reps x 3 placebo variants integrated): %.1f min\n",
            per_rep * 18 * 500 / 60))
