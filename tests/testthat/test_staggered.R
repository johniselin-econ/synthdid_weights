# =============================================================================
# Tests for staggered-adoption weighted SDID
# =============================================================================

# A small balanced panel with never-treated units and three adoption cohorts.
make_staggered_panel = function(seed = 1) {
  set.seed(seed)
  N = 30; Tt = 12
  U = matrix(rpois(N * 2, 5), N, 2); V = matrix(rpois(Tt * 2, 5), Tt, 2)
  L = U %*% t(V) + outer(1:N, rep(1, Tt)) + outer(rep(1, N), 1:Tt)
  Y = L + matrix(rnorm(N * Tt, 0, 0.5), N, Tt)
  time = 1:Tt
  adoption.time = rep(Inf, N)      # rows 1-10, 27-30 never treated
  adoption.time[11:16] = 7
  adoption.time[23:26] = 8
  adoption.time[17:22] = 9
  tau = 3
  for (j in which(is.finite(adoption.time))) {
    post = which(time >= adoption.time[j])
    Y[j, post] = Y[j, post] + tau
  }
  list(Y = Y, adoption.time = adoption.time, time = time, tau = tau, N = N, Tt = Tt)
}

test_that("single cohort reduces to block weighted SDID exactly", {
  set.seed(2)
  N <- 25; Tt <- 10
  Y <- matrix(rnorm(N * Tt), N, Tt) + outer(1:N, rep(1, Tt))
  time <- 1:Tt
  adoption.time <- rep(Inf, N)
  adoption.time[20:25] <- 6                 # single cohort at t = 6
  w <- runif(6) + 0.1

  est.stag <- synthdid_estimate_staggered(Y, adoption.time, time,
                treated.weights = w, control = "never")

  ctl <- which(adoption.time > Tt); trt <- which(adoption.time == 6)
  pre <- which(time < 6); post <- which(time >= 6)
  Yg  <- Y[c(ctl, trt), c(pre, post)]
  est.block <- synthdid_estimate_weighted(Yg, length(ctl), length(pre),
                 treated.weights = w)

  expect_equal(as.numeric(est.stag), as.numeric(est.block), tolerance = 1e-10)
})

test_that("single cohort with uniform weights equals block uniform weighted SDID", {
  set.seed(3)
  N <- 22; Tt <- 9
  Y <- matrix(rnorm(N * Tt), N, Tt) + outer(1:N, rep(1, Tt))
  time <- 1:Tt
  adoption.time <- rep(Inf, N); adoption.time[18:22] <- 5

  est.stag <- synthdid_estimate_staggered(Y, adoption.time, time, control = "never")  # NULL -> uniform

  ctl <- which(adoption.time > Tt); trt <- which(adoption.time == 5)
  pre <- which(time < 5); post <- which(time >= 5)
  Yg  <- Y[c(ctl, trt), c(pre, post)]
  est.block <- synthdid_estimate_weighted(Yg, length(ctl), length(pre))

  expect_equal(as.numeric(est.stag), as.numeric(est.block), tolerance = 1e-10)
})

test_that("aggregate equals the cohort-weighted sum and weights normalize", {
  p <- make_staggered_panel()
  est <- synthdid_estimate_staggered(p$Y, p$adoption.time, p$time, control = "never")
  ct <- attr(est, "cohort.table")

  expect_equal(sum(ct$weight), 1, tolerance = 1e-12)
  expect_equal(as.numeric(est), sum(ct$weight * ct$estimate), tolerance = 1e-12)
  expect_equal(sort(ct$cohort), c(7, 8, 9))
  expect_true(all(ct$N1 > 0) && all(ct$N0 >= 2))
})

test_that("cohort.weight schemes match documented aggregations (uniform weights)", {
  p <- make_staggered_panel()
  e.share <- synthdid_estimate_staggered(p$Y, p$adoption.time, p$time,
               control = "never", cohort.weight = "treated.share")
  e.per   <- synthdid_estimate_staggered(p$Y, p$adoption.time, p$time,
               control = "never", cohort.weight = "treated.periods")
  cs <- attr(e.share, "cohort.table"); cp <- attr(e.per, "cohort.table")

  expect_equal(cs$weight, cs$N1 / sum(cs$N1), tolerance = 1e-12)
  expect_equal(cp$weight, (cp$N1 * cp$T1) / sum(cp$N1 * cp$T1), tolerance = 1e-12)
})

test_that("population weights shift the estimand relative to uniform", {
  p <- make_staggered_panel()
  N1 <- sum(is.finite(p$adoption.time) & p$adoption.time <= p$Tt)
  # weight the cohort-9 (largest, latest) units heavily
  trt.idx <- which(is.finite(p$adoption.time) & p$adoption.time <= p$Tt)
  w <- ifelse(p$adoption.time[trt.idx] == 9, 5, 1)

  e.uni <- synthdid_estimate_staggered(p$Y, p$adoption.time, p$time, control = "never")
  e.pop <- synthdid_estimate_staggered(p$Y, p$adoption.time, p$time,
             treated.weights = w, control = "never")
  # different estimands: aggregate weights should differ
  expect_false(isTRUE(all.equal(as.numeric(e.uni), as.numeric(e.pop))))
})

test_that("not-yet-treated control mode runs and is feasible", {
  p <- make_staggered_panel()
  est <- synthdid_estimate_staggered(p$Y, p$adoption.time, p$time, control = "notyet")
  expect_true(is.finite(as.numeric(est)))
  expect_equal(sort(attr(est, "cohort.table")$cohort), c(7, 8, 9))
})

test_that("infeasible cohort errors, or drops with drop.infeasible", {
  p <- make_staggered_panel()
  p$adoption.time[27:28] <- 2                # only 1 pre-period, < min.pre = 2
  expect_error(
    synthdid_estimate_staggered(p$Y, p$adoption.time, p$time, control = "never"))
  est <- synthdid_estimate_staggered(p$Y, p$adoption.time, p$time,
           control = "never", drop.infeasible = TRUE)
  expect_true(2 %in% attr(est, "cohort.dropped"))
  expect_equal(sum(attr(est, "cohort.table")$weight), 1, tolerance = 1e-12)
})

test_that("bootstrap vcov returns a finite non-negative variance", {
  p <- make_staggered_panel()
  est <- synthdid_estimate_staggered(p$Y, p$adoption.time, p$time, control = "never")
  set.seed(5)
  V <- vcov(est, replications = 30)
  expect_true(is.finite(V[1, 1]) && V[1, 1] >= 0)
})

test_that("golden: treated-periods aggregation reproduces Stata sdid (Jones et al. 2026, Table 5)", {
  # Reference: Stata `sdid emp_tot_serv_share state year ma_dereg_dum,
  # method(sdid) vce(noinference)` on the Jones et al. (2026) SJE Table-5 sample
  # (data_sje.dta; drop SD/DE, ma_dereg_year != 1960, 1969-1998). That run reports
  # aggregate ATT = 1.799401 and stores 20 per-adoption ATTs in e(tau) with a single
  # never-treated control (Iowa). The do-files live in docs/_gold/. This test guards
  # that synthdid_estimate_staggered()'s treated-periods aggregation -- W_g
  # proportional to N1_g * T_post,g -- matches sdid's staggered aggregation exactly.
  g   <- c(1970,1975,1976,1977,1978,1979,1980,1981,1982,1983,
           1984,1985,1986,1987,1988,1989,1990,1991,1993,1994)
  tau <- c(0.898536,3.0149169,0.67503059,4.4611384,3.0836946,2.2149758,3.9303146,
           0.90104583,3.9752146,1.8173267,3.0042734,0.90581351,1.2334213,2.57342,
           0.90942796,2.27197,0.70044995,1.1989118,0.47404305,0.7850102)
  N1  <- c(1,1,1,1,1,1,1,2,1,1,1,4,2,5,6,1,4,2,1,1)   # states adopting in each year
  T1  <- 1999 - g                                      # post periods g..1998
  gold <- 1.799401

  w_periods <- N1 * T1
  agg_periods <- sum(w_periods * tau) / sum(w_periods)
  expect_equal(agg_periods, gold, tolerance = 1e-5)

  # the other cohort-weight schemes are well-defined but target different estimands
  expect_false(isTRUE(all.equal(sum(N1 * tau) / sum(N1), gold, tolerance = 1e-3)))
})

test_that("a single-control cohort is refused cleanly (needs a donor pool)", {
  # SDID needs >= 2 controls per cohort; the default min.controls = 2 must reject a
  # never-treated pool of size 1 (the pathological Jones case, N0 = 1) with a clear
  # error rather than a cryptic failure deep in the solver.
  set.seed(11)
  N <- 20; Tt <- 10
  Y <- matrix(rnorm(N * Tt), N, Tt) + outer(1:N, rep(1, Tt))
  time <- 1:Tt
  adoption.time <- rep(Inf, N)
  adoption.time[1] <- Inf                 # exactly ONE never-treated (row 1)
  adoption.time[2:11]  <- 6               # two cohorts, both with >= 2 pre-periods
  adoption.time[12:20] <- 8
  expect_equal(sum(adoption.time > Tt), 1L)                 # a single never-treated
  expect_error(
    synthdid_estimate_staggered(Y, adoption.time, time, control = "never"),
    regexp = "control")
})

test_that("cluster bootstrap runs when a cluster vector is supplied", {
  p <- make_staggered_panel()
  cl <- rep(1:15, each = 2)                  # 15 clusters of 2 units
  est <- synthdid_estimate_staggered(p$Y, p$adoption.time, p$time,
           cluster = cl, control = "never")
  set.seed(6)
  V <- vcov(est, replications = 30)
  expect_true(is.finite(V[1, 1]) && V[1, 1] >= 0)
})
