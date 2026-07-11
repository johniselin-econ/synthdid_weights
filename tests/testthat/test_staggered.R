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

test_that("aggregate recovers a known homogeneous effect (both control modes)", {
  # make_staggered_panel plants tau = 3 on every treated cell; the aggregate
  # should recover it. This is the guard a window-construction off-by-one
  # would trip that the structural identity tests cannot. Averaged over five
  # fixed seeds because the notyet mode's capped windows leave cohorts 7 and 8
  # with a single post period each (single-draw sd ~0.7 for the aggregate;
  # seeds 1-5 means: never 2.94, notyet 2.58) -- itself a live demonstration
  # of the documented notyet window-truncation caveat.
  for (ctrl in c("never", "notyet")) {
    ests <- sapply(1:5, function(s) {
      p <- make_staggered_panel(seed = s)
      as.numeric(synthdid_estimate_staggered(p$Y, p$adoption.time, p$time, control = ctrl))
    })
    expect_lt(abs(mean(ests) - 3), 0.75)
  }
})

test_that("population weights shift the estimand under heterogeneous effects", {
  # Give cohort 9 a LARGER effect (tau = 6 vs 3 elsewhere). Uniform
  # treated.share estimand = (6*3 + 4*3 + 6*6)/16 = 4.125; upweighting
  # cohort 9 by 5x gives (6*3 + 4*3 + 30*6)/40 = 5.25. The weighted estimate
  # must move toward cohort 9's effect by a gap large relative to noise.
  p <- make_staggered_panel()
  extra <- 3
  for (j in which(p$adoption.time == 9)) {
    post <- which(p$time >= 9)
    p$Y[j, post] <- p$Y[j, post] + extra
  }
  trt.idx <- which(is.finite(p$adoption.time) & p$adoption.time <= p$Tt)
  w <- ifelse(p$adoption.time[trt.idx] == 9, 5, 1)

  e.uni <- synthdid_estimate_staggered(p$Y, p$adoption.time, p$time, control = "never")
  e.pop <- synthdid_estimate_staggered(p$Y, p$adoption.time, p$time,
             treated.weights = w, control = "never")
  expect_gt(as.numeric(e.pop) - as.numeric(e.uni), 0.5)
})

test_that("interleaved treated/control rows match the reordered block fit", {
  # Treated units scattered among control rows: exercises the per-cohort
  # reordering (rows.idx = c(ctl.idx, trt.idx)) that the contiguous-block
  # tests never touch. Weights are ordered by row among treated units.
  set.seed(7)
  N <- 24; Tt <- 10; time <- 1:Tt
  Y <- matrix(rnorm(N * Tt), N, Tt) + outer(1:N, rep(1, Tt))
  adoption.time <- rep(Inf, N)
  trt <- c(3, 7, 11, 15, 19, 23)
  adoption.time[trt] <- 6
  for (j in trt) Y[j, time >= 6] <- Y[j, time >= 6] + 2
  w <- runif(length(trt)) + 0.1

  est.stag <- synthdid_estimate_staggered(Y, adoption.time, time,
                treated.weights = w, control = "never")

  ctl <- which(adoption.time > Tt)
  est.block <- synthdid_estimate_weighted(Y[c(ctl, trt), ], length(ctl), sum(time < 6),
                 treated.weights = w)
  expect_equal(as.numeric(est.stag), as.numeric(est.block), tolerance = 1e-10)
})

test_that("NA adoption.time is refused, not silently dropped", {
  p <- make_staggered_panel()
  p$adoption.time[1] <- NA
  expect_error(
    synthdid_estimate_staggered(p$Y, p$adoption.time, p$time, control = "never"),
    regexp = "anyNA")
})

test_that("panel-shaped arguments are rejected via ...", {
  p <- make_staggered_panel()
  X <- array(rnorm(p$N * p$Tt), dim = c(p$N, p$Tt, 1))
  expect_error(
    synthdid_estimate_staggered(p$Y, p$adoption.time, p$time, control = "never", X = X),
    regexp = "unsupported in the staggered estimator")
})

test_that("an all-zero-weight cohort is infeasible (share <= 0), like stratified", {
  p <- make_staggered_panel()
  trt.idx <- which(is.finite(p$adoption.time) & p$adoption.time <= p$Tt)
  w <- ifelse(p$adoption.time[trt.idx] == 8, 0, 1)   # cohort 8 carries zero weight
  expect_error(
    synthdid_estimate_staggered(p$Y, p$adoption.time, p$time,
                                treated.weights = w, control = "never"),
    regexp = "share > 0")
  est <- synthdid_estimate_staggered(p$Y, p$adoption.time, p$time,
           treated.weights = w, control = "never", drop.infeasible = TRUE)
  expect_true(8 %in% attr(est, "cohort.dropped"))
  expect_equal(sum(attr(est, "cohort.table")$weight), 1, tolerance = 1e-12)
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

test_that("golden: the treated-periods aggregation FORMULA matches Stata sdid (Jones et al. 2026, Table 5)", {
  # SCOPE: this test validates the aggregation FORMULA (W_g proportional to
  # N1_g * T_post,g) against the Stata reference, NOT the estimator end-to-end.
  # The formula is applied inline to Stata's own per-adoption ATTs; the guard
  # that synthdid_estimate_staggered() implements this same formula in code is
  # the "cohort.weight schemes match documented aggregations" test above.
  # End-to-end comparison on this panel is impossible: it has a single
  # never-treated control (Iowa), which Stata sdid runs as a degenerate
  # one-donor synthetic control and our estimator refuses (min.controls = 2;
  # see the single-control test below).
  #
  # Reference: Stata `sdid emp_tot_serv_share state year ma_dereg_dum,
  # method(sdid) vce(noinference)` on the Jones et al. (2026) SJE Table-5 sample
  # (data_sje.dta; drop SD/DE, ma_dereg_year != 1960, 1969-1998): aggregate
  # ATT = 1.799401, 20 per-adoption ATTs in e(tau).
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
  adoption.time <- rep(Inf, N)            # row 1 stays Inf: the ONE never-treated
  adoption.time[2:11]  <- 6               # two cohorts, both with >= 2 pre-periods
  adoption.time[12:20] <- 8
  expect_equal(sum(adoption.time > Tt), 1L)                 # a single never-treated
  expect_error(
    synthdid_estimate_staggered(Y, adoption.time, time, control = "never"),
    regexp = "need >= .* controls")
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
