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

test_that("cluster bootstrap runs when a cluster vector is supplied", {
  p <- make_staggered_panel()
  cl <- rep(1:15, each = 2)                  # 15 clusters of 2 units
  est <- synthdid_estimate_staggered(p$Y, p$adoption.time, p$time,
           cluster = cl, control = "never")
  set.seed(6)
  V <- vcov(est, replications = 30)
  expect_true(is.finite(V[1, 1]) && V[1, 1] >= 0)
})
