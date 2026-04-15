# =============================================================================
# Tests for weighted SDID extensions
# =============================================================================

test_that("Proposition 1: uniform weights recover unweighted estimate", {
  setup = random.low.rank()
  N1 = nrow(setup$Y) - setup$N0
  T1 = ncol(setup$Y) - setup$T0

  tau.unweighted = synthdid_estimate(setup$Y, setup$N0, setup$T0)
  tau.weighted   = synthdid_estimate_weighted(setup$Y, setup$N0, setup$T0,
                     treated.weights = rep(1/N1, N1),
                     period.weights  = rep(1/T1, T1))

  expect_equal(c(tau.unweighted), c(tau.weighted), tolerance = 1e-10)
})

test_that("Proposition 1 holds for SC and DiD variants", {
  setup = random.low.rank()
  N1 = nrow(setup$Y) - setup$N0
  T1 = ncol(setup$Y) - setup$T0

  tau.sc     = sc_estimate(setup$Y, setup$N0, setup$T0)
  tau.sc.w   = sc_estimate_weighted(setup$Y, setup$N0, setup$T0,
                 treated.weights = rep(1/N1, N1))
  tau.did    = did_estimate(setup$Y, setup$N0, setup$T0)
  tau.did.w  = did_estimate_weighted(setup$Y, setup$N0, setup$T0,
                 treated.weights = rep(1/N1, N1))

  expect_equal(c(tau.sc),  c(tau.sc.w),  tolerance = 1e-10)
  expect_equal(c(tau.did), c(tau.did.w), tolerance = 1e-10)
})

test_that("weighted estimate with NULL weights equals default", {
  setup = random.low.rank()
  tau.null = synthdid_estimate_weighted(setup$Y, setup$N0, setup$T0,
               treated.weights = NULL, period.weights = NULL)
  tau.unweighted = synthdid_estimate(setup$Y, setup$N0, setup$T0)
  expect_equal(c(tau.null), c(tau.unweighted), tolerance = 1e-10)
})

test_that("column fixed-effect invariance holds for weighted estimator", {
  setup = random.low.rank()
  N1 = nrow(setup$Y) - setup$N0
  tw = runif(N1); tw = tw / sum(tw)

  bt = 2 * matrix(1:ncol(setup$Y), nrow(setup$Y), ncol(setup$Y), byrow = TRUE)
  tau.orig   = synthdid_estimate_weighted(setup$Y, setup$N0, setup$T0, treated.weights = tw)
  tau.shifted = synthdid_estimate_weighted(setup$Y + bt, setup$N0, setup$T0, treated.weights = tw)

  expect_equal(c(tau.orig), c(tau.shifted), tolerance = 1e-10)
})

test_that("row fixed-effect invariance holds for weighted DiD", {
  setup = random.low.rank()
  N1 = nrow(setup$Y) - setup$N0
  tw = runif(N1); tw = tw / sum(tw)

  ai = 3 * matrix(1:nrow(setup$Y), nrow(setup$Y), ncol(setup$Y))
  tau.orig   = did_estimate_weighted(setup$Y, setup$N0, setup$T0, treated.weights = tw)
  tau.shifted = did_estimate_weighted(setup$Y + ai, setup$N0, setup$T0, treated.weights = tw)

  expect_equal(c(tau.orig), c(tau.shifted), tolerance = 1e-10)
})

test_that("single treated unit with weight 1.0 works", {
  setup = random.low.rank()
  # Use just 1 treated unit
  Y1 = setup$Y[c(1:setup$N0, setup$N0 + 1), ]
  tau = synthdid_estimate_weighted(Y1, setup$N0, setup$T0, treated.weights = 1)
  expect_true(is.finite(c(tau)))
})

test_that("degenerate weights (all mass on one unit) works", {
  setup = random.low.rank()
  N1 = nrow(setup$Y) - setup$N0
  tw = rep(0, N1)
  tw[1] = 1
  tau = synthdid_estimate_weighted(setup$Y, setup$N0, setup$T0, treated.weights = tw)
  expect_true(is.finite(c(tau)))
})

test_that("weight normalization works", {
  setup = random.low.rank()
  N1 = nrow(setup$Y) - setup$N0
  # Unnormalized weights that sum to 100
  tw.raw = runif(N1) * 100

  tau = synthdid_estimate_weighted(setup$Y, setup$N0, setup$T0, treated.weights = tw.raw)
  expect_true(is.finite(c(tau)))
  # Verify stored weights sum to 1
  expect_equal(sum(attr(tau, 'treated.weights')), 1)
})

test_that("variance estimators return finite positive values", {
  setup = random.low.rank()
  N1 = nrow(setup$Y) - setup$N0
  tw = runif(N1); tw = tw / sum(tw)

  tau = synthdid_estimate_weighted(setup$Y, setup$N0, setup$T0, treated.weights = tw)

  se.boot = suppressWarnings(sqrt(vcov(tau, method = 'bootstrap', replications = 20)))
  se.jk   = sqrt(vcov(tau, method = 'jackknife'))
  se.plac = sqrt(vcov(tau, method = 'placebo', replications = 20))

  expect_true(is.finite(se.boot) && se.boot > 0)
  expect_true(is.finite(se.jk)   && se.jk > 0)
  expect_true(is.finite(se.plac) && se.plac > 0)
})

test_that("placebo weight options run without error", {
  setup = random.low.rank()
  N1 = nrow(setup$Y) - setup$N0
  tw = runif(N1); tw = tw / sum(tw)
  tau = synthdid_estimate_weighted(setup$Y, setup$N0, setup$T0, treated.weights = tw)

  se.uniform    = sqrt(vcov(tau, method = 'placebo', replications = 10, placebo.weights = 'uniform'))
  se.size_match = sqrt(vcov(tau, method = 'placebo', replications = 10, placebo.weights = 'size_match'))
  se.permute    = sqrt(vcov(tau, method = 'placebo', replications = 10, placebo.weights = 'permute'))

  expect_true(is.finite(se.uniform)    && se.uniform > 0)
  expect_true(is.finite(se.size_match) && se.size_match > 0)
  expect_true(is.finite(se.permute)    && se.permute > 0)
})

test_that("effective.sample.size changes eta.omega", {
  setup = random.low.rank()
  N1 = nrow(setup$Y) - setup$N0
  # Highly concentrated weights to see a difference
  tw = rep(0.001, N1); tw[1] = 1; tw = tw / sum(tw)

  tau.default = synthdid_estimate_weighted(setup$Y, setup$N0, setup$T0,
                  treated.weights = tw, effective.sample.size = FALSE)
  tau.ess     = synthdid_estimate_weighted(setup$Y, setup$N0, setup$T0,
                  treated.weights = tw, effective.sample.size = TRUE)

  # Both should produce valid estimates
  expect_true(is.finite(c(tau.default)))
  expect_true(is.finite(c(tau.ess)))
  # The tuning parameters should differ (stored in opts)
  expect_false(isTRUE(all.equal(
    attr(tau.default, 'opts')$zeta.omega,
    attr(tau.ess, 'opts')$zeta.omega
  )))
})

test_that("event study returns correct dimensions", {
  setup = random.low.rank()
  N1 = nrow(setup$Y) - setup$N0
  tw = runif(N1); tw = tw / sum(tw)

  tau = synthdid_estimate_weighted(setup$Y, setup$N0, setup$T0, treated.weights = tw)
  es = synthdid_event_study(tau)

  expect_s3_class(es, 'synthdid_event_study')
  expect_equal(nrow(es), ncol(setup$Y))
  expect_equal(ncol(es), 3)  # relative_time, estimate, time
})

test_that("event study with SEs doesn't error", {
  setup = random.low.rank()
  tau = synthdid_estimate(setup$Y, setup$N0, setup$T0)
  es = synthdid_event_study(tau, se.method = "bootstrap", replications = 10)

  expect_true("se" %in% names(es))
  expect_true("ci_lower" %in% names(es))
  expect_true("ci_upper" %in% names(es))
  expect_true(all(is.finite(es$se)))
})

test_that("print and summary work for weighted estimates", {
  setup = random.low.rank()
  N1 = nrow(setup$Y) - setup$N0
  tw = runif(N1); tw = tw / sum(tw)

  tau = synthdid_estimate_weighted(setup$Y, setup$N0, setup$T0, treated.weights = tw)

  expect_output(print(tau), "synthdid_estimate_weighted")
  s = summary(tau, fast = TRUE)
  expect_true("estimate" %in% names(s))
  expect_true("treated.weights" %in% names(s))
  expect_true("N1.effective" %in% names(s$dimensions))
})

test_that("synthdid_controls_weighted works for all weight types", {
  setup = random.low.rank()
  N1 = nrow(setup$Y) - setup$N0
  tw = runif(N1); tw = tw / sum(tw)

  tau = synthdid_estimate_weighted(setup$Y, setup$N0, setup$T0, treated.weights = tw)

  ctrl.omega   = synthdid_controls_weighted(tau, weight.type = 'omega')
  ctrl.lambda  = synthdid_controls_weighted(tau, weight.type = 'lambda')
  ctrl.treated = synthdid_controls_weighted(tau, weight.type = 'treated')
  ctrl.period  = synthdid_controls_weighted(tau, weight.type = 'period')

  expect_true(nrow(ctrl.omega) > 0)
  expect_true(nrow(ctrl.lambda) > 0)
  expect_true(nrow(ctrl.treated) > 0)
  expect_true(nrow(ctrl.period) > 0)
})

test_that("input validation rejects bad weights", {
  setup = random.low.rank()
  N1 = nrow(setup$Y) - setup$N0

  # Wrong length
  expect_error(synthdid_estimate_weighted(setup$Y, setup$N0, setup$T0,
    treated.weights = rep(1, N1 + 5)))
  # Negative weights
  expect_error(synthdid_estimate_weighted(setup$Y, setup$N0, setup$T0,
    treated.weights = c(-1, rep(1, N1 - 1))))
  # All-zero weights
  expect_error(synthdid_estimate_weighted(setup$Y, setup$N0, setup$T0,
    treated.weights = rep(0, N1)))
})

test_that("synthdid_placebo_weighted runs and returns valid estimate", {
  setup = random.low.rank()
  N1 = nrow(setup$Y) - setup$N0
  tw = runif(N1); tw = tw / sum(tw)

  tau = synthdid_estimate_weighted(setup$Y, setup$N0, setup$T0, treated.weights = tw)
  placebo = synthdid_placebo_weighted(tau)

  expect_true(is.finite(c(placebo)))
  expect_true(inherits(placebo, 'synthdid_estimate_weighted'))
  # Placebo should use pre-treatment data only: fewer columns
  placebo.setup = attr(placebo, 'setup')
  expect_equal(ncol(placebo.setup$Y), setup$T0)
  # Placebo should pass through opts (zeta.omega should match)
  expect_equal(attr(placebo, 'opts')$zeta.omega, attr(tau, 'opts')$zeta.omega)
})

test_that("covariate path works for weighted estimator", {
  setup = random.low.rank()
  N1 = nrow(setup$Y) - setup$N0
  tw = runif(N1); tw = tw / sum(tw)

  # Create a covariate array: random noise
  X = array(rnorm(prod(dim(setup$Y))), dim = c(dim(setup$Y), 1))

  tau.nocov = synthdid_estimate_weighted(setup$Y, setup$N0, setup$T0, treated.weights = tw)
  tau.cov   = synthdid_estimate_weighted(setup$Y, setup$N0, setup$T0, treated.weights = tw, X = X)

  expect_true(is.finite(c(tau.cov)))
  # Covariate adjustment should change the estimate
  expect_false(isTRUE(all.equal(c(tau.nocov), c(tau.cov))))
  # Beta should be non-null
  expect_true(!is.null(attr(tau.cov, 'weights')$beta))
  expect_true(length(attr(tau.cov, 'weights')$beta) == 1)
})
