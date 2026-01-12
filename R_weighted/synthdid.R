#' A function mapping a numeric vector to a (presumably sparser) numeric vector of the same shape to
#' be passed onto synthdid_estimate.
#' @param v a vector
sparsify_function = function(v) { v[v <= max(v)/4] = 0; v/sum(v) }

#' Computes the synthetic diff-in-diff estimate for an average treatment effect on a treated block.
#'
#' See 'Synthetic Difference in Differences' by Arkhangelsky et al. This implements Algorithm 1.
#' @param Y the observation matrix.
#' @param N0 the number of control units (N_co in the paper). Rows 1-N0 of Y correspond to the control units.
#' @param T0 the number of pre-treatment time steps (T_pre in the paper). Columns 1-T0 of Y correspond to pre-treatment time steps.
#' @param X an optional 3-D array of time-varying covariates. Shape should be N X T X C for C covariates.
#' @param noise.level, an estimate of the noise standard deviation sigma. Defaults to the standard deviation of first differences of Y.
#' @param eta.omega  determines the tuning parameter zeta.omega = eta.omega * noise.level. Defaults to the value (N_tr T_post)^(1/4).
#' @param eta.lambda analogous for lambda.  Defaults to an 'infinitesimal' value 1e-6.
#' @param zeta.omega if passed, overrides the default zeta.omega = eta.omega * noise.level. Deprecated.
#' @param zeta.lambda analogous for lambda.
#' @param omega.intercept Binary. Use an intercept when estimating omega.
#' @param lambda.intercept Binary. Use an intercept when estimating lambda.
#' @param weights a list with fields lambda and omega. If non-null weights$lambda is passed,
#'        we use them instead of estimating lambda weights. Same for weights$omega.
#' @param update.omega If true, solve for omega using the passed value of weights$omega only as an initialization.
#'        If false, use it exactly as passed. Defaults to false if a non-null value of weights$omega is passed.
#' @param update.lambda  Analogous.
#' @param min.decrease Tunes a stopping criterion for our weight estimator. Stop after an iteration results in a decrease
#' 		        in penalized MSE smaller than min.decrease^2.
#' @param max.iter A fallback stopping criterion for our weight estimator. Stop after this number of iterations.
#' @param sparsify A function mapping a numeric vector to a (presumably sparser) numeric vector of the same shape, which must sum to one.
#'                  If not null, we try to estimate sparse weights via a second round of Frank-Wolfe optimization
#'                  initialized at sparsify( the solution to the first round ).
#' @param max.iter.pre.sparsify Analogous to max.iter, but for the pre-sparsification first-round of optimization.
#'     		                Not used if sparsify=NULL.
#' @return An average treatment effect estimate with 'weights' and 'setup' attached as attributes.
#'         'weights' contains the estimated weights lambda and omega and corresponding intercepts,
#'         as well as regression coefficients beta if X is passed.
#'         'setup' is a list describing the problem passed in: Y, N0, T0, X.
#' @export synthdid_estimate
synthdid_estimate <- function(Y, N0, T0, X = array(dim = c(dim(Y), 0)),
                              noise.level = sd(apply(Y[1:N0,1:T0], 1, diff)),
                              eta.omega = ((nrow(Y)-N0)*(ncol(Y)-T0))^(1/4), eta.lambda = 1e-6,
                              zeta.omega  = eta.omega  * noise.level,  zeta.lambda = eta.lambda * noise.level,
                              omega.intercept = TRUE, lambda.intercept = TRUE,
                              weights = list(omega = NULL, lambda = NULL),
                              update.omega = is.null(weights$omega), update.lambda = is.null(weights$lambda),
                              min.decrease = 1e-5 * noise.level, max.iter = 1e4,
			      sparsify = sparsify_function,
			      max.iter.pre.sparsify = 100) {
  stopifnot(nrow(Y) > N0, ncol(Y) > T0, length(dim(X)) %in% c(2, 3), dim(X)[1:2] == dim(Y), is.list(weights),
    is.null(weights$lambda) || length(weights$lambda) == T0, is.null(weights$omega) || length(weights$omega) == N0,
    !is.null(weights$lambda) || update.lambda, !is.null(weights$omega) || update.omega)
  if (length(dim(X)) == 2) { dim(X) = c(dim(X), 1) }
  if (is.null(sparsify)) { max.iter.pre.sparsify = max.iter }
  N1 = nrow(Y) - N0
  T1 = ncol(Y) - T0

  if (dim(X)[3] == 0) {
    weights$vals = NULL
    weights$lambda.vals = NULL
    weights$omega.vals = NULL
    if (update.lambda) {
      Yc = collapsed.form(Y, N0, T0)
      lambda.opt = sc.weight.fw(Yc[1:N0, ], zeta = zeta.lambda, intercept = lambda.intercept, lambda=weights$lambda,
				min.decrease = min.decrease, max.iter = max.iter.pre.sparsify)
      if(!is.null(sparsify)) {
	lambda.opt = sc.weight.fw(Yc[1:N0, ], zeta = zeta.lambda, intercept = lambda.intercept, lambda=sparsify(lambda.opt$lambda),
				  min.decrease = min.decrease, max.iter = max.iter)
      }
      weights$lambda = lambda.opt$lambda
      weights$lambda.vals = lambda.opt$vals
      weights$vals = lambda.opt$vals
    }
    if (update.omega) {
      Yc = collapsed.form(Y, N0, T0)
      omega.opt = sc.weight.fw(t(Yc[, 1:T0]), zeta = zeta.omega, intercept = omega.intercept, lambda=weights$omega,
			       min.decrease = min.decrease, max.iter = max.iter.pre.sparsify)
      if(!is.null(sparsify)) {
	omega.opt = sc.weight.fw(t(Yc[, 1:T0]), zeta = zeta.omega, intercept = omega.intercept, lambda=sparsify(omega.opt$lambda),
			         min.decrease = min.decrease, max.iter = max.iter)
      }
      weights$omega = omega.opt$lambda
      weights$omega.vals = omega.opt$vals
      if (is.null(weights$vals)) { weights$vals = omega.opt$vals }
      else { weights$vals = pairwise.sum.decreasing(weights$vals, omega.opt$vals) }
    }
  } else {
    Yc = collapsed.form(Y, N0, T0)
    Xc = apply(X, 3, function(Xi) { collapsed.form(Xi, N0, T0) })
    dim(Xc) = c(dim(Yc), dim(X)[3])
    weights = sc.weight.fw.covariates(Yc, Xc, zeta.lambda = zeta.lambda, zeta.omega = zeta.omega,
      lambda.intercept = lambda.intercept, omega.intercept = omega.intercept,
      min.decrease = min.decrease, max.iter = max.iter,
      lambda = weights$lambda, omega = weights$omega, update.lambda = update.lambda, update.omega = update.omega)
  }

  X.beta = contract3(X, weights$beta)
  estimate = t(c(-weights$omega, rep(1 / N1, N1))) %*% (Y - X.beta) %*% c(-weights$lambda, rep(1 / T1, T1))

  class(estimate) = 'synthdid_estimate'
  attr(estimate, 'estimator') = "synthdid_estimate"
  attr(estimate, 'weights') = weights
  attr(estimate, 'setup') = list(Y = Y, X = X, N0 = N0, T0 = T0)
  attr(estimate, 'opts') = list(zeta.omega = zeta.omega, zeta.lambda = zeta.lambda,
                                omega.intercept = omega.intercept, lambda.intercept = lambda.intercept,
                                update.omega = update.omega, update.lambda = update.lambda,
                                min.decrease = min.decrease, max.iter=max.iter)
  return(estimate)
}

# =============================================================================
# WEIGHTED VERSION: synthdid_estimate_weighted
# =============================================================================

#' Computes the synthetic diff-in-diff estimate with user-specified treated unit weights.
#'
#' This is a weighted extension of the SDID estimator that allows for population-weighted
#' or other user-specified weights when averaging across multiple treated units.
#' See 'Synthetic Difference in Differences' by Arkhangelsky et al. for the base algorithm.
#'
#' @param Y the observation matrix.
#' @param N0 the number of control units (N_co in the paper). Rows 1-N0 of Y correspond to the control units.
#' @param T0 the number of pre-treatment time steps (T_pre in the paper). Columns 1-T0 of Y correspond to pre-treatment time steps.
#' @param treated.weights A vector of weights for treated units (length N1 = nrow(Y) - N0).
#'        If NULL, uses uniform weights (1/N1 for each treated unit), equivalent to the standard estimator.
#'        Weights will be normalized to sum to 1.
#' @param period.weights A vector of weights for post-treatment periods (length T1 = ncol(Y) - T0).
#'        If NULL, uses uniform weights (1/T1 for each period).
#'        Weights will be normalized to sum to 1.
#' @param X an optional 3-D array of time-varying covariates. Shape should be N X T X C for C covariates.
#' @param noise.level an estimate of the noise standard deviation sigma. Defaults to the standard deviation of first differences of Y.
#' @param eta.omega  determines the tuning parameter zeta.omega = eta.omega * noise.level. Defaults to the value (N_tr T_post)^(1/4).
#' @param eta.lambda analogous for lambda. Defaults to an 'infinitesimal' value 1e-6.
#' @param zeta.omega if passed, overrides the default zeta.omega = eta.omega * noise.level. Deprecated.
#' @param zeta.lambda analogous for lambda.
#' @param omega.intercept Binary. Use an intercept when estimating omega.
#' @param lambda.intercept Binary. Use an intercept when estimating lambda.
#' @param weights a list with fields lambda and omega. If non-null weights$lambda is passed,
#'        we use them instead of estimating lambda weights. Same for weights$omega.
#' @param update.omega If true, solve for omega using the passed value of weights$omega only as an initialization.
#'        If false, use it exactly as passed. Defaults to false if a non-null value of weights$omega is passed.
#' @param update.lambda Analogous.
#' @param min.decrease Tunes a stopping criterion for our weight estimator. Stop after an iteration results in a decrease
#'        in penalized MSE smaller than min.decrease^2.
#' @param max.iter A fallback stopping criterion for our weight estimator. Stop after this number of iterations.
#' @param sparsify A function mapping a numeric vector to a (presumably sparser) numeric vector of the same shape, which must sum to one.
#'        If not null, we try to estimate sparse weights via a second round of Frank-Wolfe optimization
#'        initialized at sparsify( the solution to the first round ).
#' @param max.iter.pre.sparsify Analogous to max.iter, but for the pre-sparsification first-round of optimization.
#'        Not used if sparsify=NULL.
#' @return An average treatment effect estimate with 'weights', 'setup', and 'treated.weights' attached as attributes.
#'         'weights' contains the estimated weights lambda and omega and corresponding intercepts,
#'         as well as regression coefficients beta if X is passed.
#'         'setup' is a list describing the problem passed in: Y, N0, T0, X.
#'         'treated.weights' contains the normalized weights used for treated units.
#'         'period.weights' contains the normalized weights used for post-treatment periods.
#' @export synthdid_estimate_weighted
synthdid_estimate_weighted <- function(Y, N0, T0,
                                       treated.weights = NULL,
                                       period.weights = NULL,
                                       X = NULL,
                                       noise.level = NULL,
                                       eta.omega = NULL, eta.lambda = 1e-6,
                                       zeta.omega  = NULL,  zeta.lambda = NULL,
                                       omega.intercept = TRUE, lambda.intercept = TRUE,
                                       weights = list(omega = NULL, lambda = NULL),
                                       update.omega = is.null(weights$omega), update.lambda = is.null(weights$lambda),
                                       min.decrease = NULL, max.iter = 1e4,
                                       sparsify = sparsify_function,
                                       max.iter.pre.sparsify = 100) {

  # Handle default X
 if (is.null(X)) {
    X = array(dim = c(dim(Y), 0))
  }

  # Handle other defaults that depend on data
  if (is.null(noise.level)) {
    noise.level = sd(apply(Y[1:N0, 1:T0], 1, diff))
  }
  if (is.null(eta.omega)) {
    eta.omega = ((nrow(Y) - N0) * (ncol(Y) - T0))^(1/4)
  }
  if (is.null(zeta.omega)) {
    zeta.omega = eta.omega * noise.level
  }
  if (is.null(zeta.lambda)) {
    zeta.lambda = eta.lambda * noise.level
  }
  if (is.null(min.decrease)) {
    min.decrease = 1e-5 * noise.level
  }

  # Input validation
  stopifnot(nrow(Y) > N0, ncol(Y) > T0, length(dim(X)) %in% c(2, 3), dim(X)[1:2] == dim(Y), is.list(weights),
    is.null(weights$lambda) || length(weights$lambda) == T0, is.null(weights$omega) || length(weights$omega) == N0,
    !is.null(weights$lambda) || update.lambda, !is.null(weights$omega) || update.omega)

  if (length(dim(X)) == 2) { dim(X) = c(dim(X), 1) }
  if (is.null(sparsify)) { max.iter.pre.sparsify = max.iter }

  N1 = nrow(Y) - N0
  T1 = ncol(Y) - T0

  # Process treated unit weights
  if (is.null(treated.weights)) {
    treated.weights = rep(1 / N1, N1)
  } else {
    stopifnot(length(treated.weights) == N1, all(treated.weights >= 0), sum(treated.weights) > 0)
    treated.weights = treated.weights / sum(treated.weights)  # Normalize to sum to 1
  }

  # Process period weights
  if (is.null(period.weights)) {
    period.weights = rep(1 / T1, T1)
  } else {
    stopifnot(length(period.weights) == T1, all(period.weights >= 0), sum(period.weights) > 0)
    period.weights = period.weights / sum(period.weights)  # Normalize to sum to 1
  }

  # Estimate control weights (omega) and time weights (lambda)
  if (dim(X)[3] == 0) {
    weights$vals = NULL
    weights$lambda.vals = NULL
    weights$omega.vals = NULL

    if (update.lambda) {
      # Use weighted collapsed form for treated units
      Yc = collapsed.form_weighted(Y, N0, T0, treated.weights, period.weights)
      lambda.opt = sc.weight.fw(Yc[1:N0, ], zeta = zeta.lambda, intercept = lambda.intercept, lambda=weights$lambda,
                                min.decrease = min.decrease, max.iter = max.iter.pre.sparsify)
      if(!is.null(sparsify)) {
        lambda.opt = sc.weight.fw(Yc[1:N0, ], zeta = zeta.lambda, intercept = lambda.intercept, lambda=sparsify(lambda.opt$lambda),
                                  min.decrease = min.decrease, max.iter = max.iter)
      }
      weights$lambda = lambda.opt$lambda
      weights$lambda.vals = lambda.opt$vals
      weights$vals = lambda.opt$vals
    }

    if (update.omega) {
      Yc = collapsed.form_weighted(Y, N0, T0, treated.weights, period.weights)
      omega.opt = sc.weight.fw(t(Yc[, 1:T0]), zeta = zeta.omega, intercept = omega.intercept, lambda=weights$omega,
                               min.decrease = min.decrease, max.iter = max.iter.pre.sparsify)
      if(!is.null(sparsify)) {
        omega.opt = sc.weight.fw(t(Yc[, 1:T0]), zeta = zeta.omega, intercept = omega.intercept, lambda=sparsify(omega.opt$lambda),
                                 min.decrease = min.decrease, max.iter = max.iter)
      }
      weights$omega = omega.opt$lambda
      weights$omega.vals = omega.opt$vals
      if (is.null(weights$vals)) { weights$vals = omega.opt$vals }
      else { weights$vals = pairwise.sum.decreasing(weights$vals, omega.opt$vals) }
    }
  } else {
    # Handle covariates case
    Yc = collapsed.form_weighted(Y, N0, T0, treated.weights, period.weights)
    Xc = apply(X, 3, function(Xi) { collapsed.form_weighted(Xi, N0, T0, treated.weights, period.weights) })
    dim(Xc) = c(dim(Yc), dim(X)[3])
    weights = sc.weight.fw.covariates(Yc, Xc, zeta.lambda = zeta.lambda, zeta.omega = zeta.omega,
      lambda.intercept = lambda.intercept, omega.intercept = omega.intercept,
      min.decrease = min.decrease, max.iter = max.iter,
      lambda = weights$lambda, omega = weights$omega, update.lambda = update.lambda, update.omega = update.omega)
  }

  # Compute the weighted treatment effect estimate
  X.beta = contract3(X, weights$beta)
  estimate = t(c(-weights$omega, treated.weights)) %*% (Y - X.beta) %*% c(-weights$lambda, period.weights)

  class(estimate) = 'synthdid_estimate_weighted'
  attr(estimate, 'estimator') = "synthdid_estimate_weighted"
  attr(estimate, 'weights') = weights
  attr(estimate, 'treated.weights') = treated.weights
  attr(estimate, 'period.weights') = period.weights
  attr(estimate, 'setup') = list(Y = Y, X = X, N0 = N0, T0 = T0)
  attr(estimate, 'opts') = list(zeta.omega = zeta.omega, zeta.lambda = zeta.lambda,
                                omega.intercept = omega.intercept, lambda.intercept = lambda.intercept,
                                update.omega = update.omega, update.lambda = update.lambda,
                                min.decrease = min.decrease, max.iter=max.iter)
  return(estimate)
}


#' synthdid_estimate for synthetic control estimates.
#' Takes all the same parameters, but by default, passes options to use the synthetic control estimator
#' By default, this uses only 'infinitesimal' ridge regularization when estimating the weights.
#' @param Y the observation matrix.
#' @param N0 the number of control units. Rows 1-N0 of Y correspond to the control units.
#' @param T0 the number of pre-treatment time steps. Columns 1-T0 of Y correspond to pre-treatment time steps.
#' @param eta.omega determines the level of ridge regularization, zeta.omega = eta.omega * noise.level, as in synthdid_estimate.
#' @param ... additional options for synthdid_estimate
#' @return an object like that returned by synthdid_estimate
#' @export sc_estimate
sc_estimate = function(Y, N0, T0, eta.omega = 1e-6, ...) {
  estimate = synthdid_estimate(Y, N0, T0, eta.omega = eta.omega,
			       weights = list(lambda = rep(0, T0)), omega.intercept = FALSE, ...)
  attr(estimate, 'estimator') = "sc_estimate"
  estimate
}

#' Weighted synthetic control estimate.
#' Takes all the same parameters as sc_estimate, plus treated.weights and period.weights.
#' @param Y the observation matrix.
#' @param N0 the number of control units. Rows 1-N0 of Y correspond to the control units.
#' @param T0 the number of pre-treatment time steps. Columns 1-T0 of Y correspond to pre-treatment time steps.
#' @param treated.weights A vector of weights for treated units (length N1).
#' @param period.weights A vector of weights for post-treatment periods (length T1).
#' @param eta.omega determines the level of ridge regularization.
#' @param ... additional options for synthdid_estimate_weighted
#' @return an object like that returned by synthdid_estimate_weighted
#' @export sc_estimate_weighted
sc_estimate_weighted = function(Y, N0, T0, treated.weights = NULL, period.weights = NULL, eta.omega = 1e-6, ...) {
  estimate = synthdid_estimate_weighted(Y, N0, T0,
                                        treated.weights = treated.weights,
                                        period.weights = period.weights,
                                        eta.omega = eta.omega,
                                        weights = list(lambda = rep(0, T0)), omega.intercept = FALSE, ...)
  attr(estimate, 'estimator') = "sc_estimate_weighted"
  estimate
}

#' synthdid_estimate for diff-in-diff estimates.
#' Takes all the same parameters, but by default, passes options to use the diff-in-diff estimator
#' @param Y the observation matrix.
#' @param N0 the number of control units. Rows 1-N0 of Y correspond to the control units.
#' @param T0 the number of pre-treatment time steps. Columns 1-T0 of Y correspond to pre-treatment time steps.
#' @param ... additional  options for synthdid_estimate
#' @return an object like that returned by synthdid_estimate
#' @export did_estimate
did_estimate = function(Y, N0, T0, ...) {
  estimate = synthdid_estimate(Y, N0, T0, weights = list(lambda = rep(1 / T0, T0), omega = rep(1 / N0, N0)), ...)
  attr(estimate, 'estimator') = "did_estimate"
  estimate
}

#' Weighted diff-in-diff estimate.
#' Takes all the same parameters as did_estimate, plus treated.weights and period.weights.
#' @param Y the observation matrix.
#' @param N0 the number of control units. Rows 1-N0 of Y correspond to the control units.
#' @param T0 the number of pre-treatment time steps. Columns 1-T0 of Y correspond to pre-treatment time steps.
#' @param treated.weights A vector of weights for treated units (length N1).
#' @param period.weights A vector of weights for post-treatment periods (length T1).
#' @param ... additional options for synthdid_estimate_weighted
#' @return an object like that returned by synthdid_estimate_weighted
#' @export did_estimate_weighted
did_estimate_weighted = function(Y, N0, T0, treated.weights = NULL, period.weights = NULL, ...) {
  estimate = synthdid_estimate_weighted(Y, N0, T0,
                                        treated.weights = treated.weights,
                                        period.weights = period.weights,
                                        weights = list(lambda = rep(1 / T0, T0), omega = rep(1 / N0, N0)), ...)
  attr(estimate, 'estimator') = "did_estimate_weighted"
  estimate
}

#' Computes a placebo variant of our estimator using pre-treatment data only
#' @param estimate, as output by synthdid_estimate
#' @param treated.fraction, the fraction of pre-treatment data to use as a placebo treatment period
#'        Defaults to NULL, which indicates that it should be the fraction of post-treatment to pre-treatment data
#' @export synthdid_placebo
synthdid_placebo = function(estimate, treated.fraction = NULL) {
  setup = attr(estimate, 'setup')
  opts = attr(estimate, 'opts')
  weights = attr(estimate, 'weights')
  X.beta = contract3(setup$X, weights$beta)
  estimator = attr(estimate, 'estimator')

  if (is.null(treated.fraction)) { treated.fraction = 1 - setup$T0 / ncol(setup$Y) }
  placebo.T0 = floor(setup$T0 * (1 - treated.fraction))

  do.call(estimator, c(list(Y=setup$Y[, 1:setup$T0], N0=setup$N0, T0=placebo.T0, X=setup$X[, 1:setup$T0,]), opts))
}

#' Computes a placebo variant for the weighted estimator using pre-treatment data only
#' @param estimate, as output by synthdid_estimate_weighted
#' @param treated.fraction, the fraction of pre-treatment data to use as a placebo treatment period
#'        Defaults to NULL, which indicates that it should be the fraction of post-treatment to pre-treatment data
#' @export synthdid_placebo_weighted
synthdid_placebo_weighted = function(estimate, treated.fraction = NULL) {
  setup = attr(estimate, 'setup')
  opts = attr(estimate, 'opts')
  weights = attr(estimate, 'weights')
  treated.weights = attr(estimate, 'treated.weights')
  X.beta = contract3(setup$X, weights$beta)

  if (is.null(treated.fraction)) { treated.fraction = 1 - setup$T0 / ncol(setup$Y) }
  placebo.T0 = floor(setup$T0 * (1 - treated.fraction))
  placebo.T1 = setup$T0 - placebo.T0

  # For placebo, use uniform period weights for the placebo post-period
  placebo.period.weights = rep(1 / placebo.T1, placebo.T1)

  synthdid_estimate_weighted(Y = setup$Y[, 1:setup$T0], N0 = setup$N0, T0 = placebo.T0,
                             X = setup$X[, 1:setup$T0, , drop = FALSE],
                             treated.weights = treated.weights,
                             period.weights = placebo.period.weights,
                             zeta.omega = opts$zeta.omega, zeta.lambda = opts$zeta.lambda,
                             omega.intercept = opts$omega.intercept, lambda.intercept = opts$lambda.intercept,
                             min.decrease = opts$min.decrease, max.iter = opts$max.iter)
}

#' Outputs the effect curve that was averaged to produce our estimate
#' @param estimate, as output by synthdid_estimate
#' @export synthdid_effect_curve
synthdid_effect_curve = function(estimate) {
  setup = attr(estimate, 'setup')
  weights = attr(estimate, 'weights')
  X.beta = contract3(setup$X, weights$beta)
  N1 = nrow(setup$Y) - setup$N0
  T1 = ncol(setup$Y) - setup$T0

  tau.sc = t(c(-weights$omega, rep(1 / N1, N1))) %*% (setup$Y - X.beta)
  tau.curve = tau.sc[setup$T0 + (1:T1)] - c(tau.sc[1:setup$T0] %*% weights$lambda)
  tau.curve
}

#' Outputs the effect curve that was averaged to produce the weighted estimate
#' @param estimate, as output by synthdid_estimate_weighted
#' @export synthdid_effect_curve_weighted
synthdid_effect_curve_weighted = function(estimate) {
  setup = attr(estimate, 'setup')
  weights = attr(estimate, 'weights')
  treated.weights = attr(estimate, 'treated.weights')
  X.beta = contract3(setup$X, weights$beta)
  N1 = nrow(setup$Y) - setup$N0
  T1 = ncol(setup$Y) - setup$T0

  # Use treated.weights instead of uniform weights
  if (is.null(treated.weights)) {
    treated.weights = rep(1 / N1, N1)
  }

  tau.sc = t(c(-weights$omega, treated.weights)) %*% (setup$Y - X.beta)
  tau.curve = tau.sc[setup$T0 + (1:T1)] - c(tau.sc[1:setup$T0] %*% weights$lambda)
  tau.curve
}
