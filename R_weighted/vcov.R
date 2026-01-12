#' Calculate Variance-Covariance Matrix for a Fitted Model Object
#'
#' Provides variance estimates based on the following three options
#' \itemize{
#'   \item The bootstrap, Algorithm 2 in Arkhangelsky et al.
#'   \item The jackknife, Algorithm 3 in Arkhangelsky et al.
#'   \item Placebo, Algorithm 4 in Arkhangelsky et al.
#' }
#'
#' The jackknife is not recommended for SC, see section 5 in Arkhangelsky et al.
#' "placebo" is the only option that works for only one treated unit.
#'
#' @param object A synthdid model
#' @param method, the CI method. The default is bootstrap (warning: this may be slow on large
#'  data sets, the jackknife option is the fastest, with the caveat that it is not recommended
#'  for SC).
#' @param replications, the number of bootstrap replications
#' @param ... Additional arguments (currently ignored).
#'
#' @references Dmitry Arkhangelsky, Susan Athey, David A. Hirshberg, Guido W. Imbens, and Stefan Wager.
#'  "Synthetic Difference in Differences". arXiv preprint arXiv:1812.09970, 2019.
#'
#' @method vcov synthdid_estimate
#' @export
vcov.synthdid_estimate = function(object,
  method = c("bootstrap", "jackknife", "placebo"),
  replications = 200, ...) {
    method = match.arg(method)
    if(method == 'bootstrap') {
	se = bootstrap_se(object, replications)
    } else if(method == 'jackknife') {
	se = jackknife_se(object)
    } else if(method == 'placebo') {
	se = placebo_se(object, replications)
    }
    matrix(se^2)
}

#' Calculate the standard error of a synthetic diff in diff estimate. Deprecated. Use vcov.synthdid_estimate.
#' @param ... Any valid arguments for vcov.synthdid_estimate
#' @export synthdid_se
synthdid_se = function(...) { sqrt(vcov(...)) }


# The bootstrap se: Algorithm 2 of Arkhangelsky et al.
bootstrap_se = function(estimate, replications) { sqrt((replications-1)/replications) * sd(bootstrap_sample(estimate, replications)) }
bootstrap_sample = function(estimate, replications) {
    setup = attr(estimate, 'setup')
    opts = attr(estimate, 'opts')
    weights = attr(estimate, 'weights')
    if (setup$N0 == nrow(setup$Y) - 1) { return(NA) }
    theta = function(ind) {
	if(all(ind <= setup$N0) || all(ind > setup$N0)) { NA }
	else {
	    weights.boot = weights
	    weights.boot$omega = sum_normalize(weights$omega[sort(ind[ind <= setup$N0])])
	    do.call(synthdid_estimate, c(list(Y=setup$Y[sort(ind),], N0=sum(ind <= setup$N0), T0=setup$T0, X=setup$X[sort(ind), ,], weights=weights.boot), opts))
	}
    }
    bootstrap.estimates = rep(NA, replications)
    count = 0
    while(count < replications) {
	bootstrap.estimates[count+1] = theta(sample(1:nrow(setup$Y), replace=TRUE))
	if(!is.na(bootstrap.estimates[count+1])) { count = count+1 }
    }
    bootstrap.estimates
}


# The fixed-weights jackknife estimate of variance: Algorithm 3 of Arkhangelsky et al.
# if weights = NULL is passed explicitly, calculates the usual jackknife estimate of variance.
# returns NA if there is one treated unit or, for the fixed-weights jackknife, one control with nonzero weight
jackknife_se = function(estimate, weights = attr(estimate, 'weights')) {
    setup = attr(estimate, 'setup')
    opts = attr(estimate, 'opts')
    if (!is.null(weights)) {
      opts$update.omega = opts$update.lambda = FALSE
    }
    if (setup$N0 == nrow(setup$Y) - 1 || (!is.null(weights) && sum(weights$omega != 0) == 1)) { return(NA) }
    theta = function(ind) {
	weights.jk = weights
	if (!is.null(weights)) { weights.jk$omega = sum_normalize(weights$omega[ind[ind <= setup$N0]]) }
	estimate.jk = do.call(synthdid_estimate,
	    c(list(Y=setup$Y[ind, ], N0=sum(ind <= setup$N0), T0=setup$T0, X = setup$X[ind, , ], weights = weights.jk), opts))
    }
    jackknife(1:nrow(setup$Y), theta)
}

#' Jackknife standard error of function `theta` at samples `x`.
#' @param x vector of samples
#' @param theta a function which returns a scalar estimate
#' @importFrom stats var
#' @keywords internal
jackknife = function(x, theta) {
  n = length(x)
  u = rep(0, n)
  for (i in 1:n) {
    u[i] = theta(x[-i])
  }
  jack.se = sqrt(((n - 1) / n) * (n - 1) * var(u))

  jack.se
}



# The placebo se: Algorithm 4 of Arkhangelsky et al.
placebo_se = function(estimate, replications) {
    setup = attr(estimate, 'setup')
    opts = attr(estimate, 'opts')
    weights = attr(estimate, 'weights')
    N1 = nrow(setup$Y) - setup$N0
    if (setup$N0 <= N1) { stop('must have more controls than treated units to use the placebo se') }
    theta = function(ind) {
	N0 = length(ind)-N1
	weights.boot = weights
	weights.boot$omega = sum_normalize(weights$omega[ind[1:N0]])
        do.call(synthdid_estimate, c(list(Y=setup$Y[ind,], N0=N0,  T0=setup$T0,  X=setup$X[ind, ,], weights=weights.boot), opts))
    }
    sqrt((replications-1)/replications) * sd(replicate(replications, theta(sample(1:setup$N0))))
}

sum_normalize = function(x) {
    if(sum(x) != 0) { x / sum(x) }
    else { rep(1/length(x), length(x)) }
    # if given a vector of zeros, return uniform weights
    # this fine when used in bootstrap and placebo standard errors, where it is used only for initialization
    # for jackknife standard errors, where it isn't, we handle the case of a vector of zeros without calling this function.
}


# =============================================================================
# WEIGHTED VERSIONS OF VARIANCE ESTIMATION
# =============================================================================

#' Calculate Variance-Covariance Matrix for a Weighted Synthetic DID Estimate
#'
#' Provides variance estimates for weighted synthdid estimates.
#' The key difference from the unweighted version is that treated unit weights
#' must be renormalized when resampling.
#'
#' @param object A synthdid_estimate_weighted model
#' @param method, the CI method. The default is bootstrap.
#' @param replications, the number of bootstrap replications
#' @param ... Additional arguments (currently ignored).
#'
#' @method vcov synthdid_estimate_weighted
#' @export
vcov.synthdid_estimate_weighted = function(object,
  method = c("bootstrap", "jackknife", "placebo"),
  replications = 200, ...) {
    method = match.arg(method)
    if(method == 'bootstrap') {
      se = bootstrap_se_weighted(object, replications)
    } else if(method == 'jackknife') {
      se = jackknife_se_weighted(object)
    } else if(method == 'placebo') {
      se = placebo_se_weighted(object, replications)
    }
    matrix(se^2)
}

#' Calculate the standard error of a weighted synthetic diff in diff estimate
#' @param ... Any valid arguments for vcov.synthdid_estimate_weighted
#' @export synthdid_se_weighted
synthdid_se_weighted = function(...) { sqrt(vcov(...)) }


# The bootstrap se for weighted estimates: modified Algorithm 2
# Key change: renormalize treated.weights for resampled treated units
bootstrap_se_weighted = function(estimate, replications) {
  sqrt((replications-1)/replications) * sd(bootstrap_sample_weighted(estimate, replications))
}

bootstrap_sample_weighted = function(estimate, replications) {
    setup = attr(estimate, 'setup')
    opts = attr(estimate, 'opts')
    weights = attr(estimate, 'weights')
    treated.weights = attr(estimate, 'treated.weights')
    period.weights = attr(estimate, 'period.weights')
    N1 = nrow(setup$Y) - setup$N0

    if (setup$N0 == nrow(setup$Y) - 1) { return(NA) }

    theta = function(ind) {
      control.ind = sort(ind[ind <= setup$N0])
      treated.ind = sort(ind[ind > setup$N0])
      treated.ind.local = treated.ind - setup$N0  # indices within treated units

      if(length(control.ind) == 0 || length(treated.ind) == 0) { return(NA) }

      # Renormalize control weights
      weights.boot = weights
      weights.boot$omega = sum_normalize(weights$omega[control.ind])

      # Renormalize treated weights for resampled treated units
      # Count how many times each treated unit appears
      treated.counts = tabulate(treated.ind.local, nbins = N1)
      treated.weights.boot = treated.weights * treated.counts
      treated.weights.boot = sum_normalize(treated.weights.boot)

      # Reconstruct Y matrix with resampled units
      Y.boot = setup$Y[c(control.ind, treated.ind), ]
      X.boot = setup$X[c(control.ind, treated.ind), , ]
      N0.boot = length(control.ind)

      do.call(synthdid_estimate_weighted,
              c(list(Y = Y.boot, N0 = N0.boot, T0 = setup$T0, X = X.boot,
                     treated.weights = treated.weights.boot,
                     period.weights = period.weights,
                     weights = weights.boot), opts))
    }

    bootstrap.estimates = rep(NA, replications)
    count = 0
    while(count < replications) {
      bootstrap.estimates[count+1] = theta(sample(1:nrow(setup$Y), replace=TRUE))
      if(!is.na(bootstrap.estimates[count+1])) { count = count+1 }
    }
    bootstrap.estimates
}


# The fixed-weights jackknife estimate for weighted estimates: modified Algorithm 3
# Key change: renormalize treated.weights when leaving out treated units
jackknife_se_weighted = function(estimate, weights = attr(estimate, 'weights')) {
    setup = attr(estimate, 'setup')
    opts = attr(estimate, 'opts')
    treated.weights.orig = attr(estimate, 'treated.weights')
    period.weights = attr(estimate, 'period.weights')
    N1 = nrow(setup$Y) - setup$N0

    if (!is.null(weights)) {
      opts$update.omega = opts$update.lambda = FALSE
    }
    if (setup$N0 == nrow(setup$Y) - 1 || (!is.null(weights) && sum(weights$omega != 0) == 1)) {
      return(NA)
    }

    theta = function(ind) {
      control.ind = ind[ind <= setup$N0]
      treated.ind = ind[ind > setup$N0]
      treated.ind.local = treated.ind - setup$N0

      # Renormalize control weights
      weights.jk = weights
      if (!is.null(weights)) {
        weights.jk$omega = sum_normalize(weights$omega[control.ind])
      }

      # Renormalize treated weights for remaining treated units
      treated.weights.jk = treated.weights.orig[treated.ind.local]
      treated.weights.jk = sum_normalize(treated.weights.jk)

      estimate.jk = do.call(synthdid_estimate_weighted,
          c(list(Y = setup$Y[ind, ], N0 = sum(ind <= setup$N0), T0 = setup$T0,
                 X = setup$X[ind, , ],
                 treated.weights = treated.weights.jk,
                 period.weights = period.weights,
                 weights = weights.jk), opts))
    }
    jackknife(1:nrow(setup$Y), theta)
}


# The placebo se for weighted estimates: modified Algorithm 4
# For placebo, we reassign N1 control units as "treated" and estimate effect
# Key change: need to assign weights to the placebo treated units
placebo_se_weighted = function(estimate, replications) {
    setup = attr(estimate, 'setup')
    opts = attr(estimate, 'opts')
    weights = attr(estimate, 'weights')
    treated.weights.orig = attr(estimate, 'treated.weights')
    period.weights = attr(estimate, 'period.weights')
    N1 = nrow(setup$Y) - setup$N0

    if (setup$N0 <= N1) {
      stop('must have more controls than treated units to use the placebo se')
    }

    theta = function(ind) {
      N0.placebo = length(ind) - N1

      # Renormalize control weights for remaining controls
      weights.boot = weights
      weights.boot$omega = sum_normalize(weights$omega[ind[1:N0.placebo]])

      # For placebo treated units, use uniform weights
      # (original treated.weights don't apply to control units acting as placebo treated)
      treated.weights.placebo = rep(1/N1, N1)

      do.call(synthdid_estimate_weighted,
              c(list(Y = setup$Y[ind, ], N0 = N0.placebo, T0 = setup$T0,
                     X = setup$X[ind, , ],
                     treated.weights = treated.weights.placebo,
                     period.weights = period.weights,
                     weights = weights.boot), opts))
    }

    sqrt((replications-1)/replications) * sd(replicate(replications, theta(sample(1:setup$N0))))
}
