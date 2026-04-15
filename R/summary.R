#' Outputs a table of important synthetic controls and their corresponding weights, sorted by weight.
#' The table is truncated to exclude synthetic controls that do not matter for any estimate ---
#' for each estimate, the truncated controls may have total weight no larger that 1-mass.
#' @param estimates, a list of estimates output by synthdid_estimate. Or a single estimate.
#' @param sort.by, the index of the estimate to sort by. Defaults to 1.
#' @param mass, which controls the length of the table. Defaults to 0.9.
#' @param weight.type, 'omega' for units, 'lambda' for time periods
#' @export synthdid_controls
synthdid_controls = function(estimates, sort.by = 1, mass = .9, weight.type = 'omega') {
  if (class(estimates) == 'synthdid_estimate') { estimates = list(estimates) }
  if (is.null(names(estimates))) { names(estimates) = sprintf('estimate %d', 1:length(estimates)) }
  if (!weight.type %in% c('omega', 'lambda')) { stop('weight.type must be "omega" or "lambda"') } 
  weights = do.call(cbind, lapply(estimates, function(est) { attr(est, 'weights')[[weight.type]] }))
  if (is.null(dim(weights))) { dim(weights) = c(length(weights), 1) }

  Y = attr(estimates[[1]], 'setup')$Y
  o = rev(order(weights[, sort.by]))
  tab = weights[o, , drop = FALSE]
  rownames(tab) = if(weight.type == 'omega') { rownames(Y)[o] } else { colnames(Y)[o] }
  colnames(tab) = names(estimates)

  # truncate table to retain a weight sum of at least mass for each unit
  tab.len = max(apply(tab, 2, function(col) { Position(function(x) { x >= mass }, cumsum(col), nomatch=nrow(tab)) }))
  tab[1:tab.len, , drop=FALSE]
}

#' Summarize a synthdid object
#' @param object The object to summarize
#' @param weight.digits The number of digits to use when displaying weights (omega, lambda)
#' @param fast Be fast but less accurate, e.g. jackknife instead of bootstrap se estimate
#' @param ... Additional arguments (currently ignored).
#' @method summary synthdid_estimate
#' @export
summary.synthdid_estimate = function(object, weight.digits=3, fast=FALSE, ...) {
  N0 = attr(object, 'setup')$N0
  T0 = attr(object, 'setup')$T0
  list(estimate = c(object),
    se = sqrt(if(fast) { vcov(object, method = 'jackknife') } else { vcov(object) }),
    controls = round(synthdid_controls(object, weight.type='omega'),  digits=weight.digits),
    periods  = round(synthdid_controls(object, weight.type='lambda'), digits=weight.digits),
    dimensions = c( N1 = nrow(Y(object))-N0, N0 = N0, N0.effective = round(1 / sum(omega(object)^2),  weight.digits),
		    T1 = ncol(Y(object))-T0, T0 = T0, T0.effective = round(1 / sum(lambda(object)^2), weight.digits)))
}
# =============================================================================
# WEIGHTED VERSIONS
# =============================================================================

#' Outputs a table of important synthetic controls and their corresponding weights (weighted version)
#' @param estimates, a list of weighted estimates. Or a single estimate.
#' @param sort.by, the index of the estimate to sort by. Defaults to 1.
#' @param mass, which controls the length of the table. Defaults to 0.9.
#' @param weight.type, 'omega' for units, 'lambda' for time periods, 'treated' for treated unit weights
#' @export synthdid_controls_weighted
synthdid_controls_weighted = function(estimates, sort.by = 1, mass = .9, weight.type = 'omega') {
  if (class(estimates) == 'synthdid_estimate_weighted') { estimates = list(estimates) }
  if (is.null(names(estimates))) { names(estimates) = sprintf('estimate %d', 1:length(estimates)) }
  if (!weight.type %in% c('omega', 'lambda', 'treated', 'period')) {
    stop('weight.type must be "omega", "lambda", "treated", or "period"')
  }

  if (weight.type == 'treated') {
    weights = do.call(cbind, lapply(estimates, function(est) { attr(est, 'treated.weights') }))
    Y = attr(estimates[[1]], 'setup')$Y
    N0 = attr(estimates[[1]], 'setup')$N0
    if (is.null(dim(weights))) { dim(weights) = c(length(weights), 1) }
    o = rev(order(weights[, sort.by]))
    tab = weights[o, , drop = FALSE]
    rownames(tab) = rownames(Y)[(N0+1):nrow(Y)][o]
    colnames(tab) = names(estimates)
    return(tab)
  }

  if (weight.type == 'period') {
    weights = do.call(cbind, lapply(estimates, function(est) { attr(est, 'period.weights') }))
    Y = attr(estimates[[1]], 'setup')$Y
    T0 = attr(estimates[[1]], 'setup')$T0
    if (is.null(dim(weights))) { dim(weights) = c(length(weights), 1) }
    o = rev(order(weights[, sort.by]))
    tab = weights[o, , drop = FALSE]
    rownames(tab) = colnames(Y)[(T0+1):ncol(Y)][o]
    colnames(tab) = names(estimates)
    return(tab)
  }

  # omega and lambda weights - same as original
  weights = do.call(cbind, lapply(estimates, function(est) { attr(est, 'weights')[[weight.type]] }))
  if (is.null(dim(weights))) { dim(weights) = c(length(weights), 1) }

  Y = attr(estimates[[1]], 'setup')$Y
  o = rev(order(weights[, sort.by]))
  tab = weights[o, , drop = FALSE]
  rownames(tab) = if(weight.type == 'omega') { rownames(Y)[o] } else { colnames(Y)[o] }
  colnames(tab) = names(estimates)

  # truncate table to retain a weight sum of at least mass for each unit
  tab.len = max(apply(tab, 2, function(col) { Position(function(x) { x >= mass }, cumsum(col), nomatch=nrow(tab)) }))
  tab[1:tab.len, , drop=FALSE]
}

#' Summarize a weighted synthdid object
#' @param object The object to summarize
#' @param weight.digits The number of digits to use when displaying weights
#' @param fast Be fast but less accurate, e.g. jackknife instead of bootstrap se estimate
#' @param ... Additional arguments (currently ignored).
#' @method summary synthdid_estimate_weighted
#' @export
summary.synthdid_estimate_weighted = function(object, weight.digits=3, fast=FALSE, ...) {
  setup = attr(object, 'setup')
  weights = attr(object, 'weights')
  N0 = setup$N0
  T0 = setup$T0
  Y = setup$Y
  treated.weights = attr(object, 'treated.weights')
  period.weights = attr(object, 'period.weights')
  omega = weights$omega
  lambda = weights$lambda

  list(estimate = c(object),
    se = sqrt(if(fast) { vcov(object, method = 'jackknife') } else { vcov(object) }),
    controls = round(synthdid_controls_weighted(object, weight.type='omega'),  digits=weight.digits),
    periods  = round(synthdid_controls_weighted(object, weight.type='lambda'), digits=weight.digits),
    treated.weights = round(synthdid_controls_weighted(object, weight.type='treated'), digits=weight.digits),
    period.weights  = round(synthdid_controls_weighted(object, weight.type='period'),  digits=weight.digits),
    dimensions = c( N1 = nrow(Y)-N0, N0 = N0, N0.effective = round(1 / sum(omega^2),  weight.digits),
                    T1 = ncol(Y)-T0, T0 = T0, T0.effective = round(1 / sum(lambda^2), weight.digits),
                    N1.effective = round(1 / sum(treated.weights^2), weight.digits),
                    T1.effective = round(1 / sum(period.weights^2),  weight.digits)))
}
