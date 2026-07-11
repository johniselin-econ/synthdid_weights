# =============================================================================
# STAGGERED-ADOPTION WEIGHTED SDID
# Cohort-by-cohort block weighted SDID against not-yet-treated (or
# never-treated-in-window) controls, aggregated with researcher treated-unit
# weights. Each cohort fit is a call to the block estimator
# synthdid_estimate_weighted; this file only splits the panel into
# adoption-date blocks and aggregates. Structural sibling of stratified.R.
#
# Estimand: tau^(omega-tilde, stag) = sum_g W_g tau_g^(omega-tilde), where
# tau_g is the within-cohort weighted SDID effect and W_g the cohort weight.
# With treated.weights = 1/N1 and cohort.weight = "treated.periods" this
# reduces to the standard staggered SDID of Clarke et al. (2023); a single
# cohort returns synthdid_estimate_weighted exactly.
# =============================================================================

#' Staggered-adoption weighted SDID estimate.
#'
#' Splits a balanced panel into adoption-date cohorts, runs weighted SDID within
#' each cohort against a pool of not-yet-treated (or never-treated) control units,
#' and aggregates the cohort estimates with researcher-chosen cohort weights.
#'
#' @param Y the N x T observation matrix (rows = units, columns = periods),
#'        balanced, with rows and columns in a fixed order shared with
#'        `adoption.time` and `time`.
#' @param adoption.time length-N vector giving each unit's first treated period,
#'        in the same units as `time`; use `Inf` for never-treated units. NA is an
#'        error (it would otherwise silently exclude the unit from every role);
#'        recode NA to `Inf` explicitly if it means never-treated.
#' @param time length-T vector of period labels, strictly increasing (enforced).
#'        Defaults to `1:ncol(Y)`. A unit j is treated in column t iff
#'        `time[t] >= adoption.time[j]`.
#' @param treated.weights length-N1 vector of researcher weights over the treated
#'        units (those with finite `adoption.time` no later than the last period),
#'        ordered by unit (i.e., by row of Y). If NULL, uniform 1/N1. Normalized to
#'        sum to 1 across all treated units; within each cohort they are renormalized.
#' @param cluster optional length-N vector of cluster IDs (e.g. state FIPS), stored
#'        for cluster-robust `vcov()` and subset into each cohort fit.
#' @param control control pool for each cohort. "notyet" (default) uses units not
#'        yet treated over the cohort's post-window (the window is capped at the next
#'        adoption among the candidate controls so they stay untreated); "never" uses
#'        units never treated within the observed window (`adoption.time` beyond the
#'        last period), the pool the Stata `sdid` package uses and the golden-test target.
#'        CAUTION: with densely spaced adoptions (e.g. one cohort per period), the
#'        "notyet" cap truncates each cohort's post-window -- to as little as a single
#'        period -- so the estimand becomes a correspondingly short-run effect; inspect
#'        `cohort.table$T1` before interpreting.
#' @param cohort.weight aggregation weight for cohort g. "treated.share" (default)
#'        uses W_g = sum of treated.weights in g, giving the estimand
#'        tau^(omega-tilde, stag). "treated.periods" multiplies that by the number of
#'        post periods, which at uniform weights reproduces the treated-cell weighting
#'        of Clarke et al. (2023).
#' @param min.controls minimum controls for a cohort to be estimable (default 2).
#' @param min.pre minimum pre-periods for a cohort to be estimable (default 2).
#' @param drop.infeasible if TRUE, cohorts with too few pre-periods, post-periods, or
#'        controls (or that lose a side under resampling) are dropped and the remaining
#'        cohort weights renormalized, rather than raising an error. Used by the
#'        bootstrap. Default FALSE.
#' @param ... additional options passed to every within-cohort
#'        [synthdid_estimate_weighted] call (e.g. `detrend = TRUE`). `X`,
#'        `period.weights`, and `weights` are rejected with an error: cohorts have
#'        different time windows and donor pools, so panel-shaped arguments cannot
#'        be passed through.
#' @return The aggregate estimate: a scalar of class 'synthdid_estimate_staggered' with
#'         attributes 'cohort.fits' (named list of within-cohort fits), 'cohort.table'
#'         (data.frame: cohort, N1, N0, T0, T1, share, weight, estimate), 'setup'
#'         (Y, adoption.time, time), 'treated.weights', 'treated.index', 'cluster', and 'opts'.
#' @export synthdid_estimate_staggered
synthdid_estimate_staggered = function(Y, adoption.time, time = NULL,
                                       treated.weights = NULL,
                                       cluster = NULL,
                                       control = c("notyet", "never"),
                                       cohort.weight = c("treated.share", "treated.periods"),
                                       min.controls = 2, min.pre = 2,
                                       drop.infeasible = FALSE,
                                       ...) {
  control = match.arg(control)
  cohort.weight = match.arg(cohort.weight)
  N = nrow(Y); Tt = ncol(Y)
  if (is.null(time)) time = seq_len(Tt)
  # NA adoption dates would silently vanish from every role (treated, never
  # pool, not-yet candidates) via which()'s NA-dropping; refuse them instead.
  # Recode NA to Inf explicitly if it means never-treated.
  stopifnot(length(adoption.time) == N, !anyNA(adoption.time),
            length(time) == Tt, !anyNA(time), all(diff(time) > 0),
            is.null(cluster) || length(cluster) == N)
  forbidden = intersect(c("X", "period.weights", "weights"), names(list(...)))
  if (length(forbidden))
    stop("unsupported in the staggered estimator (cohorts have different time windows ",
         "and donor pools): ", paste(forbidden, collapse = ", "))

  last.time = time[Tt]
  treated.index = which(is.finite(adoption.time) & adoption.time <= last.time)
  N1 = length(treated.index)
  stopifnot(N1 >= 1)

  if (is.null(treated.weights)) {
    tw.full = rep(1 / N1, N1)
  } else {
    stopifnot(length(treated.weights) == N1, all(treated.weights >= 0), sum(treated.weights) > 0)
    tw.full = treated.weights / sum(treated.weights)
  }
  # full-length weight lookup by unit row (0 for controls); positional indexing
  w.by.unit = numeric(N)
  w.by.unit[treated.index] = tw.full

  never.pool = which(adoption.time > last.time)   # never treated in window (incl. Inf)
  cohorts = sort(unique(adoption.time[treated.index]))

  fits = list(); rows = list(); dropped = numeric(0)
  for (g in cohorts) {
    trt.idx = which(adoption.time == g)
    pre     = which(time < g)
    if (control == "never") {
      ctl.idx = never.pool
      post    = which(time >= g)
    } else {                                  # not-yet-treated, window capped so controls stay clean
      cand = which(adoption.time > g)
      cap  = if (length(cand)) min(adoption.time[cand]) else Inf
      post = which(time >= g & time < cap)
      ctl.idx = cand
    }
    share = sum(w.by.unit[trt.idx])
    infeasible = length(pre) < min.pre || length(post) < 1 ||
                 length(ctl.idx) < min.controls || share <= 0
    if (infeasible) {
      if (drop.infeasible) { dropped = c(dropped, g); next }
      stop(sprintf(paste("cohort %s infeasible: %d pre-period(s), %d post-period(s), %d control(s), weight share %g;",
                         "need >= %d pre, >= 1 post, >= %d controls, share > 0 (or set drop.infeasible = TRUE)"),
                   as.character(g), length(pre), length(post), length(ctl.idx), share, min.pre, min.controls))
    }
    rows.idx = c(ctl.idx, trt.idx)            # controls first, then treated
    cols.idx = c(pre, post)                   # pre first, then post (time order within each)
    Yg  = Y[rows.idx, cols.idx, drop = FALSE]
    N0g = length(ctl.idx); T0g = length(pre)
    tw.g = w.by.unit[trt.idx]; tw.g = tw.g / sum(tw.g)

    fit = synthdid_estimate_weighted(Yg, N0g, T0g,
                                     treated.weights = as.numeric(tw.g),
                                     cluster = if (is.null(cluster)) NULL else cluster[rows.idx],
                                     ...)
    Wg    = if (cohort.weight == "treated.share") share else share * length(post)
    key   = as.character(g)
    fits[[key]] = fit
    rows[[key]] = data.frame(cohort = g, N1 = length(trt.idx), N0 = N0g,
                             T0 = T0g, T1 = length(post), share = share,
                             weight = Wg, estimate = as.numeric(fit),
                             stringsAsFactors = FALSE)
  }
  if (length(fits) == 0) stop("no estimable cohort")

  cohort.table = do.call(rbind, rows); rownames(cohort.table) = NULL
  # renormalize if any cohort was dropped (mirrors stratified.R)
  cohort.table$share  = cohort.table$share  / sum(cohort.table$share)
  cohort.table$weight = cohort.table$weight / sum(cohort.table$weight)
  estimate = sum(cohort.table$weight * cohort.table$estimate)

  class(estimate) = 'synthdid_estimate_staggered'
  attr(estimate, 'estimator')       = 'synthdid_estimate_staggered'
  attr(estimate, 'cohort.fits')     = fits
  attr(estimate, 'cohort.table')    = cohort.table
  attr(estimate, 'cohort.dropped')  = dropped
  attr(estimate, 'setup')           = list(Y = Y, adoption.time = adoption.time, time = time)
  attr(estimate, 'treated.weights') = tw.full
  attr(estimate, 'treated.index')   = treated.index
  attr(estimate, 'cluster')         = cluster
  attr(estimate, 'opts')            = c(list(control = control, cohort.weight = cohort.weight,
                                             min.controls = min.controls, min.pre = min.pre),
                                        list(...))
  estimate
}

#' @method print synthdid_estimate_staggered
#' @export
print.synthdid_estimate_staggered = function(x, ...) {
  cat(sprintf("staggered weighted SDID: %1.5f  (%d cohorts, control = %s)\n",
              as.numeric(x), nrow(attr(x, 'cohort.table')), attr(x, 'opts')$control))
  print(attr(x, 'cohort.table'))
  invisible(x)
}

#' Variance-covariance matrix for a staggered weighted SDID estimate (bootstrap only).
#'
#' Each replication resamples units (or whole clusters when a cluster vector is
#' available), rebuilds the cohort structure from the resampled rows, and re-runs the
#' entire staggered estimator with `drop.infeasible = TRUE`, so the within-cohort
#' omega/lambda weights, any detrending, and the cohort weights are all re-fit inside
#' every draw.
#'
#' @param object a synthdid_estimate_staggered model.
#' @param method only "bootstrap" is supported.
#' @param replications number of bootstrap replications.
#' @param cluster optional cluster IDs; defaults to the cluster stored in the estimate.
#'        Non-NULL resamples clusters; pass NULL to force unit-level resampling.
#' @param ... ignored.
#' @method vcov synthdid_estimate_staggered
#' @export
vcov.synthdid_estimate_staggered = function(object,
  method = c("bootstrap"),
  replications = 200,
  cluster = attr(object, 'cluster'), ...) {
    method = match.arg(method)
    samples = staggered_bootstrap_sample(object, replications, cluster)
    samples = samples[!is.na(samples)]
    if (length(samples) < 2) {
      warning("Staggered bootstrap: fewer than 2 valid replicates; returning NA")
      return(matrix(NA_real_))
    }
    R = length(samples)
    matrix((sqrt((R - 1) / R) * sd(samples))^2)
}

# `replications` resample-and-refit draws (unit- or cluster-resampled). Exposed
# (not exported) so parallel drivers can call it with their own seeds.
staggered_bootstrap_sample = function(object, replications, cluster = attr(object, 'cluster')) {
  replicate(replications,
    staggered_boot_rep(object, if (is.null(cluster)) "unit" else "cluster", cluster))
}

# A single resample-and-refit replication. method: "unit" or "cluster".
# Returns NA if the draw leaves no treated units or the refit fails.
staggered_boot_rep = function(object, method, cluster = NULL) {
  setup = attr(object, 'setup')
  opts  = attr(object, 'opts')
  Y = setup$Y; adoption.time = setup$adoption.time; time = setup$time
  N = nrow(Y); last.time = time[length(time)]

  w.full = numeric(N)
  w.full[attr(object, 'treated.index')] = attr(object, 'treated.weights')

  if (method == "cluster") {
    stopifnot(!is.null(cluster), length(cluster) == N)
    drawn = sample(unique(cluster), replace = TRUE)
    ind = unlist(lapply(drawn, function(cl) which(cluster == cl)))
  } else {
    ind = sample(seq_len(N), replace = TRUE)
  }
  if (length(ind) == 0) return(NA_real_)

  adt.b = adoption.time[ind]
  trt.b = which(is.finite(adt.b) & adt.b <= last.time)
  if (length(trt.b) < 1) return(NA_real_)
  tw.b = w.full[ind][trt.b]
  if (sum(tw.b) <= 0) return(NA_real_)

  args = c(list(Y = Y[ind, , drop = FALSE], adoption.time = adt.b, time = time,
                treated.weights = tw.b, cluster = NULL, drop.infeasible = TRUE),
           opts)
  tryCatch(as.numeric(do.call(synthdid_estimate_staggered, args)),
           error = function(e) NA_real_)
}
