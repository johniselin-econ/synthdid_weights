#' Print a synthdid object
#' @param x The object to print
#' @param ... Additional arguments (currently ignored).
#' @method print synthdid_estimate
#' @export
print.synthdid_estimate = function(x, ...) { cat(format(x, ...), "\n") }

#' Print a weighted synthdid object
#' @param x The object to print
#' @param ... Additional arguments (currently ignored).
#' @method print synthdid_estimate_weighted
#' @export
print.synthdid_estimate_weighted = function(x, ...) { cat(format(x, ...), "\n") }

#' Format a weighted synthdid object
#' @param x The object to format
#' @param ... Additional arguments (currently ignored).
#' @method format synthdid_estimate_weighted
#' @export
format.synthdid_estimate_weighted = function(x, ...) {
  info = summary(x, fast = TRUE)
  d = as.list(info$dimensions)
  sprintf('synthdid_weighted: %1.3f +- %1.3f. Effective N0/N0 = %1.1f/%d~%1.1f. Effective T0/T0 = %1.1f/%d~%1.1f. Effective N1/N1 = %1.1f/%d~%1.1f. T1 = %d.',
          c(x), 1.96 * info$se,
          d$N0.effective, d$N0, d$N0.effective / d$N0,
          d$T0.effective, d$T0, d$T0.effective / d$T0,
          d$N1.effective, d$N1, d$N1.effective / d$N1,
          d$T1)
}
