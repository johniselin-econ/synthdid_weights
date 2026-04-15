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
  sprintf("synthdid_estimate_weighted: %.3f", c(x))
}
