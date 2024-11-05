#' near
#'
#' A safe way to compare two floating point numbers. The function
#' is based on the `near` function in the `dplyr`
#' package.
#'
#' @param x,y Numeric vectors to compare
#' @param tol Tolerance of comparison (Default: sqrt of the machine precision)
#'
#' @return `TRUE` when `x` and `y` are near, otherwise `FALSE`
#' @export
near <- function(x, y, tol = .Machine$double.eps^0.5) {
  return(abs(x - y) < tol)
}
