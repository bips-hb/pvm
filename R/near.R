#' near
#'
#' A safe way to compare two floating point numbers. The function
#' is based on the \code{near} function in the \code{dplyr}
#' package.
#'
#' @param x,y Numeric vectors to compare
#' @param tol Tolerance of comparison (Default: sqrt of the machine precision)
#'
#' @return \code{TRUE} when \code{x} and \code{y} are near, otherwise \code{FALSE}
#' @export
near <- function(x, y, tol = .Machine$double.eps^0.5) {
  return(abs(x - y) < tol)
}
