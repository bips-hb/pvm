#' Fisher's Exact Test
#'
#' Performs the (one-sided) Fisher's exact test on the table:
#' \tabular{lcc}{
#'    \tab event \tab not event\cr
#'   drug \tab \code{a} \tab \code{c}\cr
#'   not drug \tab \code{b} \tab \code{d}
#' }
#'
#' Wrapper function for the \code{Rcpp} functions \code{fishersTestGreater}
#' and \code{midPFishersTestGreater}.
#'
#' @inheritParams createTable
#' @param midpvalue The mid-p-value correction (suggested by Agresti) is applied
#'
#' @return p-value
#' @export
fisherExactTest <- function(a, b, c, d, midpvalue = FALSE) {
  if (midpvalue) {
    mapply(
      function(a, b, c, d) midPFisherTestGreater(a, b, c, d),
      a, b, c, d
    )
  } else {
    mapply(
      function(a, b, c, d) fisherTestGreater(a, b, c, d),
      a, b, c, d
    )
  }
}
