#' Relative Reporting Risk (RRR)
#'
#' Determines the proportional reporting rate for a table of the form:
#' \tabular{lcc}{
#'    \tab event \tab not event\cr
#'   drug \tab \code{a} \tab \code{c}\cr
#'   not drug \tab \code{b} \tab \code{d}
#' }
#'
#' @inheritParams createTable
#' @param E The expected count when the drug and event are (approximately) independent
#'
#' @return The RRR of the table
#' @export
RRR <- function(a, b, c, d, E = ((a + b)*(a + c)) / (a + b + c + d)) {
  a / (((a + b) * (a + c)) / (a + b + c + d))
}