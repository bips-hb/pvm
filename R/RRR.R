#' Relative Reporting Risk (RRR)
#'
#' Determines the RRR for a collection of 2 x 2 tables of the form
#' \tabular{lcc}{
#'    \tab event \tab not event\cr
#'   drug \tab `a` \tab `c`\cr
#'   not drug \tab `b` \tab `d`
#' }
#'
#' @template standardParams
#'
#' @return The RRR
#' @export
RRR <- function(a, b, c, d) {
  # The expected count when the drug and event are (approximately) independent
  E <- ((a + b)*(a + c)) / (a + b + c + d)
  a / E
}
