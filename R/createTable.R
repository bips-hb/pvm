#' Create 2 x 2 Table
#'
#' Returns a 2 x 2 contingency table of the form:
#' \tabular{lcc}{
#'    \tab event \tab not event\cr
#'   drug \tab \code{a} \tab \code{c}\cr
#'   not drug \tab \code{b} \tab \code{d}
#' }
#'
#' @param a Count in the upper left corner of the table
#' @param b Count in the lower left corner of the table
#' @param c Count in the upper right corner of the table
#' @param d Count in the lower right corner of the table
#'
#' @return A 2 x 2 matrix
#' @export
createTable <- function(a, b, c, d) {
  matrix(c(a, b, c, d), nrow = 2)
}