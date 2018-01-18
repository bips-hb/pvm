#' \eqn{\chi 2} Test
#'
#' Performs the \eqn{\chi 2} test with or without Yates's continuity
#' correction on the table
#' \tabular{lcc}{
#'    \tab event \tab not event\cr
#'   drug \tab \code{a} \tab \code{c}\cr
#'   not drug \tab \code{b} \tab \code{d}
#' }
#'
#' Warning about the Chi-squared approximation to be incorrect are
#' suppressed (this case is quite common in the type of 2x2 tables
#' one commonly observes in the field of pharmacovigilance).
#'
#' @inheritParams createTable
#' @param yates If \code{TRUE}, Yates's correction is used
#'
#' @return p-value
#' @export
chi2Test <- function(a, b, c, d, yates = FALSE) {
  mapply(
    function(a, b, c, d) {
      suppressWarnings(chisq.test(createTable(a, b, c, d), correct = yates)$p.value)
    },
    a, b, c, d
  )
}
