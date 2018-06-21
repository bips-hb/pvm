#' Chi-squared Test
#'
#' Performs the chi-squared test with or without Yates's continuity
#' correction to a collection of 2 x 2 tables of the form
#' \tabular{lcc}{
#'    \tab event \tab not event\cr
#'   drug \tab \code{a} \tab \code{c}\cr
#'   not drug \tab \code{b} \tab \code{d}
#' }
#' 
#' @note The standard warnings for when the counts are too low 
#' in the 2 x 2 tables are suppressed. Due to the sparse nature 
#' of spontaneous reporting data, this happens quite frequently.
#'
#' @inheritParams BCPNN
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
