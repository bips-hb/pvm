#' Reporting Odds Ratio (ROR)
#'
#' Determines the reporting odds ratio for a table of the form:
#' \tabular{lcc}{
#'    \tab event \tab not event\cr
#'   drug \tab \code{a} \tab \code{c}\cr
#'   not drug \tab \code{b} \tab \code{d}
#' }
#' In case the parameter \code{alpha} is set, it returns
#' the lower endpoint of the \eqn{100(1 - \alpha)} percent confidence interval.
#'
#' @inheritParams createTable
#' @param alpha Value between \eqn{(0,1)}. If set, the lower endpoint that confidence interval is returned
#'
#' @return The ROR or the lower endpoint of the confidence interval of the ROR
#' @export
ROR <- function(a, b, c, d, alpha = NULL) {
  ror <- (a * d) / (b * c)
  if (is.null(alpha)) {
    return(ror)
  } else{
    phi <- qnorm(1.0 - alpha / 2)
    return(exp(log(ror) - phi * sqrt((1 / a) + (1 / b) + (1 + c) + (1 / d))))
  }
}