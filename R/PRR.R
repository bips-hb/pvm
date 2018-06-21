#' Proportional Reporting Rate (PRR)
#'
#' Determines the proportional reporting rate to a 
#' collection of 2 x 2 tables of the form
#' \tabular{lcc}{
#'    \tab event \tab not event\cr
#'   drug \tab \code{a} \tab \code{c}\cr
#'   not drug \tab \code{b} \tab \code{d}
#' }
#' In case the parameter \code{alpha} is set, it returns
#' the lower endpoint of the \eqn{100(1 - \alpha)} percent confidence interval.
#'
#' @inheritParams BCPNN
#'
#' @return The PRR or the lower endpoint of the confidence interval
#' @export
PRR <- function(a, b, c, d, alpha = NULL) {

  # determine the marginals for the drug
  q   <- a + c
  r   <- b + d

  prr <- (a / q) * (b / r)
  if (is.null(alpha)) {
    return(prr)
  } else{
    phi <- qnorm(1.0 - alpha / 2)
    return(exp(log(prr) - phi * sqrt((1 / a) - (1 / q) + (1 / b) - (1 / r))))
  }
}