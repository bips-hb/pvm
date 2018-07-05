#' Yule's Q
#'
#' Determines Yule's Q for a collection of 2 x 2 tables 
#' of the form
#' \tabular{lcc}{
#'    \tab event \tab not event\cr
#'   drug \tab \code{a} \tab \code{c}\cr
#'   not drug \tab \code{b} \tab \code{d}
#' }
#' In case the parameter \code{alpha} is set, it returns
#' the lower endpoint of the \eqn{100(1 - \alpha)} percent confidence interval.
#'
#' @template standardParams
#' @template alphaParam
#'
#' @return Yule's Q or the lower endpoint of the confidence interval
#' @export
YulesQ <- function(a, b, c, d, alpha = NULL) {
  Q <- (a * d - b * c) / (a * d + b * c)
  if (is.null(alpha)) {
    return(Q)
  } else{
    phi <- qnorm(1.0 - alpha / 2)
    return(Q - phi * ((1 - Q^2)/2 * sqrt((1 / a) + (1 / b) + (1 + c) + (1 / d))))
  }
}