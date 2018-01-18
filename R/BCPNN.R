#' Bayesian Confidence Propagation Neural Network (BCPNN)
#'
#' Applies the BPCNN to a 2 x 2 table of the form
#' \tabular{lcc}{
#'    \tab event \tab not event\cr
#'   drug \tab \code{a} \tab \code{c}\cr
#'   not drug \tab \code{b} \tab \code{d}
#' }
#'
#' Based on the \code{BCPNN} function from the \code{PhViD}
#' package.
#'
#' @inheritParams ROR
#'
#' @return The maximum aposteriori estimate of the information component or the lower endpoint of the approximate credible interval
#' @export
BCPNN <- function(a, b, c, d, alpha = NULL) {
  n1. <- a + c
  n.1 <- a + b
  n   <- a + b + c + d
  p1  <- 1 + n1.
  p2  <- 1 + n - n1.
  q1  <- 1 + n.1
  q2  <- 1 + n - n.1
  r1  <- a + 1
  r2b <- n - a - 1 + (2 + n)^2/(q1 * p1)
  IC  <- log(2)^(-1) * (digamma(r1) - digamma(r1 + r2b) -
                          (digamma(p1) - digamma(p1 + p2) + digamma(q1) -
                             digamma(q1 + q2)))
  VICb <- log(2)^(-2) * (trigamma(r1) - trigamma(r1 +
                                                   r2b) + (trigamma(p1) - trigamma(p1 + p2) + trigamma(q1) -
                                                             trigamma(q1 + q2)))
  if (is.null(alpha)) {
    return(IC)
  } else {
    return(qnorm(alpha, IC, sqrt(VICb)))
  }
}