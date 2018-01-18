#' Gamma Poisson Shrinker (GPS)
#'
#' Applies the Gamma Poisson Shrinker (GPS) introduced by
#' DuMouchel (1999) to a 2 x 2 table of the form
#' \tabular{lcc}{
#'    \tab event \tab not event\cr
#'   drug \tab \code{a} \tab \code{c}\cr
#'   not drug \tab \code{b} \tab \code{d}
#' }
#'
#' The function is based on a function in the \code{PhVid} Packages (GPS)
#'
#' @inheritParams createTable
#' @param prior List that contains the prior parameters (see function \code{fitPriorParametersGPS})
#' @param alpha Value between \eqn{(0,1)}. If set, the lower endpoint that confidence interval is returned
#'
#' @return a tibble
#' @export
GPS <- function(a, b, c, d, prior = fitPriorParametersGPS(a, b, c, d), alpha = NULL) {
  alpha1 <- prior$alpha1
  beta1  <- prior$beta1
  alpha2 <- prior$alpha2
  beta2  <- prior$beta2
  w      <- prior$w

  expected_count <- ((a + b)*(a + c)) / (a + b + c + d)

  Q <- w * dnbinom(a, size = alpha1, prob = beta1/(beta1 + expected_count))/(w * dnbinom(a, size = alpha1, prob = beta1/(beta1 + expected_count)) + (1 - w) * dnbinom(a, size = alpha2, prob = beta2/(beta2 + expected_count)))
  EBGM <- log(2)^(-1) * (Q * (digamma(alpha1 + a) - log(beta1 + expected_count)) + (1 - Q) * (digamma(alpha2 + a) - log(beta2 + expected_count)))

  if (is.null(alpha)) {
    return(EBGM)
  } else {
    return(log(PhViD::.QuantileDuMouchel(alpha, Q, alpha1 + a, beta1 + expected_count, alpha2 + a, beta2 + expected_count), 2))
  }
}
