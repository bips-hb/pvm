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
#' @param prior List that contains the prior parameters (see function \code{\link{fitPriorParametersGPS}})
#' @param alpha Value between \eqn{(0,1)}. If set, the lower endpoint of that confidence interval is returned
#'
#' @return a tibble
#' 
#' @references DuMouchel, W. (1999). Bayesian Data Mining in Large Frequency Tables, 
#'             with an Application to the FDA Spontaneous Reporting System. 
#'             The American Statistician, 53(3), 177–190. 
#'             https://doi.org/10.1080/00031305.1999.10474456
#' @export
GPS <- function(a, b, c, d, prior = fitPriorParametersGPS(a, b, c, d), alpha = NULL) {
  alpha1 <- prior$alpha1
  beta1  <- prior$beta1
  alpha2 <- prior$alpha2
  beta2  <- prior$beta2
  w      <- prior$w

  E = ((a + b)*(a + c)) / (a + b + c + d) # expected count

  Q <- w * dnbinom(a, size = alpha1, prob = beta1 / (beta1 + E)) /
    dbinbinom(a, size1 = alpha1, prob1 = beta1 / (beta1 + E), 
                 size2 = alpha2, prob2 = beta2 / (beta2 + E), w)
              
  EBlog2 <- (Q * (digamma(alpha1 + a) - log(beta1 + E)) + 
             (1 - Q) * (digamma(alpha2 + a) - log(beta2 + E))) / log(2) 
  
  EBGM <- 2^EBlog2

  if (is.null(alpha)) {
    return(EBGM)
  } else {
    return(log2(PhViD::.QuantileDuMouchel(alpha, Q, alpha1 + a, beta1 + E, alpha2 + a, beta2 + E)))
  }
}
