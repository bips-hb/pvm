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
#'             
#'             DuMouchel, W., & Pregibon, D. (2001). Empirical bayes screening 
#'             for multi-item associations. Proceedings of the Seventh ACM 
#'             SIGKDD International Conference on Knowledge Discovery and 
#'             Data Mining - KDD ’01, (October), 67–76. 
#'             http://doi.org/10.1145/502512.502526
#' @export
GPS <- function(a, b, c, d, prior = fitPriorParametersGPS(a, b, c, d), alpha = NULL) {
  alpha1 <- prior$alpha1
  beta1  <- prior$beta1
  alpha2 <- prior$alpha2
  beta2  <- prior$beta2
  w      <- prior$w

  # to overcome possible integer overflow later
  a <- as.numeric(a)
  b <- as.numeric(b)
  c <- as.numeric(c)
  d <- as.numeric(d) 

  E = ((a + b)*(a + c)) / (a + b + c + d) # expected count

  Q <- w * dnbinom(a, size = alpha1, prob = beta1 / (beta1 + E)) /
    dbinbinom(a, size1 = alpha1, prob1 = beta1 / (beta1 + E), 
                 size2 = alpha2, prob2 = beta2 / (beta2 + E), w)
              
  EBlog <- (Q * (digamma(alpha1 + a) - log(beta1 + E)) + 
             (1 - Q) * (digamma(alpha2 + a) - log(beta2 + E))) 
  
  EBGM <- exp(EBlog)

  if (is.null(alpha)) {
    return(EBGM)
  } else {
    # estimate the lower end point of the (1 - alpha)*100 confidence interval
    n_pairs <- length(EBGM)
    EBlow <- rep(NA, n_pairs) # allocate memory
    
    # loop over all pairs
    for (p in 1:n_pairs) {
      res <- uniroot(GPSConfidenceInterval, 
                     interval = c(-1000, max(EBGM)),
                     a = a[p],
                     E = E[p],
                     alpha = alpha,
                     alpha1 = alpha1,
                     beta1 = beta1,
                     alpha2 = alpha2,
                     beta2 = beta2,
                     w = w)
      EBlow[p] <- res$root
    }
    
    return(EBlow)
  }
}
