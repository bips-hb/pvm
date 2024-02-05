#' Gamma Poisson Shrinker (GPS)
#'
#' Applies the Gamma Poisson Shrinker (GPS) introduced by
#' DuMouchel (1999) to a collection of 2 x 2 tables of the form
#' \tabular{lcc}{
#'    \tab event \tab not event\cr
#'   drug \tab `a` \tab `c`\cr
#'   not drug \tab `b` \tab `d`
#' }
#'
#' @template standardParams
#' @param E Vector with the expected values when there are no associations. By default set to 
#'          the values used by DuMouchel (1999), i.e., `((a + b)*(a + c)) / (a + b + c + d)`.
#' @param prior List that contains the prior parameters. If not specified, automatically fitted to the data, 
#'              see [fitPriorParametersGPS()]. 
#' @template alphaParam
#' 
#' @return a vector with the GPS estimates
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
#'             
#' @seealso [fitPriorParametersGPS()]
#' @export
GPS <- function(a, b, c, d, E = ((a + b)*(a + c)) / (a + b + c + d),
                prior = fitPriorParametersGPS(a, b, c, d), alpha = NULL) {
  
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

  Q <- w * dnbinom(a, size = alpha1, prob = beta1 / (beta1 + E)) /
    dbinbinom(a, size1 = alpha1, prob1 = beta1 / (beta1 + E), 
                 size2 = alpha2, prob2 = beta2 / (beta2 + E), w)
              
  EBlog <- (Q * (digamma(alpha1 + a) - log(beta1 + E)) + 
             (1 - Q) * (digamma(alpha2 + a) - log(beta2 + E))) 
  
  EBGM <- exp(EBlog)

  if (is.null(alpha)) {
    return(EBGM)
  } else {
    EBlow <- EBGM 
    
    # estimate the lower end point of the (1 - alpha)*100 confidence interval
    EBlow[!is.na(EBlow)] <- PhViD::.QuantileDuMouchel(alpha, 
                                       Q[!is.na(EBlow)], 
                                       alpha1 + a[!is.na(EBlow)], 
                                       beta1 + E[!is.na(EBlow)], 
                                       alpha2 + a[!is.na(EBlow)], 
                                       beta2 + E[!is.na(EBlow)])
    return(EBlow)
  }
}
