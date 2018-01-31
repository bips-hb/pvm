#' Confidence Interval for the GPS
#' 
#' A function that can be used by the optimization
#' function \code{uniroot}. It is used by the 
#' \code{\link{GPS}} function to determine the 
#' lower endpoint of the (1 - \code{alpha})*100 
#' confidence interval. 
#' 
#' @param x The function value
#' @param a The observed number of reports for that drug-event pair
#' @param E The expected count for the drug-event pair
#' @param alpha1,beta1 The parameters for the first Gamma distribution
#' @param alpha2,beta2 The parameters for the second Gamma distribution
#' @param w The weight of the first mode
#' 
#' @return Function value
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
#' @export
GPSConfidenceInterval <- function(x, a, E, alpha,
                                  alpha1, beta1, 
                                  alpha2, beta2, 
                                  w) { 
  alpha - w*pgamma(x, shape = alpha1 + a, rate = beta1 + E) - 
    (1 - w)*pgamma(x, shape = alpha2 + a, rate = beta2 + E)
}