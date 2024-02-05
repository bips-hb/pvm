#' Bimodal Negative Binomial
#'
#' Returns the values of the density function of a 
#' bimodal negative binomial distribution. 
#'
#' @param x The x-values
#' @param size1,shape1,prob1 The size, shape and prob parameters for the first mode
#' @param size2,shape2,prob2 The size, shape and prob parameters for the second mode
#' @param w The weight of the first mode (must lie in \eqn{[0,1]})
#'
#' @return The density for the values in x
#' 
#' @seealso \code{\link{logLikelihood2NegativeBinomial}}, \code{\link{fitPriorParametersGPS}}, \code{\link{GPS}}
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
dbinbinom <- function(x, size1, prob1, size2, prob2, w) {
  w * dnbinom(x, size = size1, prob = prob1) +
    (1 - w) * dnbinom(x, size = size2, prob = prob2)
}
