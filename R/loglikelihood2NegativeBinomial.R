#' Log-Likelihood of the Bimodal Negative Binomial 
#' 
#' Returns the log-likelihood of the bimodal negative binomial model 
#' used by the Gamma Poisson shrinker (GPS), see function 
#' [fitPriorParametersGPS()]. The function is written such 
#' that it can be used by the base function [nlminb()].
#' 
#' @param p A vector with the parameters (`alpha1`, `beta1`, 
#'          `alpha2`, `beta2` and `w`, in that order)
#' @param a A vector with the number of reports for each of the drug-event pairs 
#' @param E A vector (of the same length as `a`) with the number of reports 
#'          one would expect under the assumption of 'independence'
#' 
#' @return The negative log-likelihood (i.e., -1 * log-likelihood)
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
#' @examples 
#' alpha1 <- 0.2
#' beta1 <- 0.06
#' alpha2 <- 1.4
#' beta2 <- 1.8
#' w <- 0.1
#' 
#' a <- c(5, 1, 56, 3)
#' E <- c(3.4, 0.5, 10, 0.5) 
#' 
#' p <- c(alpha1, beta1, alpha2, beta2, w)
#' loglikelihood2NegativeBinomial(p, a, E)
#' #[1] 16.80512
#' 
#' @seealso [GPS()], [fitPriorParametersGPS()], [dbinbinom()]
#' @export
loglikelihood2NegativeBinomial <- function(p, a, E) {
  
  # unwrap p for readability's sake
  alpha1 <- p[1]
  beta1  <- p[2]
  alpha2 <- p[3]
  beta2  <- p[4]
  w      <- p[5]
  
  -sum(
    log(
      pvm::dbinbinom(a, 
                     size1 = alpha1, 
                     prob1 = beta1 / (beta1 + E), 
                     size2 = alpha2, 
                     prob2 = beta2 / (beta2 + E),
                     w = w)
    )
  )
}
