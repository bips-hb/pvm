#' Prior Parameter Fit for the GPS
#'
#' Fits the prior parameters to the data for the Gamma Poisson shrinker (GPS).
#' The initial guess for the parameter values are set the same as by DuMouchel (1999).
#'
#' @template standardParams
#' @param alpha1 Prior parameter \eqn{\alpha_1} (Default = 0.2)
#' @param beta1 Prior parameter \eqn{\beta_1} (Default = 0.1)
#' @param alpha2 Prior parameter \eqn{\alpha_2} (Default = 2.0)
#' @param beta2 Prior parameter \eqn{\beta_2} (Default = 4.0)
#' @param E Passed to `nlminb()` (Default = `((a + b)*(a + c)) / (a + b + c + d)`)
#' @param w Prior parameter \eqn{w} (Default = 1/3)
#' 
#' @return A list with the prior parameters
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
#' @examples 
#' a <- srdata$tables$a 
#' b <- srdata$tables$b 
#' c <- srdata$tables$c 
#' d <- srdata$tables$d 
#' 
#' fitPriorParametersGPS(a, b, c, d) 
#' 
#' # $alpha1
#' # [1] 98.28478
#' #
#' # $beta1
#' # [1] 16.48081
#' #
#' # $alpha2
#' # [1] 16.61439
#' #
#' # $beta2
#' # [1] 18.00642
#' #
#' # $w
#' # [1] 0.06132586
#' 
#' @seealso `loglikelihood2NegativeBinomial()`
#' @export
fitPriorParametersGPS <- function(a, b, c, d, 
                                  E = ((a + b)*(a + c)) / (a + b + c + d),
                                  alpha1 = 0.2, beta1 = 0.1,
                                  alpha2 = 2.0, beta2 = 4,
                                  w = 1/3) {
  
  # to overcome possible integer overflow later
  a <- as.numeric(a)
  b <- as.numeric(b)
  c <- as.numeric(c)
  d <- as.numeric(d) 

  # maximizing the log likelihood 
  res <- suppressWarnings(
      nlminb(c(alpha1, beta1, alpha2, beta2, w), 
             loglikelihood2NegativeBinomial, 
             a = a, 
             E = E, 
             lower = 0, 
             upper = c(Inf, Inf, Inf, Inf, 1.0))
  )
   
  # unpack the prior parameters
  list(
    alpha1 = res$par[1],
    beta1  = res$par[2],
    alpha2 = res$par[3],
    beta2  = res$par[4],
    w      = res$par[5]
  )
}
