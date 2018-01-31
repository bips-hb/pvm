#' Prior Parameter Fit for the GPS
#'
#' Fits the prior parameters to the data for the Gamma Poisson Shrinker (GPS).
#' The initial guess for the parameter values are set the same as by DuMouchel (1999).
#'
#' @param a A list that contains the counts of the upper left corners of the tables
#' @param b A list that contains the counts of the lower left corners of the tables
#' @param c A list that contains the counts of the upper right corners of the tables
#' @param d A list that contains the counts of the lower right corners of the tables
#' @param alpha1 Prior parameter \eqn{\alpha_1} (Default = 0.2)
#' @param beta1 Prior parameter \eqn{\beta_1} (Default = 0.1)
#' @param alpha2 Prior parameter \eqn{\alpha_2} (Default = 2.0)
#' @param beta2 Prior parameter \eqn{\beta_2} (Default = 4.0)
#' @param w Prior parameter \eqn{w} (Default = 1/3)
#' 
#' @return A list with the prior parameters
#' 
#' @references DuMouchel, W. (1999). Bayesian Data Mining in Large Frequency Tables, 
#'             with an Application to the FDA Spontaneous Reporting System. 
#'             The American Statistician, 53(3), 177–190. 
#'             https://doi.org/10.1080/00031305.1999.10474456
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
#' @export
fitPriorParametersGPS <- function(a, b, c, d,
                                  alpha1 = 0.2, beta1 = 0.1,
                                  alpha2 = 2.0, beta2 = 4,
                                  w = 1/3) {
  
  # to overcome possible integer overflow later
  a <- as.numeric(a)
  b <- as.numeric(b)
  c <- as.numeric(c)
  d <- as.numeric(d) 
  
  E = ((a + b)*(a + c)) / (a + b + c + d) # expected count
  
  # remove cases where E is zero
  a <- a[E != 0]
  E <- E[E != 0]
  
  # maximizing the log likelihood 
  res <- suppressWarnings(
    optim(par = c(alpha1, beta1, alpha2, beta2, w), 
                fn = pvm::loglikelihood2NegativeBinomial, 
                a = a, E = E)
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
