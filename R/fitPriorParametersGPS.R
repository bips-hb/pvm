#' Prior Parameter Fit for the GPS
#'
#' Fits the prior parameters to the data for the Gamma Poisson Shrinker (GPS).
#' The initial guess for the parameter values are set the same as by DuMouchel (1999).
#'
#' The function is based on the code from the \code{PhViD} package.
#'
#' @param a a list that contains the counts of the upper left corners of the tables
#' @param b a list that contains the counts of the lower left corners of the tables
#' @param c a list that contains the counts of the upper right corners of the tables
#' @param d a list that contains the counts of the lower right corners of the tables
#' @param alpha1 Prior parameter \eqn{\alpha_1} (Default = .2)
#' @param beta1 Prior parameter \eqn{\beta_1} (Default = .06)
#' @param alpha2 Prior parameter \eqn{\alpha_2} (Default = 1.4)
#' @param beta2 Prior parameter \eqn{\beta_2} (Default = 1.8)
#' @param w Prior parameter \eqn{w} (Default = .1)
#'
#' @return the prior parameters
#' @export
fitPriorParametersGPS <- function(a, b, c, d,
                                  alpha1 = 0.2, beta1 = 0.1,
                                  alpha2 = 2.0, beta2 = 4,
                                  w = 1/3) {
  E = ((a + b)*(a + c)) / (a + b + c + d) # expected count

  # maximizing the log likelihood 
  res <- suppressWarnings(
          optim(par = c(alpha1, beta1, alpha2, beta2, w), 
                fn = pvm::loglikelihood2NegativeBinomial, 
                a = a, E = E, 
                method="L-BFGS-B",
                lower = c(0.0, 0.0, 0.0, 0.0, 0.0), 
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
