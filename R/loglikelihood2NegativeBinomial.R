#' Log-Likelihood of the Bimodal Negative Binomial 
#' 
#' Returns the log-likelihood of the bimodal negative binomial model 
#' used by the Gamma Poisson Shrinker (GPS), see function \code{\link{GPS}}.
#' The function is written such that it can be used by the base function
#' \code{\link{nlm}}.
#' 
#' @param p A vector with the parameters (\code{alpha1}, \code{beta1}, 
#'          \code{alpha2}, \code{beta2} and \code{w}, in that order)
#' @param a A vector with the number of reports for each of the drug-event pairs 
#' @param E A vector (of the same length as \code{a}) with the number of reports 
#'          one would expect under the assumption of 'independence'
#' 
#' @return The log-likelihood 
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
      w * dnbinom(a, size = alpha1, prob = beta1 / (beta1 + E)) + 
        (1 - w) * dnbinom(a, size = alpha2, prob = beta2 / (beta2 + E))
    )
  )
}