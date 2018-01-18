#' Binomial Log-Likelihood Ratio Test
#'
#' Performs the log-likelihood ratio test assuming a binomial distribution.
#' The test was first suggested by Duan et al. (2013).
#'
#' @inheritParams createTable
#'
#' @return log-likelihood ratio
#' @export
logLikelihoodRatioBinomial <- function(a, b, c, d) {
  N <- a + b + c + d

  r  <- a / N                     # actual observed ratio
  r0 <- ((a + c) * (a + b)) / N^2 # expected ratio under H0

  # log-likelihood
  a * (log(r) - log(r0)) + (N - a) * (log(1 - r) - log(1 - r0))
}
