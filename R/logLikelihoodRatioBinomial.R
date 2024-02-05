#' Binomial Log-Likelihood Ratio Test
#'
#' Performs the log-likelihood ratio test for 
#' a collection of 2 x 2 tables of the form
#' \tabular{lcc}{
#'    \tab event \tab not event\cr
#'   drug \tab `a` \tab `c`\cr
#'   not drug \tab `b` \tab `d`
#' }
#'
#' @template standardParams
#' @return The log-likelihood ratio
#' 
#' @references Lian Duan, Khoshneshin, M., Street, W. N., 
#'             & Mei Liu. (2013). Adverse drug effect detection. 
#'             IEEE Journal of Biomedical and Health Informatics, 
#'             17(2), 305–11. https://doi.org/10.1109/TITB.2012.2227272
#' @export
logLikelihoodRatioBinomial <- function(a, b, c, d) {
  N <- a + b + c + d

  r  <- a / N                     # actual observed ratio
  r0 <- ((a + c) * (a + b)) / N^2 # expected ratio under H0

  # log-likelihood
  a * (log(r) - log(r0)) + (N - a) * (log(1 - r) - log(1 - r0))
}
