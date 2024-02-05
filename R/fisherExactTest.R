#' Fisher's Exact Test
#'
#' Performs the (one-sided) Fisher's exact test 
#' to a collection of 2 x 2 tables of the form:
#' \tabular{lcc}{
#'    \tab event \tab not event\cr
#'   drug \tab `a` \tab `c`\cr
#'   not drug \tab `b` \tab `d`
#' }
#'
#' Wrapper function for the `Rcpp` functions `fishersTestGreater`
#' and `midPFishersTestGreater`.
#'
#' @template standardParams
#' @param midpvalue The mid-p-value correction (suggested by Agresti) is applied
#'
#' @return p-value
#' @references Ahmed, I., Dalmasso, C., Haramburu, F., Thiessard, F., 
#'             Bro\"et, P., & Tubert-Bitter, P. (2010). False Discovery 
#'             Rate Estimation for Frequentist Pharmacovigilance Signal 
#'             Detection Methods. Biometrics, 66(1), 301–309. 
#'             https://doi.org/10.1111/j.1541-0420.2009.01262.x
#' @export
fisherExactTest <- function(a, b, c, d, midpvalue = FALSE) {
  if (midpvalue) {
    mapply(
      function(a, b, c, d) midPFisherTestGreater(a, b, c, d),
      a, b, c, d
    )
  } else {
    mapply(
      function(a, b, c, d) fisherTestGreater(a, b, c, d),
      a, b, c, d
    )
  }
}
