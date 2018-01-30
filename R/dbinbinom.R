#' Bimodal Negative Binomial
#'
#' @param x The x-values
#' @param size1,shape1 The size and shape parameters for the first mode
#' @param size2,shape2 The size and shape parameters for the second mode
#' @param w The weight of the first mode (must lie in \eqn{[0,1]})
#'
#' @return The density for the values in x
#' @export
dbinbinom <- function(x, size1, prob1, size2, prob2, w) {
  w * dnbinom(x, size = size1, prob = prob1) +
    (1 - w) * dnbinom(x, size = size2, prob = prob2)
}
