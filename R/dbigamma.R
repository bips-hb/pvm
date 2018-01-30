#' Bimodal Ganmma Distribution
#'
#' @param x The x-values
#' @param shape1,rate1 The shape and rate parameters for the first mode
#' @param shape2,rate2 The shape and rate parameters for the second mode
#' @param w The weight of the first mode (must lie in \eqn{[0,1]})
#'
#' @return The density for the values in x
#' @export
dbigamma<- function(x, shape1, rate1, shape2, rate2, w) {
  w * dgamma(x, shape = shape1, rate = rate1) +
    (1 - w) * dgamma(x, shape = shape2, rate = rate2) 
}
