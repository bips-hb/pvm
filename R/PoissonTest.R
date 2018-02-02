#' Test of the Poisson Mean
#'
#' Performs the test of the Poisson mean for the table
#' \tabular{lcc}{
#'    \tab event \tab not event\cr
#'   drug \tab \code{a} \tab \code{c}\cr
#'   not drug \tab \code{b} \tab \code{d}
#' }
#'
#' @inheritParams createTable
#'
#' @return p-value
#' @export
PoissonTest <- function(a, b, c, d) {
  # to overcome possible integer overflow later
  a <- as.numeric(a)
  b <- as.numeric(b)
  c <- as.numeric(c)
  d <- as.numeric(d) 

  E = ((a + b)*(a + c)) / (a + b + c + d)

  mapply(
    function(k, r) return(poisson.test(k, r = r, alternative = "greater")$p.value),
    k = a,
    r = E
  )
}
