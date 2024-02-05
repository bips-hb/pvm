#' Test of the Poisson Mean
#'
#' Performs the test of the Poisson mean 
#' to a collection of 2 x 2 tables of the form
#' \tabular{lcc}{
#'    \tab event \tab not event\cr
#'   drug \tab `a` \tab `c`\cr
#'   not drug \tab `b` \tab `d`
#' }
#'
#' @template standardParams
#'
#' @return p-value
#' @export
PoissonTest <- function(a, b, c, d) {
  # to overcome possible integer overflow 
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
