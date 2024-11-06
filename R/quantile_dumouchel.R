#' Utilities copied from PhViD
#'
#' Package was archived on CRAN and therefore broke ours,
#' but luckily the required imported functions where easy top transfer.
#' There is unfortunately no documentation, however.
#' PhViD::.FCoutQuantileDuMouchel
#' PhViD::.QuantileDuMouchel
#'
#' @noRd
#' @keywords internal
fcout_quantile_DuMouchel <- function (p, Seuil, Q, a1, b1, a2, b2)
{
  Q * pgamma(p, shape = a1, rate = b1) +
    (1 - Q) *
    pgamma(p, shape = a2, rate = b2) -
    Seuil
}

#' @noRd
#' @keywords internal
quantile_DuMouchel <- function (Seuil, Q, a1, b1, a2, b2) {
  m <- rep(-1e+05, length(Q))
  M <- rep(1e+05, length(Q))
  x <- rep(1, length(Q))
  Cout <- fcout_quantile_DuMouchel(x, Seuil, Q, a1, b1, a2,
                                  b2)
  while (max(round(Cout * 10000)) != 0) {
    S <- sign(Cout)
    xnew <- (1 + S)/2 * ((x + m)/2) + (1 - S)/2 * ((M + x)/2)
    M <- (1 + S)/2 * x + (1 - S)/2 * M
    m <- (1 + S)/2 * m + (1 - S)/2 * x
    x <- xnew
    Cout <- fcout_quantile_DuMouchel(x, Seuil, Q, a1, b1,
                                    a2, b2)
  }
  x
}
