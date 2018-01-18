#' The Alternative Bayesian Confidence Propagation Neural Network (BCPNN)
#'
#' Applies the BPCNN as proposed by Noren et al. (2006) to a 2 x 2 table of the form
#' \tabular{lcc}{
#'    \tab event \tab not event\cr
#'   drug \tab \code{a} \tab \code{c}\cr
#'   not drug \tab \code{b} \tab \code{d}
#' }
#'
#' The implementation of the algorithm is based on:
#'
#' Norén, G. N., Bate, A., Orre, R., & Edwards, I. R. (2006).
#' Extending the methods used to screen the WHO drug safety database
#' towards analysis of complex associations and improved accuracy for
#' rare events. Statistics in Medicine, 25(21), 3740–3757. http://doi.org/10.1002/sim.2473
#'
#' The notation used within the code is kept as close to the notation used in the
#' paper.
#'
#' @inheritParams ROR
#' @param mc_estimate The value is estimated using Monte Carlo runs (Default = \code{FALSE})
#' @param mc_runs The number of Monte Carlo runs used to estimate the credible interval. NOTE: might be computationally expensive. (Default: 1000)
#'
#' @return The maximum aposteriori estimate of the information component or the lower endpoint of the approximate credible interval
#' @export
BCPNN2 <- function(a, b, c, d, alpha = NULL, mc_estimate = FALSE, mc_runs = 1000) {
  if (!mc_estimate) {
    if (!is.null(alpha)) {
      if (!dplyr::near(alpha, .025)) {
        stop("ERROR BCPNN2: the lower end point of the CI can only be approximated when alpha = .025. Otherwise, use MC (set mc_estimate = TRUE)")
      }
    }
  }

  # NOTE: we could speed up the code by combining some of the expressions
  # here. We decided to keep it like this for readability's sake.

  # determine the marginals and the total
  n1. <- a + c
  n0. <- b + d
  n.1 <- a + b
  n.0 <- c + d
  n   <- a + b + c + d

  # determine the q-values (see eq. (5) in the paper by Noren et al.)
  q1. <- (n1. + .5) / (n + 1)
  q0. <- (n0. + .5) / (n + 1)
  q.1 <- (n.1 + .5) / (n + 1)
  q.0 <- (n.0 + .5) / (n + 1)

  # determine the prior parameters (see eq. (3) & (4) in the paper by Noren et al.)
  alpha.. <- 0.5 / (q1. * q.1)
  alpha11 <- q1. * q.1 * alpha..
  alpha10 <- q1. * q.0 * alpha..
  alpha01 <- q0. * q.1 * alpha..
  alpha00 <- q0. * q.0 * alpha..

  # determine the gammas on the basis of the alphas and the observed counts
  # see eq. (2) in the paper by Noren et al.
  gamma11 <- alpha11 + a
  gamma01 <- alpha01 + b
  gamma10 <- alpha10 + c
  gamma00 <- alpha00 + d
  gamma.. <- gamma11 + gamma10 + gamma01 + gamma00

  if (mc_estimate) {
    mapply(
      function(gamma11, gamma10, gamma01, gamma00) {
      # sample from the posterior distribution
      p <- rdirichlet(mc_runs, c(gamma11,
                                 gamma10,
                                 gamma01,
                                 gamma00))

      # compute the corresponding ICs for the samples from the posterior distribution
      p11 <- p[, 1]
      p1. <- p11 + p[, 2]
      p.1 <- p11 + p[, 3]

      IC <- log(p11 / (p1. * p.1), base = 2)

      if (is.null(alpha)) {
        return(mean(IC)) # return the MAP estimate of IC
      } else {
        return(sort(IC)[max(1, round(mc_runs * alpha))]) # return the estimated lower end points of the CI
      }
      },
      gamma11,
      gamma10,
      gamma01,
      gamma00)


  } else { # approximations are used

    # determine the posterior expectations of p11, p1. and p.1
    # (following directly from eq. (2) in the paper by Noren et al.)
    exp_p11 <- gamma11 / gamma..
    exp_p1. <- (gamma11 + gamma10) / gamma..
    exp_p.1 <- (gamma11 + gamma01) / gamma..

    # obtain the MAP estimate of the information component, see eq. (7)
    IC <- log(exp_p11 / (exp_p1. * exp_p.1), base = 2)

    if (is.null(alpha)) {
      return(IC)
    } else {
      # use the approximation proposed by Noren et al. See eq. () and the appendix
      # for more details

      # perform linear interpolation using a list of known values
      r_values <- seq(0.0, 1.0, by = 0.1)
      A <- c(3.09, 2.93, 2.78, 2.62, 2.45, 2.25, 2.03, 1.79, 1.61, 1.13, 0.073)
      B <- c(2.22, 2.27, 2.26, 2.25, 2.15, 2.12, 2.05, 1.93, 1.89, 1.15, -0.081)

      Delta025 <- mapply(
        function(gamma11, gamma10, gamma01, gamma00) {
          # first establish the r value
          r <- gamma11 / min(gamma10 + gamma11, gamma01 + gamma11)

          # get the lower values
          lower <- (r_values > r - 0.1) & (r_values <= r)
          r_lower <- sum(r * lower * 1)
          A_lower <- sum(A * lower * 1)
          B_lower <- sum(B * lower * 1)

          # get the upper values
          upper <- (r_values > r) & (r_values < r + 0.1)
          r_upper <- sum(r * upper * 1)
          A_upper <- sum(A * upper * 1)
          B_upper <- sum(B * upper * 1)

          # interpolate A and B
          A <- A_lower + ((r - r_lower) / 0.1) * (A_upper - A_lower)
          B <- B_lower + ((r - r_lower) / 0.1) * (B_upper - B_lower)

          # determine the difference between the IC estimate and the lower end point (eq. (9) in the paper)
          Delta025 <- A * gamma11^(-.5) + B * gamma11^(-1.5)

          # return the approximate lower end point of the credible interval (eq (8) in the paper)
          return(Delta025)
        },
        gamma11,
        gamma10,
        gamma01,
        gamma00)
      return(IC - Delta025)
    }
  }
}
