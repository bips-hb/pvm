#' Bayesian Confidence Propagation Neural Network (BCPNN)
#'
#' Applies the BCPNN to a collection of 2 x 2 tables of the form
#' \tabular{lcc}{
#'    \tab event \tab not event\cr
#'   drug \tab `a` \tab `c`\cr
#'   not drug \tab `b` \tab `d`
#' }
#' There are two versions of the BCPNN: 
#' \itemize{
#'   \item{`'original'` - The original version proposed by Bate et al. (1998)}
#'   \item{`'alternative'` - The BCPNN as proposed by Norén et al. (2006)}
#' }
#' 
#' The implementation of this function is based on the implementation in the
#' `PhViD` package. 
#' 
#' @template standardParams
#' @template alphaParam
#' @param version Version of the BCPNN used. Can either be `'original'` (Default) 
#'                for the BCPNN as proposed orignally by Bate et al. (1998), or 
#'                `'alternative'` for the BCPNN as proposed by 
#'                Norén et al. (2006).
#' @param mc_estimate The value is estimated using Monte Carlo runs (Default = `FALSE`). 
#'                    Only used when `version = 'alternative'`.
#' @param mc_runs The number of Monte Carlo runs used to estimate the credible interval. 
#'                (Default: 1000). Only used when `version = 'alternative'`.
#'
#' @return The maximum aposteriori estimate of the information component (IC) or 
#'         the lower endpoint of the approximate credible interval
#'         
#' @importFrom MCMCpack rdirichlet
#'         
#' @references Bate, A., Lindquist, M., Edwards, I. R., Olsson, S., Orre, R., 
#'             Lansner, A., & De Freitas, R. M. (1998). A Bayesian neural 
#'             network method for adverse drug reaction signal generation. 
#'             European Journal of Clinical Pharmacology, 54(4), 315–321. 
#'             http://doi.org/10.1007/s002280050466
#' 
#'             Norén, G. N., Bate, A., Orre, R., & Edwards, I. R. (2006). 
#'             Extending the methods used to screen the WHO drug safety 
#'             database towards analysis of complex associations and 
#'             improved accuracy for rare events. Statistics in Medicine, 
#'             25(21), 3740–3757. http://doi.org/10.1002/sim.2473
#' 
#' @examples 
#' # get the tables
#' a <- srdata$tables$a
#' b <- srdata$tables$b
#' c <- srdata$tables$c
#' d <- srdata$tables$d
#' 
#' # Applying the original BCPNN: 
#' BCPNN(a, b, c, d)
#' # [1]  0.349783103 -0.609077730 -0.168446711 -0.277981964 ...
#' 
#' # Getting the lower end point of the 95% confidence intervaL: 
#' BCPNN(a, b, c, d, alpha = 0.05) 
#' # [1]  0.280077253 -0.994960076 -0.293624528 -0.408661852 ...
#' 
#' # Using the alternative version: 
#' BCPNN(a, b, c, d, version = 'alternative')
#' # [1]  0.350235800 -0.595807902 -0.166901050 -0.276387348 ...
#' 
#' # Getting the lower end points of the 95% confidence interval 
#' # using the alternative version. The estimates are based on 
#' # 10,000 Monte Carlo samples:
#' BCPNN(a, b, c, d, version = 'alternative', 
#'       alpha = 0.05, mc_estimate = TRUE, mc_runs = 10^4)
#' # [1]  [1]  0.31621489 -0.92490130 -0.25601307 -0.37040303 ...
#' @export
BCPNN <- function(a, b, c, d, alpha = NULL, 
                  version = 'original', mc_estimate = FALSE, mc_runs = 1000) {
  if (version == 'original') {
    n1. <- a + c
    n.1 <- a + b
    n   <- a + b + c + d
    p1  <- 1 + n1.
    p2  <- 1 + n - n1.
    q1  <- 1 + n.1
    q2  <- 1 + n - n.1
    r1  <- a + 1
    r2b <- n - a - 1 + (2 + n) ^ 2 / (q1 * p1)
    IC  <- (digamma(r1) - digamma(r1 + r2b) -
              (digamma(p1) - digamma(p1 + p2) + digamma(q1) -
                  digamma(q1 + q2))) / log(2)
    VICb <- (trigamma(r1) - trigamma(r1 + r2b) + 
               (trigamma(p1) - trigamma(p1 + p2) + trigamma(q1) -
                    trigamma(q1 + q2))) / log(2)^2
    if (is.null(alpha)) {
      return(IC)
    } else {
      return(stats::qnorm(alpha, IC, sqrt(VICb)))
    }
  }
  
  ### The alternative version of the BCPNN ---------------
  if (!mc_estimate) {
    if (!is.null(alpha)) {
      if (!near(alpha, .025)) {
        stop("the lower end point of the CI can only be approximated when alpha = .025. Otherwise, use MC (set mc_estimate = TRUE)")
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
        p <- MCMCpack::rdirichlet(mc_runs, c(gamma11,
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
      # use the approximation proposed by Noren et al. (See the appendix)
      
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
