#' Plot the GPS Prior Fit
#'
#' Plots the fit of the prior distribution of the GPS to
#' the observed number of reports
#'
#' @inheritParams BCPNN
#' @param prior List that contains the prior parameters. If not specified, automatically fitted to the data, 
#'              see \code{\link{fitPriorParametersGPS}}
#' @param bins Number of bins used 
#' @param precision The number of sample points used for plotting (Default = 100)
#' @param xmin Minimal x-value (Default = 0)
#' @param xmax Maximal x-value (Default = 10)
#'
#' @export
plotFitPriorGPS <- function(a, b, c, d,
                            prior = fitPriorParametersGPS(a, b, c, d),
                            bins = 50,
                            precision = 100,
                            xmin = 0,
                            xmax = 10) {
  
    alpha1 <- prior$alpha1
    beta1  <- prior$beta1
    alpha2 <- prior$alpha2
    beta2  <- prior$beta2
    w      <- prior$w
    
    E <- ((a + b) * (a + c)) / (a + b + c + d) # expected count
    
    # remove cases where E is zero
    a <- a[E != 0]
    E <- E[E != 0]
    
    # relative risk
    RR <- a / E
    
    x <- seq(xmin, xmax, len = precision)
    
    p <- ggplot2::qplot(x, geom = "blank") +
      ggplot2::stat_function(
        mapping = aes(x = x, y = ..y..),
        color = "blue",
        fun = dbigamma,
        n = precision,
        args = list(
          shape1 = alpha1,
          rate1 = beta1,
          shape2 = alpha1,
          rate2 = beta2,
          w = w
        )
      ) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0.01, 0)) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
      ggplot2::labs(x = "relative risk (RR)", y = "density", title = "GPS Prior") +
      ggplot2::geom_histogram(aes(x = RR, ..density..), bins = bins) +
      xlim(xmin, xmax)
    return(p)
  }