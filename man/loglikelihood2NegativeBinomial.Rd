% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loglikelihood2NegativeBinomial.R
\name{loglikelihood2NegativeBinomial}
\alias{loglikelihood2NegativeBinomial}
\title{Log-Likelihood of the Bimodal Negative Binomial}
\usage{
loglikelihood2NegativeBinomial(p, a, E)
}
\arguments{
\item{p}{A vector with the parameters (\code{alpha1}, \code{beta1}, 
\code{alpha2}, \code{beta2} and \code{w}, in that order)}

\item{a}{A vector with the number of reports for each of the drug-event pairs}

\item{E}{A vector (of the same length as \code{a}) with the number of reports 
one would expect under the assumption of 'independence'}
}
\value{
The negative log-likelihood (i.e., -1 * log-likelihood)
}
\description{
Returns the log-likelihood of the bimodal negative binomial model 
used by the Gamma Poisson shrinker (GPS), see function 
\code{\link{fitPriorParametersGPS}}. The function is written such 
that it can be used by the base function \code{\link{nlminb}}.
}
\examples{
alpha1 <- 0.2
beta1 <- 0.06
alpha2 <- 1.4
beta2 <- 1.8
w <- 0.1

a <- c(5, 1, 56, 3)
E <- c(3.4, 0.5, 10, 0.5) 

p <- c(alpha1, beta1, alpha2, beta2, w)
loglikelihood2NegativeBinomial(p, a, E)
#[1] 16.80512

}
\references{
DuMouchel, W. (1999). Bayesian Data Mining in Large Frequency Tables, 
            with an Application to the FDA Spontaneous Reporting System. 
            The American Statistician, 53(3), 177–190. 
            https://doi.org/10.1080/00031305.1999.10474456
            
            DuMouchel, W., & Pregibon, D. (2001). Empirical bayes screening 
            for multi-item associations. Proceedings of the Seventh ACM 
            SIGKDD International Conference on Knowledge Discovery and 
            Data Mining - KDD ’01, (October), 67–76. 
            http://doi.org/10.1145/502512.502526
}
\seealso{
\code{\link{GPS}}, \code{\link{fitPriorParametersGPS}}, \code{\link{dbinbinom}}
}
