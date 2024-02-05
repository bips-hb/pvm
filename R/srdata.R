#' A Simulated Spontaneous Reporting System
#'
#' A simulated spontaneous reporting data set generated with the
#' \code{SRSim} simulator. The data set contains 10,000 reports
#' for 10 drugs and 10 adverse events (AEs). Five drug-AE pairs 
#' are associated with an odds ratio of 2. All other drug-AE pairs
#' have an odds ratio of 1. Five drugs are innocent bystanders, i.e.,
#' they are prescribed together with one other drug, but they do not 
#' cause any adverse events. 
#' \cr\cr
#' The marginal probabilities over the drugs and the AEs were drawn 
#' from a Beta distribution with parameters \eqn{\alpha = 1.0} and 
#' \eqn{\beta = 20.0}.
#' \cr\cr
#' The conditional probability of an innocent bystander given that
#' the other drug is prescribed is set to .9 (this is regulated with 
#' the argument \code{bystander_prob}). 
#' \cr\cr
#' The following commands were used for generating the data set:
#' \cr\cr
#' \code{
#' library(SRSim)
#' srdata <- SRSim::simulateSRS(n_reports = 10000,
#'                             n_drugs = 10,
#'                             n_events = 10,
#'                             n_innocent_bystanders = 5,
#'                             bystander_prob = 0.9,
#'                             n_correlated_pairs = 5,
#'                             theta = 2,
#'                             seed = 1)
#'  # create the 2x2 tables
#'  
#'  srdata$tables <- SRSim::convert2Tables(srdata)
#' }
#'
#' @format \code{srdata} contains the following elements:
#' \describe{
#'   \item{sr}{A binary data frame with 10,000 rows and 20 columns. The first
#'            10 columns represent the drugs; the latter represent the events.
#'            Each row is a report. In case of a 1, the drug/event has been reported,
#'            zero otherwise. The column names are
#'            \code{drug1} till \code{drug10} and \code{event1} till \code{event10}.
#'   }
#'   \item{dag}{The directed acycled graph as an \code{igraph} object}
#'         \item{nodes}{A tibble with all the information on each node/variate:
#'             \itemize{
#'                \item{\code{label}}{ The label for each node/variate}
#'                \item{\code{in_degree}}{ The number of edges pointing to the node}
#'                \item{\code{id}}{ The ID of each node (simple integer)}
#'                \item{\code{parent_id}}{ The ID of the parent node - if any. Otherwise equal to \code{-1}}
#'                \item{\code{margprob}}{ The marginal probability of the node/variate}
#'                \item{\code{beta0}}{ The intercept in the logistic regression model for that node}
#'                \item{\code{beta1}}{ The regression coefficient in the logistic regression model for the parent}
#'             }
#'         }
#'   \item{prob_drugs}{A vector with marginal probabilities of the drugs}
#'   \item{prob_events}{A vector with marginal probabilities of the events}
#'   \item{tables}{A data frame with 100 rows. Each row contains the data on a drug-event pair.
#'     The columns represent:
#'     \describe{
#'       \item{\code{drug_id}}{The ID of the drug}
#'       \item{\code{event_id}}{The ID of the event}
#'       \item{\code{prob_drug}}{The marginal probability of that drug}
#'       \item{\code{prob_event}}{The marginal probability of that event}
#'       \item{\code{or}}{The odds ratio}
#'       \item{\code{associated}}{\code{TRUE} is there is a non-zero correlation, \code{FALSE} otherwise}
#'       \item{\code{a}}{Number of times the drug and event appeared together in a report}
#'       \item{\code{b}}{Number of times the event appeared without the drug in a report}
#'       \item{\code{c}}{Number of times the drug appeared without the event in a report}
#'       \item{\code{d}}{Number of times the drug and event both did not appear in a report}
#'    }
#'  }
#' }
"srdata"