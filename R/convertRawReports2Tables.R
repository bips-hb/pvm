#' Convert Reports to 2 x 2 Tables
#'
#' Creates a data frame containing all 2 x 2 contingency tables
#' given a raw spontaneous reporting (SR) data set. An SR data set
#' is a binary matrix, where each row is a report. The first
#' columns represent the presence (\code{1}) or absence of a drug
#' (\code{0}), the other columns represent the presence or absence 
#' of an event.
#' \cr\cr
#' The tables are organized as follows:
#' \tabular{lccc}{
#'    \tab event \eqn{j} \tab not event \eqn{j} \tab \emph{total}\cr
#'   drug \eqn{i} \tab \code{a} \tab \code{c} \tab \code{a} + \code{c}\cr
#'   not drug \eqn{i} \tab \code{b} \tab \code{d} \tab \code{b} + \code{d}\cr
#'   \emph{total} \tab \code{a} + \code{b} \tab \code{c} + \code{d} \tab \code{n_reports}
#' }
#'
#' The code is a simplified version of the function \code{create2x2Tables}
#' in the \code{SRSim} package.
#'
#' @param reports A binary matrix. Each row is a report
#' @param n_drugs The number of drugs
#' @param n_events The number of events
#'
#' @return A data frame where each row represents a 2 x 2 table. The columns represent:
#'   \item{\code{drug_id}}{The ID of the drug}
#'   \item{\code{event_id}}{The ID of the event}
#'   \item{\code{a}}{Number of times the drug and event appeared together in a report}
#'   \item{\code{b}}{Number of times the event appeared without the drug in a report}
#'   \item{\code{c}}{Number of times the drug appeared without the event in a report}
#'   \item{\code{d}}{Number of times the drug and event both did not appear in a report}
#'
#' @seealso \code{\link{convertRawReports2Tables}}
#' @export
convertRawReports2Tables <- function(reports, n_drugs, n_events) {
  return(
    convertRawReports2TablesRcpp(as.matrix(reports), n_drugs, n_events)
  )
}
