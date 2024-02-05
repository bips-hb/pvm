#' Convert Reports to 2 x 2 Tables
#'
#' Creates a data frame containing all 2 x 2 contingency tables
#' given a raw spontaneous reporting (SR) data set. An SR data set
#' is a binary matrix, where each row is a report. The first
#' columns represent the presence (`1`) or absence of a drug
#' (`0`), the other columns represent the presence or absence 
#' of an event.
#' \cr\cr
#' The tables are organized as follows:
#' \tabular{lccc}{
#'    \tab event \eqn{j} \tab not event \eqn{j} \tab *total*\cr
#'   drug \eqn{i} \tab `a` \tab `c` \tab `a` + `c`\cr
#'   not drug \eqn{i} \tab `b` \tab `d` \tab `b` + `d`\cr
#'   *total* \tab `a` + `b` \tab `c` + `d` \tab `n_reports`
#' }
#'
#' The code is a simplified version of the function `create2x2Tables`
#' in the `SRSim` package.
#'
#' @param reports A binary matrix. Each row is a report
#' @param n_drugs The number of drugs
#' @param n_events The number of events
#'
#' @return A data frame where each row represents a 2 x 2 table. The columns represent:
#'   \item{`drug_id`}{The ID of the drug}
#'   \item{`event_id`}{The ID of the event}
#'   \item{`a`}{Number of times the drug and event appeared together in a report}
#'   \item{`b`}{Number of times the event appeared without the drug in a report}
#'   \item{`c`}{Number of times the drug appeared without the event in a report}
#'   \item{`d`}{Number of times the drug and event both did not appear in a report}
#'
#' @seealso [convertRawReports2Tables()]
#' @export
convertRawReports2Tables <- function(reports, n_drugs, n_events) {
  return(
    convertRawReports2TablesRcpp(as.matrix(reports), n_drugs, n_events)
  )
}
