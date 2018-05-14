#' Elastic Net
#'
#' Applies the Elastic Net to raw spontaneous report data.
#' Every event is regressed on all the drugs in the data set.
#' The function returns a data frame with every drug-event pair
#' and the estimated regression coefficient.
#' \cr\cr
#' In case there are not enough observations of an event (the
#' event must appear at least twice), the regression is not
#' performed. All the regression estimates for the drugs and that
#' particular event are set to \code{0}. The entries in the \code{lambda}
#' column of the data frame are set to \code{NA}.
#' \cr\cr
#' \strong{Shrinkage parameter}\cr
#' One can set the shrinkage parameter with the argument \code{lambda}
#' in a number of ways:
#' \enumerate{
#'   \item \code{lambda = NULL} (\emph{Default}). The parameter is set
#'   through cross-validation. The number of folds can be set with
#'   \code{nfolds} (Default = 10). The loss function used can be set
#'   with \code{type.measure} (Default = \code{deviance}). See for
#'   other \code{type.measure} options the function \code{glmnet::cv.glmnet}.
#'   The \code{glmnet::cv.glmnet} function returns two estimates:
#'   \code{lambda.min} and \code{lambda.1se}. To use the former, set
#'   \code{lambda.type} to \code{"min"} (default). For the latter, type
#'   \code{"1se"}.
#'   \item Set to one value, e.g., \code{lambda = 0.5}. The same shrinkage parameter
#' is used for all events.
#'   \item A vector of length \code{n_events}, e.g., \code{lambda = c(0.5, 0.8, 1)}. The
#' shrinkage parameters are specified for each event individually.
#' }
#'
#' @param reports A binary matrix, where each row represents a report
#' @param n_drugs,n_events The number of drugs and events
#' @param lambda Shrinkage parameter. Can be a list of length \code{n_events}. When not set, estimated through cross-validation
#' @param nfolds Number of folds used for cross-validation
#' @param type.measure Loss function used (Default: \code{deviance}). See for more options \code{glmnet::cv.glmnet}
#' @param lambda.type Type of estimate that is used (either \code{"min"} - default - or \code{"1se"})
#' @param alpha The elastic net mixing parameter (Default: 0.5)
#' @param event_ids IDs of the events evaluated (Default: all)
#' @param verbose Verbosity (Default: \code{FALSE})
#'
#' @return A data frame with the columns
#'         \item{drug_id}{ID for the drug (simply numbered 1,2,3,...etc.)}
#'         \item{event_id}{ID for the event (simply numbered 1,2,3,...etc.)}
#'         \item{lambda}{The shrinkage parameter \eqn{\lambda} that was used for this pair.
#'                       In case the regression was not performed (because the event was not observed
#'                       or only observed once), the entry is \code{NA}}
#'         \item{LASSO}{The regression parameter after regressing all drugs to the event in question}
#'         \item{alpha}{The alpha value used for each case}
#' @export
ElasticNet <- function(reports,
                       n_drugs,
                       n_events,
                       lambda = NULL,
                       nfolds = 10,
                       type.measure = "deviance",
                       lambda.type = "min",
                       alpha = 0.5,
                       event_ids = 1:n_events,
                       verbose = FALSE) {
  
  
  res <- LASSO(reports = reports, 
               n_drugs = n_drugs, 
               n_events = n_events, 
               lambda = lambda, 
               nfolds = nfolds,
               type.measure = type.measure,
               lambda.type = lambda.type, 
               alpha = alpha, 
               event_ids = event_ids, 
               verbose = verbose)
  
  res$alpha = alpha 
  
  return(res)
}
