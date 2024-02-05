#' LASSO
#'
#' Applies the LASSO to raw spontaneous report data.
#' Every event is regressed on all the drugs in the data set.
#' The function returns a data frame with every drug-event pair
#' and the estimated regression coefficient.
#' \cr\cr
#' In case there are not enough observations of an event (the
#' event must appear at least twice), the regression is not
#' performed. All the regression estimates for the drugs and that
#' particular event are set to `0`. The entries in the `lambda`
#' column of the data frame are set to `NA`.
#' \cr\cr
#' **Shrinkage parameter**\cr
#' One can set the shrinkage parameter with the argument `lambda`
#' in a number of ways:
#' \enumerate{
#'   \item `lambda = NULL` (*Default*). The parameter is set
#'   through cross-validation. The number of folds can be set with
#'   `nfolds` (Default = 10). The loss function used can be set
#'   with `type.measure` (Default = `deviance`). See for
#'   other `type.measure` options the function `glmnet::cv.glmnet`.
#'   The `glmnet::cv.glmnet` function returns two estimates:
#'   `lambda.min` and `lambda.1se`. To use the former, set
#'   `lambda.type` to `"min"` (default). For the latter, type
#'   `"1se"`.
#'   \item Set to one value, e.g., `lambda = 0.5`. The same shrinkage parameter
#' is used for all events.
#'   \item A vector of length `n_events`, e.g., `lambda = c(0.5, 0.8, 1)`. The
#' shrinkage parameters are specified for each event individually.
#' }
#'
#' @param reports A binary matrix, where each row represents a report
#' @param n_drugs,n_events The number of drugs and events
#' @param lambda Shrinkage parameter. Can be a list of length `n_events`. When not set, estimated through cross-validation
#' @param nfolds Number of folds used for cross-validation
#' @param type.measure Loss function used (Default: `deviance`). See for more options `glmnet::cv.glmnet`
#' @param lambda.type Type of estimate that is used (either `"min"` - default - or `"1se"`)
#' @param alpha The elastic net mixing parameter (Default: 1.0 - LASSO)
#' @param event_ids IDs of the events evaluated (Default: all)
#' @param verbose Verbosity (Default: `FALSE`)
#'
#' @return A data frame with the columns
#'         \item{drug_id}{ID for the drug (simply numbered 1,2,3,...etc.)}
#'         \item{event_id}{ID for the event (simply numbered 1,2,3,...etc.)}
#'         \item{lambda}{The shrinkage parameter \eqn{\lambda} that was used for this pair.
#'                       In case the regression was not performed (because the event was not observed
#'                       or only observed once), the entry is `NA`}
#'         \item{LASSO}{The regression parameter after regressing all drugs to the event in question}
#' @export
LASSO <- function(reports,
                  n_drugs,
                  n_events,
                  lambda = NULL,
                  nfolds = 10,
                  type.measure = "deviance",
                  lambda.type = "min",
                  alpha = 1.0,
                  event_ids = 1:n_events,
                  verbose = FALSE) {

  # process the arguments
  if (!is.null(lambda)) { # in case lambda is set
    lambda_set <- TRUE
    if (length(lambda) == 1) {
      lambda <- rep(lambda, n_events)
    }
    if (length(lambda) != n_events) {
      stop("The length of the lambda vector must equal the number of events\n")
    }
  } else {
    lambda_set <- FALSE
  }

  n_reports <- nrow(reports)

  res <- as.data.frame(
    expand.grid(
      drug_id = 1:n_drugs,
      event_id = 1:n_events
    )
  )

  # initialize the LASSO and lambda column
  res$lambda <- NA
  res$LASSO <- 0

  # get only the drugs
  x <- as.matrix( reports[, 1:n_drugs] )

  # walk over all events
  for (event_id in event_ids) {
    y <- unlist( reports[, event_id + n_drugs])  # get the vector for the event

    if (verbose) {
      cat("Fitting event no.", event_id, '\n')
    }

    # check whether there enough observations
    if (sum(y) > 1 & sum(y) < n_reports - 1) {
      if (!lambda_set) { # estimate lambda in case it is not given
        suppressWarnings(
          est <- glmnet::cv.glmnet(x = x,
                                   y = y,
                                   family = "binomial",
                                   nfolds = nfolds,
                                   alpha = alpha,
                                   type.measure = type.measure)
          )
        if (lambda.type == "min") {
          shrinkage_parameter <- est$lambda.min
        } else {
          shrinkage_parameter <- est$lambda.1se
        }
      } else {
        shrinkage_parameter <- lambda[event_id]
      }

      suppressWarnings(
        fit <- glmnet::glmnet(x, y, lambda = shrinkage_parameter, alpha = alpha, family = "binomial")
      )
      beta <- coef(fit)[-1] # get all the parameters (but discard the intercept)

      # add the betas to the data frame
      res[res$event_id == event_id, ]$LASSO = beta
      res[res$event_id == event_id, ]$lambda = shrinkage_parameter
    }
  }

  return(res)
}
