#' Calculate Concordance Index (C-index) for Survival Predictions
#'
#' Computes the concordance index between predicted risk scores and observed survival times
#'
#' @param newdata Dataframe containing 'time' and 'status' columns
#' @param model Trained survival model object
#' @param model_name Model type from supported list:
#'        "Lasso", "Ridge", "Enet", "RFRSF", "GBM", "CoxBoost", "plsRcox",
#'        "XGBoost", "BlackBoost", "DeepHit", "DeepSurv", "SurvivalSVM"
#' @return survcomp::concordance.index object containing:
#'         - c.index: concordance index value
#'         - se: standard error
#'         - lower: lower CI
#'         - upper: upper CI
#'         - p.value: significance
#' @export
#' @examples
#' \dontrun{
#' data(train_data)
#' cindex_result <- cal_cindex(
#'   newdata = train_data,
#'   model = trained_model,
#'   model_name = "CoxBoost"
#' )
#' }
cal_cindex <- function(newdata, model, model_name) {
  # Validate input data
  if (is.null(newdata) || !all(c("time", "status") %in% colnames(newdata))) {
    stop("Input data must contain 'time' and 'status' columns")
  }

  # Check for required packages
  if (!requireNamespace("survcomp", quietly = TRUE)) {
    stop("Package 'survcomp' required. Install with: install.packages('survcomp')")
  }

  # Generate predictions
  pred_df <- tryCatch(
    cal_pred(newdata, model, model_name),
    error = function(e) {
      stop("Prediction failed: ", conditionMessage(e))
    }
  )

  # Calculate concordance index
  survcomp::concordance.index(
    x = pred_df$pred,
    surv.time = newdata$time,
    surv.event = newdata$status,
    method = "noether"
  )
}
