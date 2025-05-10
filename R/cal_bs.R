#' Calculate Integrated Brier Score (IBS) for New Dataset
#'
#' Computes Integrated Brier Score using trained survival analysis models
#'
#' @param newdata Dataframe containing time-to-event data with 'time' and 'status' columns
#' @param model Trained survival analysis model object
#' @param model_name Model specification. Supported models:
#'   "Lasso", "Ridge", "Enet", "RFRSF", "GBM", "CoxBoost", "plsRcox",
#'   "XGBoost", "BlackBoost", "DeepHit", "DeepSurv", "SurvivalSVM"
#' @return Numeric Integrated Brier Score (IBS), where lower values indicate better performance
#' @export
#' @examples
#' \donttest{
#' # Requires pre-trained model
#' library(survival)
#' data(example_data)
#' cox_model <- coxph(Surv(time, status) ~ ., data = example_data)
#' ibs_score <- calculate_ibs(
#'   newdata = example_data,
#'   model = cox_model,
#'   model_name = "CoxBoost"
#' )
#' }
cal_bs <- function(newdata, model, model_name) {
  # Parameter validation ----------------------------------------------------
  required_cols <- c("time", "status")
  if (!all(required_cols %in% colnames(newdata))) {
    stop("Input data must contain 'time' and 'status' columns")
  }

  valid_models <- c("Lasso", "Ridge", "Enet", "RFRSF", "GBM", "CoxBoost",
                    "plsRcox", "XGBoost", "BlackBoost", "DeepHit",
                    "DeepSurv", "SurvivalSVM","GLMBoost","SuperPC")
  model_name <- match.arg(model_name, choices = valid_models)

  if (is.null(model)) stop("Must provide valid trained model")

  # Generate predictions ----------------------------------------------------
  pred_df <- tryCatch(
    expr = cal_pred(newdata, model, model_name),
    error = function(e) {
      message("Prediction failed: ", conditionMessage(e))
      return(NULL)
    }
  )

  if (is.null(pred_df)) return(NA_real_)

  # Prepare evaluation data -------------------------------------------------
  eval_data <- data.frame(
    time = as.numeric(newdata$time),
    event = as.numeric(newdata$status),
    score = as.numeric(pred_df$pred)
  )

  # Calculate IBS -----------------------------------------------------------
  if (!requireNamespace("survcomp", quietly = TRUE)) {
    stop("Package survcomp required: install.packages('survcomp')")
  }

  survcomp::sbrier.score2proba(
    data.tr = eval_data,
    data.ts = eval_data,
    method = "cox"
  )$bsc.integrated
}
