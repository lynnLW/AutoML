#' Calculate Survival Risk Scores and Metrics (C-index, IBS, AUC)
#'
#' Computes comprehensive survival analysis metrics including risk predictions,
#' concordance index, integrated Brier score, and time-dependent AUC values.
#'
#' @param data Dataframe containing 'time' and 'status' columns
#' @param model Trained survival model object
#' @param model_name Model type from supported list:
#'        "Lasso", "Ridge", "Enet", "RFRSF", "GBM", "CoxBoost", "plsRcox",
#'        "XGBoost", "BlackBoost", "DeepHit", "DeepSurv", "SurvivalSVM"
#' @return List containing:
#'         - pred_df: Full prediction dataframe
#'         - pred_coxb: Risk scores
#'         - cindex_coxb: C-index object
#'         - cindex: Numeric C-index value
#'         - bs: Integrated Brier Score
#'         - km_auc_X: AUC values for X years (1,2,3,5,7,10)
#' @export
#' @examples
#' \dontrun{
#' data(example_data)
#' results <- cal_metrics(
#'   data = example_data,
#'   model = trained_model,
#'   model_name = "CoxBoost"
#' )
#' }
cal_metrics <- function(data, model, model_name) {
  # Early return for invalid data
  if (is.null(data) || nrow(data) == 0 || sum(data$status) == 0) {
    km_auc_names <- paste0("km_auc_", c(1, 2, 3, 5, 7, 10))
    km_auc_values <- setNames(as.list(rep(NA, 6)), km_auc_names)
    return(c(
      list(
        pred_df = NA,
        pred_coxb = NA,
        cindex_coxb = NA,
        cindex = NA,
        bs = NA
      ),
      km_auc_values
    ))
  }

  # Calculate all metrics
  pred_df <- cal_pred(data, model, model_name)
  cindex_obj <- cal_cindex(data, model, model_name)
  bs_score <- cal_bs(data, model, model_name)
  auc_values <- cal_multi_auc(data, model, model_name)

  # Return combined results
  c(
    list(
      pred_df = pred_df,
      pred_coxb = pred_df$pred,
      cindex_coxb = cindex_obj,
      cindex = cindex_obj$c.index,
      bs = bs_score
    ),
    auc_values
  )
}
