#' Calculate Time-Dependent AUC for Multiple Time Points
#'
#' Computes area under the ROC curve (AUC) at specified time points for survival models
#'
#' @param newdata Data frame containing test data with 'time' and 'status' columns
#' @param model Trained survival model object
#' @param model_name Model type from supported list:
#'   c("RSF", "Enet", "StepCox","CoxBoost","plsRcox","superpc","GBM","survivalsvm","Ridge","Lasso")
#' @param cuts Numeric vector of time points to evaluate (default: c(1,2,3,5,7,10))
#' @param unit Time unit: "m" (months), "d" (days), or "y" (years)
#' @return List containing AUC values, false positive rates, and true positive rates
#' @export
#' @examples
#' \donttest{
#' data(example_data)
#' model <- ML.survival.model(example_data)
#' auc_results <- calculate_multi_auc(
#'   newdata = example_data,
#'   model = model,
#'   model_name = "CoxBoost",
#'   cuts = c(2,5),
#'   unit = "y"
#' )
#' }
cal_multi_auc <- function(newdata, model, model_name, cuts = NULL, unit = "m") {
  # Parameter Validation ----------------------------------------------------
  ## Check required columns
  required_cols <- c("time", "status")
  if (!all(required_cols %in% colnames(newdata))) {
    stop("Input data must contain 'time' and 'status' columns")
  }

  ## Validate model name
  valid_models <- c("Lasso", "Ridge", "Enet", "RFRSF", "GBM", "CoxBoost",
                    "plsRcox", "XGBoost", "BlackBoost", "DeepHit",
                    "DeepSurv", "SurvivalSVM","GLMBoost","SuperPC")
  model_name <- match.arg(model_name, valid_models)

  ## Validate time unit
  unit <- match.arg(unit, c("m", "d", "y"))

  ## Set default time points
  if (is.null(cuts)) {
    cuts <- c(1, 2, 3, 5, 7, 10)
  } else {
    stopifnot(is.numeric(cuts))
  }

  # Dependency Management ---------------------------------------------------
  if (!requireNamespace("survivalROC", quietly = TRUE)) {
    stop("Please install survivalROC: install.packages('survivalROC')")
  }

  # Generate Predictions ----------------------------------------------------
  pred_df <- tryCatch(
    cal_pred(newdata, model, model_name),
    error = function(e) {
      message("Prediction failed: ", conditionMessage(e))
      return(NULL)
    }
  )

  if (is.null(pred_df)) return(list())

  # Convert Time Units ------------------------------------------------------
  time_conversion <- switch(unit,
                            "m" = 12,   # months to years
                            "d" = 365,  # days to years
                            "y" = 1     # years
  )

  # Core Calculation Function -----------------------------------------------
  calculate_time_auc <- function(cutoff_time) {
    ## Convert to years
    converted_cutoff <- cutoff_time*time_conversion

    ## Check data validity
    if (max(pred_df$time) < converted_cutoff) {
      return(list(
        fp = NA_real_,
        tp = NA_real_,
        auc = NA_real_
      ))
    }

    ## Calculate ROC
    roc_res <- survivalROC::survivalROC(
      Stime = pred_df$time,
      status = pred_df$status,
      marker = pred_df$pred,
      predict.time = converted_cutoff,
      method = "KM"
    )

    list(
      fp = roc_res$FP,
      tp = roc_res$TP,
      auc = roc_res$AUC
    )
  }

  # Main Calculation --------------------------------------------------------
  results<-list()
  for (i in cuts){
    res <- calculate_time_auc(i)
    results[[paste0("km_fp_", i)]] = res$fp
    results[[paste0("km_tp_", i)]]= res$tp
    results[[paste0("km_auc_", i)]] = res$auc
  }

  # Format Output -----------------------------------------------------------
  results
}
