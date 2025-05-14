#' Determine Optimal Cutoff for Survival Analysis
#'
#' @description Identifies optimal risk stratification cutoff using maximally selected rank statistics
#' @param data Dataframe containing survival data
#' @param time Column name specifying observed time (quoted string)
#' @param status Column name specifying event status (quoted string, 0/1 format)
#' @param pred Column name containing prediction scores (quoted string)
#' @return Optimal cutoff value (numeric)
#' @export
optimal_cutoff <- function(data, time, status, pred) {

  required_cols <- c(time, status, pred)
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!all(data[[status]] %in% c(0, 1))) {
    stop("Status variable must be binary (0/1)")
  }

  # package check----------------------------------------------------------------
  if (!requireNamespace("survminer", quietly = TRUE)) {
    stop("Package survminer required. Install with install.packages('survminer')")
  }

  # main function ------------------------------------------------------------
  result <- tryCatch({
    survminer::surv_cutpoint(
      data = data,
      time = time,
      event = status,
      variables = pred
    )
  }, error = function(e) {
    stop("Cutpoint calculation failed: ", conditionMessage(e))
  })

  # result ----------------------------------------------------------------
  if (is.null(result$cutpoint)) {
    warning("No optimal cutoff found. Check input data distribution.")
    return(NA_real_)
  }

  result$cutpoint[[1]]$cutpoint
}

