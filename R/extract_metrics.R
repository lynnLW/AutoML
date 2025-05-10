#' Calculate Average Performance Metrics for Cross-Validated Models
#'
#' Computes mean performance metrics across training, validation and test sets
#' from multiple cross-validated models
#'
#' @param model_list List of models containing train/valid/test performance metrics
#' @param metric_names Vector of metric names to extract (default includes C-index,
#'        IBS and time-dependent AUCs)
#' @return List containing three data frames:
#'         - train: Training set metrics
#'         - valid: Validation set metrics
#'         - test: Test set metrics
#' @export
#' @examples
#' \dontrun{
#' cv_results <- extract_metrics(model_list = cv_models)
#' }
extract_metrics <- function(model_list, metric_names = NULL) {
  # Set default metrics if not specified
  if (is.null(metric_names)) {
    metric_names <- c("cindex", "bs", paste0("km_auc_", c(1, 2, 3, 5, 7, 10)))
  }

  # Enhanced safe extraction with type checking
  safe_extract <- function(model, metric) {
    tryCatch({
      val <- model[[metric]]
      if (is.null(val)) return(NA_real_)
      if (!is.numeric(val)) return(NA_real_)
      val
    }, error = function(e) NA_real_)
  }

  # Optimized metric extraction using vectorization
  extract_metrics_df <- function(data_list) {
    metric_matrix <- vapply(
      data_list,
      function(x) {
        vapply(metric_names, function(m) safe_extract(x, m), numeric(1))
      },
      numeric(length(metric_names))
    )

    stats::setNames(
      as.data.frame(t(metric_matrix)),
      metric_names
    )
  }

  # Extract metrics for each dataset type
  list(
    train = extract_metrics_df(lapply(model_list, `[[`, "train")),
    valid = extract_metrics_df(lapply(model_list, `[[`, "valid")),
    test = extract_metrics_df(lapply(model_list, `[[`, "test"))
  )
}
