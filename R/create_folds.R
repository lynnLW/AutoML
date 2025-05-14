#' Create Repeated Stratified Folds for Cross-Validation
#'
#' Generates repeated stratified cross-validation folds while preserving class distribution.
#'
#' @param data A data frame containing the dataset to be split.
#' @param fold Integer, number of folds (k) to create.
#' @param nrepeats Integer, number of times to repeat the k-fold splitting.
#' @param strata Character, name of the column in `data` used for stratification.
#' @param seed Integer, random seed for reproducibility.
#'
#' @return A list of length `nrepeats`, where each element is a list of fold indices.
#' @export
create_folds <- function(data, fold, nrepeats, strata, seed) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.")
  }
  if (!strata %in% colnames(data)) {
    stop("`strata` column not found in `data`.")
  }
  if (fold < 2) {
    stop("`fold` must be >= 2.")
  }
  if (nrepeats < 1) {
    stop("`nrepeats` must be >= 1.")
  }


  y <- data[[strata]]

  if (length(unique(y)) < 2) {
    stop("Stratification requires at least 2 classes in `strata`.")
  }

  set.seed(seed)
  folds_list <- lapply(
    1:nrepeats,
    function(x) caret::createFolds(as.factor(y), k = fold)
  )

  return(folds_list)
}

