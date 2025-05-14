#' Display Progress Updates During Iterative Operations
#'
#' Prints progress percentage to the console at specified intervals.
#' Useful for long-running loops to provide user feedback.
#'
#' @param index Integer, current iteration index (1-based).
#' @param totalN Integer, total number of iterations.
#' @param breakN Integer, number of progress updates to display (default: 20).
#'
#' @return Invisibly returns `NULL`. Prints progress to console as a side effect.
#'
#' @export
display_progress <- function(index, totalN, breakN = 20) {
  # Input validation
  if (!is.numeric(index) || length(index) != 1 || index < 1) {
    stop("`index` must be a positive integer.")
  }
  if (!is.numeric(totalN) || length(totalN) != 1 || totalN < 1) {
    stop("`totalN` must be a positive integer.")
  }
  if (!is.numeric(breakN) || length(breakN) != 1 || breakN < 1) {
    stop("`breakN` must be a positive integer.")
  }
  if (index > totalN) {
    stop("`index` cannot exceed `totalN`.")
  }

  # Calculate and print progress
  if (index %% ceiling(totalN / breakN) == 0 || index == totalN) {
    percentage <- round(index * 100 / totalN)
    cat(sprintf("%d%% ", percentage))
    utils::flush.console()  # Explicit package reference
  }

  invisible(NULL)
}
