# R/internal_data.R

#' @title Load Internal Gene Sets
#' @description Makes internal_genesets available within package environment
#' @keywords internal
.load_internal_genesets <- function() {
  utils::data("internal_genesets", envir = environment())
  get("internal_genesets", envir = environment())
}