#' Gene Set Collection
#'
#' Predefined gene sets for pathway enrichment analysis.
#'
#' @format A list where each element is a character vector representing a gene set:
#' \describe{
#'   \item{geneset_name}{Name of the gene set (character)}
#'   \item{genes}{List of gene symbols (character vector)}
#' }
#' @source Curated from MSigDB (https://www.gsea-msigdb.org).
"internal_genesets"

#' Training Dataset for Demonstration
#'
#' A sample dataset containing training data and gene lists for testing pathway
#' analysis functions. This dataset is provided as an example for users to
#' familiarize themselves with the package functionality.
#'
#' @format A list containing two components:
#' \describe{
#'   \item{training_data}{A data frame with example gene expression values}
#'   \item{gene_list}{Character vector of gene symbols for testing purposes}
#'   \item{list_train_vali_Data}{A list with external validata data frame}
#'   \item{list_train_vali_meta}{A list with external validata meta information}
#' }
#' @source Derived from GEO \url{https://www.ncbi.nlm.nih.gov/geo/}
"demo_list"
