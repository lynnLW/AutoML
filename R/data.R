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

#' Results Example
#'
#' Precomputed  results demonstrating the expected output format
#' of the analysis functions. This dataset serves as a reference for users to
#' understand and validate their own analysis results.
#'
#' @format A list containing analysis results:
#' \describe{
#'   \item{feature_select}{a data frame containing results from feature selection}
#'   \item{candidate_genes}{a vector of gene list for modeling}
#'   \item{10_5_model_list}{a list containing results from ML.survival.model}
#'   \item{test_index}{List of each model containing performace metrics in external validation cohorts}
#' }
#' @source Generated from example results using data in demo list
"demo_result"
