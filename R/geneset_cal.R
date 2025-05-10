#' Optimized Gene Set Enrichment Analysis Function
#'
#' Performs gene set enrichment analysis with enhanced functionality
#'
#' @param expr_matrix Expression matrix with genes as rows and samples as columns
#' @param category Gene set collection (MSigDB categories or built-in collections)
#' @param subcategory MSigDB subcategory (optional)
#' @param geneset Input gene sets, supports following formats:
#'   - Character vector (e.g. c("Gene1", "Gene2")) - will be converted to list named "custom_set"
#'   - Named list (e.g. list(pathway1 = c("Gene1", "Gene2")))
#' @param species Species, defaults to "Homo sapiens"
#' @param method GSVA method, defaults to "ssgsea"
#' @param prefix Output filename prefix
#' @param output_dir Output directory, defaults to "./GSVA_Results"
#' @return GSVA enrichment score matrix
#' @export
#'
#' @examples
#' \dontrun{
#' # Using named list
#' custom_sets <- list(DNA_repair = c("BRCA1", "BRCA2", "TP53"),
#'                    Cell_cycle = c("CDK1", "CDK2", "CCND1"))
#' gsva_results <- geneset_cal(expr_matrix = expr_data,
#'                               geneset = custom_sets)
#' }

geneset_cal <- function(expr_matrix,
                        category = NULL,
                        subcategory = NULL,
                        geneset = NULL,
                        species = "Homo sapiens",
                        method = "ssgsea",
                        prefix = NULL,
                        output_dir = "./GSVA_Results") {

  # parameter ----------------------------------------------------------------
  suppressPackageStartupMessages({
    require(GSVA, quietly = TRUE)
    require(msigdbr, quietly = TRUE)
    require(dplyr, quietly = TRUE)
  })

  # expr matrix
  if (!is.matrix(expr_matrix) && !is.data.frame(expr_matrix)) {
    stop("expr_matrix must be matrix or dataframe")
  }

  # the output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # gene set or category
  if (!is.null(category) && !is.null(geneset)) {
    warning("When setting category and geneset,using category in priority")
    geneset <- NULL
  }

  # main function----------------------------------------------------------------

  # getting gene set
  get_genesets <- function(){
    # loading internal gene sets
    if(F){
      internal_genesets <- get(load("data/internal_genesets.rdata"))
      msigdb_categories <- c("C3","C2","C8","C6","C7","C4","C5","H","C1")
    }


    # Msigdb
    if (!is.null(category)) {
      if (category %in% msigdb_categories) {

        gs_data <- msigdbr(species = species,
                           category = category,
                           subcategory = subcategory) %>%
          group_by(gs_name) %>%
          distinct(gene_symbol, .keep_all = TRUE) %>%
          ungroup()

        gs_list <- split(gs_data$gene_symbol, gs_data$gs_name)
        prefix <- if (is.null(prefix)) {
          paste(c(category, subcategory), collapse = "_")
        } else {
          prefix
        }
        return(list(gs = gs_list, prefix = prefix))
      }

      # internal gene sets
      if (category %in% names(internal_genesets)) {
        return(list(gs = internal_genesets[[category]],
                    prefix = ifnull(prefix, category)))
      }

      stop("unknown category:", category)
    }

    # undermined gene sets
    if (!is.null(geneset)) {
      gs=validate_custom_geneset(geneset)
      return(list(gs = gs,
                  prefix = ifnull(prefix, "custom_geneset")))
    }

    stop("Must input category or geneset")
  }

  # undermined gene sets
  validate_custom_geneset <- function(gs) {

    if (!is.list(gs)) {

      set_name <- if (!is.null(names(gs))) names(gs)[1] else "custom_set"
      gs <- list(set_name = gs)
      names(gs) <- set_name
    }


    if (!all(sapply(gs, is.character))) {
      invalid_types <- unique(sapply(gs, class))
      stop("Gene set must be character:",
           paste(invalid_types, collapse = ", "))
    }

    # null gene sets
    if (any(sapply(gs, length) == 0)) {
      empty_sets <- names(gs)[sapply(gs, length) == 0]
      stop("Null gene sets", paste(empty_sets, collapse = ", "))
    }

    # Auto name the input gene set
    if (is.null(names(gs))) {
      names(gs) <- paste0("GeneSet_", seq_along(gs))
      warning("Undefined gene set:",
              paste(names(gs), collapse = ", "))
    }

    return(gs)
  }

  ### Other function --------------------------------------------------------------
  ifnull <- function(x, y){if(is.null(x)) y else x}

  format_expr<-function(expr){
    ##
    names(expr)[1]<-"symbol"
    expr<-expr[which(expr$symbol!="NA"),]
    # Calculates row means for all duplicate entries
    index=order(rowMeans(expr[,!colnames(expr) %in% c("symbol")]),decreasing = T)
    expr_ordered=expr[index,]
    # Identifies the row with the highest mean expression, retains only this first occurring, highest-expressed instance
    keep=!duplicated(expr_ordered$symbol)
    # final exper
    expr_max=expr_ordered[keep,]
    row.names(expr_max)<-expr_max$symbol
    expr_max<-expr_max[,!colnames(expr_max) %in% c("symbol")]
    ##
    return(expr_max)
  }

  # main function --------------------------------------------------------------
  tryCatch({
    # getting gene sets
    gs_info <- get_genesets()
    gs_list <- gs_info$gs
    prefix <- gs_info$prefix

    # GSVA
    gsva_scores <- gsva(
      expr = as.matrix(expr_matrix),
      gset.idx.list = gs_list,
      method = method,
      kcdf = "Gaussian",
      parallel.sz = 4,
      min.sz = 1,
      verbose = T
    )

    # MPI transform
    if(category=="MPI"){
      mapping<-internal_genesets[['MPI_mapping']]
      intersect_genes<-intersect(row.names(gsva_scores),unique(mapping$aids))
      mapping<-mapping[intersect_genes,]
      gsva_scores<-gsva_scores[intersect_genes,]
      gsva_scores<-data.frame(symbol=mapping$a,gsva_scores)
      gsva_scores<-format_expr(gsva_scores)
    }

    # result
    result_df <- as.data.frame(t(gsva_scores))
    output_file <- file.path(output_dir,
                             sprintf("GSVA_%s_%s.csv", prefix, method))

    # saving
    write.csv(result_df, file = output_file, row.names = TRUE)
    message("Finshed! Results saving in:", output_file)

    return(result_df)

  }, error = function(e) {
    stop("Failed:", e$message)
  })
}
