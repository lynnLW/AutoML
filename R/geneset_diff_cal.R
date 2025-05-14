#' Train Survival Models Using Multiple Machine Learning Algorithms
#'
#' Performs differential analysis on gene sets and trains survival models using
#' various machine learning approaches. Supports both limma and basic statistical methods.
#'
#' @param gs.exp Gene expression matrix (genes as rows, samples as columns)
#' @param grouping Dataframe with sample grouping information (must contain 'accession' and 'group' columns)
#' @param prefix Prefix string for output filenames
#' @param geneset_name Name of the gene set being analyzed
#' @param method Analysis method: "limma" for linear models or alternative statistical method (default: "limma")
#' @param normalized Logical indicating whether to normalize expression data (default: TRUE)
#' @param outdir Output directory path (default: "./")
#' @importFrom dplyr %>% group_by summarise across
#'
#' @return Invisibly returns a list containing:
#' \itemize{
#'   \item For limma method: List of differential expression results for each contrast
#'   \item For other methods: Dataframe with mean expression and p-values
#' }
#'
#' @export
#'
geneset_diff_cal <- function(gs.exp,
                             grouping,
                             prefix,
                             geneset_name,
                             method=NULL,
                             normalized=T,
                             outdir="./"){

  ## Set default method parameter
  if (is.null(method)){
    method="limma"
  } else {
    method=method
  }

  ## Create output directory if it doesn't exist
  if(!dir.exists(outdir)){
    dir.create(outdir)
  }

  ## Format input data
  names(grouping)[1] <- "accession"
  group_info <- data.frame(row.names=grouping$accession, group=grouping$group)
  gs.exp <- gs.exp[row.names(group_info), , drop=FALSE]

  # Standardize dataframe column names to comply with R naming conventions
  clean_column_names <- function(df, replace_char = "_") {

    # Get original column names
    original_names <- colnames(df)

    # Step 1: Replace all illegal characters
    clean_names <- gsub(
      pattern = "[^a-zA-Z0-9]",  # Match non-alphanumeric characters
      replacement = replace_char,
      x = original_names
    )

    # Step 2: Handle consecutive replacement characters
    clean_names <- gsub(
      pattern = paste0(replace_char, "+"),  # Match consecutive replacement chars
      replacement = replace_char,
      x = clean_names
    )

    # Step 3: Remove leading/trailing replacement characters
    clean_names <- gsub(
      pattern = paste0("^", replace_char, "|", replace_char, "$"),
      replacement = "",
      x = clean_names
    )

    # Step 4: Handle names starting with numbers
    clean_names <- ifelse(
      grepl("^[0-9]", clean_names),
      paste0("X_", clean_names),  # Add prefix
      clean_names
    )

    # Step 5: Ensure uniqueness
    colnames(df) <- make.unique(clean_names)

    # Return modified dataframe
    return(df)
  }

  # Execute column name cleaning
  cleaned_df <- clean_column_names(gs.exp)

  ## Perform differential analysis between groups
  if (method=="limma"){
    groups <- group_info$group
    exp2 <- as.data.frame(t(cleaned_df))

    if(normalized==T){
      exp2 <- scale(exp2)
    } else {
      exp2 <- exp2
    }

    design <- stats::model.matrix(~0+groups)
    colnames(design) = levels(factor(groups))
    rownames(design) = colnames(exp2)
    head(design)

    # Generate contrasts for all pairwise comparisons
    pairwise_comparisons <- utils::combn(unique(groups), 2, simplify = FALSE)
    for (i in 1:length(pairwise_comparisons)){
      contrasts <- paste0(pairwise_comparisons[[i]][1],"-",pairwise_comparisons[[i]][2])
      compare <- limma::makeContrasts(contrasts=contrasts, levels=design)
      fit <- limma::lmFit(exp2, design)
      fit2 <- limma::contrasts.fit(fit, compare)
      fit3 <- limma::eBayes(fit2)
      Diff <- limma::topTable(fit3, number=Inf)

      ## Save differential expression results
      write.table(Diff, file=paste0(prefix,"_",contrasts,"_diff_DEG.csv"), quote=F, sep=",")
      sig <- Diff[Diff$adj.P.Val<0.05,]
      head(sig)
      write.table(sig, file=paste0(prefix,"_",contrasts,"_sig_DEG.csv"), quote=F, sep=",")
    }

  } else {
    ## Alternative analysis method (not limma)
    genelist = colnames(cleaned_df)
    exp2 <- cbind(cleaned_df, group=group_info$group)

    ## Calculate average expression by group
    mean_values <- exp2 %>%
      group_by(.data$group) %>%
      summarise(across(all_of(genelist), mean))
    mean_values <- mean_values %>% t() %>% as.data.frame()
    names(mean_values) <- mean_values[1,]
    mean_values <- mean_values[-1,]

    ## Generate all possible group combinations
    group_combinations <- utils::combn(unique(exp2$group), 2, simplify = FALSE)

    ## Function to calculate p-values for group comparisons
    calculate_pvalue_for_combination <- function(gene, groups) {
      subset_data <- exp2 %>% filter(.data$group %in% groups)
      result <- ggpubr::compare_means(stats::as.formula(paste(gene, "~ group")),
                              data = subset_data,
                              method = method)
      as.data.frame(result)
    }

    # Calculate p-values for all genes across all group combinations
    pvalues_list <- list()
    pvalues_list <- lapply(genelist, function(gene) {
      pvalue_pairwise_result <- lapply(group_combinations, function(groups) {
        calculate_pvalue_for_combination(gene, groups)
      })
      pvalue_result <- data.frame()
      if(length(pvalue_pairwise_result)>1){
        for(i in 1:length(pvalue_pairwise_result)){
          if (nrow(pvalue_result) > 0) {
            pvalue_result <- cbind(pvalue_result, pvalue_pairwise_result[[i]])
          } else {
            pvalue_result <- pvalue_pairwise_result[[i]]  # Initial assignment
          }
        }
      } else {
        if (nrow(pvalue_result) > 0) {
          pvalue_result <- cbind(pvalue_result, pvalue_pairwise_result[[1]])
        } else {
          pvalue_result <- pvalue_pairwise_result[[1]]  # Initial assignment
        }
      }
      return(pvalue_result)
    })

    # Combine results into dataframe
    pvalues <- do.call(rbind, pvalues_list)
    row.names(pvalues) <- pvalues$.y.

    ## Merge and save results
    static_table <- merge(mean_values, pvalues, by="row.names")
    utils::write.table(static_table,
                file=paste0(outdir,"/",prefix,"_",geneset_name,"_mean_diff.csv"),
                quote=F, sep=",", row.names = F)
  }
}
