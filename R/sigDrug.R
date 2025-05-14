#' Perform Differential Expression Analysis for Drug Targets
#'
#' This function calculates differential expression between groups for a list of genes,
#' computes p-values using Wilcoxon test, and saves results.
#'
#' @param merge_expr A data frame containing expression data with a 'group' column for sample grouping.
#' @param genelist A character vector of gene names to analyze.
#' @param outdir Output directory path for saving results.
#'
#' @return A data frame with mean expression values, p-values, fold change, and filtered results (p < 0.05).
#' @export
#'
sigDrug <- function(merge_expr, genelist, outdir) {
  # Check input validity
  if (!"group" %in% colnames(merge_expr)) {
    stop("Input data must contain a 'group' column")
  }

  # Create output directory if missing
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
    message("Created output directory: ", outdir)
  }

  ### Diff expression calculation
  df_mean <- merge_expr %>%
    dplyr::group_by(.data$group) %>%
    dplyr::summarise(dplyr::across(
      tidyselect::all_of(genelist),
      \(x) mean(x, na.rm = TRUE)
    )) %>%
      t() %>%
      as.data.frame()

    # Format mean matrix
    colnames(df_mean) <- df_mean[1, ]
    df_mean <- df_mean[-1, ]
    df_mean <- stats::na.omit(df_mean)

    ### Function to calculate p-values
    cal_pvalue <- function(input, gene) {
      stats::pairwise.wilcox.test(
        input[, gene],
        input[, "group"],
        p.adjust.method = "bonf"
      )$p.value
    }

    ### Calculate p-values for all genes
    pvalue <- vector("numeric", length = nrow(df_mean))
    for (i in seq_len(nrow(df_mean))) {
      gene <- rownames(df_mean)[i]
      pvalue[i] <- cal_pvalue(merge_expr, gene)
    }

    ### Add results to data frame
    df_mean$pvalue <- pvalue
    df_mean[, 1] <- as.numeric(df_mean[, 1])
    df_mean[, 2] <- as.numeric(df_mean[, 2])
    df_mean$fc <- df_mean[, 1] / df_mean[, 2]
    df_mean <- df_mean[df_mean$pvalue < 0.05, ]

    # Save results
    utils::write.csv(
      df_mean,
      file = file.path(outdir, "1.diff_drug_0.05.csv"),
      row.names = TRUE
    )

    return(df_mean)
}
