#' Calculate and visualize differences between two groups
#'
#' @param expr_matrix Expression matrix (features in columns, samples in rows)
#' @param grouping Factor vector with group assignments
#' @param top_n Number of top differentially expressed features to show (default: 10)
#' @param width Plot width in cm (default: 24)
#' @param height Plot height in cm (default: 16)
#' @param outdir Output directory path (default: "./")
#' @param p_adjust_method Method for p-value adjustment (default: "none")
#' @param plot_type Type of plot ("boxplot" or "violin", default: "boxplot")
#' @return Invisibly returns a data frame with statistical results
#' @importFrom stats wilcox.test p.adjust
#' @importFrom dplyr bind_cols select all_of %>%
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_violin scale_fill_manual
#' @importFrom ggplot2 facet_wrap theme_bw labs theme element_text element_blank
#' @importFrom ggplot2 ggsave after_stat
#' @importFrom ggpubr stat_compare_means
#' @importFrom ggsci scale_fill_d3
#' @importFrom utils head
#' @export
cal_diff <- function(expr_matrix,
                     grouping,
                     top_n = 10,
                     width = 24,
                     height = 16,
                     outdir = "./",
                     p_adjust_method = "none",
                     plot_type = "boxplot") {
  # output directory
  if (!dir.exists(outdir)){
    dir.create(outdir,recursive = T)
  }

  if (!is.factor(grouping)) {
    grouping <- as.factor(grouping)
  }

  if (nlevels(grouping) != 2) {
    stop("Grouping factor must have exactly 2 levels")
  }

  # Calculate Wilcoxon p-values for all features
  p_values <- apply(expr_matrix, 2, function(x) {
    tryCatch(
      wilcox.test(x ~ grouping)$p.value,
      error = function(e) NA_real_
    )
  })

  # Adjust p-values if requested
  if (p_adjust_method != "none") {
    p_values <- p.adjust(p_values,method = p_adjust_method)
  }

  # Get top N features with smallest p-values
  top_features <- names(sort(p_values, na.last = TRUE))[1:min(top_n, length(p_values))]
  top_features <- top_features[!is.na(top_features)]

  if (length(top_features) == 0) {
    warning("No significant features found")
    return(invisible(NULL))
  }

  # Prepare data for plotting
  plot_data <- expr_matrix %>%
    as.data.frame() %>%
    dplyr::select(all_of(top_features)) %>%
    dplyr::bind_cols(group = grouping) %>%
    tidyr::pivot_longer(
      cols = -.data$group,
      names_to = "Feature",
      values_to = "Expression"
    )

  # Create plot base
  p <-ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[["group"]],
                                              y = .data[["Expression"]],
                                              fill = .data[["group"]]))

  # Add the appropriate geometry based on plot_type
  if (plot_type == "violin") {
    p <- p +
      ggplot2::geom_violin(trim = FALSE, alpha =1) +
      ggplot2::geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA)
  } else {
    p <- p +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha =1)
  }

  p<-p+
    ggsci::scale_fill_d3() +
    ggplot2::facet_wrap(~Feature, scales = "free_y", nrow = 2) +
    ggpubr::stat_compare_means(
      ggplot2::aes(label = ggplot2::after_stat(.data[['p.signif']])),
      method = "wilcox.test",
      label.x = 1.5,
      size = 3
    ) +
    ggplot2::labs(
      y = "Values",
      title = paste("Top", length(top_features), "Differentially Expressed Features")
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = element_text(colour = "black", size = 10),
      axis.title.x = element_blank(),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(colour = "black", size = 10),
      title = element_text(colour = "black", size = 10)
    )
  print(p)
  file_path <- file.path(outdir, paste0("top_diff_combined.jpg"))
  # Save plot
  ggpubr::ggexport(
    filename = file_path,
    width = width*300,
    height = height*300,
    device = "jpg",
    res=600
  )
}

