#' Generate Random Forest/Survival Forest Feature Importance Plot and Results
#'
#' This function generates a feature importance plot and CSV file for random forest
#' or survival forest models. It visualizes variable importance rankings and saves
#' both the plot and data.
#'
#' @param fit Model fit object containing importance scores (requires `fit$importance`).
#' @param prefix Prefix string for output filenames.
#' @param outdir Output directory path for saving results.
#'
#' @return A data frame containing feature importance scores sorted in descending order.
#'
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_bar geom_text theme_bw theme element_text
#' @importFrom ggplot2 labs ggsave
#' @importFrom utils write.csv
#' @importFrom stats reorder
#'
generat_rfrsf_importance <- function(fit, prefix, outdir) {
  # Extract importance scores
  importance <- as.data.frame(fit$importance)
  names(importance)[1] <- "Importance"
  importance$gene <- row.names(importance)

  # Order by importance
  importance <- importance[order(importance$Importance), ]

  # Convert to factor for plotting
  importance$gene <- factor(importance$gene, levels = importance$gene)

  # Create importance plot
  p <- ggplot2::ggplot(importance,
                       ggplot2::aes(x = stats::reorder(.data$gene, .data$Importance),
                                    y = .data$Importance)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.3f", .data$Importance)),
      position = ggplot2::position_dodge(width = 0.9),
      hjust = -0.1,
      size = 3,
      color = "black"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 10, color = "black"),
      axis.text.x = ggplot2::element_text(size = 10, color = "black"),
      axis.title.x = ggplot2::element_text(size = 12, color = "black"),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = "",
      x = "Importance",
      y = ""
    )

  # Print and save plot
  print(p)
  ggplot2::ggsave(
    file.path(outdir, paste0(prefix, "_importance_plot.jpg")),
    plot = p,
    width = 5,
    height = 5,
    dpi = 600
  )

  # Save sorted results
  importance <- importance[order(importance$Importance, decreasing = TRUE), ]
  utils::write.csv(
    importance,
    file = file.path(outdir, paste0(prefix, "_importance.csv")),
    row.names = FALSE
  )

  return(importance)
}
