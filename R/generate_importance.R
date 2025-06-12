#' Generate Random Forest/Survival Forest Feature Importance Plot
#'
#' This function generates a horizontal feature importance plot and CSV file for
#' random forest or survival forest models. It visualizes variable importance rankings
#' and saves both the plot and data.
#'
#' @param fit Model fit object containing importance scores (requires `fit$importance`).
#' @param prefix Prefix string for output filenames.
#' @param outdir Output directory path for saving results.
#' @param top_n Number of top features to display (NULL for all).
#' @param bar_color Color for importance bars (default: "steelblue").
#' @param text_size Size of importance value labels (default: 3).
#' @param width Plot width in inches (default: 8).
#' @param height Plot height in inches (default: 6).
#'
#' @return A data frame containing feature importance scores sorted in descending order.
#'
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_bar geom_text theme_bw theme element_text
#' @importFrom ggplot2 labs ggsave coord_flip scale_y_continuous
#' @importFrom utils write.csv
#' @importFrom stats reorder
generat_rf_importance <- function(fit,
                                     prefix,
                                     outdir,
                                     top_n = NULL,
                                     bar_color = "steelblue",
                                     text_size = 3,
                                     width = 8,
                                     height = 6) {

  # Extract importance scores
  importance <- as.data.frame(fit$importance)
  names(importance)[1] <- "Importance"
  importance$gene <- row.names(importance)

  # Order by importance (descending)
  importance <- importance[order(importance$Importance, decreasing = TRUE), ]

  # Subset top N features if specified
  if (!is.null(top_n)) {
    importance <- head(importance, top_n)
  }

  # Convert to factor for plotting (ordered by importance)
  importance$gene <- factor(importance$gene, levels = rev(importance$gene))

  # Create horizontal importance plot
  p <- ggplot2::ggplot(
    importance,
    ggplot2::aes(
      x = stats::reorder(.data$gene, .data$Importance),
      y = .data$Importance
    )
  ) +
    ggplot2::geom_bar(
      stat = "identity",
      fill = bar_color,
      width = 0.7
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = sprintf("%.3f", .data$Importance)
      ),
      hjust = -0.1,
      size = text_size,
      color = "black"
    ) +
    ggplot2::coord_flip() +  # Flip coordinates for horizontal bars
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1))) + # Add space for labels
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 10, color = "black"),
      axis.text.x = ggplot2::element_text(size = 10, color = "black"),
      axis.title.y = ggplot2::element_blank(),  # Remove y-axis title
      axis.title.x = ggplot2::element_text(size = 12, color = "black"),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"),
      panel.grid.major.y = ggplot2::element_blank(),  # Remove horizontal grid lines
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(color = "grey90"),  # Light vertical grid
      panel.grid.minor.x = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      x = "Importance Score",
      y = "",
      title = "Feature Importance Ranking"
    )

  # Print and save plot
  print(p)
  ggpubr::ggexport(
    file.path(outdir, paste0(prefix, "_importance_plot.jpg")),
    plot = p,
    width = width*250,
    height = height*250,
    res=600
  )

  # Save sorted results
  utils::write.csv(
    importance,
    file = file.path(outdir, paste0(prefix, "_importance.csv")),
    row.names = FALSE
  )

  return(importance)
}
