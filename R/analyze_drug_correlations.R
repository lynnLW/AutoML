#' Drug-Risk Score Correlation Analysis
#'
#' Identifies drugs significantly correlated with risk scores and generates diagnostic plots.
#'
#' @param data A dataframe containing drug sensitivity data with risk scores
#' @param risk_score_col Name of column containing risk scores (default: "risk_score")
#' @param cor_method Correlation method ("pearson" or "spearman", default: "pearson")
#' @param r_threshold Absolute correlation cutoff (default: 0.5)
#' @param p_threshold Significance threshold (default: 0.05)
#' @param top_n Number of top hits to plot (default: 10)
#' @param point_alpha Point transparency (default: 0.6)
#' @param point_color Point color (default: "#1E88E5")
#' @param line_color Trendline color (default: "#D81B60")
#' @param ncol Number of columns in combined plot (default: 2)
#' @param outdir the output directory
#'
#' @return A list containing:
#' \itemize{
#'   \item{cor_results: Full correlation results (tibble)}
#'   \item{top_associations: Filtered significant results (tibble)}
#'   \item{plots: List of ggplot objects for top associations}
#' }
#'
#' @importFrom dplyr filter arrange mutate %>% bind_rows
#' @import ggplot2
#' @importFrom utils head
#' @export
analyze_drug_correlations <- function(data,
                                      risk_score_col = "risk_score",
                                      cor_method = "pearson",
                                      r_threshold = 0.35,
                                      p_threshold = 0.05,
                                      top_n = 10,
                                      point_alpha = 0.6,
                                      point_color = "#1E88E5",
                                      line_color = "#D81B60",
                                      ncol=2,
                                      outdir="./") {

  # Input validation
  if (!risk_score_col %in% colnames(data)) {
    stop("Risk score column '", risk_score_col, "' not found in input data")
  }

  if (!cor_method %in% c("pearson", "spearman")) {
    stop("cor_method must be either 'pearson' or 'spearman'")
  }

  # output directory
  if (!dir.exists(outdir)){
    dir.create(outdir,recursive = T)
  }

  drug_cols <- setdiff(colnames(data), risk_score_col)

  # Compute correlations
  cor_results <- lapply(drug_cols, function(drug) {
    ct <- suppressWarnings(
      stats::cor.test(
        x = data[[drug]],
        y = data[[risk_score_col]],
        method = cor_method,
        exact = FALSE
      )
    )

    tibble::tibble(
      drug = drug,
      r = ct$estimate,
      p.value = ct$p.value,
      n = sum(stats::complete.cases(data[, c(drug, risk_score_col)])))
  }) %>% dplyr::bind_rows()

  utils::write.csv(cor_results,file=paste0(outdir,"/cor_results.csv"))

  # Apply significance filters and FDR correction
  sig_results <- cor_results %>%
    dplyr::filter(abs(.data$r) > r_threshold &
                      .data$p.value < p_threshold) %>%
      dplyr::arrange(.data$p.value) %>%
      dplyr::mutate(p.adj = stats::p.adjust(.data$p.value, method = "BH"))

  # Initialize empty result list
  result <- list(
    cor_results = cor_results,
    top_associations = sig_results,
    plots = list()
  )

  if (nrow(sig_results)>0){
    # Determine actual number of plots to show
    actual_n <- min(top_n, nrow(sig_results))
    top_drugs <- head(sig_results, top_n)
    plot_list <- list()

    for (i in seq_len(nrow(top_drugs))) {
      drug <- top_drugs$drug[i]
      r_val <- round(top_drugs$r[i], 3)
      p_val <- format.pval(top_drugs$p.value[i], digits = 3, eps = 0.001)

      # Create annotation text
      annot_text <- paste0("r = ", r_val, "\np = ", p_val)

      plot_list[[drug]] <- ggplot2::ggplot(
        data = data,
        mapping = aes(
          x = .data[[drug]],
          y = .data[[risk_score_col]]
        )
      ) +
        ggplot2::geom_point(
          alpha = point_alpha,
          color = point_color,
          na.rm = TRUE
        ) +
        ggplot2::geom_smooth(
          method = "lm",
          color = line_color,
          se = FALSE,
          na.rm = TRUE
        ) +
        # Add correlation annotation
        ggplot2::annotate(
          "text",
          x = Inf, y = Inf,
          label = annot_text,
          hjust = 1.1, vjust = 1.1,
          size = 3,
          color = "black"
        ) +
        ggplot2::labs(
          title = "",
          x = paste(drug),
          y = "Risk Score"
        ) +
        ggplot2::theme_bw()+
        ggplot2::theme(
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.title = element_text(colour = "black",size=8),
          axis.text = element_text(colour = "black",size=6),
          title = element_text(colour = "black",size=8)
        )
    }

    result$plots <- plot_list

    # Create combined plot with dynamic height
    n_plots <- length(plot_list)
    nrow <- ceiling(n_plots / ncol)

    # Calculate dynamic height (base height per row is 8cm)
    plot_height <- max(8 * nrow, 12)  # Minimum height of 12cm

    # Create combined plot
    combined_plot <- patchwork::wrap_plots(plot_list, ncol = ncol) +
      patchwork::plot_annotation(
        title = "Top Drug-Risk Score Correlations",
        subtitle = paste("Showing", n_plots, "significant associations"),
        theme = theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 10)
        )
      )

    print(combined_plot)

    # Define output path
    file_path <- file.path(outdir, paste0("top_cor_combined.jpg"))

    # Save plot with dynamic dimensions
    ggplot2::ggsave(
      filename = file_path,
      plot = combined_plot,
      width = 12 * min(ncol, 2),  # Adjust width based on number of columns
      height = plot_height,
      device = "jpg",
      units = "cm",
      dpi = 600
    )
  }

}

