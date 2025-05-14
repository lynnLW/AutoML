#' Select core prognostic features from the candidate genes
#'
#' Identify Key Genes Through Feature Frequency Analysis and Visualization
#'
#' @param selected.feature Feature selection results (output from feature selection function)
#' @param featurelist Alternative feature list (can be used instead of selected.feature)
#' @param outdir Output directory path (default: current directory)
#' @param sets Methods to display (default: all available methods)
#' @param top Number of top frequent features to display (default: 20)
#' @param top_select Logical indicating whether to select features by:
#'        - TRUE: top N most frequent features
#'        - FALSE: features meeting nmethod threshold
#' @param nmethod min frequency selected by 11 ML algorithms (Default 5)
#' @param width plot width (cm)
#' @param height plot height (cm)
#'
#' @return The final selected genes and saving plots
#' @export
#'
top_feature_select <- function(
    selected.feature,
    featurelist = NULL,
    sets = NULL,
    top = 20,
    top_select = FALSE,
    nmethod = 5,
    width = 7.5,
    height = 9,
    outdir = '1.feature_select/') {

  # parameters ----------------------------------------------------------------
  if (is.null(selected.feature) && is.null(featurelist)) {
    stop("Must input selected.feature or featurelist")
  }

  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
    message("Creating the output directory: ", outdir)
  }

  # parameters --------------------------------------------------------------
  sets <- if (is.null(sets)) {
    if (!is.null(selected.feature)) unique(selected.feature$method) else names(featurelist)
  } else {
    intersect(sets, if (!is.null(selected.feature)) unique(selected.feature$method) else names(featurelist))
  }

  # preprocessing--------------------------------------------------------------
  core_feature_list <- if (!is.null(selected.feature)) {
    lapply(sets, function(m) {
      selected.feature$selected.fea[selected.feature$method == m]
    }) %>% stats::setNames(sets)
  } else {
    featurelist[sets]
  }

  # Frequency ------------------------------------------------------------
  feature_freq <- if (!is.null(selected.feature)) {
    dplyr::count(selected.feature, .data$selected.fea, name = "Frequence") %>%
      dplyr::arrange(dplyr::desc(.data$Frequence)) %>%
      dplyr::mutate(selected.fea = gsub("\\.", "-", .data$selected.fea))
  } else {
    purrr::map_df(featurelist[sets], ~ tibble::tibble(gene = .x), .id = "method") %>%
      dplyr::count(.data$gene, name = "Frequence") %>%
      dplyr::arrange(dplyr::desc(.data$Frequence)) %>%
      dplyr::mutate(gene = gsub("\\.", "-", .data$gene))
  }

  # saving results ----------------------------------------------------------------
  utils::write.csv(feature_freq, file.path(outdir, "feature_frequency.csv"), row.names = FALSE)

  # view plot --------------------------------------------------------------
  # Helper 1: UpSet Plot --------------------------------------------------------
  plot_upset <- function(feature_list, outdir) {

    grDevices::jpeg(
      filename = file.path(outdir, "feature_upset_plot.jpg"),
      res = 600,
      width = 12,
      height = 12,
      units = "cm"
    )
    on.exit(grDevices::dev.off())

    print(
      UpSetR::upset(
        UpSetR::fromList(feature_list),
        sets = names(feature_list),
        order.by = "freq",
        mb.ratio = c(0.6, 0.4),
        mainbar.y.label = "Intersection Size",
        sets.x.label = "Set Size",
        text.scale = 0.8,
        sets.bar.color = "#E41A1C",
        main.bar.color = "#377EB8"
      )
    )
  }

  # Helper 2: Frequency Plot ----------------------------------------------------
  plot_frequency <- function(freq_data, top, outdir, width, height) {
    top <- min(top, nrow(freq_data))
    plot_data <- utils::head(freq_data, top)

    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = stats::reorder(.data[[names(plot_data)[1]]], .data$Frequence),
        y = .data$Frequence
      )
    ) +
      ggplot2::geom_segment(
        ggplot2::aes(xend = .data[[names(plot_data)[1]]], yend = 0),
        color = "#377EB8"
      ) +
      ggplot2::geom_point(
        ggplot2::aes(size = .data$Frequence),
        color = "#4DAF4A",
        alpha = 0.8
      ) +
      ggplot2::labs(
        title = paste("Top", top, "Frequent Features"),
        x = "",
        y = "Frequency"
      ) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.3),
        axis.text.x = ggplot2::element_text(size = 8, colour = "black"),
        axis.text.y = ggplot2::element_text(size = 8, colour = "black"),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 10),
        legend.title = ggplot2::element_text(size = 8, colour = "black"),
        panel.background = ggplot2::element_rect(fill = "white")
      ) +
      ggplot2::coord_flip()

    print(p)

    # saving plots
    ggplot2::ggsave(
      filename = file.path(outdir, "top_feature_selection.jpg"),
      plot = p,
      width = width,
      height = height,
      dpi = 600,
      units = "cm"
    )
  }

  ## UpSet Plot ---------------------------
  plot_upset(core_feature_list, outdir = outdir)

  ## Frequency Lollipop Plot --------------
  plot_frequency(feature_freq, top = top, outdir = outdir, width = width, height = height)

  # feature selection ----------------------------------------------------------------
  final.feature <- if (top_select) {
    utils::head(feature_freq, top)[[1]]
  } else {
    feature_freq[feature_freq$Frequence >= nmethod, ][[1]]
  }

  save(final.feature, file = file.path(outdir, "final_selected_feature.Rdata"))
  message("The final selected genes: ", length(final.feature))

  return(final.feature)
}
