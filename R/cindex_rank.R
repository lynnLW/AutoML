#' C-index Distribution Visualization
#'
#' Creates ranked performance plots for survival models
#'
#' @param vali_auc_list List of models containing performance metrics of all cohorts
#' @param index #cindex,bs,km_auc_1, km_auc_3, km_auc_5, km_auc_10
#' @param train # the training cohort default "TCGA"
#' @param plot_type Visualization type ("boxplot", "forestplot", or "barplot")
#' @param outdir Output directory for saving plots
#' @importFrom ggpubr ggexport
#' @importFrom dplyr mutate arrange desc left_join %>%
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradient2
#' @importFrom ggplot2 theme element_blank element_text element_rect labs
#' @importFrom ggplot2 scale_fill_manual coord_flip
#' @importFrom aplot insert_top insert_right
#' @importFrom grDevices dev.off
#' @return ggplot object or list of plots
#' @export
cindex_rank <- function(vali_auc_list,
                        index, #cindex,bs,km_auc_1, km_auc_3, km_auc_5, km_auc_10
                        train="TCGA",
                        plot_type="barplot",#barplot, forestplot, horiz_forestplot
                        outdir="./"
){

  # loading package---------------------------------------------------------
  if (!requireNamespace(c("dplyr","aplot"), quietly = TRUE)) {
    stop("Package dplyr aplot required: install.packages('dplyr','aplot')")
  }

  # the output directory
  if (!dir.exists(outdir)){
      dir.create(outdir,recursive = T)
  }

  dataset_col <- c(RColorBrewer::brewer.pal(9, "Set1"),
      "#3182BDFF", "#E6550DFF", "#31A354FF", "#756BB1FF", "#636363FF", "#6BAED6FF", "#FD8D3CFF", "#74C476FF",
      "#9E9AC8FF", "#969696FF", "#9ECAE1FF", "#FDAE6BFF", "#A1D99BFF", "#BCBDDCFF", "#BDBDBDFF", "#C6DBEFFF",
      "#FDD0A2FF", "#C7E9C0FF", "#DADAEBFF", "#D9D9D9FF") ## default 20 color values


  ##heatmap show all dataset cindex and mean cindex in each model
  ##show cindex in each validated data
  index_cohort_df<-data.frame()
  for (i in 1:length(vali_auc_list)){
        model<-vali_auc_list[[i]]
        model_name<-names(vali_auc_list)[i]
        index_cohort<-data.frame()
        for (j in seq_along(model)){
          cohort<-model[[j]]
          cohort_name<-names(model)[j]
          print(cohort_name)

          # Extract and process index values
          index_df<-lapply(cohort,function(df)df[[index]])
          index_df<- do.call(rbind, index_df)
          index_df<-as.data.frame(index_df)
          names(index_df)<-"index"
          index_df$Rep<-paste0("rep",1:nrow(index_df))
          index_df$Dataset<-cohort_name
          index_cohort<-rbind(index_cohort,index_df)
        }
        index_cohort$Model<-model_name
        index_cohort_df<-rbind(index_cohort_df,index_cohort)
  }
  ##mean index table
  mean_cohort_df <- stats::aggregate(
      x = index_cohort_df$index,
      by = list(index_cohort_df$Dataset,index_cohort_df$Model),
      FUN = mean)
  colnames(mean_cohort_df) <- c("Dataset","Model","index")
  mean_cohort_df$index <- as.numeric(sprintf("%.2f", mean_cohort_df$index))
  ##mean index in all cohorts
  all_mean_index <- stats::aggregate(
      x = index_cohort_df$index,
      by = list(index_cohort_df$Model),
      FUN = mean)
  colnames(all_mean_index) <- c("Model", "mean")
  all_mean_index$mean <- as.numeric(sprintf("%.2f", all_mean_index$mean))
  all_mean_index$Value <- "Mean value in all cohorts"
  ##mean index in test cohorts
  test_cohort_df<-index_cohort_df[index_cohort_df$Dataset != train,]
  test_mean_index <- stats::aggregate(
    x = test_cohort_df$index,
    by = list(test_cohort_df$Model),
    FUN = mean)
  colnames(test_mean_index) <- c("Model", "mean")
  test_mean_index$mean <- as.numeric(sprintf("%.2f", test_mean_index$mean))
  test_mean_index$Value <- "Mean value in validation cohorts"
  ##
  mean_index <- merge(
    all_mean_index,
    test_mean_index[, c("Model", "mean")],
    by = "Model",
    suffixes = c("", "_test"),
    all.x = TRUE)
  ##
  mean_index <- mean_index[
    order(
      mean_index$mean,
      mean_index$mean_test),]
  ##order the cohort data-frame
  index_cohort_df$Model <- factor(index_cohort_df$Model, levels = mean_index$Model)
  mean_cohort_df$Model <- factor(mean_cohort_df$Model, levels = mean_index$Model)
  ##index label name
  base_vars <- c("cindex", "bs")
  base_labels <- c("C-index", "Integrated Brier Score")
  km_years <- seq(1:10)
  km_vars <- paste0("km_auc_", km_years)
  km_labels <- paste0(km_years, "-year AUC")
  vars <- c(base_vars, km_vars)
  labels <- c(base_labels, km_labels)
  label_df <- data.frame(index = vars, index_label = labels)


  extract_index_label <- function(index) {
    if (index %in% label_df$index) {
      return(label_df[label_df$index == index, ]$index_label)
    }

    if (grepl("^km_auc_\\d+", index)) {
      year <- gsub("km_auc_", "", index)
      return(paste0(year, "-year AUC"))
    }

    warning("Unknown index: ", index)
    return(index)
  }

  ## heatmap plot
  generate_mean_heatmap <- function(mean_cohort_df, all_mean_index, index_label, outdir, color) {

    # Main heatmap
    p1 <- ggplot2::ggplot(
      mean_cohort_df,
      ggplot2::aes(
        x = .data$Dataset,
        y = .data$Model
      )
    ) +
      ggplot2::geom_tile(
        ggplot2::aes(fill = .data$index),
        color = "white",
        linewidth = 0.5
      ) +
      ggplot2::geom_text(
        ggplot2::aes(label = .data$index),
        vjust = 0.5,
        color = "black",
        size = 5
      ) +
      ggplot2::scale_fill_gradient2(
        low = "#0084A7",
        mid = "#F5FACD",
        high = "#E05D00",
        midpoint = median(mean_cohort_df$index),
        name = index_label
      ) +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(size = 14, color = "black"),
        panel.grid = ggplot2::element_blank(),
        panel.background = ggplot2::element_rect(fill = "white"),
        axis.text.x = ggplot2::element_blank()
      )

    # Cohort labels
    cohort_labels <- data.frame(
      Cohort = unique(mean_cohort_df$Dataset)
    )

    p2 <- ggplot2::ggplot(
      cohort_labels,
      ggplot2::aes(.data$Cohort, y = 1)
    ) +
      ggplot2::geom_tile(
        ggplot2::aes(fill = .data$Cohort),
        color = "white",
        size = 0.5
      ) +
      ggplot2::scale_fill_manual(values = color, name = "Cohort") +
      ggplot2::theme_void()

    # Side bar plots
    p3 <- ggplot2::ggplot(
      all_mean_index,
      ggplot2::aes(
        x = .data$Model,
        y = .data$mean,
        fill = .data$Value
      )
    ) +
      ggplot2::geom_bar(
        position = "dodge",
        stat = "identity"
      ) +
      ggplot2::scale_fill_manual(values = "#79AF97") +
      ggplot2::geom_text(
        ggplot2::aes(label = .data$mean),
        position = ggplot2::position_dodge(width = 0.9),
        vjust = 0.5,
        hjust = 1.2,
        size = 5
      ) +
      ggplot2::theme(
        axis.title = ggplot2::element_text(size = 14, color = "black"),
        axis.ticks.x = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5),
        panel.background = ggplot2::element_rect(fill = "white")
      ) +
      ggplot2::labs(y = "") +
      ggplot2::coord_flip()

    # combine plot
    p <- p1 %>%
      aplot::insert_top(p2, height = 0.01) %>%
      aplot::insert_right(p3, width = 0.25)

    # Save plot
    ggplot2::ggsave(
      filename = file.path(outdir, paste0("all_test_", index_label, "_heatmap.jpg")),
      plot = p,
      width = 4400,
      height = 4800,
      units = "px",
      dpi = 600
    )
    return(p)
  }
  ## bar plot
  ## barplot
  generate_bar_plot <- function(index_cohort,model,index,index_label,outdir){
      p<-ggbarplot(
        index_cohort,
        x = "Dataset",
        y = "index",
        width = 0.8,
        fill = "Dataset",
        color = "white",
        add="mean_se",
        add.params = list(color="black",size=0.5),
        palette = index_cohort$col,
        x.text.angle = 45
      )  +
        labs(
          x = NULL,
          y = index_label,
          title = paste0(model)
        ) + theme(
          legend.position = "none",
          axis.text.x = element_text(size = 10, color="black",angle = 45, hjust = 1),
          axis.text.y = element_text(size = 10,color="black")
        )+
        geom_hline(aes(yintercept = mean(index, na.rm = TRUE)), color = "black",
                   linetype = "dashed", size =0.4) +
        scale_x_discrete(
          limits = rev(index_cohort[order(index_cohort$index),]$Dataset))

      print(p)
      ggexport(p,filename=paste0(outdir,"/",model,".",index,".barplot.jpg"),width = 1700,height = 1700,res=600)

    }

  ## calculate statics for forestplot
  calculate_summary_stats <- function(data) {
    data %>%
      group_by(.data$Dataset) %>%
      summarise(
        mean_index = mean(.data$index, na.rm = TRUE),
        sd_index = ifelse(dplyr::n() == 1, 0, sd(.data$index, na.rm = TRUE)),
        n = dplyr::n(),
        se_index = .data$sd_index / ifelse(dplyr::n() == 1, 1, sqrt(.data$n))
      )
  }
  ## forestplot
  generate_forest_plot <- function(index_cohort_df,model,index,index_label,outdir) {
      df_summary <- calculate_summary_stats(index_cohort_df)
      mean_value <- mean(df_summary$mean_index)

      df_summary <- df_summary %>%
        dplyr::arrange(.data$mean_index) %>%
        dplyr::mutate(Dataset = factor(.data$Dataset, levels = unique(.data$Dataset)))

      colors <- rev(index_cohort$col)

      p <- ggplot2::ggplot(
        df_summary,
        ggplot2::aes(
          x = .data$mean_index,
          y = .data$Dataset,
          color = .data$Dataset
        )
      ) +
        ggplot2::geom_errorbarh(
          ggplot2::aes(
            xmin = .data$mean_index - .data$se_index,
            xmax = .data$mean_index + .data$se_index
          ),
          height = 0.2,
          linewidth = 0.4,
          color = "darkgrey"
        ) +
        ggplot2::geom_vline(
          xintercept = mean_value,
          linetype = "dashed",
          color = "black",
          linewidth = 0.5
        ) +
        ggplot2::geom_point(size = 3, color = colors) +
        ggplot2::theme_bw() +
        ggplot2::labs(x = index_label, y = "") +
        ggplot2::theme(
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_text(
            angle = 0,
            hjust = 1,
            color = "black",
            size = 10
          ),
          axis.text.x = ggplot2::element_text(
            angle = 0,
            hjust = 1,
            color = "black",
            size = 10
          )
        )

      # Save plot
      ggplot2::ggsave(
        plot = p,
        filename = file.path(outdir, paste0(model, ".", index, ".forestplot.jpg")),
        width = 9,
        height = 9,
        units = "cm",
        dpi = 600
      )

    }

  generate_horiz_forest_plot <-function(index_cohort_df,model,index,index_label,outdir){
      df_summary <- calculate_summary_stats(index_cohort_df)
      mean_value <- mean(df_summary$mean_index)

      df_summary <- df_summary %>%
        arrange(desc(.data$mean_index))%>%
        mutate(Dataset = factor(.data$Dataset, levels = unique(.data$Dataset)))
      colors <- index_cohort$col

      p<-ggplot(df_summary, aes(x = .data$Dataset,
                                y = .data$mean_index,
                                color=.data$Dataset))+
        geom_errorbar(aes(ymin = .data$mean_index - .data$se_index,
                           ymax = .data$mean_index + .data$se_index),
                       width = 0.2,
                       linewidth = 0.4,
                       color="darkgrey"
        ) +
        geom_hline(yintercept = mean_value, linetype = "dashed", color = "black", linewidth = 0.5) +
        geom_point(size = 3,color=colors)+
        theme_bw()+
        ggtitle(model)+
        labs(x ="",y=index_label)+
        theme( panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               axis.text.y = element_text(angle = 0, hjust = 1,color="black",size=10),
               axis.text.x = element_text(angle = 45, hjust = 1,color="black",size=10),
               plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
      print(p)
      ggsave(p,filename=paste0(outdir,"/",model,".",index,".horiz.forestplot.jpg"),width = 9,height = 9,units = "cm",dpi=600)

    }

  ###heatmap
  index_label=extract_index_label(index)
  generate_mean_heatmap(mean_cohort_df,all_mean_index,index_label,outdir,dataset_col)

  ###combine all model
  model_name<-unique(index_cohort_df$Model)

  ###
  dataset_col_df<-data.frame(Dataset=unique(index_cohort_df$Dataset),
                             col=dataset_col[1:length(unique(index_cohort_df$Dataset))])

  summary<-data.frame()
  for(i in 1:length(model_name)){
      model<-model_name[i]
      index_cohort=index_cohort_df[index_cohort_df$Model==model,]

      index_cohort<-index_cohort %>%
        left_join(dataset_col_df,by="Dataset")

      index_cohort$label<-index_label

      summary<-rbind(summary,index_cohort)
      ###output
      subdir=paste0(outdir,"/",model)
      if(!dir.exists(subdir)){
        dir.create(subdir,recursive = T)
      }

      ###
      if(plot_type=="barplot"){
        generate_bar_plot(index_cohort,model,index,index_label,subdir)
      } else if(plot_type=="forestplot"){
        generate_forest_plot(index_cohort,model,index,index_label,subdir)
      } else if(plot_type=="horiz_forestplot"){
        generate_horiz_forest_plot(index_cohort,model,index,index_label,subdir)
      }
  }

  utils::write.table(summary,file = paste0(outdir,"/",index,".summary.csv"),sep=",")
}
