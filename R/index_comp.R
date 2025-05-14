#' Compare Model Performance with Published Results
#'
#' Creates a visualization comparing AUC distributions between models
#'
#' @param own_auc_list Output from cal_vali_index (list of AUC values)
#' @param published_auc_list Optional output from cal_vali_index (published results)
#' @param model_name Specific model name to compare from published_auc_list
#' @param dataset_col Optional color vector for datasets (default uses ggplot colors)
#' @param dataset Character vector of dataset names
#' @param index Metrics to compare ("cindex", "AUC_1", "AUC_3", etc.)
#' @import ggplot2
#'
#' @return ggplot2 object showing performance comparison
#' @export
#'
index_comp <- function(own_auc_list,
                     published_auc_list=NULL,
                     model_name, ## input specific model name
                     dataset_col = NULL, # color value for cohort
                     dataset, # input datasets name
                     index #Cindex, AUC_1,3,5,7
) {

  if (is.null(dataset_col) == T) {
    dataset_col <- c(
      "#3182BDFF", "#E6550DFF", "#31A354FF", "#756BB1FF", "#636363FF", "#6BAED6FF", "#FD8D3CFF", "#74C476FF",
      "#9E9AC8FF", "#969696FF", "#9ECAE1FF", "#FDAE6BFF", "#A1D99BFF", "#BCBDDCFF", "#BDBDBDFF", "#C6DBEFFF",
      "#FDD0A2FF", "#C7E9C0FF", "#DADAEBFF", "#D9D9D9FF"
    ) ## default 20 color values
  } else {
    dataset_col <- dataset_col
  }


  index_df <- data.frame()
  for (i in names(own_auc_list)) {
    tmp <- own_auc_list[[i]]
    tmp4<-data.frame()
    for (n in names(tmp)) {
      tmp2 <- tmp[[n]]
      tmp3<-data.frame()
      for (j in 1:length(tmp2)) {
        data<-tmp2[[j]][[1]]
        tmp3 <- rbind(tmp3, data.frame(
          "Cindex" = unique(data[["cindex"]]),
          "AUC_1" = unique(data[["km_auc_1"]]),
          "AUC_2" = unique(data[["km_auc_2"]]),
          "AUC_3" = unique(data[["km_auc_3"]]),
          "AUC_5" = unique(data[["km_auc_5"]]),
          "AUC_7" = unique(data[["km_auc_7"]]),
          "AUC_10" = unique(data[["km_auc_10"]]),
          "ID"=names(tmp2)[j]))
      }
      tmp3$Model <- n
      tmp4 <- rbind(tmp4, tmp3)
    }
    tmp4$Signature <- i
    index_df <- rbind(index_df, tmp4)
  }

  index_df$text <- "red"

  if(!is.null(published_auc_list)){
    index_df2 <- data.frame()
    for (i in names(published_auc_list)) {
      tmp <- published_auc_list[[i]]
      tmp4<-data.frame()
      for (n in names(tmp)) {
        tmp2 <- tmp[[n]]
        tmp3<-data.frame()
        for (j in 1:length(tmp2)) {
          data<-tmp2[[j]][[1]]
          tmp3 <- rbind(tmp3, data.frame(
            "Cindex" = unique(data[["cindex"]]),
            "AUC_1" = unique(data[["km_auc_1"]]),
            "AUC_2" = unique(data[["km_auc_2"]]),
            "AUC_3" = unique(data[["km_auc_3"]]),
            "AUC_5" = unique(data[["km_auc_5"]]),
            "AUC_7" = unique(data[["km_auc_7"]]),
            "AUC_10" = unique(data[["km_auc_10"]]),
            "ID"=names(tmp2)[j]))
        }
        tmp3$Model <- n
        tmp4 <- rbind(tmp4, tmp3)
      }
      tmp4$Signature <- i
      index_df2 <- rbind(index_df2, tmp4)
    }
    index_df2$text <- "black"
  }

  if(!is.null(published_auc_list)){
    index_df <- rbind(index_df, index_df2)
  }
  ###
  generate_p<-function(index_df_select,index){
    ggplot(index_df_select, aes(x =.data$Signature, y = .data[[index]])) +
      geom_segment(aes(x = .data$Signature, xend = .data$Signature, y = 0, yend = .data[[index]], color = .data$ID)) +
      geom_point(aes(color = .data$ID)) +
      scale_color_manual(values = dataset_col[t], name = "Cohort") +
      theme(
        panel.grid = element_blank(),
        axis.text.y.left = (element_text(color = index_df_select$text,size = 8)),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
        legend.position = "",
        axis.text.x=(element_text(color = "black",size = 8)),
        plot.title = element_text(hjust = 0.5,color = "black",size = 10),
        panel.background = element_rect(fill = "white")
      ) +
      # scale_y_continuous(position = "right")+
      labs(y = index, x = "", title = "") +
      coord_flip()
  }

  index_df<-index_df[index_df$Model==model_name,]
  plot_list <- list()
  for (t in 1:length(dataset)) {
    index_df_select <- index_df[index_df$ID == dataset[t], ]
    index_df_select <- index_df_select[order(index_df_select[[index]], decreasing = F), ]
    index_df_select$Signature <- factor(index_df_select$Signature, levels = unique(index_df_select$Signature))

    plot_list[[dataset[t]]] <-generate_p(index_df_select,index=index)
  }

  p1 <- aplot::plot_list(gglist = plot_list, ncol = length(dataset))
  print(p1)


}
