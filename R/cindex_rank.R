#' C-index Distribution Visualization
#'
#' Creates ranked performance plots for survival models
#'
#' @param vali_auc_list List of models containing performance metrics of all cohorts
#' @param index #cindex,bs,km_auc_1, km_auc_3, km_auc_5, km_auc_10
#' @param train # the training cohort default "TCGA"
#' @param plot_type Visualization type ("boxplot", "forestplot", or "barplot")
#' @param outdir Output directory for saving plots
#' @return ggplot object or list of plots
#' @export
#' @examples
#' \dontrun{
#' cindex_rank(
#'   vali_auc_list = test.index,
#'   index="cindex",
#'   train="TCGA",
#'   plot_type = "boxplot",
#'   outdir="./"
#' )
#' }
cindex_rank <- function(vali_auc_list,
                        index, #cindex,bs,km_auc_1, km_auc_3, km_auc_5, km_auc_10
                        train="TCGA",
                        plot_type="barplot",#barplot, forestplot, horiz_forestplot
                        outdir="./"
){
  library(ggplot2)
  library(aplot)
  library(ggpubr)
  library(RColorBrewer)
  library(ggplot2)
  library(dplyr)
  library(forcats)
  library(survminer)
  library(survival)
  library(gridExtra)

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
  mean_cohort_df <- aggregate(
      x = index_cohort_df$index,
      by = list(index_cohort_df$Dataset,index_cohort_df$Model),
      FUN = mean)
  colnames(mean_cohort_df) <- c("Dataset","Model","index")
  mean_cohort_df$index <- as.numeric(sprintf("%.2f", mean_cohort_df$index))
  ##mean index in all cohorts
  all_mean_index <- aggregate(
      x = index_cohort_df$index,
      by = list(index_cohort_df$Model),
      FUN = mean)
  colnames(all_mean_index) <- c("Model", "mean")
  all_mean_index$mean <- as.numeric(sprintf("%.2f", all_mean_index$mean))
  all_mean_index$Value <- "Mean value in all cohorts"
  ##mean index in test cohorts
  test_cohort_df<-index_cohort_df[index_cohort_df$Dataset != train,]
  test_mean_index <- aggregate(
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
  # 定义基础变量和标签
  base_vars <- c("cindex", "bs")
  base_labels <- c("C-index", "Integrated Brier Score")

  # 动态生成km_auc相关变量名和标签
  km_years <- seq(1:10)  # 可扩展添加年份
  km_vars <- paste0("km_auc_", km_years)
  km_labels <- paste0(km_years, "-year AUC")

  # 合并所有变量和标签
  vars <- c(base_vars, km_vars)
  labels <- c(base_labels, km_labels)
  label_df <- data.frame(index = vars, index_label = labels)

  #增强版标签提取函数
  extract_index_label <- function(index) {
    # 先检查预定义标签
    if (index %in% label_df$index) {
      return(label_df[label_df$index == index, ]$index_label)
    }

    # 动态处理未预定义的km_auc_*指标
    if (grepl("^km_auc_\\d+", index)) {
      year <- gsub("km_auc_", "", index)
      return(paste0(year, "-year AUC"))
    }

    # 未知指标返回原名称
    warning("Unknown index: ", index)
    return(index)
  }

  ## heatmap plot
  generate_mean_heatmap<-function(mean_cohort_df,all_mean_index,index_label,outdir,color){
    p1 <- ggplot(mean_cohort_df, aes(x = Dataset, y = Model)) +
      geom_tile(aes(fill = index), color = "white", linewidth = 0.5) +
      geom_text(aes(label = index), vjust = 0.5, color = "black", size = 5) +
      scale_fill_gradient2(low = "#0084A7", mid = "#F5FACD", high = "#E05D00",
                           midpoint = median(mean_cohort_df$index),
                           name = index_label) +
      theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size =14,color="black"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_blank())
    ##cohort dataframe
    cohort_labels <- as.data.frame(unique(mean_cohort_df$Dataset))
    colnames(cohort_labels) <- "Cohort"

    p2 <- ggplot(cohort_labels, aes(Cohort, y = 1)) +
      geom_tile(aes(fill = Cohort), color = "white", size = 0.5) +
      scale_fill_manual(values = dataset_col, name = "Cohort") +
      theme_void()

    p3 <- ggplot(data = all_mean_index, aes(x = Model, y = mean, fill = Value)) +
      geom_bar(position = "dodge", stat = "identity") +
      scale_fill_manual(values =  "#79AF97") +
      geom_text(aes(label = mean), position = position_dodge(width = 0.9), vjust = 0.5, hjust = 1.2, size = 5) +
      theme(
        axis.title = element_text(size = 14,color="black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white")) + labs(y = "") + coord_flip()

    p4 <- ggplot(data = test_mean_index, aes(x = Model, y = mean, fill = Value)) +
      geom_bar(position = "dodge", stat = "identity") +
      scale_fill_manual(values =  "#8491B4") +
      geom_text(aes(label = mean), position = position_dodge(width = 0.9), vjust = 0.5, hjust = 1.2, size = 5) +
      theme(
        axis.title = element_text(size = 14,color="black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white")) + labs(y = "") + coord_flip()

    p<-print(p1 %>%
            insert_top(p2, height = 0.01) %>%
            insert_right(p3, width = 0.25) %>%
              insert_right(p4, width = 0.25)
              )
    ggexport(p,filename =paste0(outdir,"/all_test_",index,"_heatmap.jpg"),width = 4400,height =4800,res = 600)
    }

  ## barplot
  generate_bar_plot <- function(index_cohort,model,index,index_label,outdir){
      # 绘制柱状图
      p<-ggbarplot(
        index_cohort,
        x = "Dataset",
        y = "index",
        width = 0.8,
        fill = "Dataset", # 按模型分配颜色
        color = "white", # 边框颜色
        add="mean_se",
        add.params = list(color="black",size=0.5),
        palette = index_cohort$col,
        x.text.angle = 45 # X轴标签倾斜
      )  +
        labs(
          x = NULL,
          y = index_label,
          title = paste0(model)
        ) + theme(
          legend.position = "none", # 去掉图例
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
  calculate_summary_stats <- function(data, index) {
    data %>%
      group_by(Dataset) %>%
      summarise(
        mean_index = mean(index, na.rm = TRUE),
        sd_index = ifelse(n() == 1, 0, sd(index, na.rm = TRUE)),
        n = n(),
        se_index = sd_index / ifelse(n() == 1, 1, sqrt(n()))  # 当n=1时，分母设为1保持se_index=0
      )
  }
  ## forestplot
  generate_forest_plot <- function(index_cohort,model,index,index_label,outdir) {
      df_summary <- calculate_summary_stats(index_cohort,index)
      mean_value <- mean(df_summary$mean_index)

      df_summary <- df_summary %>%
        arrange(mean_index) %>%
        mutate(Dataset = factor(Dataset, levels = unique(Dataset)))

      colors <- rev(index_cohort$col)

      p<-ggplot(df_summary, aes(x = mean_index, y = Dataset, color = Dataset)) +
        geom_errorbarh(aes(xmin = mean_index - se_index, xmax = mean_index + se_index),
                       height = 0.2, linewidth = 0.4, color = "darkgrey") +
        geom_vline(xintercept = mean_value, linetype = "dashed", color = "black", linewidth = 0.5) +
        geom_point(size = 3, color = colors) +
        theme_bw() +
        labs(x = index_label, y = "") +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text.y = element_text(angle = 0, hjust = 1, color = "black", size = 10),
              axis.text.x = element_text(angle = 0, hjust = 1, color = "black", size = 10))
      print(p)
      ggsave(p,filename=paste0(outdir,"/",model,".",index,".forestplot.jpg"),width = 9,height = 9,units = "cm",dpi=600)

    }

  generate_horiz_forest_plot <-function(index_cohort,model,index,index_label,outdir){
      df_summary <- calculate_summary_stats(index_cohort,index)
      mean_value <- mean(df_summary$mean_index)

      df_summary <- df_summary %>%
        arrange(desc(mean_index))%>%
        mutate(Dataset = factor(Dataset, levels = unique(Dataset)))
      colors <- index_cohort$col

      p<-ggplot(df_summary, aes(x = Dataset, y = mean_index,color=Dataset))+
        geom_errorbar(aes(ymin = mean_index - se_index,
                           ymax = mean_index + se_index),  # 水平误差条
                       width = 0.2,      # 水平误差条的高度
                       linewidth = 0.4,
                       color="darkgrey"# 水平误差条的线条宽度
        ) +
        geom_hline(yintercept = mean_value, linetype = "dashed", color = "black", linewidth = 0.5) +
        geom_point(size = 3,color=colors)+
        theme_bw()+ # 自定义颜色
        ggtitle(model)+
        labs(x ="",y=index_label)+
        theme( panel.grid.major = element_blank(),  # 去除主要网格线
               panel.grid.minor = element_blank(),  # 去除次要网格线
               axis.text.y = element_text(angle = 0, hjust = 1,color="black",size=10),
               axis.text.x = element_text(angle = 45, hjust = 1,color="black",size=10),
               plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
      print(p)
      ggsave(p,filename=paste0(outdir,"/",model,".",index,".horiz.forestplot.jpg"),width = 9,height = 9,units = "cm",dpi=600)

    }

  ###heatmap
  index_label=extract_index_label(index)
  generate_mean_heatmap(mean_cohort_df,all_mean_index,index_label,outdir)

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

  write.table(summary,file = paste0(outdir,"/",index,".summary.csv"),sep=",")
}
