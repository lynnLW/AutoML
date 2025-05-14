#' C-index Distribution Visualization
#'
#' Creates ranked performance plots for survival models
#'
#' @param index_list Optional pre-computed performance metrics list
#' @param model_list List of models containing performance metrics
#' @param order Metric to use for ordering ("train", "valid", or "test")
#' @param index Metrics to include ("all" or specific metrics)
#' @param outdir Output directory for saving plots
#' @param plot_type Visualization type ("boxplot", "forestplot", or "barplot")
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @importFrom dplyr mutate arrange desc left_join %>% summarise
#' @importFrom stats sd median
#' @importFrom ggpubr ggboxplot ggbarplot ggexport
#' @import ggplot2
#' @return ggplot object or list of plots
#' @export
cindex_rank2 <- function(index_list=NULL,
                         model_list=NULL,
                         order="valid", # train, valid, test
                         index="all", #all,cindex, bs, km_auc_1, km_auc_3,km_auc_5,km_auc_7,km_auc_10,
                         outdir="./",
                         plot_type="boxplot", # boxplot; forestplot; barplot
                         width = 1900, # width of plot
                         height = 1700 # height of plot
){

  if (!dir.exists(outdir)){
      dir.create(outdir,recursive = T)
  }

  ## default 20 color values
  model_col <- c(RColorBrewer::brewer.pal(9, "Set1"),
      "#3182BDFF", "#E6550DFF", "#31A354FF", "#756BB1FF", "#636363FF", "#6BAED6FF", "#FD8D3CFF", "#74C476FF",
      "#9E9AC8FF", "#969696FF", "#9ECAE1FF", "#FDAE6BFF", "#A1D99BFF", "#BCBDDCFF", "#BDBDBDFF", "#C6DBEFFF",
      "#FDD0A2FF", "#C7E9C0FF", "#DADAEBFF", "#D9D9D9FF")

  if (is.null(index_list)==F){

      ##show cindex in train data to selected best model
      model_cindex_df<-data.frame()

      ##convert all index data into one dataframe
      for(i in 1:length(index_list)){
          name=names(index_list)[i]
          df<-index_list[[i]]
          df<-data.frame(df,model=name)
          model_cindex_df<-rbind(model_cindex_df,df)
      }
      ##calculate the mean
      utils::write.table(model_cindex_df,file=paste0(outdir,"/rep10_index.csv"),sep=",")
      mean_table <- model_cindex_df %>%
        tidytable::pivot_longer(cols = -.data$model, names_to = "variable", values_to = "value") %>% # Convert to long format
        dplyr::group_by(.data$model, .data$variable) %>% # Group by model and variable
        dplyr::summarise(mean_value = mean(.data$value, na.rm = TRUE), .groups = "drop") %>% # Calculate mean
        tidytable::pivot_wider(names_from = "variable", values_from = "mean_value") # Convert back to wide format
      utils::write.table(mean_table ,file=paste0(outdir,"/rep10_mean_index.csv"),sep=",")
      ###order
        if (order=="train"){
          mean_table<-mean_table[order(mean_table$train.cindex,decreasing = T),]
          model_cindex_df$model <- factor(model_cindex_df$model,
                                          levels = mean_table$model)
        } else if (order=="valid") {
          mean_table<-mean_table[order(mean_table$valid.cindex,decreasing = T),]
          model_cindex_df$model <- factor(model_cindex_df$model,
                                          levels = mean_table$model)
        } else if (order=="test") {
          mean_table<-mean_table[order(mean_table$test.cindex,decreasing = T),]
          model_cindex_df$model <- factor(model_cindex_df$model,
                                          levels = mean_table$model)
        }

      ## calculate statics for forestplot
      calculate_summary_stats <- function(data, index) {
        index_sym <- rlang::sym(index)
        data %>%
          group_by(.data$model) %>%
          summarise(
            mean_index = mean(!!index_sym, na.rm = TRUE),
            sd_index = ifelse(dplyr::n() == 1, 0, sd(!!index_sym, na.rm = TRUE)),
            n = dplyr::n(),
            se_index = .data$sd_index / ifelse(dplyr::n() == 1, 1, sqrt(.data$n))
          )
      }

      ###index plot
      generate_box_plot <- function(model_cindex_df,var,label,width,height,outdir) {
        p <- ggboxplot(model_cindex_df,
                         x = "model",
                         y = var,
                         palette = model_col[1:length(unique(model_cindex_df[["model"]]))],
                         color = "model",
                         legend = "none",
                         ylab = label,
                         xlab = "",
                         width = 0.6,
                         add = "jitter",
                         add.params = list(size = 0.5),
                         font.label = list(size = 10, color = "black"),
                         label.rectangle = TRUE) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        print(p)
        ggexport(p,filename=paste0(outdir,"/",var,".boxplot.jpg"),width = width,height =height,res=600)
      }

      ###index plot
      generate_bar_plot <- function(model_cindex_df,var,label,width,height,outdir) {
        # bar plot
        p<-ggbarplot(
          model_cindex_df,
          x = "model",
          y = var,
          width = 0.8,
          fill = "model",
          color = "white",
          add="mean_sd",
          add.params = list(color="dimgray",
                            size=0.2),
          palette =  model_col[1:length(unique(model_cindex_df[["model"]]))],
          x.text.angle = 45
        )  +
          labs(
            x = NULL,
            y = label
          ) + theme(
            legend.position = "none",
            axis.text.x = element_text(size = 10, color="black",angle = 45, hjust = 1),
            axis.text.y = element_text(size = 10,color="black")
          )

        print(p)
        ggexport(p,filename=paste0(outdir,"/",var,".barplot.jpg"),width = width,height = height,res=600)
      }

      # Function to generate AUC plot
      generate_forest_plot <- function(model_cindex_df, var, label,width,height,outdir) {
        df_summary <- calculate_summary_stats(model_cindex_df, var)
        median_value <- median(df_summary$mean_index)
        df_summary$model <- forcats::fct_rev(df_summary$model)

        colors <- model_col[1:length(unique(model_cindex_df$model))]

        p<-ggplot(df_summary, aes(x = .data$mean_index,
                                  y = .data$model,
                                  color = .data$model)) +
          geom_errorbarh(aes(xmin = .data$mean_index - .data$se_index,
                             xmax = .data$mean_index + .data$se_index),
                         height = 0.15, linewidth = 0.3, color = "dimgrey") +
          geom_vline(xintercept = median_value, linetype = "dashed", color = "red", linewidth = 0.5) +
          geom_point(size =3.5, color = colors) +
          theme_bw() +
          labs(x = label, y = "") +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text.y = element_text(angle = 0, hjust = 1, color = "black", size = 10),
                axis.text.x = element_text(angle = 0, hjust = 1, color = "black", size = 10))
        print(p)
        ggexport(p,filename=paste0(outdir,"/",var,".forestplot.jpg"),width = width,height =height,res=600)
      }

      # index label
      base_vars <- c("cindex", "bs")
      base_labels <- c("C-index", "Integrated Brier Score")
      km_years <- seq(1:10)
      km_vars <- paste0("km_auc_", km_years)
      km_labels <- paste0(km_years, "-year AUC")
      vars <- c(base_vars, km_vars)
      labels <- c(base_labels, km_labels)
      label_df <- data.frame(index = vars, index_label = labels)

      # function to extract index label
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

      vars<-colnames(model_cindex_df)[1:(ncol(model_cindex_df)-1)]
      # generat for the index
      if(index=="all" && plot_type=="boxplot"){
        for (i in 1:(ncol(model_cindex_df)-1)){
          var=vars[i]
          var.1=stringr::str_split(vars[i],"\\.")[[1]][2]
          label=extract_index_label(var.1)
          generate_box_plot(model_cindex_df,var,label,width,height,outdir)
          generate_forest_plot(model_cindex_df,var,label,width,height,outdir)
        }
      } else if(index=="bs"){
          ids=vars[grep("bs",vars)]
          for(var in ids){
            if(plot_type=="boxplot"){
              generate_box_plot(model_cindex_df,var,"Integrated Brier Score",width,height,outdir)
            } else if (plot_type=="forestplot"){
              generate_forest_plot(model_cindex_df,var,"Integrated Brier Score",width,height,outdir)
            } else if (plot_type=="barplot"){
              generate_bar_plot(model_cindex_df,var,"Integrated Brier Score",width,height,outdir)
            }
          }
      } else if(index=="auc"){
        ids=vars[grep("auc",vars)]
        for(var in ids){
          y=stringr::str_split(var,"_",simplify=T)[3]
          if(plot_type=="boxplot"){
            generate_box_plot(model_cindex_df,var,paste0(y,"-year AUC"),width,height,outdir)
          } else if (plot_type=="forestplot"){
            generate_forest_plot(model_cindex_df,var,paste0(y,"-year AUC"),width,height,outdir)
          } else if (plot_type=="barplot"){
            generate_bar_plot(model_cindex_df,var,paste0(y,"-year AUC"),width,height,outdir)
          }
        }
      } else if(index=="cindex"){
        ids=vars[grep("cindex",vars)]
        for(var in ids){
          if(plot_type=="boxplot"){
            generate_box_plot(model_cindex_df,var,"C-index",width,height,outdir)
          } else if (plot_type=="forestplot"){
            generate_forest_plot(model_cindex_df,var,"C-index",width,height,outdir)
          } else if (plot_type=="barplot"){
            generate_bar_plot(model_cindex_df,var,"C-index",width,height,outdir)
          }
        }
      }
  }

  ##show cindex in train data to selected best model
  if (is.null(model_list)==F){
     model_cindex_list<-lapply(model_list,function(model){
      return(model$valid$cindex)
    })
     model_bs_list<-lapply(model_list,function(model){
      return(model$valid$bs)
    })
    model_cindex_df<-as.data.frame(t(as.data.frame(model_cindex_list)))
    names(model_cindex_df)[1]<-"cindex"
    model_cindex_df$model<-row.names(model_cindex_df)
    model_cindex_df<-model_cindex_df[order(model_cindex_df$cindex),]

    p1<-ggbarplot(model_cindex_df,
                  x = "model",
                  y = "cindex",
                  fill="model",
                  palette = rev(model_col[1:length(unique(model_cindex_df[["model"]]))]),
                  color = "model",
                  legend="none",
                  ylab="C-index",
                  xlab="",
                  orientation = "horiz")
    print(p1)
    ggexport(p1,filename = paste0(outdir,"/final_all_models_cindex.jpg"),width = 2000,height = 2200,res=600)

    model_bs_df<-as.data.frame(t(as.data.frame(model_bs_list)))
    names(model_bs_df)[1]<-"bs"
    model_bs_df$model<-row.names(model_bs_df)
    model_bs_df<-model_bs_df[row.names(model_cindex_df),]
    p2<-ggbarplot(model_bs_df,
                  x = "model",
                  y = "bs",
                  fill="model",
                  palette = model_col[1:length(unique(model_cindex_df[["model"]]))],
                  color = "model",
                  legend="none",
                  order = rev(model_cindex_df$model),
                  ylab="Integrated Brier Score",
                  xlab="",
                  lab.vjust = 45) + theme(
      legend.position = "none",
      axis.text.x = element_text(size = 10, color="black",angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10,color="black")
    )
    print(p2)
    ggexport(p2,filename =paste0(outdir,"/final_all_models_bs.jpg"),width = 2400,height = 2200,res=600)

    p<-ggarrange(p1, p2,
              ncol = 2,
              nrow = 1)
    print(p)
    ggexport(p,filename =paste0(outdir,"/final_all_models_cindex_bs.jpg"),width = 4200,height = 2600,res=600)
  }

}
