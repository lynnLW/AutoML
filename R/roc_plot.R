#' Generate ROC Curve Plots for Model Validation
#'
#' @description Creates ROC curve visualizations for survival model validation results
#' @param vali_auc_list List containing validation AUC results from model training
#' @param model Specific model name to plot or "all" for all models (default: "all")
#' @param cohort Optional cohort identifier for multi-cohort analysis
#' @param auc_time Specific time point for time-dependent ROC analysis
#' @param outdir Output directory path for saving plots
#' @param colors Custom color palette for cohorts (length should match number of cohorts)
#' @return ggplot object or list of ggplot objects
#' @export
#' @examples
#' \dontrun{
#' roc_plot(vali_auc_list = validation_results,
#'          model = "CoxBoost",
#'          outdir = "./plots")
#' }
roc_plot <- function(vali_auc_list,
                     model="all", ## modelname or all
                     cohort=NULL,
                     auc_time=NULL,
                     outdir=NULL,
                     colors = NULL # color value for cohort
){

  if (is.null(colors) == T) {
    colors <- c(RColorBrewer::brewer.pal(9, "Set1"),
                     "#3182BDFF", "#E6550DFF", "#31A354FF", "#756BB1FF", "#636363FF", "#6BAED6FF", "#FD8D3CFF", "#74C476FF",
                     "#9E9AC8FF", "#969696FF", "#9ECAE1FF", "#FDAE6BFF", "#A1D99BFF", "#BCBDDCFF", "#BDBDBDFF", "#C6DBEFFF",
                     "#FDD0A2FF", "#C7E9C0FF", "#DADAEBFF", "#D9D9D9FF") ## default 20 color values
  } else {
    colors <- colors
  }

  ###time dependent plot
  generate_time_roc_plot<-function(index_df,model_name,cohort_name,subdir){
    ###
    roc<-timeROC::timeROC(T=index_df$time,
                          delta=index_df$status,
                          marker = index_df$pred,
                          cause = 1,
                          weighting = "marginal",
                          times=c(1,3,5)*12,
                          iid=T)

    jpeg(filename=paste0(subdir,"/",cohort_name,"_time_roc_plot.jpg"),width = 10, height = 10, units = "cm", res = 600)
    plot(roc,time=1*12,col="#31A354FF",lty=1,lwd=2,title="")
    plot(roc,time=3*12,col="#3182BDFF",add=T,lty=1,lwd=2)
    plot(roc,time=5*12,col="#E6550DFF",add=T,lty=1,lwd=2)
    ##
    legend("bottomright",
           c(paste0("1-year AUC: ",round(roc[["AUC"]][1],2)),
             paste0("3-year AUC: ",round(roc[["AUC"]][2],2)),
             paste0("5-year AUC: ",round(roc[["AUC"]][3],2))),
           col=c("#31A354FF","#3182BDFF","#E6550DFF"),
           lty=1,lwd=2,bty="n")
    title(paste0(cohort_name))
    dev.off()
  }

  ###cohort dependet plot
  generate_cohort_roc_plot<-function(dataset,dataset_name,auc_time,subdir,dataset_col_df){
    roc_data <- data.frame()
    auc_df<-data.frame()
    for(n in 1:length(dataset)) {
      feature_data <- dataset[[n]][[1]]
      feature_name<-names(dataset)[n]

      temp_df <- data.frame(
        fpr = feature_data[[paste0('km_fp_',auc_time)]],
        tpr = feature_data[[paste0('km_tp_',auc_time)]],
        feature = feature_name
      )
      roc_data <- rbind(roc_data, temp_df)
      auc_df<- rbind(auc_df,data.frame(feature=feature_name,
                                       auc=round(feature_data[[paste0('km_auc_',auc_time)]],3)))
    }

    feature_order <- auc_df %>%
      dplyr::distinct(feature, auc) %>%
      dplyr::arrange(desc(auc)) %>%
      dplyr::mutate(feature_with_auc = sprintf("%s (AUC=%.2f)", feature, auc)) %>%
      dplyr::left_join(dataset_col_df, by = c("feature" = "Feature"))


    roc_data <- roc_data %>%
      dplyr::left_join(feature_order, by = c("feature")) %>%
      dplyr::mutate(feature = factor(feature, levels = feature_order$feature))


    p <-ggplot2::ggplot(roc_data, aes(x = fpr, y = tpr,
                              color = feature)) +
      geom_line(size = 1) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
      labs(
        title = paste0(auc_time,"-year ROC Curve for ", dataset_name),
        x = "False Positive Rate",
        y = "True Positive Rate",
        color = "Features (AUC)"
      ) +
      theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black", size = 8),
        axis.title = element_text(color = "black", size = 10),
        title = element_text(color = "black", size = 10),
        legend.position = c(0.98, 0.02),
        legend.justification = c(1, 0),
        legend.box.background = element_rect(color = "black", fill = "white"),
        legend.margin = margin(1, 1, 1, 1),
        legend.box.spacing = unit(0, "mm"),
        legend.key.size = unit(2, "mm"),
        legend.spacing.x = unit(0.1, "mm"),
        legend.spacing.y = unit(0.1, "mm"),
        legend.text = element_text(color = "black", size =8),
        legend.title = element_text(color = "black", size = 8)
      ) +
      scale_color_manual(
        values = feature_order$col,
        breaks = feature_order$feature,
        labels=feature_order$feature_with_auc,
        guide = guide_legend(
          direction = "vertical",
          title.position = "top",
          ncol = 1,
          keyheight = unit(4, "mm")
        )
      )

    print(p)
    ggplot2::ggsave(p,filename=paste0(subdir,"/",auc_time, "_feature_ROC.jpg"), width = 9, height = 9, dpi = 600,units = "cm")
  }

  ###
  if(!is.null(model)){
    if (model=="all"){
      for (i in 1:length(vali_auc_list)){
        models<-vali_auc_list[[i]]
        model_name<-names(vali_auc_list)[i]
        print(model_name)
        ###output
        subdir=paste0(outdir,"/",model_name)
        if(!dir.exists(subdir)){
          dir.create(subdir,recursive = T)
        }
        for (j in 1:length(models)){
          dataset<-models[[j]]
          dataset_name<-names(models)[j]
          index_df<-dataset[[1]]$pred_df
          ###time dependent roc curve
          generate_time_roc_plot(index_df,model_name,dataset_name,subdir)
        }
      }
  } else {
    models<-vali_auc_list[[model]]
    for (j in 1:length(models)){
      dataset<-models[[j]]
      dataset_name<-names(models)[j]
      index_df<-dataset[[1]]$pred_df
      ###time dependent roc curve
      generate_time_roc_plot(index_df,model_name,dataset_name,outdir)
    }
  }}

  ###
  if(!is.null(cohort)){
    if (cohort=="all"){
      ###setting features colors
      all_features <- unique(unlist(lapply(vali_auc_list, function(dataset) names(dataset))))
      dataset_col_df<-data.frame(Feature=all_features,
                                 col=colors[1:length(all_features)])

      for (i in 1:length(vali_auc_list)){
        dataset<-vali_auc_list[[i]]
        dataset_name<-names(vali_auc_list)[i]
        print(dataset_name)
        ###output
        subdir=paste0(outdir,"/",dataset_name)
        if(!dir.exists(subdir)){
          dir.create(subdir,recursive = T)
        }

        ###cohort dependent roc curve
        generate_cohort_roc_plot(dataset,dataset_name,auc_time,subdir,dataset_col_df)
      }
    } else {
      dataset<-vali_auc_list[[cohort]]
      dataset_name<-names(vali_auc_list)[i]
      print(dataset_name)
      ###output
      subdir=paste0(outdir,"/",dataset_name)
      if(!dir.exists(subdir)){
        dir.create(subdir,recursive = T)
      }
      generate_cohort_roc_plot(dataset,dataset_name,auc_time,subdir)

    }

  }

}
