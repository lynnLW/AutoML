#' Survival Analysis Visualization
#'
#' @description Generates Kaplan-Meier survival curves for model validation results
#' @param vali_auc_list List object containing validation results with survival data
#' @param model Specific model name to plot or "all" for all models (default: "all")
#' @param optim Logical indicating whether to show optimal risk cutpoint (default: TRUE)
#' @param outdir Output directory path (optional)
#' @param color Custom color vector for risk groups (length should match number of groups)
#' @return ggplot object or list of ggplot objects
#' @export
surv_plot <- function(vali_auc_list,
                      model="all", ## modelname or all
                      optim=T,
                      outdir=NULL,
                      color = NULL # color value for cohort
){

  if (is.null(color) == T) {
    color <- c("#0084A7", "#F5FACD", "#E05D00", "#79AF97", "#8491B4") ## default color value
  } else {
    color <- color
  }

  if (is.null(model) == T) {
    model <- "all"
  } else {
    model <- model
  }

  ###time dependent plot
  generate_surv_plot<-function(index_df,model_name,cohort_name,subdir,best_cutoff){
    if(is.null(best_cutoff)){
      index_df$group<-ifelse(index_df$pred>median(index_df$pred),"High","Low")
    } else {
      index_df$group<-ifelse(index_df$pred>best_cutoff,"High","Low")
    }
    ###output

    fit<-survival::survfit(Surv(time,status)~group, data=index_df)
    ###time dependent roc curve

    p<-survminer::ggsurvplot(fit,data=index_df,
                  title = paste0(cohort_name),
                  font.title = c(10,"black"),
                  break.y.by=0.2,
                  axes.offset=FALSE,
                  legend.title = 'Risk Score',
                  legend="right",
                  pval=TRUE,
                  pval.size = 2,
                  pval.color = "black",
                  pval.coord = c(4, 0.06),
                  risk.table = F,
                  palette=c("#925E9fff","#42B540ff"),
                  xlab="Time(months)",
                  ggtheme = theme_bw()+
                    theme(  panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(),
                            legend.background = element_rect(fill = "transparent"),
                            legend.key.size = unit(1, "mm"),
                            legend.spacing.y = unit(1, "mm"),
                            legend.margin = margin(1, 2, 2, 2),
                            legend.box.spacing = unit(3, "mm"),
                            legend.text = element_text(color = "black", size = 8),
                            legend.title = element_text(color = "black", size = 8),
                            axis.title.y = element_text(color = "black", size = 10),
                            axis.title.x = element_text(color = "black", size = 10),
                            axis.text.y = element_text(angle = 0, hjust = 1,color="black",size=10),
                            axis.text.x = element_text(angle = 0, hjust = 1,color="black",size=10),
                            panel.border = element_rect(color = "black", size = 0.8),
                            plot.title = element_text(hjust = 0, vjust = 1)))
    print(p)
    ggpubr::ggexport(p,filename=paste0(subdir,"/",cohort_name,"_surv_plot.jpg"),width = 2000,height = 1400,res = 600)
  }

  optimal_cutoff <- function(data, time, status, pred) {

    cutpoint <- survminer::surv_cutpoint(
      data,
      time = time,
      event = status,
      variables = pred
    )
    return(cutpoint$cutpoint$cutpoint)
  }

  ###
  if (model=="all"){
    for (i in 1:length(vali_auc_list)){
      models<-vali_auc_list[[i]]
      model_name<-names(vali_auc_list)[i]
      print(model_name)

      ###
      subdir=paste0(outdir,"/",model_name)
      if(!dir.exists(subdir)){
        dir.create(subdir,recursive = T)
      }

      for (j in 1:length(models)){
        cohort<-models[[j]]
        cohort_name<-names(models)[j]
        print(cohort_name)
        index_df<-cohort[[1]]$pred_df

        if (optim){

          best_cutoff <- optimal_cutoff(
            data = index_df,
            time = "time",
            status = "status",
            pred = "pred"
          )
          print(best_cutoff)
        } else {
          best_cutoff=NULL
        }
        generate_surv_plot(index_df,model_name,cohort_name,subdir,best_cutoff)
      }
    }
  } else {
    models<-vali_auc_list[[model_name]]
    ###
    subdir=paste0(outdir,"/",model_name)
    if(!dir.exists(subdir)){
      dir.create(subdir,recursive = T)
    }

    for (j in 1:length(models)){
      cohort<-models[[j]]
      cohort_name<-names(models)[j]
      index_df<-cohort[[1]]$pred_df
      ###time dependent roc curve
      if (optim){

        best_cutoff <- optimal_cutoff(
          data = index_df,
          time = "time",
          status = "status",
          pred = "pred"
        )
        print(best_cutoff)
      } else {
        best_cutoff=NULL
      }
      generate_surv_plot(index_df,model_name,cohort_name,subdir,best_cutoff)
    }
  }
}
