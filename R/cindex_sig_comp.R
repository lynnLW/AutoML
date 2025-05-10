#' Compare the AUC of our model with previously published models
#'
#' Creates a distribution plot of AUC among different mdoels
#'
#' @param own_auc_list Output of function cal_vali_index
#' @param published_auc_list Output of function cal_vali_index
#' @param model_name Model name used to compare
#' @param dataset A vector of names for all datasets
#' @return a list
#'
indexC_comp <- function(own_auc_list,
                        published_auc_list=NULL,
                        model_name,
                        dataset
) {
  library(compareC)
  df_index <- list()
  ###
  own_auc_list<-lapply(own_auc_list, function(x){x[[model_name]]})
  for(i in 1:length(own_auc_list)){
    feature_name=names(own_auc_list)[i]
    tmp<-own_auc_list[[i]]
    tmp<-lapply(tmp, function(x){x[[1]]$pred_df})
    for(n in 1:length(dataset)){
      dataset_name=dataset[n]
      df_index[[dataset_name]][[feature_name]]<-tmp[[dataset_name]]
    }
  }
  ###
  if(!is.null(published_auc_list)){
    published_auc_list<-lapply(published_auc_list, function(x){x[[model_name]]})
    for(i in 1:length(published_auc_list)){
      feature_name=names(published_auc_list)[i]
      tmp<-published_auc_list[[i]]
      tmp<-lapply(tmp, function(x){x[[1]]$pred_df})
      for(n in 1:length(dataset)){
        dataset_name=dataset[n]
        df_index[[dataset_name]][[feature_name]]<-tmp[[dataset_name]]
      }
    }
  }

  ###
  result_list <-data.frame()
  for(n in 1:length(dataset)){
    dataset_name=dataset[n]
    tmp<-df_index[[dataset_name]]
    df<-data.frame(row.names = row.names(tmp[[1]]),
                   time=tmp[[1]]$time,
                   status=tmp[[1]]$status)
    for(i in 1:length(tmp)){
      feature_name=names(tmp)[i]
      d<-tmp[[i]]
      d<-d[row.names(df),]
      df[[feature_name]]=d$pred
    }

    # saving result
    features<-colnames(df)[3:ncol(df)]
    feature_combinations <- combn(features, 2, simplify = FALSE)

    for (comb in feature_combinations) {
        x <- df[[comb[1]]]
        y <- df[[comb[2]]]

        # compareC to calculate difference
        cindex_result <- compareC(df$time,df$status,x, y)

        # saving results
        sig<-data.frame(comparison=paste(comb, collapse = " vs "),pvalue=cindex_result$pval,dataset=dataset_name)
        result_list<-rbind(result_list,sig)
      }
  }

  return(result_list)
}
