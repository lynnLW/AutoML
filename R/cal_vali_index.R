#' Calculate AUC scores of Machine Learning Models in all data
#'
#' @param list_train_vali_Data A list contain the dataframe (colnames:ID,time,status,feature)
#' @param candidate_genes  features that used to build the model
#' @param model_list A list contain trained model
#' @param filter_OS_time if keep survival time >30 days
#' @param meta_time follow up time used by days or months or years c("m","y","d") default m
#' @param rep repeats
#' @param outdir the output directory
#' @return test.auc list containing:
#'         - model:
#'           - test data:
#'            - rep 1:
#'              - pred: risk score
#'              - cindex: c-index
#'              - bs : integrated brier score
#'              - auc : 1,2,3,5,7,10- year auc values
#' @export
cal_vali_index <- function(list_train_vali_Data, # A list contain the dataframe for testing
                    candidate_genes, # features used to train model
                    model_list, #model list
                    filter_OS_time=F,
                    meta_time="m",
                    rep=1,
                    outdir="4.test/"
) {

  if(!dir.exists(outdir)){
      dir.create(outdir,recursive = T)
  }

  ####  Data preprocessing #####
  print("Checking data feasibility")
  list_train_vali_Data <- lapply(list_train_vali_Data,function(x){
    colnames(x) = gsub('-','.',colnames(x))
    colnames(x)[2]<-"time"
    colnames(x)[3]<-"status"
    return(x)})

  # Matching candidate genes to genes in each cohort
  print("Matching candidate genes to genes in each cohort")
  common_feature = c('time','status',candidate_genes)
  common_feature = intersect(common_feature,colnames(list_train_vali_Data[[1]]))
  for (i in names(list_train_vali_Data)) {
    common_feature = intersect(common_feature, colnames(list_train_vali_Data[[i]]))
  }

  message(paste0('---the number of the raw candidate genes is ', length(candidate_genes),' ---'))
  message(paste0('---the number of the common feature across all dataset is ', length(common_feature)-2,' ---'))

  # Matching common feature in each cohort
  for (i in names(list_train_vali_Data)) {
    list_train_vali_Data[[i]] =list_train_vali_Data[[i]][,common_feature]
  }

  ## If keep follow up days more than 30 days
  if (filter_OS_time){
    if (meta_time=="m") {
      list_train_vali_Data<-lapply(list_train_vali_Data,function(x){
        x=x[x$time>1,]
        return(x)
      })
    } else if (meta_time=="y"){
      list_train_vali_Data<-lapply(list_train_vali_Data,function(x){
        x=x[x$time>0.083,]
        return(x)
      })
    } else {
      list_train_vali_Data<-lapply(list_train_vali_Data,function(x){
        x=x[x$time>30,]
        return(x)
      })
    }
  } else {
    list_train_vali_Data<-list_train_vali_Data
  }

  ## calculate cindex
  model_auc_list<-list()
  for (i in 1:length(model_list)){
    models=model_list[[i]]
    model_name=names(model_list)[i]
    auc_list<-list()
    print(model_name)
    for (n in 1:length(list_train_vali_Data)){
      test_data<-list_train_vali_Data[[n]]
      data_name<-names(list_train_vali_Data)[n]
      print(data_name)
      result<-c()
      if (rep==1){
        model=models$model
        result[[1]]<-cal_metrics(test_data,model,model_name)
      } else {
        for ( j in 1:rep){
          print(j)
          model<-models[[j]]
          result[[j]]<-cal_metrics(test_data,model,model_name)
        }
      }
      auc_list[[n]]<-result
      names(auc_list)[n]<-data_name
    }
    model_auc_list[[i]]<-auc_list
    names(model_auc_list)[i]<-model_name
  }
  save(model_auc_list,file=paste0(outdir,"/test_index.Rdata"))
  return(model_auc_list)
}
