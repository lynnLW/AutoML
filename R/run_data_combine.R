#' Data preprocessing for model performance evaluation
#'
#' @param list_data_surv #column 1 ID # 2 survival time # 3 status 1 or 0 # 4 gene variables
#' @param sets  c(1,2,3,4)
#' @param correct_method sva,scale,min_max,none
#' @param outdir output path
#' @return data list
#' @export
#'
run_data_combine<-function(list_data_surv,
                  sets=NULL,
                  correct_method='none',
                  outdir="dataset/"){
  ##loading packages
  library(tinyarray)
  library(FactoMineR)
  library(factoextra)
  library(sva)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(cowplot)
  ##output directory
  if(!dir.exists(outdir)){
      dir.create(outdir,recursive = T)
  }
  ##sets to combat
  if(is.null(sets)==T){
    sets=seq_along(list_data_surv)
  } else {
    sets=sets
  }

  ##list all files that need to be combined
  print(paste0(length(list_data_surv)," datasets"))
  ##format each datasets
  list_data_surv<-lapply(list_data_surv,function(x){
    names(x)[1]<-"ID"
    names(x)[2]<-"OS_time"
    names(x)[3]<-"OS_status"
    return(x)
  })

  ##sets
  list_vali_data<-c()
  list_vali_data<-list_data_surv[sets]
  print(paste0(length(list_vali_data)," datasets is preparing to be scaled and combined"))

  ##dataset meta info
  meta_info<-data.frame()
  for(i in 1:length(list_vali_data)){
    df<-list_vali_data[[i]]
    dataset_name<-names(list_vali_data)[i]
    data<-data.frame(df[,1:3],batch=dataset_name)
    meta_info<-rbind(meta_info,data)
  }

  ##common feature
  common.feature<-colnames(list_vali_data[[1]])[-c(1:3)]
  for (i in names(list_vali_data)){
    common.feature<-intersect(common.feature,colnames(list_vali_data[[i]]))
    print(length(common.feature))
  }
  print(paste0(length(common.feature)," common features across all datasets"))

  ##filter data.frame
  list_vali_data<-lapply(list_vali_data,function(x){
    x<-x[,common.feature]
    return(x)
  })

  ##combine all datasets to pca
  expr<-data.frame()
  for (i in names(list_vali_data)){
    expr<-rbind(expr,list_vali_data[[i]])
  }

  ##format batch info
  expr<-expr[meta_info$ID,]
  ##pca before combat
  pre.pca <- PCA(expr,graph = FALSE)
  fviz_pca_ind(pre.pca,
               geom= "point",
               col.ind = meta_info$batch,
               addEllipses = TRUE,
               legend.title="Cohort",
               ellipse.level=0.95)
  ##saving result
  ggsave(filename =paste0(outdir,"/pca/pre_pca.jpg"),width = 4,height = 4,dpi=600)

  if(correct_method=="sva"){
    print("Performing the combat process")
    mod <- model.matrix(~as.factor(OS_status),data=meta_info)
    combat_expr <- ComBat(dat = t(expr),batch = meta_info$batch,mod=mod)
    combat.pca <- PCA(t(combat_expr),graph = FALSE)
    fviz_pca_ind(combat.pca,
                 geom= "point",
                 col.ind = meta_info$batch,
                 addEllipses = TRUE,
                 ellipse.level=0.95,
                 legend.title="Batch")
    ggsave(filename =paste0(outdir,"/pca/post_pca_sva.jpg"),
           width = 4,height = 4,dpi=600)

    ##sperate the dataset
    datasets<-unique(meta_info$batch)
    for (i in 1:length(datasets)){
      name=datasets[i]
      sample=meta_info[(meta_info$batch %in% name),]$ID
      data<-as.data.frame(t(combat_expr[,sample]))
      list_vali_data[[name]]<-data
    }

  } else if (correct_method=="scale"){
    print("Performing the scale process")
    list_vali_data<-lapply(list_vali_data, function(df) {
      scaled_df<-scale(df,center = TRUE, scale =T) %>% as.data.frame()
      return(scaled_df)
    })

    ##combine all datasets to pca
    combat_expr<-data.frame()
    for (i in names(list_vali_data)){
      combat_expr<-rbind(combat_expr,list_vali_data[[i]])
    }

    ##format batch info
    combat_expr<-combat_expr[meta_info$ID,]
    ##pca before combat
    post.pca <- PCA(combat_expr,graph = FALSE)
    fviz_pca_ind(post.pca,
                 geom= "point",
                 col.ind = meta_info$batch,
                 addEllipses = TRUE,
                 legend.title="Batch",
                 ellipse.level=0.95)
    ##saving result
    ggsave(filename =paste0(outdir,"/pca/post_pca_scale.jpg"),width = 4,height = 4,dpi=600)
  } else if (correct_method=="min_max"){
    print("Performing the min-max scale process")
    min_max_scale<- function(x) {
      return((x - min(x)) / (max(x) - min(x)) * 2 - 1)
    }
    # scale features to 1
    list_vali_data<-lapply(list_vali_data, function(df) {
      scaled_df <- as.data.frame(lapply(df, min_max_scale))
      row.names(scaled_df)<-row.names(df)
      return(scaled_df)
    })

    ##combine all datasets to pca
    combat_expr<-data.frame()
    for (i in names(list_vali_data)){
      combat_expr<-rbind(combat_expr,list_vali_data[[i]])
    }

    ##format batch info
    combat_expr<-combat_expr[meta_info$ID,]
    ##pca before combat
    post.pca <- PCA(combat_expr,graph = FALSE)
    fviz_pca_ind(post.pca,
                 geom= "point",
                 col.ind = meta_info$batch,
                 addEllipses = TRUE,
                 legend.title="Batch",
                 ellipse.level=0.95)
    ##saving result
    ggsave(filename =paste0(outdir,"/pca/post_pca_min_max.jpg"),width = 4,height = 4,dpi=600)

  } else if (correct_method=="none") {
    print("Not performing the correct process")
    list_vali_data<-list_vali_data
  }

  datasets<-names(list_vali_data)
  ##combine the metainfo
  for (i in 1:length(list_vali_data)){
    name=datasets[i]
    sample=meta_info[(meta_info$batch %in% name),]$ID
    data<-list_vali_data[[i]]
    meta<-meta_info[(meta_info$batch %in% name),]
    meta<-meta[row.names(data),]
    data<-cbind(meta[,1:3],data)
    list_vali_data[[name]]<-data
  }
  save(list_vali_data,file=paste0(outdir,"/list_vali_data.Rdata"))

  print("Finish the combine process")
  return(list_vali_data)

}
