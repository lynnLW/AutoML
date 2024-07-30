list_dataset_generate<-function(filepath,metapath,dname=NULL,output){
  ##loading packages
  library(dplyr)
  library(stringr)
  ##list all files that need to be combined
  files<-list.files(filepath,full.names = T)
  files
  print(paste0(length(files)," datasets need to be combined"))
  ##combine
  list_data<-c()
  for(i in 1:length(files)){
    df<-read.table(files[i],sep=",",header=T,row.names = 1)
    list_data[[i]]<-df
  }
  ##names each dataset
  if(is.null(dname)==T){
    filename=str_split(basename(files),"_",simplify = T)[,1]
  } else {
    filename=dname
  }
  names(list_data)<-filename
  ###
  list_data <- lapply(list_data, function(df) {
    df[rowSums(df != 0) > 0, ]
  })
  ###log2 transformation for RNA-seq data
  for(i in names(list_data)){
    #boxplot(list_data[[i]][,1:20])
    if (max(list_data[[i]][,1])>20){
    list_data[[i]]<-log2(list_data[[i]]+1)
    }
    #boxplot(list_data[[i]][,1:20])
  }
  ###view data
  boxplot(c(list_data[[1]][,1:20]))
  ###
  if(!dir.exists(paste0(output,"/nocombat"))){
    dir.create(paste0(output,"/nocombat"))
  }
  save(list_data,file=paste0(output,"/nocombat/list_data.Rdata"))
  ###scale all data
  list_data_scale<-lapply(list_data, function(df) {
    as.data.frame(scale(df,center = TRUE, scale =T))
  })
  ###
  boxplot(list_data_scale[[1]][,1:20])
  ###
  save(list_data_scale,file=paste0(output,"/nocombat/list_data_scale.Rdata"))
  ###BCR meta info
  BCR_files<-list.files(metapath,full.names = T)
  BCR_files
  BCR_data<-c()
  for(i in 1:length(filename)){
    file=BCR_files[grep(filename[i],BCR_files)]
    df<-read.table(file,sep=",",header=T,row.names = 1)
    BCR_data[[i]]<-df
  }
  ##transform
  list_data_t<-lapply(list_data,function(df){
    df<-df %>% t() %>% as.data.frame()
  })
  ##combine with metainfo
  list_data_surv<-c()
  for(i in 1:length(filename)){
    meta<-BCR_data[[i]][row.names(list_data_t[[i]]),]
    list_data_surv[[i]]<-data.frame(ID=row.names(list_data_t[[i]]),meta,list_data_t[[i]])
  }
  ##transform
  list_data_scale<-lapply(list_data_scale,function(df){
    df<-df %>% t() %>% as.data.frame()
  })
  ##combine with metainfo
  list_data_scale_surv<-c()
  for(i in 1:length(filename)){
    meta<-BCR_data[[i]][row.names(list_data_scale[[i]]),]
    list_data_scale_surv[[i]]<-data.frame(ID=row.names(list_data_scale[[i]]),meta,list_data_scale[[i]])
  }
  names(list_data_surv)<-filename
  names(list_data_scale_surv)<-filename
  save(list_data_surv,file=paste0(output,"/nocombat/list_data_surv.Rdata"))
  save(list_data_scale_surv,file=paste0(output,"/nocombat/list_data_scale_surv.Rdata"))
  ##
  common.feature<-row.names(list_data[[1]])
  for (i in names(list_data)){
    common.feature<-intersect(common.feature,row.names(list_data[[i]]))
  }
  ##
  list_data_filter<-lapply(list_data,function(x){
    x<-x[common.feature,]
    return(x)
  })
  if(!dir.exists(paste0(output,"/intersected_genes"))){
    dir.create(paste0(output,"/intersected_genes"))
  }
  write.table(common.feature,file=paste0(output,"/intersected_genes/common_feature.csv"),sep=",",row.names = F)
  ###scale all data
  list_data_filter_scale<-lapply(list_data_filter, function(df) {
    as.data.frame(scale(df))
  })
  ###
  boxplot(c(list_data_filter_scale[[1]][,1:20],list_data_filter_scale[[2]][,1:20],
            list_data_filter_scale[[3]][,1:20],list_data_filter_scale[[4]][,1:20],
            list_data_filter_scale[[5]][,1:20],list_data_filter_scale[[6]][,1:20],
            list_data_filter_scale[[7]][,1:20],list_data_filter_scale[[8]][,1:20],
            list_data_filter_scale[[9]][,1:20],list_data_filter_scale[[10]][,1:20]))
  ##transform
  list_data_filter<-lapply(list_data_filter,function(df){
    df<-df %>% t() %>% as.data.frame()
  })
  ##combine with metainfo
  list_data_filter_surv<-c()
  for(i in 1:length(filename)){
    meta<-BCR_data[[i]][row.names(list_data_filter[[i]]),]
    list_data_filter_surv[[i]]<-data.frame(ID=row.names(list_data_filter[[i]]),meta,list_data_filter[[i]])
  }
  ##
  list_data_filter_scale<-lapply(list_data_filter_scale,function(df){
    df<-df %>% t() %>% as.data.frame()
  })
  ##combine with metainfo
  list_data_filter_scale_surv<-c()
  for(i in 1:length(filename)){
    meta<-BCR_data[[i]][row.names(list_data_filter_scale[[i]]),]
    list_data_filter_scale_surv[[i]]<-data.frame(ID=row.names(list_data_filter_scale[[i]]),meta,list_data_filter_scale[[i]])
  }
  names(list_data_filter_surv)<-filename
  names(list_data_filter_scale_surv)<-filename
  save(list_data_filter_surv,file=paste0(output,"/nocombat/list_data_filter_surv.Rdata"))
  save(list_data_filter_scale_surv,file=paste0(output,"/nocombat/list_data_filter_scale_surv.Rdata"))
  ##
  if(!dir.exists(paste0(output,"/combat/"))){
    dir.create(paste0(output,"/combat/"))
  }
  ###
  meta<-c()
  bcr_files<-list.files("../BCR_meta/",full.names = T)
  for (i in 1:length(list)){
    bcr=read.table(bcr_files[i],sep=",",header = T)
    sample=colnames(list[[i]])
    batch=rep(paste0("batch",i),length(sample))
    d<-data.frame(batch=batch,dataset=names[i])
    d<-cbind(bcr,d)
    meta<-rbind(meta,d)
  }
  write.table(meta,file="combat/all_dataset_type.csv",sep=",")
  ##32,12,11分别是提取出的每个样本的数量，样本1命名为batch1，样本2命名为batch2，样本3命名为batch3

}
