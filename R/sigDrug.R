sigDrug<-function(merge_expr,genelist,outdir){

  ###diff expression
  df_mean<- merge_expr %>%
    group_by(group)%>%
    summarise(across(all_of(genelist), \(x) mean(x, na.rm = TRUE )))%>%
    t() %>% as.data.frame()
  names(df_mean)<-df_mean[1,]
  df_mean<-df_mean[-1,]
  df_mean<-na.omit(df_mean)
  ###function calculate pvalue
  cal_pvalue<-function(input,gene){
    pvalue<-pairwise.wilcox.test(input[,gene],
                                 input[,'group'],
                                 p.adjust.method = "bonf")$p.value
    return(pvalue)
  }
  ###cal pvalue
  pvalue<-c()
  for(i in 1:length(row.names(df_mean))){
    gene=row.names(df_mean)[i]
    p<-cal_pvalue(merge_expr,gene)
    pvalue<-append(pvalue,p)
  }
  ###
  df_mean$pvalue<-pvalue
  df_mean[,1] <- as.numeric(df_mean[,1] )
  df_mean[,2]  <- as.numeric(df_mean[,2] )
  df_mean$fc <- df_mean[,1] / df_mean[,2]
  df_mean<-df_mean[df_mean$pvalue<0.05,]

  write.csv(df_mean,file=paste0(outdir,"/1.diff_drug_0.05.csv"))
  return(df_mean)

}
