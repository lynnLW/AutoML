boruta_feature_select<-function(input,seed){
  set.seed(seed)
  data<-input[,-1]
  data$OS<-as.factor(data$OS)
  boruta<-Boruta(OS ~.,data = data,doTrace=2,maxRuns = 1000,ntree=1000)
  boruta
  ##
  print("plot the importance of feature")
  boruta.variable.imp <- boruta.imp(boruta)
  head(boruta.variable.imp)
  library("YSX")
  p<-sp_boxplot(boruta.variable.imp, melted=T, xvariable = "Variable", yvariable = "Importance",
                legend_variable = "finalDecision",x_label = "",coordinate_flip=F,
                title="Feature importance",
                legend_variable_order = c("shadowMax", "shadowMean", "shadowMin", "Tentative","Confirmed"),
                xtics_angle = 90)

  print(p)
  ggsave(p,filename = "feature_select/boruta_feature_importance.jpg",dpi=600,units="cm",width=10,height =8,scale = 1.5)
  write.table(boruta$finalDecision,file="feature_select/boruta_finalDecision_feature.csv",sep=",")
  result<-as.data.frame(boruta$finalDecision)
  result<-data.frame(gene=row.names(result),decision=result[,1])
  selected.feature<-result[result$decision=="Confirmed",]$gene
  return(selected.feature)
}

boruta_surv_feature_select<-function(input,seed){
  set.seed(seed)
  data<-input
  Y=Surv(data$OS.time,data$OS)
  boruta<-Boruta(y=Y,x = data[,-c(1:2)],maxRuns = 1000,ntree=1000)
  ##
  print("plot the importance of feature")
  boruta.variable.imp <- boruta.imp(boruta)
  head(boruta.variable.imp)
  library("YSX")
  p<-sp_boxplot(boruta.variable.imp, melted=T, xvariable = "Variable", yvariable = "Importance",
                legend_variable = "finalDecision",x_label = "",coordinate_flip=F,
                title="Feature importance",
                legend_variable_order = c("shadowMax", "shadowMean", "shadowMin", "Tentative","Confirmed"),
                xtics_angle = 90)
  print(p)
  ggsave(p,filename = "feature_select/boruta_surv_feature_importance.jpg",dpi=600,units="cm",width=10,height =8,scale = 1.5)
  write.table(boruta$finalDecision,file="feature_select/boruta_surv_finalDecision_feature.csv",sep=",")
  result<-as.data.frame(boruta$finalDecision)
  result<-data.frame(gene=row.names(result),decision=result[,1])
  selected.feature<-result[result$decision=="Confirmed",]$gene
  return(selected.feature)
}
