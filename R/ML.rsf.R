ML.rsf<-function(est_dd,seed=123456){
  ##### 6.RSF ###########
  message("--- 6.RSF  ---")
  set.seed(seed)
  ns_res<-tune.nodesize(Surv(OS.time, OS) ~ ., est_dd)
  ###
  nodesize=ns_res$nsize.opt
  print(paste0("The optimal nodesize is ",nodesize))
  set.seed(seed)
  fit <- rfsrc(Surv(OS.time, OS) ~ .,
               data = est_dd,
               ntree = 1000, nodesize = nodesize, # 该值建议多调整
               splitrule = "logrank",
               importance = T,
               proximity = T,
               forest = T,
               seed = seed
  )
  rid <- var.select(object = fit, conservative = "high")
  rid <- rid$topvars

  result <- data.frame(
    method = c(rep("RSF", length(rid))),
    selected.fea = rid
  )

  selected.feature <- rbind(selected.feature, result)
  write.table(result,file="6.RSF_select_features.csv",sep=",",row.names = F)
  return(selected.feature)
}

