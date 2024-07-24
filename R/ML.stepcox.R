ML.stepCox<-function(est_dd,seed=123456){
  ##### 8.StepCox ###########
  message("--- 8.StepCox ---")
  all_result<-c()
  for (direction in c("both", "backward", "forward")) {
    fit <- stats::step(coxph(Surv(OS.time, OS) ~ ., est_dd), direction = direction)
    rid <- names(coef(fit)) # 这里不用卡P值，迭代的结果就是可以纳入的基因

    result <- data.frame(
      method = c(rep(paste0("StepCox", "+", direction), length(rid))),
      selected.fea = rid
    )
    all_result <- rbind(all_result, result)

  }
  selected.feature <- rbind(selected.feature, all_result)
  write.table(result,file="8.stepcox_select_features.csv",sep=",",row.names = F)
  return(selected.feature)
}
