ML.coxboost<-function(est_dd,seed=123456){
  ##### 7.CoxBoost ###########
  message("--- 7.CoxBoost  ---")

  set.seed(seed)
  pen <- optimCoxBoostPenalty(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]),
                              trace = TRUE, start.penalty = 500, parallel = T
  )

  cv.res <- cv.CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]),
                        maxstepno = 500, K = 10, type = "verweij", penalty = pen$penalty
  )
  fit <- CoxBoost(est_dd[, "OS.time"], est_dd[, "OS"], as.matrix(est_dd[, -c(1, 2)]),
                  stepno = cv.res$optimal.step, penalty = pen$penalty
  )
  rid <- as.data.frame(coef(fit))
  rid$id <- rownames(rid)
  rid <- rid[which(rid$`coef(fit)` != 0), "id"]
  result <- data.frame(
    method = c(rep("CoxBoost", length(rid))),
    selected.fea = rid
  )

  selected.feature <- rbind(selected.feature, result)
  write.table(result,file="7.Coxboost_select_features.csv",sep=",",row.names = F)
  return(selected.feature)
}


