library(CoxBoost)
library(survival)
library(survivalROC)
library(Hmisc)
library(survcomp)

cox_boost <- function(x.train,fold){
  set.seed(12345621)
  pen <- optimCoxBoostPenalty(x.train[,'time'],x.train[,'status'],as.matrix(x.train[,-c(1,2)]),
                              trace=TRUE,start.penalty=500,parallel = T)
  #K>>number of folds to be used for cross-validation
  cv.res <- cv.CoxBoost(x.train[,'time'],x.train[,'status'],as.matrix(x.train[,-c(1,2)]),
                        maxstepno=500,K=fold,type="verweij",penalty=pen$penalty)

  fit <- CoxBoost(x.train[,'time'],x.train[,'status'],as.matrix(x.train[,-c(1,2)]),
                  stepno=cv.res$optimal.step,penalty=pen$penalty)

  pred_coxb <- predict(fit,newdata=x.test[,-c(1,2)],
                       newtime=x.test[,'time'],
                       newstatus=x.test[,'status'],
                       type="lp")
  pred_coxb <- as.numeric(pred_coxb)
  cindex_coxb = 1-rcorr.cens(pred_coxb,Surv(x.test$time, x.test$status))[[1]]
  print(sprintf("cox_%d",fold))
  print(cindex_coxb)

  coxboost_model <- list()
  coxboost_model$best_param =list(stepno=cv.res$optimal.step,penalty=pen$penalty)
  coxboost_model$model = fit
  coxboost_model$cindex = cindex_coxb
  coxboost_model$pred = pred_coxb

  ### auc km ###
  cutoff = 12*1
  if ( min(x.test$time) < cutoff )
  {
    y <- survivalROC(Stime = x.test$time, status = x.test$status, marker = pred_coxb,
                     predict.time = cutoff, method = "KM")
    coxboost_model$km_fp_1 = y$FP
    coxboost_model$km_tp_1 = y$TP
    coxboost_model$km_auc_1 = y$AUC
  }else {
    coxboost_model$km_fp_1 = NA
    coxboost_model$km_tp_1 = NA
    coxboost_model$km_auc_1 = NA
  }

  cutoff=12*3
  y <- survivalROC(Stime = x.test$time, status = x.test$status, marker = pred_coxb,
                   predict.time = cutoff,method = "KM")
  coxboost_model$km_fp_3 = y$FP
  coxboost_model$km_tp_3 = y$TP
  coxboost_model$km_auc_3 = y$AUC

  cutoff=12*5
  y <- survivalROC(Stime = x.test$time, status = x.test$status, marker = pred_coxb,
                   predict.time = cutoff,method = "KM")
  coxboost_model$km_fp_5 = y$FP
  coxboost_model$km_tp_5 = y$TP
  coxboost_model$km_auc_5 = y$AUC


  dd_ext <- data.frame("time"=x.test$time, "event"=x.test$status, "score"= pred_coxb)
  Brier_score <- sbrier.score2proba(data.tr=dd_ext, data.ts=dd_ext, method="cox")
  coxboost_model$brier_Score <- Brier_score
  save("coxboost_model", file = sprintf("%d_coxboost_result.RData", fold))

}
