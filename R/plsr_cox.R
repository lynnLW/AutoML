
library(plsRcox)
library(survival)
library(survivalROC)
library(Hmisc)
library(survcomp)

plsR_cox <- function(x.train,fold){
  set.seed(12345621)
  cv.plsRcox.res=cv.plsRcox(list(x=x.train[,3:ncol(x.train)],
                                 time=x.train$time,
                                 status=x.train$status),
                            nfold = fold,
                            nt=10,
                            verbose = FALSE)


  fit <- plsRcox(x.train[,3:ncol(x.train)],
                 time=x.train$time,
                 event=x.train$status,
                 nt=as.numeric(cv.plsRcox.res[5]))

  pred_pls <- predict(fit,type="lp",newdata=x.test[,-c(1,2)])
  pred_pls <- as.numeric(pred_pls)
  cindex_pls = 1-rcorr.cens(pred_pls,Surv(x.test$time, x.test$status))[[1]]
  print(sprintf("plscox_%d",fold))
  print(cindex_pls)

  plsRcox_model <- list()
  plsRcox_model$best_param =list(nt=as.numeric(cv.plsRcox.res[5]))
  plsRcox_model$model = fit
  plsRcox_model$cindex = cindex_pls
  plsRcox_model$pred = pred_pls


  ### auc km ###

  cutoff = 12*1
  if ( min(x.test$time) < cutoff ) {
    y <- survivalROC(Stime = x.test$time, status = x.test$status, marker = pred_pls,
                     predict.time = cutoff, method = "KM")
    plsRcox_model$km_fp_1 = y$FP
    plsRcox_model$km_tp_1 = y$TP
    plsRcox_model$km_auc_1 = y$AUC
  }else {
    plsRcox_model$km_fp_1 = NA
    plsRcox_model$km_tp_1 = NA
    plsRcox_model$km_auc_1 = NA
  }

  cutoff=12*3
  y <- survivalROC(Stime = x.test$time, status = x.test$status, marker = pred_pls,
                   predict.time = cutoff,method = "KM")
  plsRcox_model$km_fp_3 = y$FP
  plsRcox_model$km_tp_3 = y$TP
  plsRcox_model$km_auc_3 = y$AUC

  cutoff=12*5
  y <- survivalROC(Stime = x.test$time, status = x.test$status, marker = pred_pls,
                   predict.time = cutoff,method = "KM")
  plsRcox_model$km_fp_5 = y$FP
  plsRcox_model$km_tp_5 = y$TP
  plsRcox_model$km_auc_5 = y$AUC


  dd_ext <- data.frame("time"=x.test$time, "event"=x.test$status, "score"= pred_pls)
  Brier_score <- survcomp::sbrier.score2proba(data.tr=dd_ext, data.ts=dd_ext, method="cox")
  plsRcox_model$brier_Score <- Brier_score

  save("plsRcox_model", file = sprintf("%d_plsRcox_result.RData", fold))


}
