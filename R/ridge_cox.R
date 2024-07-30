library(glmnet)
library(Hmisc)
library(survival)
library(survcomp)
library(survivalROC)
ridge_cox <- function(x.train,fold){

  time = as.double(x.train$time)
  status = as.double(x.train$status)

  set.seed(12345621)
  fit = cv.glmnet(as.matrix(x.train[,3:ncol(x.train)]),as.matrix(x.train[,1:2]),
                  family = "cox",alpha=0,nfolds = fold)

  pred_cox = predict(fit,type='link',
                     newx=as.matrix(x.test[,3:ncol(x.test)]),
                     s = fit$lambda.min)

  pred_cox <- as.numeric(pred_cox)

  cindex_cox = 1-rcorr.cens(pred_cox,Surv(x.test$time, x.test$status))[[1]]
  print(sprintf("cox_%d",fold))
  print(cindex_cox)


  ridge_model <- list()
  ridge_model$model = fit
  ridge_model$best_param <- list(s = fit$lambda.min)
  ridge_model$cindex = 1-rcorr.cens(pred_cox,Surv(t.test, s.test))[[1]]
  ridge_model$pred = pred_cox
  ### auc km ###

  cutoff = 12*1
  if ( min(x.test$time) < cutoff ) {
    y <- survivalROC(Stime = x.test$time, status = x.test$status, marker = pred_cox,
                     predict.time = cutoff, method = "KM")
    ridge_model$km_fp_1 = y$FP
    ridge_model$km_tp_1 = y$TP
    ridge_model$km_auc_1 = y$AUC
  } else {
    ridge_model$km_fp_1 = NA
    ridge_model$km_tp_1 = NA
    ridge_model$km_auc_1 = NA
  }

  cutoff=12*3
  y <- survivalROC(Stime = x.test$time, status = x.test$status, marker = pred_cox,
                   predict.time = cutoff,method = "KM")
  ridge_model$km_fp_3 = y$FP
  ridge_model$km_tp_3 = y$TP
  ridge_model$km_auc_3 = y$AUC

  cutoff=12*5
  y <- survivalROC(Stime = x.test$time, status = x.test$status, marker = pred_cox,
                   predict.time = cutoff,method = "KM")
  ridge_model$km_fp_5 = y$FP
  ridge_model$km_tp_5 = y$TP
  ridge_model$km_auc_5 = y$AUC


  dd_ext <- data.frame("time"=x.test$time, "event"=x.test$status, "score"= pred_cox)
  Brier_score <- sbrier.score2proba(data.tr=dd_ext, data.ts=dd_ext, method="cox")
  ridge_model$brier_Score <- Brier_score

  save("ridge_model", file = sprintf("%d_ridge_cox_result.RData", fold))

}
