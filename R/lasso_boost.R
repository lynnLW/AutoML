library(glmnet)
library(Hmisc)
library(survcomp)
library(survivalROC)
library(survival)
lasso_cox <- function(x.train,fold){

  set.seed(12345621)
  fit = cv.glmnet(as.matrix(x.train[,3:ncol(x.train)]),as.matrix(x.train[,1:2]) ,
                  family = "cox",alpha=1,nfolds = fold)

  pred_cox = predict(fit,type='link',
                     newx=as.matrix(x.test[,3:ncol(x.test)]),
                     s = fit$lambda.min)
  pred_cox <- as.numeric(pred_cox)
  cindex_cox = 1-rcorr.cens(pred_cox,Surv(x.test$time, x.test$status))[[1]]
  print(sprintf("cox_%d",fold))
  print(cindex_cox)


  lasso_model <- list()
  lasso_model$best_param <- list( s = fit$lambda.min)
  lasso_model$model = fit
  lasso_model$cindex = 1-rcorr.cens(pred_cox,Surv(t.test, s.test))[[1]]
  lasso_model$pred = pred_cox

  ### auc km ###

  cutoff = 12*1
  if ( min(x.test$time) < cutoff )
  {
    y <- survivalROC(Stime = x.test$time, status =x.test$status, marker = pred_cox,
                     predict.time = cutoff, method = "KM")
    lasso_model$km_fp_1 = y$FP
    lasso_model$km_tp_1 = y$TP
    lasso_model$km_auc_1 = y$AUC
  }else
  {
    lasso_model$km_fp_1 = NA
    lasso_model$km_tp_1 = NA
    lasso_model$km_auc_1 = NA
  }

  cutoff=12*3
  y <- survivalROC(Stime = x.test$time, status =x.test$status, marker = pred_cox,
                   predict.time = cutoff,method = "KM")
  lasso_model$km_fp_3 = y$FP
  lasso_model$km_tp_3 = y$TP
  lasso_model$km_auc_3 = y$AUC

  cutoff=12*5
  y <- survivalROC(Stime = x.test$time, status =x.test$status, marker = pred_cox,
                   predict.time = cutoff,method = "KM")
  lasso_model$km_fp_5 = y$FP
  lasso_model$km_tp_5 = y$TP
  lasso_model$km_auc_5 = y$AUC

  dd_ext <- data.frame("time"=x.test$time, "event"=x.test$status, "score"= pred_cox)
  Brier_score <- sbrier.score2proba(data.tr=dd_ext, data.ts=dd_ext, method="cox")
  lasso_model$brier_Score <- Brier_score


  save("lasso_model", file = sprintf("%d_lasso_cox_result.RData", fold))

}
