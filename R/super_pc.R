library(superpc)
library(survival)
library(survivalROC)
library(Hmisc)
library(survcomp)

superpc <- function(x.train,fold){
  set.seed(12345621)
  data <- list(x=t(x.train[,-c(1,2)]),
               y=x.train$time,
               censoring.status=x.train$status,
               featurenames=colnames(x.train)[-c(1,2)])

  fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default

  cv.fit <- superpc.cv(fit,data,
                       n.threshold = 8,#default ,Number of thresholds to consider
                       n.fold = fold, #Number of cross-validation folds
                       n.components=3,
                       min.features=5,#default
                       max.features=nrow(data$x),
                       compute.fullcv= TRUE,
                       compute.preval=TRUE)

  test <- list(x=t(x.test[,-c(1,2)]),
               y=x.test$time,
               censoring.status=x.test$status,
               featurenames=colnames(x.test)[-c(1,2)])

  pred_superpc <- superpc.predict(fit,data,
                                  test,
                                  threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],
                                  n.components = 1)$v.pred

  pred_superpc <- as.numeric(pred_superpc)
  cindex_superpc = 1-rcorr.cens(pred_superpc ,Surv(x.test$time, x.test$status))[[1]]
  print(sprintf("superpc_%d",fold))
  print(cindex_superpc)

  superpc_model <- list()
  superpc_model$best_param =list(threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])])
  superpc_model$model = cv.fit
  superpc_model$cindex = cindex_superpc
  superpc_model$pred = pred_superpc

  ### auc km ###

  cutoff = 12*1
  if ( min(x.test$time) < cutoff )
  {
    y <- survivalROC(Stime = x.test$time, status = x.test$status, marker = pred_superpc,
                     predict.time = cutoff, method = "KM")
    superpc_model$km_fp_1 = y$FP
    superpc_model$km_tp_1 = y$TP
    superpc_model$km_auc_1 = y$AUC
  }else
  {
    superpc_model$km_fp_1 = NA
    superpc_model$km_tp_1 = NA
    superpc_model$km_auc_1 = NA
  }

  cutoff=12*3
  y <- survivalROC(Stime = x.test$time, status = x.test$status, marker = pred_superpc,
                   predict.time = cutoff,method = "KM")
  superpc_model$km_fp_3 = y$FP
  superpc_model$km_tp_3 = y$TP
  superpc_model$km_auc_3 = y$AUC

  cutoff=12*5
  y <- survivalROC(Stime = x.test$time, status = x.test$status, marker = pred_superpc,
                   predict.time = cutoff,method = "KM")
  superpc_model$km_fp_5 = y$FP
  superpc_model$km_tp_5 = y$TP
  superpc_model$km_auc_5 = y$AUC



  dd_ext <- data.frame("time"=x.test$time, "event"=x.test$status, "score"= pred_superpc)
  Brier_score <- sbrier.score2proba(data.tr=dd_ext, data.ts=dd_ext, method="cox")
  superpc_model$brier_Score <- Brier_score

  save("superpc_model", file = sprintf("%d_superpc_result.RData", fold))
}
