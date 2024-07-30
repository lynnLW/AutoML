library(gbm)
library(survival)
library(caret)
library(Hmisc)
library(survivalROC)
library(survcomp)

gbm_model<-function(x.train,fold){
  set.seed(12345621)
  # 定义参数网格
  param_grid <- expand.grid(interaction.depth = c(1,3,5,7),
                            n.trees = c(100,500,1000,2000),
                            shrinkage = c(0.005,0.001,0.005,0.01,0.1),n.minobsinnode = c(5,10,15,20))

  # 定义保存结果的向量
  results <- data.frame(interaction.depth = integer(),
                        n.trees = integer(),
                        best_iter=integer(),
                        shrinkage = double(),
                        n.minobsinnode = integer(),
                        cindex = double())

  # 开始网格搜索
  for (i in 81:nrow(param_grid)) {
    # 获取当前参数
    params <- param_grid[i, ]
    # 训练模型
    fit <- gbm(Surv(time, status) ~ .,
               data = x.train,
               distribution = "coxph",
               n.trees = params$n.trees,
               interaction.depth = params$interaction.depth,
               shrinkage = params$shrinkage,
               n.minobsinnode = params$n.minobsinnode,
               cv.folds = fold,
               keep.data = FALSE,
               verbose = FALSE)
    best.iter <- which.min(fit$cv.error)
    best.iter
    ####
    set.seed(12345621)
    fit <- gbm(formula = Surv(time,status)~.,
               data = x.train,
               distribution = 'coxph',
               n.trees = best.iter,
               interaction.depth = params$interaction.depth,
               shrinkage = params$shrinkage,
               n.minobsinnode = params$n.minobsinnode,
               cv.folds = fold)

    # 计算C-index
    pred_coxb=as.numeric(predict(fit,x.test,
                                 n.trees = best.iter,
                                 type = 'link'))
    cindex_coxb = 1-rcorr.cens(pred_coxb,Surv(x.test$time, x.test$status))[[1]]
    print(sprintf("cox_%d",fold))
    print(cindex_coxb)

    # 保存结果
    results <- rbind(results, data.frame(
      interaction.depth = params$interaction.depth,
      n.trees = params$n.trees,
      best_iter=best.iter,
      shrinkage = params$shrinkage,
      n.minobsinnode = params$n.minobsinnode,
      cindex = cindex_coxb))
  }
  # 查看最佳参数组合
  best_params <- results[which.max(results$cindex), ]
  print(best_params)

  # 使用最佳参数训练最终模型
  set.seed(12345621)
  fit <- gbm(Surv(time, status) ~ .,
             data = x.train,
             distribution = "coxph",
             n.trees = best_params$best_iter,
             interaction.depth = best_params$interaction.depth,
             shrinkage = best_params$shrinkage,
             n.minobsinnode = best_params$n.minobsinnode,
             cv.folds = fold,
             keep.data = TRUE,
             verbose = FALSE)


  # 在测试集上评估模型
  pred_coxb=as.numeric(predict(fit,x.test,
                               n.trees = best_params$best_iter,
                               type = 'link'))
  cindex_coxb = 1-rcorr.cens(pred_coxb,Surv(x.test$time, x.test$status))[[1]]
  print(sprintf("cox_%d",fold))
  print(cindex_coxb)

  gbm_model <- list()
  gbm_model$best_param =list(best_params)
  gbm_model$model = fit
  gbm_model$cindex = cindex_coxb
  gbm_model$pred = pred_coxb

  ### auc km ###
  cutoff = 12*1
  if ( min(x.test$time) < cutoff )
  {
    y <- survivalROC(Stime = x.test$time, status = x.test$status, marker = pred_coxb,
                     predict.time = cutoff, method = "KM")
    gbm_model$km_fp_1 = y$FP
    gbm_model$km_tp_1 = y$TP
    gbm_model$km_auc_1 = y$AUC
  }else {
    gbm_model$km_fp_1 = NA
    gbm_model$km_tp_1 = NA
    gbm_model$km_auc_1 = NA
  }

  cutoff=12*3
  y <- survivalROC(Stime = x.test$time, status = x.test$status, marker = pred_coxb,
                   predict.time = cutoff,method = "KM")
  gbm_model$km_fp_3 = y$FP
  gbm_model$km_tp_3 = y$TP
  gbm_model$km_auc_3 = y$AUC

  cutoff=12*5
  y <- survivalROC(Stime = x.test$time, status = x.test$status, marker = pred_coxb,
                   predict.time = cutoff,method = "KM")
  gbm_model$km_fp_5 = y$FP
  gbm_model$km_tp_5 = y$TP
  gbm_model$km_auc_5 = y$AUC

  dd_ext <- data.frame("time"=x.test$time, "event"=x.test$status, "score"= pred_coxb)
  Brier_score <- sbrier.score2proba(data.tr=dd_ext, data.ts=dd_ext, method="cox")
  gbm_model$brier_Score <- Brier_score
  save("gbm_model", file = sprintf("%d_gbm_result.RData", fold))

}
