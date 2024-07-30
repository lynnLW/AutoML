library(randomForestSRC)
library(caret)

rfsrc_model<-function(x.train,fold){
  # 定义参数网格
  set.seed(12345621)
  param_grid <- expand.grid(ntree = c(100, 500, 1000),
                            nodesize = c(5, 10, 15),
                            mtry = c(1,2,round(sqrt(ncol(x.train) - 2)),(ncol(x.train) - 2),5,7,9))

  # 定义保存结果的向量
  results <- data.frame(ntree = integer(),
                        nodesize = integer(),
                        mtry = integer(),
                        cindex = double())

  # 定义交叉验证折数
  folds <- fold
  set.seed(12345621)
  folds_indices <- createFolds(x.train$status, k = folds)

  # 开始网格搜索
  set.seed(12345621)  # 保持结果可重复
  for (i in 1:nrow(param_grid)) {
    # 获取当前参数
    params <- param_grid[i, ]
    cindex_fold <- numeric(folds)

    for (j in 1:folds) {
      # 划分训练和验证集
      train_fold <- x.train[-folds_indices[[j]], ]
      valid_fold <- x.train[folds_indices[[j]], ]

      # 交叉验证
      fit <-rfsrc(Surv(time, status) ~ .,
                  data = train_fold,
                   ntree = params$ntree,
                   nodesize = params$nodesize,
                   mtry = params$mtry,
                   splitrule = 'logrank',
                   importance = TRUE,
                   proximity = TRUE,
                   forest = TRUE,
                   seed = 12345621)

    pred_coxb <- as.numeric(predict(fit,valid_fold)$predicted)
    cindex_fold[j] = 1-rcorr.cens(pred_coxb,Surv(valid_fold$time, valid_fold$status))[[1]]
    print(sprintf("cox_%d",fold))
    print(cindex_fold[j])
    }

    # 保存结果
    results <- rbind(results, data.frame(ntree = params$ntree,
                                         nodesize = params$nodesize,
                                         mtry = params$mtry,
                                         cindex = mean(cindex_fold)))
  }

  # 查看最佳参数组合
  best_params <- results[which.max(results$cindex), ]
  print(best_params)
  # 使用最佳参数训练最终模型
  fit <- rfsrc(Surv(time, status) ~ .,
               data = x.train,
               ntree = best_params$ntree,
               nodesize = best_params$nodesize,
               mtry = best_params$mtry,
               splitrule = 'logrank',
               importance = TRUE,
               proximity = TRUE,
               forest = TRUE,
               seed = 12345621)

  # 查看模型
  pred_coxb <- as.numeric(predict(fit,x.test)$predicted)
  cindex_coxb= 1-rcorr.cens(pred_coxb,Surv(x.test$time, x.test$status))[[1]]
  print(sprintf("cox"))
  print(cindex_coxb)

  rfsrc_model <- list()
  rfsrc_model$best_param =list(best_params)
  rfsrc_model$model = fit
  rfsrc_model$cindex = cindex_coxb
  rfsrc_model$pred = pred_coxb

  ### auc km ###
  cutoff = 12*1
  if ( min(x.test$time) < cutoff )
  {
    y <- survivalROC(Stime = x.test$time, status = x.test$status, marker = pred_coxb,
                     predict.time = cutoff, method = "KM")
    rfsrc_model$km_fp_1 = y$FP
    rfsrc_model$km_tp_1 = y$TP
    rfsrc_model$km_auc_1 = y$AUC
  }else {
    rfsrc_model$km_fp_1 = NA
    rfsrc_model$km_tp_1 = NA
    rfsrc_model$km_auc_1 = NA
  }

  cutoff=12*3
  y <- survivalROC(Stime = x.test$time, status = x.test$status, marker = pred_coxb,
                   predict.time = cutoff,method = "KM")
  rfsrc_model$km_fp_3 = y$FP
  rfsrc_model$km_tp_3 = y$TP
  rfsrc_model$km_auc_3 = y$AUC

  cutoff=12*5
  y <- survivalROC(Stime = x.test$time, status = x.test$status, marker = pred_coxb,
                   predict.time = cutoff,method = "KM")
  rfsrc_model$km_fp_5 = y$FP
  rfsrc_model$km_tp_5 = y$TP
  rfsrc_model$km_auc_5 = y$AUC


  dd_ext <- data.frame("time"=x.test$time, "event"=x.test$status, "score"= pred_coxb)
  Brier_score <- sbrier.score2proba(data.tr=dd_ext, data.ts=dd_ext, method="cox")
  rfsrc_model$brier_Score <- Brier_score
  save("rfsrc_model", file = sprintf("%d_rfsrc_result.RData", fold))

}

