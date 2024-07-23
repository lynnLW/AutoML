ML.xgboost<-function(est_dd,
                     seed=123456){
          ##### 5.Xgboost ###########
          message("--- 5.Xgboost  ---")
          set.seed(seed)
          train <- apply(est_dd[, -c(1)], 2, as.numeric) %>% as.data.frame()
          train_matrix <- sparse.model.matrix(OS ~ . - 1, data = train)
          train_label <- as.numeric(train$OS)
          train_fin <- list(data = train_matrix, label = train_label)
          dtrain <- xgb.DMatrix(data = train_fin$data, label = train_fin$label)
          # 模型训练
          xgb <- xgboost(
            data = dtrain, max_depth = 6, eta = 0.5,
            objective = "binary:logistic", nround = 25
          )
          # 重要重要性排序
          importance <- xgb.importance(train_matrix@Dimnames[[2]], model = xgb)
          head(importance)
          importance$rel.imp <- importance$Gain / max(importance$Gain)

          # cutoff0.05
          xgboost.variable.imp <- importance[!importance$rel.imp < 0.05, ]
          xgboost.finalVars <- xgboost.variable.imp$Feature

          result <- data.frame(
            method = c(rep("Xgboost", length(xgboost.finalVars))),
            selected.fea = xgboost.finalVars
          )

          selected.feature <- rbind(selected.feature, result)
          write.table(result,file="5.xgboost_select_features.csv",sep=",",row.names = F)
          return(selected.feature)
}
