install.packages("randomForestSRC")
install.packages("caret")
library(randomForestSRC)
library(caret)

# 加载数据
# 假设您的数据已经存储在 est_dd 中

# 定义参数网格
param_grid <- expand.grid(ntree = c(100, 500, 1000),
                          nodesize = c(5, 10, 15),
                          mtry = c(2, sqrt(ncol(est_dd) - 2), ncol(est_dd) - 2))

# 定义交叉验证折数
folds <- 10

# 定义保存结果的向量
results <- data.frame(ntree = integer(), nodesize = integer(), mtry = integer(), cindex = double())

# 开始网格搜索
set.seed(123456)  # 保持结果可重复
for (i in 1:nrow(param_grid)) {
  # 获取当前参数
  params <- param_grid[i, ]

  # 交叉验证
  cv_results <- train(Surv(time, status) ~ ., data = est_dd,
                      method = "rf",
                      trControl = trainControl(method = "cv", number = folds),
                      metric = "C",
                      importance = TRUE)

  # 保存结果
  results <- rbind(results, data.frame(ntree = params$ntree,
                                       nodesize = params$nodesize,
                                       mtry = params$mtry,
                                       cindex = cv_results$results$cindex[1]))
}

# 查看最佳参数组合
best_params <- results[which.max(results$cindex), ]
print(best_params)

# 使用最佳参数训练最终模型
fit <- rfsrc(Surv(time, status) ~ ., data = est_dd,
             ntree = best_params$ntree,
             nodesize = best_params$nodesize,
             mtry = best_params$mtry,
             splitrule = 'logrank',
             importance = TRUE,
             proximity = TRUE,
             forest = TRUE,
             seed = 42)

# 查看模型
print(fit)
