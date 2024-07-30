ML.enet<-function(est_dd,pre_var,iter.times,selected.feature,seed=12345621){
  ##### setting the pamameters ######

  seed <- seed
  final.iter.times <- iter.times
  test.iter.times<-50
  ####
  message("--- 2.Enet  ---")
  library(glmnet)
  library(survival)
  library(pbapply)

  x1 <- as.matrix(est_dd[, pre_var])
  x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))

  print("This step will select best alpha value")
  ####
  alpha_values <- seq(0.1, 0.9, 0.1)
  mean_cv_errors_list <- list()  # cv error results
  fea_list<- list()
  ###

  for (i in 1:length(alpha_values)) {
     alpha <- alpha_values[i]
     message(paste0("--- 2.Enet ", alpha, "---"))

     set.seed(seed)  # set seed before run
     # test 50 time to select best alpha
     list.of.seed <- 1:test.iter.times
     res_list  <- pblapply(list.of.seed, function(x) { # 大概运行2天
       set.seed(x)  # 只需要设置一次随机种子
       cvfit <- cv.glmnet(
         x = x1,
         y = x2,
         nfolds = 10,  # 10-fold交叉验证选取最优lambda
         alpha = alpha,
         family = "cox",  # 依赖cox模型
         maxit = 1000
       )
       # 取出最优lambda
       mean_cv_error <- mean(cvfit$cvm)
       return(list(cv_error = mean_cv_error))
     })

     # extract all cv errors
     mean_cv_errors_list[[i]] <- sapply(res_list, function(res) res$cv_error)
  }

  # extract the mean cv errors of each alpha value
  mean_cv_errors_df <- do.call(rbind, lapply(1:length(mean_cv_errors_list), function(i) {
    data.frame(
      alpha = factor(rep(alpha_values[i], length(mean_cv_errors_list[[i]]))),
      mean_cv_error = mean_cv_errors_list[[i]]
   )}))
  ###
  jpeg(filename="2.enet_mean_cross_error_different_alpha.jpg",res=600,width = 12,height = 10,units = "cm")
  # 绘图
  p <- ggplot(mean_cv_errors_df, aes(x = alpha, y = mean_cv_error,color = factor(alpha))) +
     geom_boxplot() +
     labs(title = "Mean Cross-Validation Errors for Different Alpha Values",
          x = "Alpha",
          y = "Mean Cross-Validation Error") +
     theme_classic()+
     theme(
       plot.title = element_text(size = 10, face = "bold"),
       axis.title.x = element_text(size = 8,color="black"),
       axis.title.y = element_text(size = 8,color="black"),
       axis.text.x = element_text(size = 8,color="black"),
       axis.text.y = element_text(size = 8,color="black"),
       legend.position = "none"  # 隐藏图例
     ) +
     scale_color_brewer(palette = "Set1")
    # 打印图表
   print(p)
   dev.off()

   ##select optimal alpha
   # 找到最优的alpha值
   mean_cv_errors_sum_df <- do.call(rbind, lapply(1:length(alpha_values), function(i) {
     data.frame(
       alpha = alpha_values[i],
       mean_cv_error = mean(mean_cv_errors_list[[i]])
     )}))
   best_alpha_index <- which.min(mean_cv_errors_sum_df$mean_cv_error)
   best_alpha <- alpha_values[best_alpha_index]

   # 输出最优的alpha值
   print(paste("最优的alpha值是:", best_alpha))

   ### best alpha
   alpha <- best_alpha
   message(paste0("--- 2.Enet ", alpha, "---"))

   set.seed(seed)  # 你需要在合适的位置定义 seed
   # 1000 time
   list.of.seed <- 1:final.iter.times
   res_list <- pblapply(list.of.seed, function(x) { # 大概运行2天
       set.seed(x)  # 只需要设置一次随机种子
       cvfit <- cv.glmnet(
         x = x1,
         y = x2,
         nfolds = 10,  # 10-fold交叉验证选取最优lambda
         alpha = alpha,
         family = "cox",  # 依赖cox模型
         maxit = 1000
    )
    # 取出最优lambda
       fea <- rownames(coef(cvfit, s = "lambda.min"))[coef(cvfit, s = "lambda.min")[, 1] != 0]
       if (is.element("(Intercept)", fea)) {
         lasso_fea <- sort(fea[-1])  # 去掉截距项并排序
       } else {
         lasso_fea <- sort(fea)
       }
       return(list(fea = lasso_fea))
     })

     # 提取所有种子对应的均方误差并存储
     fea_list<-lapply(list.of.seed, function(x) res_list[[x]]$fea)

     ####
     genes <- sort(table(unlist(fea_list)), decreasing = T) # 根据基因出现的频次排序
     freq.cutoff <- iter.times * 0.05
     genes <- names(genes[genes > freq.cutoff]) # 这里选择出现频次大于50的基因，认为是多次lasso的共识基因

     result <- data.frame(
       method = c(rep(paste0("Enet", "[α=", alpha, "]"), length(genes))),
       selected.fea = genes
     )
     write.table(result,file="2.enet_select_features.csv",sep=",",row.names = F)
     selected.feature <- rbind(selected.feature, result)
     return(selected.feature)
}
