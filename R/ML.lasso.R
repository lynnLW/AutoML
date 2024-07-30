ML.lasso<-function(est_dd,pre_var,iter.times,selected.feature,seed=12345621){
  ### 1. Repeated Lasso  #############
  message("--- 1.Repeated lasso ---")
  x1 <- as.matrix(est_dd[, pre_var])
  x2 <- as.matrix(Surv(est_dd$OS.time, est_dd$OS))
  print("1000 time lasso penalty")
  # 1000 time lasso penalty
  lasso_fea_list <- list()
  list.of.seed <- 1:iter.times
  set.seed(seed)
  print("This step will probably take several hours")

  lasso_fea_list <- pbapply::pblapply(list.of.seed, function(x) { # about 2 days
    set.seed(list.of.seed[x])
    cvfit <- cv.glmnet(
      x = x1,
      y = x2,
      nfolds = 10, # 10-fold交叉验证选取最优lambda
      alpha = 1, # alpha = 1 意味着 lasso
      family = "cox", # 依赖cox模型
      maxit = 1000
    )

    # optimal lambda
    fea <- rownames(coef(cvfit, s = "lambda.min"))[coef(cvfit, s = "lambda.min")[, 1] != 0]
    if (is.element("(Intercept)", fea)) {
      lasso_fea <- sort(fea[-1]) # 去掉截距项并排序
    } else {
      lasso_fea <- sort(fea)
    }
    return(lasso_fea)
  })

  # 输出每次运行的基因集合
  lasso_res <- NULL
  for (i in 1:iter.times) {
    lasso_res <- rbind.data.frame(lasso_res,
                                  data.frame(
                                    iteration = i,
                                    n.gene = length(lasso_fea_list[[i]]),
                                    genelist = paste0(lasso_fea_list[[i]], collapse = " | "),
                                    stringsAsFactors = F
                                  ),
                                  stringsAsFactors = F
    )
  }

  genes <- sort(table(unlist(lasso_fea_list)), decreasing = T) # 根据基因出现的频次排序
  freq.cutoff <- iter.times * 0.05
  genes <- names(genes[genes > freq.cutoff]) # 这里选择出现频次大于50的基因，认为是多次lasso的共识基因. 95%


  result <- data.frame(
    method = c(rep("Lasso", length(genes))),
    selected.fea = genes
  )
  write.table(result,file="1.lasso_select_features.csv",sep=",",row.names = F)
  selected.feature <- rbind(selected.feature,result)
  return(selected.feature)
}
