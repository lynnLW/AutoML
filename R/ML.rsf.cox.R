ML.rsf.cox<-function(est_dd,val_dd_list,seed=123456){
    message("---1-1 RSF ---")
    set.seed(seed)
    trainIndex <- createDataPartition(est_dd$status, p = .7,
                                      list = FALSE,
                                      times = 1)
    x.train <- est_dd[trainIndex,]
    x.test <- est_dd[-trainIndex,]

    table(x.train$status)
    table(x.test$status)

    # 定义参数网格
    tune_grid <- expand.grid(
      ntree = c(100, 200, 500,1000),
      mtry = c(3, 5, 7,9),
      nodesize = c(3, 5, 10))

    # 训练模型
    cindex_list<-data.frame()
    best_param<-c()
    best_cindex<-0.1
    set.seed(123456)
    for (i in 1:nrow(tune_grid)){
      params<-tune_grid[i,]
      fit <- rfsrc(Surv(time,status)~., data = x.train,
                   nodesize = 6,  #
                   splitrule = 'logrank',
                   importance = T,
                   proximity = T,
                   forest = T,
                   seed = seed)


      # 预测测试集
      predictions <- predict(fit, x.test)
      probs <- predict(fit, x.test, type = "prob")

      # 计算C-index（生存模型中常用的评估指标）
      library(Hmisc)
      c_index <- 1-rcorr.cens(predict(fit, x.test)$predicted, Surv(x.test$time, x.test$status))[1]
      ##
      #rs<-cbind(x.test[,c(1,2)],RS=predict(fit, x.test)$predicted)
      #c<-as.numeric(summary(coxph(Surv(time, status) ~ RS, rs))$concordance[1])
      #print(c)
      # 输出结果
      print(paste0("C-index: ",round(c_index,3)))
      rs<-data.frame(params,c_index)
      cindex_list<-rbind(cindex_list,rs)
      #
      if (c_index > best_cindex) {
        best_cindex = c_index
        best_param = params
      }
    }
    print("best_cindex_now")
    print(best_cindex)
    print("best_params_now")
    print(best_param)

    ###
    set.seed(123456)
    seed<-123456
    best_fit <- rfsrc(Surv(time,status)~., data = x.train,
                 ntree = best_param$ntree,
                 mtry = best_param$mtry,
                 nodesize = best_param$nodesize,  #
                 splitrule = 'logrank',
                 importance = T,
                 proximity = T,
                 forest = T,
                 seed = seed)

    rs <- lapply(val_dd_list, function(x){cbind(x[,1:2], RS=predict(fit, newdata = x)$predicted)})
    rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)
    ###
    cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])})) %>%
      rownames_to_column('ID')
    cc$Model <- 'RSF'
    result <- rbind(result, cc)
    ml.res[['RSF']] =best_fit
    riskscore[['RSF']] = rs

}
