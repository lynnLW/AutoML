library(xgboost)
library(survival)
library(survivalROC)
library(Hmisc)
library(survcomp)
library(caret)

xgboost_model<- function(x.train,fold){
  ########
  ########	a,b
  ########	max_depth #3
  ########	min_child_weight#18
  ########
  y<-ifelse(x.train$status==1,x.train$time,-x.train$time)
  d.train=xgb.DMatrix(data=data.matrix(x.train[,-c(1:2)]),label=y)
  y2<-ifelse(x.test$status==1,x.test$time,-x.test$time)
  d.test=xgb.DMatrix(data=data.matrix(x.test[,-c(1:2)]),label=y2)
  params<-expand.grid(max_depth = c(3:11),
                      eta = seq(0.01,0.13,0.01),
                      gamma = seq(0,0.15,0.01),
                      subsample = seq(0.5,0.9,0.05),
                      colsample_bytree = seq(0.5,0.9,0.05),
                      min_child_weight =c(3:18),
                      lambda=seq(0,1.5,0.1),
                      alpha=seq(0,1.5,0.1))

  result=data.frame(max_depth = c(),
                    eta =  c(),
                    gamma =  c(),
                    subsample =  c(),
                    colsample_bytree =  c(),
                    min_child_weight =  c(),
                    lambda= c(),
                    alpha= c(),
                    min_loss=c())
  for (i in 1:nrow(params)){
    param <-params[i,]
    param.list <- list(objective = "survival:cox",
                  eval_metric = "cox-nloglik",
                  max_depth = param$max_depth,
                  eta = param$eta,
                  gamma = param$gamma,
                  subsample = param$subsample,
                  colsample_bytree = param$colsample_bytree,
                  min_child_weight = param$min_child_weight,
                  lambda=param$lambda,
                  alpha=param$alpha)
    cv.nround = 5000
    cv.nfold = 5
    nthreads=48
    set.seed(12345621)
    mdcv <- xgb.cv(d.train,
                   params = param.list,
                   nthread=nthreads,
                   nfold=cv.nfold,
                   nrounds=cv.nround,
                   verbose = F,
                   watchlist=list(d.train),
                   early_stopping_rounds=30,
                   maximize=FALSE )

    min_loss = min(mdcv$evaluation_log[,'test_cox_nloglik_mean'])
    result=rbind(result,data.frame(
      max_depth = param$max_depth,
      eta = param$eta,
      gamma = param$gamma,
      subsample = param$subsample,
      colsample_bytree = param$colsample_bytree,
      min_child_weight = param$min_child_weight,
      lambda=param$lambda,
      alpha=param$alpha,
      min_loss=min_loss
    ))

  }

  best_param<-result[which.min(result$min_loss),]
}
