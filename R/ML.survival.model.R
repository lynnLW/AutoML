#' Training Survival Models Using 13 Machine Learning Algorithms
#'
#' @param train_data Data.frame for training (colnames: ID, time, status, features)
#' @param candidate_genes Features used to build the model
#' @param filter_OS_time Logical, filter samples with survival time >30 days
#' @param meta_time Time unit for follow-up ("d", "m", "y"), default "m"
#' @param cor Logical, remove highly correlated features
#' @param cor_threshold Correlation threshold (0-1), default 0.85
#' @param fold Number of cross-validation folds
#' @param rep Number of repeats for cross-validation
#' @param p Proportion of data for training (0-1), default 1
#' @param deep Logical, use deep learning models (deepsurv/deephit)
#' @param outdir Output directory path
#' @param seed Random seed
#' @param ncore Number of CPU cores for parallelization
#' @importFrom dplyr %>% group_by summarise
#' @importFrom ggplot2 ggplot aes geom_point
#' @importFrom stats predict median model.weights
#' @importFrom survival Surv
#'
#' @return List containing trained models and evaluation metrics
#' @export
ML.survival.model = function(train_data,
                       candidate_genes,
                       filter_OS_time=F,
                       meta_time="m",
                       cor=F,
                       cor_threshold=0.85,
                       fold=5,
                       rep=10,
                       p=1,
                       deep=F,
                       outdir="2.train/",
                       seed=123,
                       ncore=4
){
  # package check---------------------------------------------------------------
  required_pkgs <- c(
    "survival", "glmnet", "randomForestSRC", "plsRcox",
    "gbm", "CoxBoost", "xgboost", "survcomp","mboost","superpc"
  )

  if (deep) {
    required_pkgs <- c(required_pkgs, "keras", "tensorflow", "survivalmodels")
  }

  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(
      "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
      "\nInstall with: install.packages(c(",
      paste0("'", missing_pkgs, "'", collapse = ", "), "))"
    )
  }

  ##
  if (!dir.exists(outdir)){
      dir.create(outdir,recursive = T)
  }

  # Replace '-' in column names with '.'
  print("Checking data feasibility")
  candidate_genes = gsub('-','.',candidate_genes)
  colnames(train_data) = gsub('-','.',colnames(train_data))
  colnames(train_data)[2]="time"
  colnames(train_data)[3]="status"

  # Matching candidate genes to genes in each cohort
  print("Matching candidate genes to genes in each cohort")
  common_feature = c('ID','time','status',candidate_genes)
  common_feature = intersect(common_feature,colnames(train_data))

  # Matching common feature in each cohort
  train_data<-train_data[,common_feature]

  # higly correlated genes
  if(cor){
      ###correlation
      data<-train_data
      cor_matrix <- cor(data[,-c(1:3)])
      cor_threshold <- cor_threshold

      highly_correlated_features <- caret::findCorrelation(cor_matrix, cutoff = cor_threshold)

      if ( length(highly_correlated_features)>0){
        print("Highly correlation genes:")
        print(colnames(data[,-c(1:3)])[highly_correlated_features])

        common_feature<-c('ID','time','status',colnames(data[,-c(1:3)])[-highly_correlated_features])
        train_data<-train_data[, common_feature]
      } else {
        train_data<-train_data
      }

   }

  message(paste0('---the number of the raw candidate genes is ', length(candidate_genes),' ---'))
  message(paste0('---the number of common genes across the training data is ', length(common_feature)-3,' ---'))
  print(common_feature)
  save(common_feature,file=paste0(outdir,"/train_features.Rdata"))

  ## If keep follow up days more than 30 days
  if (filter_OS_time){
    if (meta_time=="m") {
      train_data <- train_data[train_data$time >= 1, ] # days more than 30
      nsample<-nrow(train_data)
      print(paste0(nsample," with more than 30 follow-up days for the next step"))
    } else if (meta_time=="y"){
      train_data <- train_data[train_data$time >= 0.083, ]
      nsample<-nrow(train_data)
      print(paste0(nsample," with more than 30 follow-up days for the next step"))
    } else {
      train_data <- train_data[train_data$time >= 30, ]
      nsample<-nrow(train_data)
      print(paste0(nsample," with more than 30 follow-up days for the next step"))
    }
  } else {
    train_data <- train_data # days more than 30
  }

  # split into 75% training and 25% testing data
  set.seed(seed)
  est_dd <- as.data.frame(train_data)[, common_feature[-1]]
  train.index=caret::createDataPartition(est_dd$status,p=p,list=F)
  d.train=est_dd[train.index,]
  d.test=est_dd[-train.index,]
  print("Data preprocessing completed")

  ## the output data
  model_list<-list()

  t1<-Sys.time()

  ## the supporting function--------------------------------------------------
  rfrsf_model<-function(d.train,d.test,fold,rep,outdir,seed,ncore){

    # Record start time
    t1<-Sys.time()

    # Set up parallel backend using `future`
    if (.Platform$OS.type == "unix") {
      # Unix: Mac/Linux
      future::plan(future::multicore, workers = ncore)
    } else {
      # Windows
      future::plan(future::multisession, workers = ncore)
    }

    param_grid <- expand.grid(
      ntree = c(100,500,1000,2000),
      nodesize = c(1,3,5,7,10,15),
      mtry = c(1, 2, round(sqrt(ncol(d.train) - 2)), (ncol(d.train) - 2), 5, 7, 9)
    )

    # Create k-folds for cross-validation
    set.seed(seed)
    folds_list <- create_folds(d.train, fold = fold, nrepeats = rep, strata = "status", seed = seed)

    # the supporting function-----------------------------
    rfrsf_cox<-function(train_fold,valid_fold,test_fold,params,seed){
      ## modeling
      set.seed(seed)
      fit <-randomForestSRC::rfsrc(Surv(time, status) ~ .,
                  data = train_fold,
                  ntree = params$ntree,
                  nodesize = params$nodesize,
                  mtry = params$mtry,
                  splitrule = 'logrank',
                  importance = TRUE,
                  proximity = TRUE,
                  forest = TRUE,
                  seed = seed)

      model_name="RFRSF"
      # saving result
      model <- list(
        best_param = list(params),
        model = fit,
        train = cal_metrics(train_fold, fit, model_name),
        valid = cal_metrics(valid_fold, fit, model_name),
        test = cal_metrics(test_fold, fit, model_name)
      )
      return(model)
    }


    # the main function-----------------------------------
    # Parallel grid search over hyperparameters
    results <- future.apply::future_lapply(1:nrow(param_grid), function(i) {

      # Get current parameters
      params <- param_grid[i, ]
      output_list <- list()  # Collect results for each repetition
      display_progress(index = i, totalN = nrow(param_grid))

      # Parallelize the repetitions loop
      output_list <- future.apply::future_lapply(1:rep, function(j) {
        tryCatch({
          folds <- folds_list[[j]]
          model_list <- vector("list", length = fold)

          for (k in 1:fold) {
            # Get train and valid datasets
            train_fold <- d.train[-folds[[k]], ]
            valid_fold <- d.train[folds[[k]], ]
            model_list[[k]] <- rfrsf_cox(train_fold, valid_fold, test_fold = d.test, params, seed)

          }

          # Calculate cindex and bs score
          index_list <- extract_metrics(model_list)
          return(data.frame(index_list$train, index_list$valid, index_list$test))
        }, error = function(e) {
          # Catch errors and log them without stopping the loop
          paste("Error in parameter set", i, "repeat", j, ":", e$message)
          return(data.frame(NA, NA, NA, NA, NA, NA, NA, NA))  # Return NA for error cases
        })
      }, future.seed = seed)

      # Aggregate results from all repetitions
      output <- do.call(rbind, output_list)

      # Return results for this parameter set
      return(data.frame(
        ntree = params$ntree,
        nodesize = params$nodesize,
        mtry = params$mtry,
        cindex_tr = mean(output[, 1]),
        bs_tr = mean(output[, 2]),
        cindex_tv = mean(output[, 9]),
        bs_tv = mean(output[, 10])))
    }, future.seed = seed)

    # Aggregate results from all repetitions
    results_df<-data.frame()
    results_df <- do.call(rbind, results)

    # View the best parameter combination based on validation cindex
    best_params <- results_df[which.max(results_df$cindex_tv), ]
    print(best_params[,1:3])

    # Saving final model with the best parameters
    folds_list<- create_folds(d.train,fold=fold,nrepeats = 2,strata="status",seed=seed)
    folds_list<-c(folds_list[[1]],folds_list[[2]])
    ##
    totalN=fold*2
    model_list<-list()
    for (j in 1:totalN) {
      print(j)
      test_index<-folds_list[[j]]
      train_index <- setdiff(seq_len(nrow(d.train)), test_index)
      train_fold <- d.train[train_index, ]
      valid_fold <- d.train[test_index, ]

      model_list[[j]]<-rfrsf_cox(train_fold,valid_fold,d.test,best_params,seed)
    }

    ###total data model
    final_model<-rfrsf_cox(d.train,d.test,NULL,best_params,seed)
    metrics_list<-extract_metrics(model_list)

    save("model_list", file = paste0(outdir,"/",sprintf("%d_%d_RFRSF_result.RData",rep,fold)))
    save("final_model", file = paste0(outdir,"/",sprintf("%d_%d_final_RFRSF_result.RData",rep,fold)))
    save("metrics_list", file = paste0(outdir,"/",sprintf("%d_%d_RFRSF_cindex_result.RData",rep,fold)))

    # record time
    t2 <- Sys.time()
    run_time<-t2-t1
    cat("Total elapsed time:", run_time, "seconds\n")

    return(list(final_model=final_model,metrics_list=metrics_list))

  }

  lasso_cox<-function(train_fold,valid_fold,test_fold,alpha,fold,model_name,seed){
    ## set different random seeds
    set.seed(seed)
    cv_fit <- glmnet::cv.glmnet(as.matrix(train_fold[,3:ncol(train_fold)]),
                        as.matrix(train_fold[,1:2]),
                        family = "cox",
                        alpha = alpha,
                        nfolds = fold)

    final_fit <- glmnet::glmnet(as.matrix(train_fold[,3:ncol(train_fold)]),
                        as.matrix(train_fold[,1:2]),
                        family = "cox",
                        alpha=alpha,
                        lambda = cv_fit$lambda.min)

    # saving result
    model <- list(
      best_param = list(s = cv_fit$lambda.min, alpha = alpha),
      model = final_fit,
      train = cal_metrics(train_fold, final_fit, model_name),
      valid = cal_metrics(valid_fold, final_fit, model_name),
      test = cal_metrics(test_fold, final_fit, model_name)
    )
    return(model)
  }

  lasso_model <- function(d.train,d.test,fold,rep,model_name,outdir,seed=123){

    if (is.null(d.train)) stop("Training data cannot be NULL")

    ##alpha=0 ridge; alpha=1, lasso; alpha=0-1, enet;
    set.seed(seed)
    folds_list<- create_folds(d.train,fold=fold,nrepeats = 2,strata="status",seed=seed)
    folds_list<-c(folds_list[[1]],folds_list[[2]])
    totalN=fold*2
    ## saving result
    model_list<-list()

    ## alpha
    if (model_name=="Lasso"){
      alpha=1
    } else if (model_name=="Ridge"){
      alpha=0
    }

    ## --------------
    for (j in 1:totalN) {
      print(j)
      display_progress(j,totalN)

      test_index<-folds_list[[j]]
      train_index <- setdiff(seq_len(nrow(d.train)),test_index)

      train_fold <- d.train[train_index, ]
      valid_fold <- d.train[test_index, ]

      ## train data
      model_list[[j]]<-lasso_cox(train_fold,valid_fold,d.test,alpha,fold,model_name,seed)
    }

    metrics_list<-extract_metrics(model_list)

    final_model<-lasso_cox(d.train,d.test,NULL,alpha,fold,model_name,seed)

    save("model_list", file = paste0(outdir,"/",sprintf("%d_%d_%s_result.RData",rep,fold,model_name)))
    save("final_model", file = paste0(outdir,"/",sprintf("%d_%d_final_%s_result.RData",rep,fold,model_name)))
    save("metrics_list", file = paste0(outdir,"/",sprintf("%d_%d_%s_cindex_result.RData",rep,fold,model_name)))
    return(list(final_model=final_model,metrics_list=metrics_list))
  }

  enet_model <- function(d.train,d.test,fold,rep,outdir,seed){

    ##alpha=0 ridge; alpha=1, lasso; alpha=0-1, enet;
    set.seed(seed)
    folds_list<- create_folds(d.train,fold=fold,nrepeats = 1,strata="status",seed)

    ### searching the grid
    print("Searching the grid")
    model_list<-list()
    metrics_list<-data.frame()
    alpha_output<-data.frame(alpha=double(),
                             cindex_tr=double(),
                             bs_tr=double(),
                             cindex_ts=double(),
                             bs_ts=double())


    for (alpha in seq(0.1,0.9,0.1)) {
      # current parameter
      print("alpha")
      print(alpha)

      for (j in 1:fold) {

        test_index<-folds_list[[1]][[j]]
        train_index <- setdiff(seq_len(nrow(d.train)), test_index)

        train_fold <- d.train[train_index, ]
        valid_fold <- d.train[test_index, ]

        ## train data
        model_list[[j]]<-lasso_cox(train_fold,valid_fold,d.test,alpha,fold,"Enet",seed)

      }

      ###calculate mean cindex and bs score
      metrics_list<-extract_metrics(model_list)
      alpha_output<-rbind(alpha_output,data.frame(alpha=alpha,
                                                  cindex_tr=mean(metrics_list$train$cindex),
                                                  bs_tr=mean(metrics_list$train$bs),
                                                  cindex_ts=mean(metrics_list$valid$cindex),
                                                  bs_ts=mean(metrics_list$valid$bs)))

    }

    # best parameter
    best_alpha <- alpha_output[which.max(alpha_output$cindex_ts),]$alpha
    print(best_alpha)

    # best model
    set.seed(seed)
    folds_list<- create_folds(d.train,fold=fold,nrepeats = 2,strata="status",seed=seed)
    folds_list<-c(folds_list[[1]],folds_list[[2]])
    totalN=fold*2
    model_list<-list()
    metrics_list<-data.frame()

    for (j in 1:totalN) {
      print(j)
      test_index<-folds_list[[j]]
      train_index <- setdiff(seq_len(nrow(d.train)), test_index)

      train_fold <- d.train[train_index, ]
      valid_fold <- d.train[test_index, ]

      ## train data
      model_list[[j]]<-lasso_cox(train_fold,valid_fold,d.test,alpha=best_alpha,fold,"Enet",seed)

    }

    ###total data model
    final_model<-lasso_cox(d.train,d.test,NULL,alpha=best_alpha,fold,"Enet",seed)
    ###calculate mean cindex and bs score
    metrics_list<-extract_metrics(model_list)

    save("model_list", file = paste0(outdir,"/",sprintf("%d_%d_Enet_result.RData",rep,fold)))
    save("final_model", file = paste0(outdir,"/",sprintf("%d_%d_final_Enet_result.RData",rep,fold)))
    save("metrics_list", file = paste0(outdir,"/",sprintf("%d_%d_Enet_cindex_result.RData",rep,fold)))
    return(list(final_model=final_model,metrics_list=metrics_list))
  }

  cox_boost <- function(d.train,d.test,fold,rep,outdir,seed){
    set.seed(seed)
    folds_list<- create_folds(d.train,fold=fold,nrepeats = 2,strata="status",seed=seed)
    folds_list<-c(folds_list[[1]],folds_list[[2]])
    # the main function-----------------------------
    coxboost_cox<-function(train_fold,valid_fold,test_fold,fold,seed){
      ##
      set.seed(seed)
      pen <- CoxBoost::optimCoxBoostPenalty(train_fold[,'time'],train_fold[,'status'],
                                  as.matrix(train_fold[,-c(1,2)]),
                                  trace=TRUE,start.penalty=500,parallel = T)

      #K>>number of folds to be used for cross-validation
      set.seed(seed)
      cv.res <- CoxBoost::cv.CoxBoost(train_fold[,'time'],
                            train_fold[,'status'],
                            as.matrix(train_fold[,-c(1,2)]),
                            maxstepno=500,
                            K=fold,
                            type="verweij",
                            penalty=pen$penalty)

      set.seed(seed)
      fit <- CoxBoost::CoxBoost(train_fold[,'time'],train_fold[,'status'],
                      as.matrix(train_fold[,-c(1,2)]),
                      stepno=cv.res$optimal.step,
                      penalty=pen$penalty)

      # saving result
      model <- list(
        best_param = list(stepno=cv.res$optimal.step, penalty=pen$penalty),
        model = fit,
        train = cal_metrics(train_fold, fit, "CoxBoost"),
        valid = cal_metrics(valid_fold, fit, "CoxBoost"),
        test = cal_metrics(test_fold, fit, "CoxBoost")
      )
      return(model)

    }
    ##
    totalN=fold*2
    model_list<-list()
    for (j in 1:totalN) {
      print(j)
      display_progress(index = j, totalN = totalN)
      test_index<-folds_list[[j]]
      train_index <- setdiff(seq_len(nrow(d.train)), test_index)

      train_fold <- d.train[train_index, ]
      valid_fold <- d.train[test_index, ]

      ## train data
      model_list[[j]]<-coxboost_cox(train_fold,valid_fold,d.test,fold,seed)
    }

    ###total data model
    final_model<-coxboost_cox(d.train,d.test,NULL,fold,seed)

    metrics_list<-extract_metrics(model_list)

    save("model_list", file = paste0(outdir,"/",sprintf("%d_%d_CoxBoost_result.RData",rep,fold)))
    save("final_model", file = paste0(outdir,"/",sprintf("%d_%d_final_CoxBoost_result.RData",rep,fold)))
    save("metrics_list", file = paste0(outdir,"/",sprintf("%d_%d_CoxBoost_cindex_result.RData",rep,fold)))
    return(list(final_model=final_model,metrics_list=metrics_list))
  }

  plsR_model <- function(d.train,d.test,fold,rep,outdir,seed){

    set.seed(seed)
    folds_list<- create_folds(d.train,fold=fold,nrepeats = 2,strata="status",seed=seed)
    folds_list<-c(folds_list[[1]],folds_list[[2]])
    ## modeling -------------------------------------------
    plsr_cox<-function(train_fold,valid_fold,test_fold,fold,seed){
      ##
      set.seed(seed)
      cv.plsRcox.res=plsRcox::cv.plsRcox(list(x=train_fold[,3:ncol(train_fold)],
                                     time=train_fold$time,
                                     status=train_fold$status),
                                nfold = fold,
                                nt=10,
                                verbose = FALSE)
      fit <- plsRcox::plsRcox(train_fold[,3:ncol(train_fold)],
                     time=train_fold$time,
                     event=train_fold$status,
                     nt=as.numeric(cv.plsRcox.res[5]))

      model_name="plsRcox"
      # saving
      model <- list(
        best_param = list(nt=as.numeric(cv.plsRcox.res[5])),
        model = fit,
        train = cal_metrics(train_fold, fit, model_name),
        valid = cal_metrics(valid_fold, fit, model_name),
        test = cal_metrics(test_fold, fit, model_name)
      )
      return(model)
    }

    ## repetition------------------------------------------
    totalN=fold*2
    model_list<-list()
    for (j in 1:totalN) {
      print(j)
      test_index<-folds_list[[j]]
      train_index <- setdiff(seq_len(nrow(d.train)), test_index)

      train_fold <- d.train[train_index, ]
      valid_fold <- d.train[test_index, ]

      ## train data
      model_list[[j]]<-plsr_cox(train_fold,valid_fold,d.test,fold,seed)
    }

    ###total data model
    final_model<-plsr_cox(d.train,d.test,NULL,fold,seed)
    metrics_list<-extract_metrics(model_list)

    ###
    save("model_list", file = paste0(outdir,"/",sprintf("%d_%d_plsRcox_result.RData",rep,fold)))
    save("final_model", file = paste0(outdir,"/",sprintf("%d_%d_final_plsRcox_result.RData",rep,fold)))
    save("metrics_list", file = paste0(outdir,"/",sprintf("%d_%d_plsRcox_cindex_result.RData",rep,fold)))
    return(list(final_model=final_model,metrics_list=metrics_list))

  }

  superpc_model <- function(d.train,d.test,fold,rep,outdir,seed){

    set.seed(seed)
    folds_list<- create_folds(d.train,fold=fold,nrepeats = 2,strata="status",seed=seed)
    folds_list<-c(folds_list[[1]],folds_list[[2]])
    ## modeling-----------------------------------------------
    superpc_cox<-function(train_fold,valid_fold,test_fold,fold,seed){
      ##
      set.seed(seed)
      data <- list(x=t(train_fold[,-c(1,2)]),
                   y=train_fold$time,
                   censoring.status=train_fold$status,
                   featurenames=colnames(train_fold)[-c(1,2)])
      valid <- list(x=t(valid_fold[,-c(1,2)]),
                    y=valid_fold$time,
                    censoring.status=valid_fold$status,
                    featurenames=colnames(valid_fold)[-c(1,2)])

      ####
      set.seed(seed)
      fit <-superpc::superpc.train(data = data,
                           type = 'survival',s0.perc = 0.5) #default
      cv.fit <- superpc::superpc.cv(fit,data,
                           n.threshold = 20,
                           n.fold = fold, #Number of cross-validation folds
                           n.components=3,
                           min.features=2,#default
                           max.features=nrow(data$x),
                           compute.fullcv= TRUE,
                           compute.preval=TRUE)

      model_name="SuperPC"
      fit.list<-list(fit=fit,cv.fit=cv.fit,train_data=data)
      # saving
      result <- list(
        best_param =list(threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])]),
        model = fit.list,
        train = cal_metrics(train_fold,fit.list, model_name),
        valid = cal_metrics(valid_fold, fit.list, model_name),
        test = cal_metrics(test_fold,fit.list, model_name)
      )
      return(result)
    }

    ## -------------------------------------------------------
    totalN=fold*2
    model_list<-list()
    for (j in 1:totalN) {
      print(j)

      test_index<-folds_list[[j]]
      train_index <- setdiff(seq_len(nrow(d.train)), test_index)

      train_fold <- d.train[train_index, ]
      valid_fold <- d.train[test_index, ]

      ## train data
      model_list[[j]]<-superpc_cox(train_fold,valid_fold,d.test,fold,seed)
    }

    ###total data model
    final_model<-superpc_cox(d.train,d.test,NULL,fold,seed)
    metrics_list<-extract_metrics(model_list)

    save("model_list", file = paste0(outdir,"/",sprintf("%d_%d_SuperPC_result.RData",rep,fold)))
    save("final_model", file = paste0(outdir,"/",sprintf("%d_%d_final_SuperPC_result.RData",rep,fold)))
    save("metrics_list", file = paste0(outdir,"/",sprintf("%d_%d_SuperPC_cindex_result.RData",rep,fold)))
    return(list(final_model=final_model,metrics_list=metrics_list))
  }

  xgboost_model<- function(d.train,d.test,fold,rep,outdir,seed){

    # the main function-----------------------------------

    xgboost_cox<-function(train_fold,valid_fold,test_fold,best_param,nround,seed){
      ###
      y<-ifelse(train_fold$status==1,train_fold$time,-train_fold$time)
      train=xgboost::xgb.DMatrix(data=data.matrix(train_fold[,-c(1:2)]),label=y)
      ###
      set.seed(seed)
      fit <- xgboost::xgboost(data=train,
                     params=best_param,
                     nrounds=nround)

      model_name="XGBoost"
      model <- list(
        best_param = list(s=best_param,nround=nround),
        model = fit,
        train = cal_metrics(train_fold, fit, model_name),
        valid = cal_metrics(valid_fold, fit, model_name),
        test = cal_metrics(test_fold, fit, model_name)
      )
      return(model)
    }



    ###searching parameter--------------------------------
    set.seed(seed)
    train.index=caret::createDataPartition(d.train$status,p=0.7,list=F)
    x.train=d.train[train.index,]
    x.test=d.train[-train.index,]

    ########	a,b
    ########	max_depth #
    ########	min_child_weight#
    max_depth. <- c(2:15)
    min_child_weight. <- c(2:20)
    nthreads=10
    best_loss=10
    y<-ifelse(x.train$status==1,x.train$time,-x.train$time)
    train=xgboost::xgb.DMatrix(data=data.matrix(x.train[,-c(1:2)]),label=y)
    y2<-ifelse(x.test$status==1,x.test$time,-x.test$time)
    test=xgboost::xgb.DMatrix(data=data.matrix(x.test[,-c(1:2)]),label=y2)

    for (a in max_depth.){
      for (b in min_child_weight.){
        print("max_depth, min_child_weight")
        print(a)
        print(b)
        param <- list(objective = "survival:cox",
                      eval_metric = "cox-nloglik",
                      max_depth = a,
                      eta = 0.01,
                      gamma = 0,
                      subsample = 0.8,
                      colsample_bytree = 0.8,
                      min_child_weight = b,
                      lambda=1,
                      alpha=0
        )
        cv.nround = 5000
        cv.nfold = fold
        set.seed(seed)
        mdcv <- xgboost::xgb.cv(train,
                       params = param,
                       nthread=nthreads,
                       nfold=cv.nfold,
                       nrounds=cv.nround,
                       verbose = F,
                       early_stopping_rounds=30,
                       maximize=FALSE )

        if("try-error" %in% class(mdcv)) {
          next
        }
        min_loss = min(mdcv$evaluation_log[,'test_cox_nloglik_mean'])
        min_loss_index = which.min(as.numeric(unlist(mdcv$evaluation_log[,'test_cox_nloglik_mean'])))
        if (min_loss < best_loss) {
          best_loss<-min_loss
          best_param = param
        }

      }}

    max_depth. = best_param$max_depth  ##3
    min_child_weight. = best_param$min_child_weight  ##18

    ########	c
    ########	gamma
    gamma. <- seq(0,0.15,0.01)
    for (c in gamma.){
      print("gamma")
      print(c)
      param <- list(objective = "survival:cox",
                    eval_metric = "cox-nloglik",
                    max_depth = max_depth.,
                    eta = 0.01,
                    gamma = c,
                    subsample = 0.8,
                    colsample_bytree = 0.8,
                    min_child_weight = min_child_weight.,
                    lambda=1,
                    alpha=0
      )
      cv.nround = 5000
      cv.nfold = fold
      set.seed(seed)
      mdcv <- xgboost::xgb.cv(train,
                     params = param,
                     nthread=nthreads,
                     nfold=cv.nfold,
                     nrounds=cv.nround,
                     verbose = F,
                     early_stopping_rounds=30,
                     maximize=FALSE )

      if("try-error" %in% class(mdcv))
      {
        next
      }

      min_loss = min(mdcv$evaluation_log[,'test_cox_nloglik_mean'])
      if (min_loss < best_loss) {
        best_loss = min_loss
        best_param = param
      }

    }

    ########	d,e
    ########	subsample
    ########	colsample_bytree

    gamma. = best_param$gamma  ##0.14

    subsample. <- seq(0.60,0.95,0.05)
    colsample_bytree. <- seq(0.60,0.95,0.05)

    for (d in subsample.){
      for (e in colsample_bytree.){
        print("subsample, colsample_bytree")
        print(d)
        print(e)
        param <- list(objective = "survival:cox",
                      eval_metric = "cox-nloglik",
                      max_depth = max_depth.,
                      eta = 0.01,
                      gamma = gamma.,
                      subsample = d,
                      colsample_bytree = e,
                      min_child_weight = min_child_weight.,
                      lambda=1,
                      alpha=0
        )
        cv.nround = 5000
        cv.nfold = fold
        set.seed(seed)
        mdcv <- xgboost::xgb.cv(train,
                       params = param,
                       nthread=nthreads,
                       nfold=cv.nfold,
                       nrounds=cv.nround,
                       verbose = F,
                       early_stopping_rounds=30,
                       maximize=FALSE )

        if("try-error" %in% class(mdcv))
        {
          next
        }

        min_loss = min(mdcv$evaluation_log[,'test_cox_nloglik_mean'])
        print(min_loss)

        if (min_loss < best_loss) {
          best_loss = min_loss
          best_param = param
        }
      }}

    ########	f,g
    ########	alpha
    ########	lambda

    subsample. = best_param$subsample
    colsample_bytree. = best_param$colsample_bytree

    alpha. <- seq(0,1.5,0.1)
    lambda. <- seq(0,1.5,0.1)

    for (f in alpha.){
      for (g in lambda.){
        print("alpha, lambda")
        print(f)
        print(g)
        param <- list(objective = "survival:cox",
                      eval_metric = "cox-nloglik",
                      max_depth = max_depth.,
                      eta = 0.01,
                      gamma = gamma.,
                      subsample = subsample.,
                      colsample_bytree = colsample_bytree.,
                      min_child_weight = min_child_weight.,
                      lambda=g,
                      alpha=f)
        cv.nround = 5000
        cv.nfold = fold
        set.seed(seed)
        mdcv <- xgboost::xgb.cv(train,
                       params = param,
                       nthread=nthreads,
                       nfold=cv.nfold,
                       nrounds=cv.nround,
                       verbose = F,
                       early_stopping_rounds=30,
                       maximize=FALSE )

        if("try-error" %in% class(mdcv))
        {
          next
        }

        min_loss = min(mdcv$evaluation_log[,'test_cox_nloglik_mean'])

        print(min_loss)
        if (min_loss < best_loss) {
          best_loss = min_loss
          best_param = param
        }
      }}

    ########	h
    ########	eta

    alpha. = best_param$alpha   ##0
    lambda. = best_param$lambda   ##1

    eta. <- seq(0.01,0.5,0.01)

    for (h in eta.){
      print("eta")
      print(h)
      param <- list(objective = "survival:cox",
                    eval_metric = "cox-nloglik",
                    max_depth = max_depth.,
                    eta = h,
                    gamma = gamma.,
                    subsample = subsample.,
                    colsample_bytree = colsample_bytree.,
                    min_child_weight = min_child_weight.,
                    lambda=lambda.,
                    alpha=alpha.
      )
      cv.nround = 5000
      cv.nfold = fold
      set.seed(seed)
      mdcv <- xgboost::xgb.cv(train,
                     params = param,
                     nthread=nthreads,
                     nfold=cv.nfold,
                     nrounds=cv.nround,
                     verbose = F,
                     early_stopping_rounds=30,
                     maximize=FALSE )

      if("try-error" %in% class(mdcv))
      {
        next
      }

      min_loss = min(mdcv$evaluation_log[,'test_cox_nloglik_mean'])

      if (min_loss < best_loss) {
        best_loss = min_loss
        best_param = param
      }

    }

    cv.nfold = fold
    cv.nround <- seq(10,30,1)
    best_cindex=0
    for (j in cv.nround){
      print("cv.nround")
      print(j)
      set.seed(seed)
      fit <- xgboost::xgboost(data = train,verbose = F,
                     params = best_param, nrounds = j)

      pred_coxb = predict(fit,data.matrix(x.test[,-c(1:2)]))
      pred_coxb <- as.numeric(pred_coxb)

      cindex_xgb = 1-Hmisc::rcorr.cens(pred_coxb,Surv(x.test$time, x.test$status))[[1]]
      print(cindex_xgb)

      if (cindex_xgb > best_cindex) {
        best_cindex = cindex_xgb
        best_param = best_param
        best_nround=j
      }
    }

    nround = best_nround
    set.seed(seed)
    folds_list<- create_folds(d.train,fold=fold,nrepeats = 2,strata="status",seed=seed)
    folds_list<-c(folds_list[[1]],folds_list[[2]])
    ##
    totalN=fold*2
    model_list<-list()
    for (j in 1:totalN) {
      print(j)

      test_index<-folds_list[[j]]
      train_index <- setdiff(seq_len(nrow(d.train)), test_index)

      train_fold <- d.train[train_index,]
      valid_fold <- d.train[test_index,]

      ## train data
      model_list[[j]]<-xgboost_cox(train_fold,valid_fold,d.test,
                                   best_param,nround,seed)
    }

    ###total data model
    final_model<-xgboost_cox(d.train,d.test,NULL,
                             best_param,nround,seed)
    metrics_list<-extract_metrics(model_list)

    save("model_list", file = paste0(outdir,"/",sprintf("%d_%d_XGBoost_result.RData",rep,fold)))
    save("final_model", file = paste0(outdir,"/",sprintf("%d_%d_final_XGBoost_result.RData",rep,fold)))
    save("metrics_list", file = paste0(outdir,"/",sprintf("%d_%d_XGBoost_cindex_result.RData",rep,fold)))
    return(list(final_model=final_model,metrics_list=metrics_list))

  }

  black_model<-function(d.train,d.test,fold,rep,outdir,seed){

    folds_list<- create_folds(d.train,fold=fold,nrepeats = 2,strata="status",seed=seed)
    folds_list<-c(folds_list[[1]],folds_list[[2]])
    ##the main function-------------------------
    black_cox<-function(train_fold,valid_fold,test_fold,seed){
      set.seed(seed)
      fit <-mboost::blackboost(Surv(time, status) ~ .,
                        data = train_fold,
                        family = mboost::CoxPH(),
                        control = mboost::boost_control(mstop = 500))

      cv <- mboost::cv(stats::model.weights(fit), type = "kfold")
      cvm <- mboost::cvrisk(fit, folds = cv, papply = lapply)
      n=mboost::mstop(cvm)
      set.seed(seed)
      final.fit <- mboost::blackboost(Surv(time, status) ~ .,
                              data=train_fold,
                              family = mboost::CoxPH(),
                              control=mboost::boost_control(mstop = n))

      model_name="BlackBoost"
      model <- list(
        best_param = list(mstop=n),
        model = final.fit,
        train = cal_metrics(train_fold, final.fit, model_name),
        valid = cal_metrics(valid_fold, final.fit, model_name),
        test = cal_metrics(test_fold, final.fit, model_name)
      )
      return(model)
    }


    ##------------------------------------------
    totalN=fold*2
    model_list<-list()
    for (j in 1:totalN) {
      display_progress(index = j, totalN = totalN)
      test_index<-folds_list[[j]]
      train_index <- setdiff(seq_len(nrow(d.train)), test_index)

      train_fold <- d.train[train_index, ]
      valid_fold <- d.train[test_index, ]

      ## train data
      model_list[[j]]<-black_cox(train_fold,valid_fold,d.test,seed)
    }

    ###total data model
    final_model<-black_cox(d.train,d.test,NULL,seed)

    metrics_list<-extract_metrics(model_list)

    save("model_list", file = paste0(outdir,"/",sprintf("%d_%d_BlackBoost_result.RData",rep,fold)))
    save("final_model", file = paste0(outdir,"/",sprintf("%d_%d_final_BlackBoost_result.RData",rep,fold)))
    save("metrics_list", file = paste0(outdir,"/",sprintf("%d_%d_BlackBoost_cindex_result.RData",rep,fold)))
    return(list(final_model=final_model,metrics_list=metrics_list))

  }

  DeepHit_model<-function(d.train,d.test,fold,rep,outdir,seed){

    t1<-Sys.time()

    ### created nrepeats k-fold
    set.seed(seed)
    folds_list<- create_folds(d.train,fold=fold,nrepeats = 1,strata="status",seed)

    # param grid
    generate_param_grid <- function(layer) {
      param_list <- list(lr = c(0.001, 0.01, 0.1),
                         nd1 = c(8, 16, 32, 64, 128),
                         epochs = c(50, 100, 200),
                         batch = c(32, 64, 128, 256))

      if (layer > 1) param_list$nd2 <- c(8, 16, 32, 64, 128)
      if (layer > 2) param_list$nd3 <- c(8, 16, 32, 64, 128)
      if (layer > 3) param_list$nd4 <- c(8, 16, 32, 64, 128)

      param_list$layer=layer

      expand.grid(param_list)
    }

    # modeling function---------------------------------------
    DeepHit_cox<-function(train_fold,valid_fold,test_fold,params,seed){
      ###
      if(params$layer==1){
        ##1 layer
        survivalmodels::set_seed(seed)
        fit =survivalmodels::deephit(Surv(time, status) ~ .,
                      data = train_fold,
                      frac = 0.2,
                      learning_rate=params$lr,
                      activation = "relu",
                      num_nodes = params$nd1,
                      early_stopping = TRUE,
                      epochs = params$epochs,
                      batch_size = params$batch)
      } else if(params$layer==2){
        ##2 layer
        survivalmodels::set_seed(seed)
        fit = survivalmodels::deephit(Surv(time, status) ~ .,
                      data = train_fold,
                      frac = 0.2,
                      learning_rate=params$lr,
                      activation = "relu",
                      num_nodes = c(params$nd1,params$nd2),
                      early_stopping = TRUE,
                      epochs = params$epochs,
                      batch_size = params$batch)
      } else if(params$layer==3){
        ##3 layer
        survivalmodels::set_seed(seed)
        fit = survivalmodels::deephit(Surv(time, status) ~ .,
                      data = train_fold,
                      frac = 0.2,
                      learning_rate=params$lr,
                      activation = "relu",
                      num_nodes = c(params$nd1,params$nd2,params$nd3),
                      early_stopping = TRUE,
                      epochs = params$epochs,
                      batch_size = params$batch)
      } else if(params$layer==4){
        ##4 layer
        survivalmodels::set_seed(seed)
        fit = survivalmodels::deephit(Surv(time, status) ~ .,
                      data = train_fold,
                      frac = 0.2,
                      learning_rate=params$lr,
                      activation = "relu",
                      num_nodes = c(params$nd1,params$nd2,params$nd3,params$nd4),
                      early_stopping = TRUE,
                      epochs = params$epochs,
                      batch_size = params$batch)
      }

      model_name="DeepHit"
      # saving result
      model <- list(
        best_param = list(params),
        model = fit,
        train = cal_metrics(train_fold, fit, model_name),
        valid = cal_metrics(valid_fold, fit, model_name),
        test = cal_metrics(test_fold, fit, model_name)
      )
      return(model)
    }


    # param grid search and  modeling-------------------------------------------
    search_grid <- function(layer) {
      param_grid <- generate_param_grid(layer)
      best_cindex<-0
      best_param<-c()

      for(i in 1:nrow(param_grid)){
        display_progress(index = i, totalN = nrow(param_grid))

        params <- param_grid[i, ]
        model_list<-list()

        for (j in 1:fold){

          test_index <- folds_list[[1]][[j]]
          train_index <- setdiff(seq_len(nrow(d.train)), test_index)

          train_fold <- d.train[train_index, ]
          valid_fold <- d.train[test_index, ]

          model_list[[j]]<-DeepHit_cox(train_fold,valid_fold,d.test,params,seed)

        }

        # calculate c-index
        cindex_list <- extract_metrics(model_list)

        cindex_tv <- mean(cindex_list$valid$cindex,na.rm=T)

        # best parameter
        if (cindex_tv > best_cindex) {
          best_cindex <- cindex_tv
          best_param <- params
        }
      }
      return(list(best_cindex=best_cindex,best_param=best_param))
    }

    final_best_cindex=0
    final_best_param<-c()
    # searching layer
    for (layer in 1:4) {
      print(paste("Layer:", layer))
      output<-search_grid(layer)
      if(output$best_cindex>final_best_cindex){
        final_best_cindex<-output$best_cindex
        final_best_param<-output$best_param
        print(final_best_cindex)
        print(final_best_param)
      }
    }

    print(final_best_cindex)
    print(final_best_param)
    best_param<-final_best_param

    set.seed(seed)
    folds_list<- create_folds(d.train,fold=fold,nrepeats = 2,strata="status",seed=seed)
    folds_list<-c(folds_list[[1]],folds_list[[2]])
    ##
    totalN=fold*2
    model_list<-list()
    for (j in 1:totalN) {
      print(j)
      ## splitting into training and validation data
      test_index<-folds_list[[j]]
      train_index <- setdiff(seq_len(nrow(d.train)), test_index)

      train_fold <- d.train[train_index, ]
      valid_fold <- d.train[test_index, ]

      ## train data
      model_list[[j]]<-DeepHit_cox(train_fold,valid_fold,d.test,best_param,seed)

    }

    ## total data model
    final_model<-DeepHit_cox(d.train,d.test,NULL,best_param,seed)

    ## calculate mean cindex and bs score
    metrics_list<-extract_metrics(model_list)

    save("model_list", file = paste0(outdir,"/",sprintf("%d_%d_DeepHit_result.RData",rep,fold)))
    save("final_model", file = paste0(outdir,"/",sprintf("%d_%d_final_DeepHit_result.RData",rep,fold)))
    save("metrics_list", file = paste0(outdir,"/",sprintf("%d_%d_DeepHit_cindex_result.RData",rep,fold)))

    t2<-Sys.time()
    run_time <- t2 - t1
    print(run_time)

    return(list(final_model=final_model,metrics_list=metrics_list))

  }

  DeepSurv_model<-function(d.train,d.test,fold,rep,outdir,seed){
    t1<-Sys.time()

    ### created nrepeats k-fold
    set.seed(seed)
    folds_list<- create_folds(d.train,fold=fold,nrepeats = 1,
                              strata="status",seed)

    # param grid
    generate_param_grid <- function(layer) {
      param_list <- list(lr = c(0.001, 0.01, 0.1),
                         nd1 = c(8, 16, 32, 64, 128),
                         epochs = c(50, 100, 200),
                         batch = c(32, 64, 128, 256))

      if (layer > 1) param_list$nd2 <- c(8, 16, 32, 64, 128)
      if (layer > 2) param_list$nd3 <- c(8, 16, 32, 64, 128)
      if (layer > 3) param_list$nd4 <- c(8, 16, 32, 64, 128)

      param_list$layer=layer

      expand.grid(param_list)
    }

    # the main function-------------------------------------------
    DeepSurv_cox<-function(train_fold,valid_fold,test_fold,params,seed){
      ###
      if(params$layer==1){
        ##1 layer
        survivalmodels::set_seed(seed)
        fit =survivalmodels::deepsurv(Surv(time, status) ~ .,
                      data = train_fold,
                      frac = 0.2,
                      learning_rate=params$lr,
                      activation = "relu",
                      num_nodes = params$nd1,
                      early_stopping = TRUE,
                      epochs = params$epochs,
                      batch_size = params$batch)
      } else if(params$layer==2){
        ##2 layer
        survivalmodels::set_seed(seed)
        fit = survivalmodels::deepsurv(Surv(time, status) ~ .,
                      data = train_fold,
                      frac = 0.2,
                      learning_rate=params$lr,
                      activation = "relu",
                      num_nodes = c(params$nd1,params$nd2),
                      early_stopping = TRUE,
                      epochs = params$epochs,
                      batch_size = params$batch)
      } else if(params$layer==3){
        ##3 layer
        survivalmodels::set_seed(seed)
        fit = survivalmodels::deepsurv(Surv(time, status) ~ .,
                      data = train_fold,
                      frac = 0.2,
                      learning_rate=params$lr,
                      activation = "relu",
                      num_nodes = c(params$nd1,params$nd2,params$nd3),
                      early_stopping = TRUE,
                      epochs = params$epochs,
                      batch_size = params$batch)
      } else if(params$layer==4){
        ##4 layer
        survivalmodels::set_seed(seed)
        fit = survivalmodels::deepsurv(Surv(time, status) ~ .,
                      data = train_fold,
                      frac = 0.2,
                      learning_rate=params$lr,
                      activation = "relu",
                      num_nodes = c(params$nd1,params$nd2,params$nd3,params$nd4),
                      early_stopping = TRUE,
                      epochs = params$epochs,
                      batch_size = params$batch)
      }

      model_name="DeepSurv"
      # saving result
      model <- list(
        best_param = list(params),
        model = fit,
        train = cal_metrics(train_fold, fit, model_name),
        valid = cal_metrics(valid_fold, fit, model_name),
        test = cal_metrics(test_fold, fit, model_name)
      )
      return(model)
    }

    # param grid search and modeling-------------------------------
    search_grid <- function(layer) {
      param_grid <- generate_param_grid(layer)
      best_cindex<-0
      best_param<-c()

      for(i in 1:nrow(param_grid)){
        display_progress(index = i, totalN = nrow(param_grid))
        params <- param_grid[i, ]
        model_list<-list()

        for(j in 1:fold){
          test_index <- folds_list[[1]][[j]]
          train_index <- setdiff(seq_len(nrow(d.train)), test_index)

          train_fold <- d.train[train_index, ]
          valid_fold <- d.train[test_index, ]

          model_list[[j]]<-DeepSurv_cox(train_fold,valid_fold,d.test,params,seed)
        }

        # Calculate cindex and bs score
        index_list <- extract_metrics(model_list)

        # Return results for this parameter set
        cindex_tv = mean(index_list$valid$cindex,na.rm=T)
        if(cindex_tv > best_cindex){
          best_cindex<-cindex_tv
          best_param<-params
        }
      }
      return(list(best_cindex=best_cindex,best_param=best_param))
    }

    final_best_cindex=0
    final_best_param<-c()
    # search layer
    for (layer in 1:4) {
      print(paste("Layer:", layer))
      output<-search_grid(layer)

      if(output$best_cindex>final_best_cindex){
        final_best_cindex<-output$best_cindex
        final_best_param<-output$best_param
        print(final_best_cindex)
        print(final_best_param)
      }
    }

    print(final_best_cindex)
    print(final_best_param)
    best_param<-final_best_param
    set.seed(seed)
    folds_list<- create_folds(d.train,fold=fold,nrepeats = 2,strata="status",seed=seed)
    folds_list<-c(folds_list[[1]],folds_list[[2]])
    ##
    totalN=fold*2
    model_list<-list()
    for (j in 1:totalN) {
      print(j)
      ## training and testing data
      test_index<-folds_list[[j]]
      train_index <- setdiff(seq_len(nrow(d.train)), test_index)

      train_fold <- d.train[train_index, ]
      valid_fold <- d.train[test_index, ]

      ## train data
      model_list[[j]]<-DeepSurv_cox(train_fold,valid_fold,d.test,best_param,seed)

    }

    ###total data model
    final_model<-DeepSurv_cox(d.train,d.test,NULL,best_param,seed)

    ###calculate mean cindex and bs score
    metrics_list<-extract_metrics(model_list)

    save("model_list", file = paste0(outdir,"/",sprintf("%d_%d_DeepSurv_result.RData",rep,fold)))
    save("final_model", file = paste0(outdir,"/",sprintf("%d_%d_final_DeepSurv_result.RData",rep,fold)))
    save("metrics_list", file = paste0(outdir,"/",sprintf("%d_%d_DeepSurv_cindex_result.RData",rep,fold)))

    t2<-Sys.time()
    run_time <- t2 - t1
    print(run_time)

    return(list(final_model=final_model,metrics_list=metrics_list))

  }

  gbm_model<-function(d.train,d.test,fold,rep,outdir,seed,ncore){

    # Record start time
    start_time <- Sys.time()

    # Set up parallel backend using `future`
    if (.Platform$OS.type == "unix") {
      # Unix: Mac/Linux
      future::plan(future::multicore, workers = ncore)
    } else {
      # Windows
      future::plan(future::multisession, workers = ncore)
    }

    ## create folds
    folds_list<- create_folds(d.train,fold=fold,nrepeats = 1,strata="status",seed)

    # Save results in a list for efficiency
    hyper_grid <- expand.grid(
      learning_rate = c(0.3, 0.1, 0.05, 0.01, 0.005),
      RMSE=NA,
      cindex = NA,
      bs=NA
    )

    # execute grid search in learning rate
    for(i in seq_len(nrow(hyper_grid))) {
      results<-data.frame()
      lr=hyper_grid$learning_rate[i]
      result<-future.apply::future_lapply(
        1:fold,function(j,d.train,d.test,fold,lr,seed,ncore){

          test_index <- folds_list[[1]][[j]]
          train_index <- setdiff(seq_len(nrow(d.train)), test_index)

          train_fold=d.train[train_index, ]
          valid_fold=d.train[test_index, ]

          # fit gbm
          set.seed(seed)  # for reproducibility
          fit <- gbm::gbm(
            formula =Surv(time, status) ~ .,
            data = train_fold,
            distribution = "coxph",
            n.trees = 5000,
            interaction.depth = 3,
            shrinkage =lr,
            n.minobsinnode = 10,
            cv.folds = fold,
            keep.data = FALSE,
            verbose = FALSE,
            n.cores = ncore)

          best.iter <- which.min(fit$cv.error)

          valid = cal_metrics(valid_fold, fit,"GBM")

          list(
            trees=  best.iter,
            RMSE = sqrt(min(fit$cv.error)),
            cindex = valid$cindex,
            bs = valid$bs)

        },
        future.seed = seed,
        d.train = d.train,
        d.test = d.test,
        fold = fold,
        lr=lr,
        seed = seed,
        ncore=ncore)

      results<-do.call(rbind,lapply(result,function(x){
        data.frame(lr=hyper_grid$learning_rate[i],
                   trees=x$trees,
                   RMSE=x$RMSE,
                   cindex=x$cindex,
                   bs=x$bs)
      }))
      # add SSE, trees, and training time to results
      hyper_grid$RMSE[i]  <- mean(results$RMSE)
      hyper_grid$bs[i]  <- mean(results$bs)
      hyper_grid$cindex[i]<-mean(results$cindex)
    }

    best_lr<-hyper_grid[which.max(hyper_grid$cindex),]$learning_rate
    print(paste0("The best learning rate is ",best_lr))

    # search grid in depth and n.minobsinnode
    hyper_grid <- expand.grid(
      j=seq_len(fold),
      n.trees = 5000,
      shrinkage = best_lr,
      interaction.depth = c(1,3, 5, 7),
      n.minobsinnode = c(5, 10, 15)
    )

    ## the modeling function----------------------------------------------
    model_fit<-function(j,n.trees,shrinkage,interaction.depth,n.minobsinnode){

      # split raw data into training and validation data
      test_index <- folds_list[[1]][[j]]
      train_index <- setdiff(seq_len(nrow(d.train)), test_index)
      train_fold <- d.train[train_index, ]
      valid_fold <- d.train[test_index, ]

      # gbm modeling
      set.seed(seed) # for reproducibility
      fit <-gbm::gbm(Surv(time, status) ~ .,
                     data = train_fold,
                     distribution = "coxph",
                     n.trees = n.trees,
                     interaction.depth = interaction.depth,
                     shrinkage = shrinkage,
                     n.minobsinnode = n.minobsinnode,
                     cv.folds = fold,
                     keep.data = FALSE,
                     verbose = FALSE,
                     n.cores = ncore)
      best.iter <- which.min(fit$cv.error)
      valid = cal_metrics(valid_fold,fit,"GBM")

      # result
      cindex = valid$cindex
    }

    # perform search grid with functional programming
    hyper_grid$cindex <- purrr::pmap_dbl(
      hyper_grid,
      ~ model_fit(
        j= ..1,
        n.trees = ..2,
        shrinkage = ..3,
        interaction.depth = ..4,
        n.minobsinnode = ..5
      )
    )

    result <- hyper_grid %>%
      dplyr::group_by(.data$interaction.depth, .data$n.minobsinnode) %>%
      dplyr::summarise(mean_cindex = mean(.data$cindex), .groups = 'drop') %>%
      dplyr::arrange(desc(.data$mean_cindex))

    best_param <-as.data.frame(result[which.max(result$mean_cindex), ])

    # best parameter
    message("\nOptimal Parameters Identified:")
    print(best_param)

    # performance evaluation of best params in repetation
    set.seed(seed)
    folds_list<- create_folds(d.train,fold=fold,nrepeats = 2,strata="status",seed=seed)
    folds_list<-c(folds_list[[1]],folds_list[[2]])
    totalN=fold*2
    ## saving result
    model_list<-list()

    ## the modeling function----------------------------------------------
    gbm_cox<-function(train_fold,valid_fold,test_fold,fold,params,seed,ncore){
      set.seed(seed)
      fit <-gbm::gbm(Surv(time, status) ~ .,
                     data = train_fold,
                     distribution = "coxph",
                     n.trees = 5000,
                     interaction.depth = params$interaction.depth,
                     shrinkage = best_lr,
                     n.minobsinnode = params$n.minobsinnode,
                     cv.folds = fold,
                     keep.data = FALSE,
                     verbose = FALSE,
                     n.cores = ncore)
      best.iter <- which.min(fit$cv.error)

      # result
      model <- list(
        best_param =  list(interaction.depth=params$interaction.depth,
                           ntrees=best.iter,
                           shrinkage = best_lr,
                           n.minobsinnode = params$n.minobsinnode),
        model = fit,
        train = cal_metrics(train_fold, fit, "GBM"),
        valid = cal_metrics(valid_fold, fit,"GBM"),
        test = cal_metrics(test_fold, fit, "GBM")
      )
      return(model)
    }
    for (j in 1:totalN) {
      print(j)
      test_index<-folds_list[[j]]
      train_index <- setdiff(seq_len(nrow(d.train)), test_index)

      train_fold <- d.train[train_index, ]
      valid_fold <- d.train[test_index, ]

      ## train data
      model_list[[j]]<-gbm_cox(train_fold,valid_fold,d.test,fold,best_param,seed=seed,ncore)

    }

    ###total data model
    final_model<-gbm_cox(d.train,d.test,NULL,fold,best_param,seed,ncore)

    metrics_list<-extract_metrics(model_list)

    end_time<-Sys.time()
    run_time<-end_time-start_time
    print(run_time)

    save("model_list", file = paste0(outdir,"/",sprintf("%d_%d_GBM_result.RData",rep,fold)))
    save("final_model", file = paste0(outdir,"/",sprintf("%d_%d_final_GBM_result.RData",rep,fold)))
    save("metrics_list", file = paste0(outdir,"/",sprintf("%d_%d_GBM_cindex_result.RData",rep,fold)))
    return(list(final_model=final_model,metrics_list=metrics_list))
  }

  glm_model<-function(d.train,d.test,fold,rep,outdir,seed,ncore){

    # Set up parallel backend using `future`
    if (.Platform$OS.type == "unix") {
      # Unix: Mac/Linux
      future::plan(future::multicore, workers = ncore)
    } else {
      # Windows
      future::plan(future::multisession, workers = ncore)
    }

    ## define param grid
    param_grid <- expand.grid(nu=seq(0.1,1,0.1))

    ## modeling
    glm_cox<-function(train_fold,valid_fold,test_fold,params,seed){
      set.seed(seed)
      fit <-mboost::glmboost(Surv(time, status) ~.,
                      data = train_fold,
                      family = mboost::CoxPH(),
                      control = mboost::boost_control(mstop=2000,
                                              nu=params),
                      center = FALSE)
      cv10f <- mboost::cv(model.weights(fit), type = "kfold")
      cvm <- mboost::cvrisk(fit, folds = cv10f, papply = lapply)
      n=mboost::mstop(cvm)
      final.fit <- mboost::glmboost(Surv(time, status) ~.,
                            data = train_fold,
                            family = mboost::CoxPH(),
                            control=mboost::boost_control(mstop =n,
                                                  nu=params),
                            center = FALSE)

      model_name="GLMBoost"
      # saving result
      model <- list(
        best_param = list(params),
        model = final.fit,
        train = cal_metrics(train_fold, final.fit, model_name),
        valid = cal_metrics(valid_fold, final.fit, model_name),
        test = cal_metrics(test_fold, final.fit, model_name)
      )
      return(model)
    }

    ## dataframe to save results
    results <- data.frame()

    ####
    folds_list<- create_folds(d.train,fold=fold,nrepeats = 1,strata="status",seed)

    # searching the grid
    print("Searching the grid")
    for (i in 1:nrow(param_grid)){
      params<-param_grid[i,]
      display_progress(index = i, totalN = nrow(param_grid))
      model_list<-list()
      for(j in 1:fold){
        tryCatch({
          test_index<-folds_list[[1]][[j]]
          train_index <- setdiff(seq_len(nrow(d.train)), test_index)

          train_fold <- d.train[train_index,]
          valid_fold <- d.train[test_index,]

          ## train data
          model_list[[j]]<-glm_cox(train_fold,valid_fold,d.test,params,seed)
        }, error = function(e) {
          # Catch errors and log them without stopping the loop
          paste("Error in parameter set repeat", j, ":", e$message)
          return(data.frame(NA, NA, NA, NA, NA, NA, NA, NA))  # Return NA for error cases
        })
      }

      index_list<-extract_metrics(model_list)
      ###
      results<-rbind(results,data.frame(
        nu = params,
        cindex_tr = mean(index_list$train$cindex,na.rm=T),
        bs_tr = mean(index_list$train$bs,na.rm=T),
        cindex_tv = mean(index_list$valid$cindex,na.rm=T),
        bs_tv = mean(index_list$valid$bs,na.rm=T)))
    }

    # best param
    best_param <- results[which.max(results$cindex_tv),][1]
    print(best_param)

    set.seed(seed)
    folds_list<- create_folds(d.train,fold=fold,nrepeats = 2,strata="status",seed=seed)
    folds_list<-c(folds_list[[1]],folds_list[[2]])
    totalN=fold*2
    model_list<-list()
    for (j in 1:totalN) {

      test_index<-folds_list[[j]]
      train_index <- setdiff(seq_len(nrow(d.train)), test_index)

      train_fold <- d.train[train_index, ]
      valid_fold <- d.train[test_index, ]

      ## train data
      model_list[[j]]<-glm_cox(train_fold,valid_fold,d.test,best_param$nu,seed)
    }

    ###total data model
    final_model<-glm_cox(d.train,d.test,NULL,best_param$nu,seed)

    metrics_list<-extract_metrics(model_list)

    save("model_list", file = paste0(outdir,"/",sprintf("%d_%d_GLMBoost_result.RData",rep,fold)))
    save("final_model", file = paste0(outdir,"/",sprintf("%d_%d_final_GLMBoost_result.RData",rep,fold)))
    save("metrics_list", file = paste0(outdir,"/",sprintf("%d_%d_GLMBoost_cindex_result.RData",rep,fold)))
    return(list(final_model=final_model,metrics_list=metrics_list))

  }


  ## the main function--------------------------------------------------
  # 1.RSF --------
  message("---1 RSF ---")
  set.seed(seed)
  model_list[[1]]<-rfrsf_model(d.train,d.test,fold=fold,rep=rep,outdir,seed=seed,ncore)
  names(model_list)[[1]]<-"RFRSF"

  # 2.lasso -------
  message("---2 lasso ---")
  set.seed(seed)
  model_list[[2]]<-lasso_model(d.train,d.test,fold,rep,"Lasso",outdir,seed)
  names(model_list)[[2]]<-"Lasso"

  # 3.enet -------
  message("---3 enet ---")
  set.seed(seed)
  model_list[[3]]<-enet_model(d.train,d.test,fold=fold,rep=rep,outdir,seed)
  names(model_list)[[3]]<-"Enet"

  # 4.ridge -------
  message("---4 ridge ---")
  set.seed(seed)
  model_list[[4]]<-lasso_model(d.train,d.test,fold,rep,"Ridge",outdir,seed)
  names(model_list)[[4]]<-"Ridge"


  # 5.cox boost -------
  message("---5 cox boost ---")
  set.seed(seed)
  model_list[[5]]<-cox_boost(d.train,d.test,fold=fold,rep=rep,outdir,seed)
  names(model_list)[[5]]<-"CoxBoost"


  # 6.superpc -------
  message("---6 superpc ---")
  set.seed(seed)
  model_list[[6]]<-superpc_model(d.train,d.test,fold=fold,rep=rep,outdir,seed)
  names(model_list)[[6]]<-"SuperPC"


  # 7.plsr cox -------
  message("---7 plsr cox ---")
  set.seed(seed)
  model_list[[7]]<-plsR_model(d.train,d.test,fold=fold,rep=rep,outdir,seed)
  names(model_list)[[7]]<-"plsRcox"

  # 8.xgboost -------
  message("--8 xgboost ---")
  set.seed(seed)
  model_list[[8]]<-xgboost_model(d.train,d.test,fold=fold,rep=rep,outdir,seed)
  names(model_list)[[8]]<-"XGBoost"

  # 9.gbm -------
  message("---9 GBM ---")
  set.seed(seed)
  model_list[[9]]<-gbm_model(d.train,d.test,fold=fold,rep=rep,outdir,seed,ncore)
  names(model_list)[[9]]<-"GBM"

  # 10.glm -------
  message("---10 Glm ---")
  set.seed(seed)
  model_list[[10]]<-glm_model(d.train,d.test,fold=fold,rep=rep,outdir,seed,ncore)
  names(model_list)[[10]]<-"GLMBoost"

  # 11.black -------
  message("---11 blackboost ---")
  set.seed(seed)
  model_list[[11]]<-black_model(d.train,d.test,fold=fold,rep=rep,outdir,seed)
  names(model_list)[[11]]<-"BlackBoost"

  if (deep==T){
    # 12.deepsurv -------
    message("---12 deepsurv ---")
    set.seed(seed)

    model_list[[12]]<-DeepSurv_model(d.train,d.test,fold,rep,outdir,seed)
    names(model_list)[[12]]<-"DeepSurv"

    # 13.deephit -------
    message("---13 deephit ---")
    set.seed(seed)
    model_list[[13]]<-DeepHit_model(d.train,d.test,fold,rep,outdir,seed)
    names(model_list)[[13]]<-"DeepHit"
  }

  t2<-Sys.time()
  run_time <- t2 - t1
  print(run_time)

  save(model_list,file=paste0(outdir,"/",sprintf("%d_%d_model_list.RData",rep,fold)))
  print("Model training is finished")
  return(model_list)
}
