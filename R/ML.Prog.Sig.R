ML.Prog.Sig = function(train_data, # cohort data used for training, the colnames of which inlcuding ID, OS.time, OS, and the other candidate genes
                           list_train_vali_Data, # a list of the validation data and the training data. The cohort data is the same to the cohort data used for training
                           # 要求队列的测序深度不能太低，将genelist和队列列名取交集之后尽量保证丢失的信息不超过20%
                           candidate_genes = NULL,
                           unicox.filter.for.candi = NULL, # 是否使用unicox 对基因进行筛选
                           unicox_p_cutoff = NULL, # 默认为0.05 unicox 筛选阈值
                           mode = NULL, # all, single, double
                           single_ml = NULL,# c("RSF", "Enet", "StepCox","CoxBoost","plsRcox","superpc","GBM","survivalsvm","Ridge","Lasso")
                           alpha_for_Enet = NULL , # 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9
                           direction_for_stepcox = NULL, #  c("both", "backward", "forward")
                           double_ml1 = NULL,#c('RSF', "StepCox","CoxBoost","Lasso")
                           double_ml2 = NULL,# c("RSF", "Enet", "StepCox","CoxBoost","plsRcox","superpc","GBM","survivalsvm","Ridge","Lasso")
                           nodesize = NULL, # reference 5-10
                           seed = NULL, # 5201314
                       d="m", # days, months, or years for survival time
                           cores_for_parallel = NULL  #cores for gbm
){
  if (is.null(alpha_for_Enet) == T){
    alpha_for_Enet<-0.1 ## default 0.35
  } else {
    alpha_for_Enet<-alpha_for_Enet
  }

  if (is.null(cores_for_parallel) == T){
    cores_for_parallel<-6 ## default 6
  } else {
    cores_for_parallel<-cores_for_parallel
  }

  if (is.null(direction_for_stepcox) == T){
    direction_for_stepcox<-'both' ## default 0.35
  } else {
    direction_for_stepcox<-direction_for_stepcox
  }

  #loading the packages
  if(T) {
    library(Matrix)
    library(survival)
    library(randomForestSRC)
    library(glmnet)
    library(plsRcox)
    library(superpc)
    library(gbm)
    library(CoxBoost)
    library(survivalsvm)
    library(dplyr)
    library(tibble)
    library(BART)
    library(miscTools)
    library(compareC)
    library(ggplot2)
    library(ggsci)
    library(tidyr)
    library(ggbreak)
    library(mixOmics)
    library(data.table)

    Sys.setenv(LANGUAGE = "en") #显示英文报错信息
    options(stringsAsFactors = FALSE) #禁止chr转成factor
  }

  # Checking data feasibility
  # Replace '-' in column names with '.'
  print("Checking data feasibility")
  candidate_genes = gsub('-','.',candidate_genes)
  colnames(train_data) = gsub('-','.',colnames(train_data))
  colnames(train_data)[2]="time"
  colnames(train_data)[3]="status"

  list_train_vali_Data <- lapply(list_train_vali_Data,function(x){
    colnames(x) = gsub('-','.',colnames(x))
    colnames(x)[2]<-"time"
    colnames(x)[3]<-"status"
    return(x)})

  # Matching candidate genes to genes in each cohort
  print("Matching candidate genes to genes in each cohort")
  common_feature = c('ID', 'time', 'status',candidate_genes)
  common_feature = intersect(common_feature,colnames(list_train_vali_Data[[1]]))
  for (i in names(list_train_vali_Data)) {
    common_feature = intersect(common_feature, colnames(list_train_vali_Data[[i]]))
  }

  message(paste0('---the number of the raw candidate genes is ', length(candidate_genes),' ---'))
  message(paste0('---the number of the common feature across all dataset is ', length(common_feature)-3,' ---'))
  write.table(common_feature[4:length(common_feature)],file="1.common_features_across_datasets.csv",sep=",",row.names = F, quote = F)

  # Matching common feature in each cohort
  train_data<-train_data[,common_feature]

  for (i in names(list_train_vali_Data)) {
    list_train_vali_Data[[i]] =list_train_vali_Data[[i]][,common_feature]
  }

  # keep follow up days more than 30 days
  if (d=="m") {
    train_data <- train_data[train_data$time > 1, ] # days more than 30
    list_train_vali_Data<-lapply(list_train_vali_Data,function(x){
      x=x[x$time>1,]
      return(x)
    })
  } else if (d=="y"){
    train_data <- train_data[train_data$time > 0.083, ]
    list_train_vali_Data<-lapply(list_train_vali_Data,function(x){
      x=x[x$time>0.083,]
      return(x)
    })
  } else {
    train_data <- train_data[train_data$time > 30, ]
    list_train_vali_Data<-lapply(list_train_vali_Data,function(x){
      x=x[x$time>30,]
      return(x)
    })
  }
  ##
  nsample<-nrow(train_data)
  print(paste0(nsample," with more than 30 follow-up days for the next step"))
  print("Data preprocessing completed")

  #unicox filter for candidate genes
  if(unicox.filter.for.candi==T){
    #unicox and km selection
    source("R/sigUnicox.R")
    source("R/sigKMcox.R")
    genelist.1 <- SigUnicox(gene_list = common_feature[-c(1:3)], inputSet = train_data, unicox_pcutoff = 0.05)
    genelist.2 <- SigKMcox(gene_list = genelist.1, inputSet = train_data, KM_pcutoff = 0.05)
    common_feature<-c("ID","time","status",genelist.2)
    message(paste0('---the number of the final unicox filtered candidate genes is ', length(genelist.2),' ---'))
    write.table(genelist.2,paste("3.KM_unicox_select_genes.csv",sep = ""),row.names = F, quote = F,sep=",")
    print("----- finish the preprocess of the unicox and km analysis-----")
  } else {
    message(paste0('---the number of the final not unicox filtered candidate genes is ', length(common_feature)-3,' ---'))
  }

  est_dd <- as.data.frame(train_data)[, common_feature[-1]]
  val_dd_list <- lapply(list_train_vali_Data, function(x){x[, common_feature[-1]]})
  pre_var = common_feature[-c(1:3)]

  seed=seed

  returnIDtoRS = function(rs.table.list, rawtableID){

    for (i in names(rs.table.list)) {
      rs.table.list[[i]] $ID = rawtableID[[i]]$ID
      rs.table.list[[i]] = rs.table.list[[i]] %>% dplyr::select('ID', everything())
    }

    return(rs.table.list)
  }

  if(mode == 'all') {

    # return a dataframe  containing the C-index of models with different
    # machine learning algorithms in different queues
    result <- data.frame()
    ml.res = list()
    riskscore = list()

    # 1-1.RSF --------------------------------------------------------------
    message("---1-1 RSF ---")
    set.seed(seed)
    ns_res<-tune.nodesize(Surv(time,status) ~ ., est_dd)
    rf_nodesize<-ns_res$nsize.opt
    fit <- rfsrc(Surv(time,status)~., data = est_dd,
                 ntree = 1000, nodesize = rf_nodesize,  #
                 splitrule = 'logrank',
                 importance = T,
                 proximity = T,
                 forest = T,
                 seed = seed)
    rs <- lapply(val_dd_list, function(x){cbind(x[, 1:2], RS  = predict(fit, newdata = x)$predicted)})
    rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

    cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(time, status) ~ RS, x))$concordance[1])})) %>%
      rownames_to_column('ID')
    cc$Model <- 'RSF'
    result <- rbind(result, cc)
    ml.res[[ 'RSF']] =fit
    riskscore[[ 'RSF']] = rs

    # 1-2.RSF + CoxBoost #########
    message("---1-2.RSF + CoxBoost ---")

    set.seed(seed)
    ns_res<-tune.nodesize(Surv(time,status) ~ ., est_dd)
    rf_nodesize<-ns_res$nsize.opt
    fit <- rfsrc(Surv(time, status)~., data = est_dd,
                   ntree = 1000,
                   nodesize = rf_nodesize,
                   splitrule = 'logrank',
                   importance = T,
                   proximity = T,
                   forest = T,
                   seed = seed)
      rid <- var.select(object = fit, conservative = "high")

      rid <- rid$topvars

      if(length(rid)>1) {

        est_dd2 <- train_data[, c('time', 'status', rid)]
        val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('time', 'status', rid)]})
        set.seed(seed)
        pen <- optimCoxBoostPenalty(est_dd2[, 'time'], est_dd2[, 'status'], as.matrix(est_dd2[, -c(1,2)]),
                                    trace=TRUE, start.penalty = 500, parallel = T)


        cv.res <- cv.CoxBoost(est_dd2[, 'time'], est_dd2[, 'status'], as.matrix(est_dd2[, -c(1, 2)]),
                              maxstepno = 500, K = 10, type = "verweij",  penalty = pen$penalty)
        fit <- CoxBoost(est_dd2[, 'OS.time'], est_dd2[, 'OS'], as.matrix(est_dd2[, -c(1, 2)]),
                        stepno = cv.res$optimal.step, penalty = pen$penalty)
        rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, newdata = x[, -c(1, 2)], newtime = x[, 1],  newstatus = x[, 2], type = "lp")))})
        cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
          rownames_to_column('ID')
        cc$Model <- paste0('RSF + ','CoxBoost')
        result <- rbind(result, cc)
        ml.res[[ paste0('RSF + ','CoxBoost')]] =fit

        rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

        riskscore[[ paste0('RSF + ','CoxBoost')]] = rs

      } else {
        warning('The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2')
      }


      # 1-3.RSF + Enet####

      message("---1-3.RSF + Enet ---")

      set.seed(seed)
      fit <- rfsrc(Surv(OS.time, OS)~., data = est_dd,
                   ntree = 1000, nodesize = rf_nodesize, #
                   splitrule = 'logrank',
                   importance = T,
                   proximity = T,
                   forest = T,
                   seed = seed)
      rid <- var.select(object = fit, conservative = "high")


      rid <- rid$topvars


      if(length(rid)>1) {

        est_dd2 <- train_data[, c('OS.time', 'OS', rid)]
        val_dd_list2 <- lapply(list_train_vali_Data, function(x){x[, c('OS.time', 'OS', rid)]})
        x1 <- as.matrix(est_dd2[, rid])
        x2 <- as.matrix(Surv(est_dd2$OS.time, est_dd2$OS))
        for (alpha in seq(0.1, 0.9, 0.1)) {
          set.seed(seed)
          fit = cv.glmnet(x1, x2, family = "cox", alpha = alpha, nfolds = 10)
          rs <- lapply(val_dd_list2, function(x){cbind(x[, 1:2], RS = as.numeric(predict(fit, type = 'link', newx = as.matrix(x[, -c(1, 2)]), s = fit$lambda.min)))})
          cc <- data.frame(Cindex = sapply(rs, function(x){as.numeric(summary(coxph(Surv(OS.time, OS) ~ RS, x))$concordance[1])})) %>%
            rownames_to_column('ID')
          cc$Model <- paste0('RSF + ', 'Enet', '[α=', alpha, ']')
          result <- rbind(result, cc)
          ml.res[[paste0('RSF + ', 'Enet', '[α=', alpha, ']')]] =fit
          rs =returnIDtoRS(rs.table.list = rs,rawtableID = list_train_vali_Data)

          riskscore[[paste0('RSF + ', 'Enet', '[α=', alpha, ']')]] = rs

        }

      } else {
        warning('The number of seleted candidate gene by RSF, the first machine learning algorithm, is less than 2')
      }
