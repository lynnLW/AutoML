#' Screening out the core variables for the prognosis with the machine learning algorithms
#'
#' A function can be used for screening out the core features from the given candidate genes with eight machine learning algorithms.
#'
#' @param InputMatrix A gene expression dataframe after log2(x+1) scaled. The first three of the column names are, in order, ID,OS.time, OS. Columns starting with the fourth are gene symbols. OS.time is a numeric variable in days. OS is a numeric variable containing 0, 1. 0: Alive, 1: Dead.
#' @param candidate_genes The input genes, that you want to screen out from, for identifying the core features.
#' @param mode  We provide three modes including 'all', 'single', and 'all_without_SVM'. The 'all' mode means using all eight methods for selecting. The 'single' mode means using only one method for running. Since SVM takes so much time, we're singling him out. The 'all_without_SVM' mode means the other seven methods used for selecting.
#' @param seed  The seed. You can set it as any number. For example, 5201314.
#' @param single_ml The one method from the eight methods including "RSF", "Enet", "Boruta", "Xgboost", "SVM-REF", "Lasso", "CoxBoost', 'StepCox'.
#' @param nodesize The node size parameter for the RSF method. The default is 5. You can try another positive integer. For example, 10,15,20, etc.

Corefeature.Prog.Screen <- function(InputMatrix, ### first column ID,second OS.time, third OS, (0/1), feature list...
                                       genelist,
                                       mode = NULL, # all, single,all_without_SVM
                                       seed = NULL,
                                       d="m", ##days, months, or years
                                       single_ml = NULL, # c("RSF", "Enet", "Boruta","Xgboost","SVM-REF","Lasso","CoxBoost','StepCox')
                                       unicox_km=T) {
  ### Screen out the core features via the multiple machine leaning algorithms
  ### loading the packages ####

  if (T) {
    Biocductor_packages <- c(
      "tidyverse",
      "scales",
      "Hmisc",
      "survival",
      "randomForestSRC",
      "glmnet",
      "plsRcox",
      "CoxBoost",
      "survivalsvm",
      "dplyr",
      "tibble",
      "BART",
      "miscTools",
      "compareC",
      "tidyr",
      "mixOmics",
      "data.table",
      "pbapply",
      "e1071",
      "Boruta",
      "xgboost",
      "Ckmeans.1d.dp",
      "Matrix"
    )

    lapply(Biocductor_packages, function(x) {
      library(x,
              character.only = T
              # ,  lib.loc = "/export/bioinfo-team/home/xiongzj/R/x86_64-pc-linux-gnu-library/4.1"
      )
    })
  }

  print("-----format the inputmatrix and genelist-----")
  ## 将genelist和表达矩阵的基因名称格式统一
  genelist <- gsub("-", ".", genelist)
  genelist <- gsub("_", ".", genelist)
  colnames(InputMatrix)[4:ncol(InputMatrix)] <- gsub("-", ".", colnames(InputMatrix)[4:ncol(InputMatrix)])
  colnames(InputMatrix)[4:ncol(InputMatrix)] <- gsub("_", ".", colnames(InputMatrix)[4:ncol(InputMatrix)])
  ##
  names(InputMatrix)[2]<-"OS.time"
  names(InputMatrix)[3]<-"OS"
  ##
  print("Starting the data preprocess")
  ############### data preprocess#######
  # null value replace with zero
  InputMatrix[is.na(InputMatrix)] <- 0
  # table(is.na(inputSet))
  InputMatrix <- InputMatrix %>% as.data.frame()
  InputMatrix$OS.time <- as.numeric(InputMatrix $OS.time)

  # keep follow up days more than 30 days
  if (d=="m") {
    InputMatrix <- InputMatrix[InputMatrix$OS.time > 1, ] # days more than 30
  } else if (d=="y"){
      InputMatrix <- InputMatrix[InputMatrix$OS.time > 0.083, ]
  } else {
      InputMatrix <- InputMatrix[InputMatrix$OS.time > 30, ]
  }
  ##
  nsample<-nrow(InputMatrix)
  print(paste0(nsample," with more than 30 follow-up days for the next step"))

  print("Gets the intersection of genelist and expression profile")
  # 获取genelist和表达谱的交集
  comsa1 <- intersect(colnames(InputMatrix)[4:ncol(InputMatrix)], genelist)
  # write.table(comsa1,"2.intersection_genelist_exprSet_gene.txt", row.names = F, quote = F)

  print("Processing the  input representation matrix")
  # 对输入的表达矩阵进行处理
  InputMatrix <- InputMatrix[, c("ID", "OS.time", "OS", comsa1)]

  InputMatrix[, c(1:2)] <- apply(InputMatrix[, c(1:2)], 2, as.factor)
  InputMatrix[, c(2:ncol(InputMatrix))] <- apply(InputMatrix[, c(2:ncol(InputMatrix))], 2, as.numeric)
  InputMatrix <- as.data.frame(InputMatrix)
  # rownames(inputSet) <- inputSet$ID
  print("Data preprocessing completed")

  if(unicox_km){
    #unicox and km selection
    #source("R/sigUnicox.R")
    #source("R/sigKMcox.R")
    genelist.1 <- SigUnicox(gene_list = genelist, inputSet = InputMatrix, unicox_pcutoff = 0.05)
    genelist.2 <- SigKMcox(gene_list = genelist.1, inputSet = InputMatrix, KM_pcutoff = 0.05)

    candidate_genes <- genelist.2
    write.table(candidate_genes,paste("3.KM_unicox_select_genes.csv",sep = ""),row.names = F, quote = F,sep=",")
    print("----- finish the preprocess of the unicox and km analysis-----")
    ##### setting the pamameters ######
  }

    seed <- seed
    iter.times <- 1000
    # Checking data feasibility
    message("--- check data feasibility ---")

    # Matching candidate genes to genes in each cohort
    common_feature <- c("ID", "OS.time", "OS", candidate_genes)
    common_feature <- intersect(common_feature, colnames(InputMatrix))

    message(paste0("---the number of candidate genes is ", length(candidate_genes), " ---"))

    ######### the main of the function ##########

    if (!is.na(seed) &
        mode %in% c("all", "single", "all_without_SVM") &
        identical(c("ID", "OS.time", "OS"), colnames(InputMatrix)[1:3]) &
        length(candidate_genes) > 0 &
        identical(c("ID", "OS.time", "OS"), common_feature[1:3]) &
        length(common_feature) > 3) {
        message("--- Data preprocessing ---")
      # Data preprocessing

      # Matching candidate genes to genes in each cohort

      InputMatrix <- InputMatrix[, common_feature]
      InputMatrix[, -c(1:3)] <- apply(InputMatrix[, -c(1:3)], 2, as.numeric)
      InputMatrix[, c(1:2)] <- apply(InputMatrix[, c(1:2)], 2, as.factor)
      InputMatrix[, c(2:3)] <- apply(InputMatrix[, c(2:3)], 2, as.numeric)
      InputMatrix <- InputMatrix[!is.na(InputMatrix$OS.time) & !is.na(InputMatrix$OS), ]

      InputMatrix[, -c(1:3)] <- apply(InputMatrix[, -c(1:3)], 2, function(x) {
        x[is.na(x)] <- mean(x, na.rm = T)
        return(x)
      })

      est_dd <- as.data.frame(InputMatrix)[, common_feature[-1]]
      pre_var <- common_feature[-c(1:3)]
      selected.feature <- data.frame()

      if (mode == "all") {

        ### 1. Repeated Lasso  #############
        message("--- 1.Repeated lasso ---")
        #source("R/ML.lasso.R")
        selected.feature<-ML.lasso(est_dd,pre_var,iter.times,seed=123456)

        ##### 2.Enet ###########
        message("--- 2.Enet  ---")
        #source("R/ML.enet.r")
        selected.feature<-ML.enet(est_dd,pre_var,iter.times,seed=123456)

        ##### 3.Boruta ###########
        message("--- 3.Boruta  ---")
        #source("R/ML.boruta.r")
        selected.feature<-ML.boruta(est_dd,seed=123456)

        ##### 4.SVM-REF ##########
        message("--- 4.SVM-REF  ---")
        print("This step will probably take several hours")
        input <- est_dd[, -1]
        #source("R/ML.svm.R")
        # 10CV (k-fold crossValidation）
        svmRFE(input, k = 10, halve.above = 100) # 分割数据，分配随机数
        nfold <- 10
        nrows <- nrow(input)
        folds <- rep(1:nfold, len = nrows)[sample(nrows)]
        folds <- lapply(1:nfold, function(x) which(folds == x))
        results <- lapply(folds, svmRFE.wrap, input, k = 10, halve.above = 100) # 特征选择
        top.features <- WriteFeatures(results, input, save = F) # 查看主要变量
        n.features <- nrow(top.features)
        if (n.features > 300) {
          n.svm <- 300
        } else {
          n.svm <- n.features
        }
        featsweep <- base::lapply(1:n.svm, FeatSweep.wrap, results, input)
        no.info <- min(prop.table(table(input[, 1])))
        errors <- sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))
        fea <- top.features[1:which.min(errors), "FeatureName"]

        result <- data.frame(
          method = c(rep("SVM-REF", length(fea))),
          selected.fea = fea
        )
        write.table(result,file="4.svm_select_features.csv",sep=",",row.names = F)
        selected.feature <- rbind(selected.feature, result)

        ##### 5.xgboost ##########
        message("--- 5.xgboost  ---")
        #source("R/ML.xgboost.R")
        selected.feature<-ML.xgboost(est_dd,seed=123456)

        ##### 6.rsf ##########
        message("--- 6.rsf  ---")
        #source("R/ML.rsf.R")
        selected.feature<-ML.rsf(est_dd,seed=123456)

        ##### 7.coxboost ##########
        message("--- 7.coxboost  ---")
        #source("R/ML.coxboost.R")
        selected.feature<-ML.coxboost(est_dd,seed=123456)

        ##### 8.xgboost ##########
        message("--- 8.stepcox  ---")
        #source("R/ML.stepcox.R")
        selected.feature<-ML.stepCox(est_dd,seed=123456)

        return(selected.feature)
      }
    }

}
