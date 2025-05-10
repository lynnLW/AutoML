#' Screening out the core prognostic genes with the machine learning algorithms
#'
#' A function can be used for screening out the core features from the given candidate genes with 9 machine learning algorithms.
#'
#' @param InputMatrix ### 1: ID, 2: OS_time, 3: OS_status (0/1), feature list...
#' @param genelist The candidate genes
#' @param fold default 5
#' @param filter_OS_time default F If keep follow up time > 30 days
#' @param meta_time When filter_OS_time=T, set follow up time months(m) or days(d) or years(y)
#' @param unicox_km unicox and km filtering
#' @param deg differential expression filtering
#' @param up when up=T, HR>1 in unicox analysis and log2FC>0 in differential analysis
#' @param method  "all" using 9 ML algorithms or boruta or xgboost or RSF or rfe or lasso
#' @param svm_method when method="all", setting if using svm algorithms, which is time-consuming
#' @param outdir the outpt directory
#' @param seed  The seed
#' @param ncore multisession workers
#' @return feature list
#' @export
#' @examples
#' \donttest{
#' # Requires trained data and genelist
#' data(train_data)
#' selected.feature<-feature_selection(InputMatrix,genelist=genelist)
#' }
feature_selection<- function(InputMatrix,
                              genelist=NULL,
                              fold=5,
                              filter_OS_time=F,
                              meta_time="m",
                              unicox_km=T,
                              deg=T,
                              up=F,
                              method="all",
                              svm_method=F,
                              outdir='1.feature_select/',
                              seed=123,
                              ncore=4
                              ){
  ### loading the packages ####

  if (T) {
    Biocductor_packages <- c(
      "tidyverse",
      "scales",
      "future.apply",
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
      "caret",
      "xgboost",
      "Ckmeans.1d.dp",
      "Matrix",
      "YSX"
    )

    lapply(Biocductor_packages, function(x) {
      library(x,character.only = T)
    })
  }

  ## setting the parameters
  if (T) {
    if(is.null(genelist)==T){
      genelist<-colnames(InputMatrix)[-c(1:3)]
    } else {
      genelist<-genelist
    }

    if(!dir.exists(outdir)){
      dir.create(outdir,recursive = T)
    }
  }

  ## format the inputmatrix and genelist
  print("-----format the inputmatrix and genelist-----")
  genelist <- gsub("-", ".", genelist)
  genelist <- gsub("_", ".", genelist)
  colnames(InputMatrix)[4:ncol(InputMatrix)] <- gsub("-", ".", colnames(InputMatrix)[4:ncol(InputMatrix)])
  colnames(InputMatrix)[4:ncol(InputMatrix)] <- gsub("_", ".", colnames(InputMatrix)[4:ncol(InputMatrix)])
  names(InputMatrix)[1]<-"ID"
  names(InputMatrix)[2]<-"OS_time"
  names(InputMatrix)[3]<-"OS_status"

  ## preprocessing
  print("Starting the data preprocessing")
  ## null value replace with zero
  InputMatrix[,is.na(InputMatrix[,4:ncol(InputMatrix)])] <- 0
  InputMatrix <- InputMatrix %>% as.data.frame()
  InputMatrix$OS_time <- as.numeric(InputMatrix$OS_time)

  ## If keep follow up days more than 30 days
  if (filter_OS_time){
    if (meta_time=="m") {
      InputMatrix <- InputMatrix[InputMatrix$OS_time >= 1, ] # days more than 30
      nsample<-nrow(InputMatrix)
      print(paste0(nsample," with more than 30 follow-up days for the next step"))
    } else if (meta_time=="y"){
      InputMatrix <- InputMatrix[InputMatrix$OS_time >= 0.083, ]
      nsample<-nrow(InputMatrix)
      print(paste0(nsample," with more than 30 follow-up days for the next step"))
    } else {
      InputMatrix <- InputMatrix[InputMatrix$OS_time >= 30, ]
      nsample<-nrow(InputMatrix)
      print(paste0(nsample," with more than 30 follow-up days for the next step"))
    }
  } else {
    InputMatrix <- InputMatrix # days more than 30
  }

  ## get the common genes in the training cohort
  print("Gets the intersection of genelist and inputmatrix")
  comsa1 <- intersect(colnames(InputMatrix)[4:ncol(InputMatrix)], genelist)
  print(paste0("Gets ",length(comsa1)," genes"))
  InputMatrix <- InputMatrix[, c("ID", "OS_time", "OS_status", comsa1)]
  print("Data preprocessing completed")

  ## record the start time
  t1<-Sys.time()
  ## differential expression filtering
  if (deg==T){
    print("Processing differential expression selection")
    ##########function
    run_diff_gene<-function(input,genes,diff_pcutoff = 0.05,outdir){
      ###diff expression
      df_mean<- input %>%
        group_by(OS_status)%>%
        summarise(across(everything(),\(x) mean(x, na.rm = TRUE)))%>%
        t() %>% as.data.frame()
      names(df_mean)<-df_mean[1,]
      df_mean<-df_mean[-1,]
      ###function calculate pvalue
      cal_pvalue<-function(input,gene){
        pvalue<-pairwise.wilcox.test(input[,gene],
                                     input[,'OS_status'],
                                     p.adjust.method = "bonf")$p.value
        return(pvalue)
      }
      ###cal pvalue
      pvalue<-c()
      for(i in 1:length(genes)){
        gene=genes[i]
        p<-cal_pvalue(input,gene)
        pvalue<-append(pvalue,p)
      }
      ###
      df_mean$pvalue<-pvalue
      df_mean$Yes <- as.numeric(df_mean$Yes)
      df_mean$No <- as.numeric(df_mean$No)
      df_mean$fc <- df_mean[['Yes']] / df_mean[['No']]
      df_mean<-df_mean[df_mean$pvalue<diff_pcutoff,]

      return(df_mean)
    }

    #########
    expr=as.data.frame(InputMatrix[,-c(1,2)])
    expr[,1]<-ifelse(expr[,1]==1,"Yes","No")
    expr[,1]<-as.factor(expr[,1])
    diff_gene <- run_diff_gene(input = expr,
                               genes = colnames(expr)[-1],
                               diff_pcutoff = 0.05,
                               outdir = outdir)
    ########
    diff_gene$log2fc<-log2((diff_gene$Yes+1)/(diff_gene$No+1))
    if (up==T){
      diff_gene<-diff_gene[diff_gene$log2fc>0,]
      print(paste0("Gets ",length(row.names(diff_gene))," differential genes with a pvalue <0.05 and log2Fc>0"))
      write.csv(diff_gene,file=paste0(outdir,"/1.diff_gene_0.05.csv"))
      geneset.1<-row.names(diff_gene)
      candidate_genes <- geneset.1
    } else {
      diff_gene<-diff_gene
      print(paste0("Gets ",length(row.names(diff_gene))," differential genes with a pvalue <0.05"))
      write.csv(diff_gene,file=paste0(outdir,"/1.diff_gene_0.05.csv"))
      geneset.1<-row.names(diff_gene)
      candidate_genes <- geneset.1
    }

  } else {
    candidate_genes <- genelist
  }

  ## unicox and km analysis filtering
  if(unicox_km==T){
    print("Processing unicox and km selection")

    ## unicox function
    SigUnicox <- function(gene_list,
                          inputSet,
                          unicox_pcutoff=0.05,
                          outdir){
      ##display progress
      display.progress <- function(index, totalN, breakN = 20) {
        if (index %% ceiling(totalN / breakN) == 0) {
          cat(paste(round(index * 100 / totalN), "% ", sep = ""))
        }
      }

      ######## unicox#######
      print("Stating the univariable cox regression")
      unicox <- data.frame()
      for (i in 1:ncol(inputSet[, 4:ncol(inputSet)])) {
        display.progress(index = i, totalN = ncol(inputSet[, 4:ncol(inputSet)]))
        gene <- colnames(inputSet[, 4:ncol(inputSet)])[i]
        tmp <- data.frame(
          expr = as.numeric(inputSet[, 4:ncol(inputSet)][, i]),
          futime = inputSet$OS_time,
          fustat = inputSet$OS_status,
          stringsAsFactors = F
        )
        cox <- coxph(Surv(futime, fustat) ~ expr, data = tmp)
        coxSummary <- summary(cox)
        unicox <- rbind.data.frame(unicox,
                                   data.frame(
                                     gene = gene,
                                     HR = as.numeric(coxSummary$coefficients[, "exp(coef)"])[1],
                                     z = as.numeric(coxSummary$coefficients[, "z"])[1],
                                     pvalue = as.numeric(coxSummary$coefficients[, "Pr(>|z|)"])[1],
                                     lower = as.numeric(coxSummary$conf.int[, 3][1]),
                                     upper = as.numeric(coxSummary$conf.int[, 4][1]),
                                     stringsAsFactors = F
                                   ),
                                   stringsAsFactors = F
        )
      }

      sunicox <- unicox[which(unicox$pvalue < unicox_pcutoff),]

      if(up==T){
        sunicox<-sunicox[sunicox$HR>1,]
      } else{
        sunicox<-sunicox
      }
      write.table(sunicox,file=paste0(outdir,"/2.unicox_",unicox_pcutoff,"_result.csv"),row.names = F, quote = F,sep=",")

      print("Finished the univariable cox regression")
      selgene <- sunicox$gene
      return(selgene)
    }

    geneset.2 <- SigUnicox(gene_list = genelist, inputSet = InputMatrix,unicox_pcutoff = 0.05,outdir = outdir)
    print(paste0("Gets ",length(geneset.2)," unicox gene with a pvalue <0.05"))

    ## km selection
    SigKMcox <- function(gene_list,
                         inputSet,
                         KM_pcutoff=0.05,
                         outdir
    ) {
      display.progress <- function(index, totalN, breakN = 20) {
        if (index %% ceiling(totalN / breakN) == 0) {
          cat(paste(round(index * 100 / totalN), "% ", sep = ""))
        }
      }
      ############### km selection#######
      print("Stating the KM survival")

      kmoutput <- NULL
      for (i in 1:ncol(inputSet[, 4:ncol(inputSet)])) {
        display.progress(index = i, totalN = ncol(inputSet), breakN = 20)
        g <- colnames(inputSet[, 4:ncol(inputSet)])[i]
        tmp <- inputSet[, c("OS_time", "OS_status", g)]
        tmp$group <- ifelse(tmp[, 3] > median(tmp[, 3]), "High", "Low")
        fitd <- survdiff(Surv(OS_time, OS_status) ~ group, data = tmp, na.action = na.exclude)
        p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
        kmoutput <- rbind(kmoutput, data.frame(
          gene = g,
          pvalue = p.val,
          stringsAsFactors = F
        ))
      }

      print("Finished the KM selection")
      skmoutput <- kmoutput[which(kmoutput$pvalue < KM_pcutoff),]
      write.table(skmoutput,paste0(outdir,"/3.KM_",KM_pcutoff,"_result.csv"),row.names = F, quote = F,sep=",")
      selgene <- kmoutput[which(kmoutput$pvalue < KM_pcutoff), "gene"]
      return(selgene)
    }
    geneset.3 <- SigKMcox(gene_list = genelist, inputSet = InputMatrix,KM_pcutoff = 0.05,outdir = outdir)
    print(paste0("Gets ",length(geneset.2)," KM gene with a pvalue <0.05"))

    ####
    print("Performing the intersection between km and unicox method")
    geneset.op<-intersect(geneset.2,geneset.3)
    print(paste0("Gets ",length(geneset.op)," unicox and km gene with a pvalue <0.05"))
    write.table(geneset.op,paste0(outdir,"/4.op_KM_unicox_gene.csv"),row.names = F, quote = F,sep=",")

    ####
    if(deg==T){
      print("Performing the intersection between km_unicox and diff method")
      candidate_genes<-intersect(geneset.op,candidate_genes)
      print(paste0("Gets ",length(candidate_genes)," km_unicox and diff gene with a pvalue <0.05"))
      write.table(candidate_genes,paste0(outdir,"/5.op_unicox_KM_diff_gene.csv"),row.names = F, quote = F,sep=",")
    } else {
      candidate_genes<-geneset.op
    }
    print("----- finish the preprocess of the unicox and km analysis-----")
  } else {
    candidate_genes<-candidate_genes
  }

  # Matching candidate genes to genes in each cohort
  common_feature <- c("ID", "OS_time", "OS_status", candidate_genes)
  common_feature <- intersect(common_feature, colnames(InputMatrix))
  message(paste0("---the number of candidate genes is ", length(common_feature)-3, " ---"))

  ######### function of each method ##########

  ##1.Lasso
  FS.lasso<-function(est_dd,
                     pre_var,
                     iter.times,
                     fold,
                     selected.feature,
                     outdir,
                     seed=seed){
    ### 1. Repeated Lasso  #############
    message("--- 1.Repeated lasso ---")
    x1 <- as.matrix(est_dd[, pre_var])
    x2 <- as.matrix(Surv(est_dd$OS_time, est_dd$OS_status))
    print("1000 time lasso penalty")
    # 1000 time lasso penalty
    lasso_fea_list <- list()
    list.of.seed <- 1:iter.times
    set.seed(seed)
    print("This step will probably take several minutes")

    lasso_fea_list <- future.apply::future_lapply(list.of.seed, function(x) {
      set.seed(list.of.seed[x])
      cvfit <- cv.glmnet(
        x = x1,
        y = x2,
        nfolds = fold,
        alpha = 1,
        family = "cox",
        maxit = 1000)

      # optimal lambda
      fea <- rownames(coef(cvfit, s = "lambda.min"))[coef(cvfit, s = "lambda.min")[, 1] != 0]
      if (is.element("(Intercept)", fea)) {
        lasso_fea <- sort(fea[-1])
      } else {
        lasso_fea <- sort(fea)
      }
      return(lasso_fea)
    },future.seed=TRUE)

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

    genes <- sort(table(unlist(lasso_fea_list)), decreasing = T)
    freq.cutoff <- iter.times * 0.05
    genes <- names(genes[genes > freq.cutoff])


    result <- data.frame(
      method = c(rep("Lasso", length(genes))),
      selected.fea = genes
    )
    write.table(result,file=paste0(outdir,"/1.lasso_select_features.csv"),sep=",",row.names = F)
    selected.feature <- rbind(selected.feature,result)
    return(selected.feature)
  }

  ##2.Enet
  FS.enet<-function(est_dd,
                    pre_var,
                    iter.times,
                    fold,
                    selected.feature,
                    outdir,
                    seed=seed){
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
    x2 <- as.matrix(Surv(est_dd$OS_time, est_dd$OS_status))

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
      res_list  <- future.apply::future_lapply(list.of.seed, function(x) {
        set.seed(x)
        cvfit <- cv.glmnet(
          x = x1,
          y = x2,
          nfolds = 10,
          alpha = alpha,
          family = "cox",
          maxit = 1000)
        mean_cv_error <- mean(cvfit$cvm)
        return(list(cv_error = mean_cv_error))
      },future.seed=TRUE)

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
    jpeg(filename=paste0(outdir,"/2.enet_mean_cross_error_different_alpha.jpg"),
         res=600,width = 12,height = 10,units = "cm")
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
        legend.position = "none"
      ) +
      scale_color_brewer(palette = "Set1")
    print(p)
    dev.off()

    ##select optimal alpha
    mean_cv_errors_sum_df <- do.call(rbind, lapply(1:length(alpha_values), function(i) {
      data.frame(
        alpha = alpha_values[i],
        mean_cv_error = mean(mean_cv_errors_list[[i]])
      )}))
    best_alpha_index <- which.min(mean_cv_errors_sum_df$mean_cv_error)
    best_alpha <- alpha_values[best_alpha_index]

    print(paste("Best alpha is:", best_alpha))

    ### best alpha
    alpha <- best_alpha
    message(paste0("--- 2.Enet ", alpha, "---"))

    set.seed(seed)
    # 1000 time
    list.of.seed <- 1:final.iter.times
    res_list <- future.apply::future_lapply(list.of.seed, function(x) {
      set.seed(x)
      cvfit <- cv.glmnet(
        x = x1,
        y = x2,
        nfolds = fold,
        alpha = alpha,
        family = "cox",
        maxit = 1000)
      fea <- rownames(coef(cvfit, s = "lambda.min"))[coef(cvfit, s = "lambda.min")[, 1] != 0]
      if (is.element("(Intercept)", fea)) {
        lasso_fea <- sort(fea[-1])
      } else {
        lasso_fea <- sort(fea)
      }
      return(list(fea = lasso_fea))
    },future.seed = TRUE)

    fea_list<-lapply(list.of.seed, function(x) res_list[[x]]$fea)

    ####
    genes <- sort(table(unlist(fea_list)), decreasing = T)
    freq.cutoff <- iter.times * 0.05
    genes <- names(genes[genes > freq.cutoff])

    result <- data.frame(
      method = c(rep(paste0("Enet","[alpha=",alpha,"]"), length(genes))),
      selected.fea = genes
    )
    write.table(result,file=paste0(outdir,"/2.Enet_select_features.csv"),sep=",",row.names = F)
    selected.feature <- rbind(selected.feature, result)
    return(selected.feature)
  }

  ##3.Ridge
  FS.ridge<-function(est_dd,
                     pre_var,
                     iter.times,
                     fold,
                     selected.feature,
                     outdir,
                     seed=seed){
    ### 1. Repeated Ridge  #############
    message("--- 3.Repeated Ridge ---")
    x1 <- as.matrix(est_dd[, pre_var])
    x2 <- as.matrix(Surv(est_dd$OS_time, est_dd$OS_status))
    print("1000 time ridge penalty")
    # 1000 time lasso penalty
    lasso_fea_list <- list()
    list.of.seed <- 1:iter.times
    set.seed(seed)
    print("This step will probably take several minutes")

    lasso_fea_list <- future.apply::future_lapply(list.of.seed, function(x) {
      set.seed(list.of.seed[x])
      cvfit <- cv.glmnet(
        x = x1,
        y = x2,
        nfolds = fold,
        alpha = 0,
        family = "cox",
        maxit = 1000)

      # optimal lambda
      fea <- rownames(coef(cvfit, s = "lambda.min"))[coef(cvfit, s = "lambda.min")[, 1] != 0]
      if (is.element("(Intercept)", fea)) {
        lasso_fea <- sort(fea[-1])
      } else {
        lasso_fea <- sort(fea)
      }
      return(lasso_fea)
    },future.seed = TRUE)

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

    genes <- sort(table(unlist(lasso_fea_list)), decreasing = T)
    freq.cutoff <- iter.times * 0.05
    genes <- names(genes[genes > freq.cutoff])


    result <- data.frame(
      method = c(rep("Ridge", length(genes))),
      selected.fea = genes
    )
    write.table(result,file=paste0(outdir,"/3.Ridge_select_features.csv"),sep=",",row.names = F)
    selected.feature <- rbind(selected.feature,result)
    return(selected.feature)
  }

  ##4.Boruta
  boruta.imp <- function(x) {
    imp <- reshape2::melt(x$ImpHistory, na.rm = T)[, -1]
    colnames(imp) <- c("Variable", "Importance")
    imp <- imp[is.finite(imp$Importance), ]

    variableGrp <- data.frame(
      Variable = names(x$finalDecision),
      finalDecision = x$finalDecision
    )

    showGrp <- data.frame(
      Variable = c("shadowMax", "shadowMean", "shadowMin"),
      finalDecision = c("shadowMax", "shadowMean", "shadowMin")
    )

    variableGrp <- rbind(variableGrp, showGrp)

    boruta.variable.imp <- merge(imp, variableGrp, all.x = T)

    sortedVariable <- boruta.variable.imp %>%
      group_by(Variable) %>%
      summarise(median = median(Importance)) %>%
      arrange(median)
    sortedVariable <- as.vector(sortedVariable$Variable)


    boruta.variable.imp$Variable <- factor(boruta.variable.imp$Variable, levels = sortedVariable)

    invisible(boruta.variable.imp)
  }

  FS.boruta<-function(est_dd,
                      selected.feature,
                      outdir,
                      seed=seed){
    ##### 4.Boruta ###########
    set.seed(seed)
    message("--- 4.Boruta  ---")
    boruta <- Boruta(
      x = as.matrix(est_dd[,-c(1,2)]),
      y = as.factor(est_dd[, c(2)]),
      pValue = 0.01,
      mcAdj = T,
      maxRuns = 1000
    )

    # head(boruta.variable.imp)
    boruta.finalVars <- data.frame(Item = getSelectedAttributes(boruta, withTentative = T), Type = "Boruta")

    result <- data.frame(
      method = c(rep("Boruta", length(boruta.finalVars$Item))),
      selected.fea = boruta.finalVars$Item
    )

    selected.feature <- rbind(selected.feature, result)
    write.table(result,file=paste0(outdir,"/4.Boruta_select_features.csv"),sep=",",row.names = F)
    return(selected.feature)
  }

  ##5.CoxBoost
  FS.coxboost<-function(est_dd,
                        selected.feature,
                        fold,
                        outdir,
                        seed=seed){
    ##### 5.CoxBoost ###########
    message("--- 5.CoxBoost  ---")

    set.seed(seed)
    pen <- optimCoxBoostPenalty(est_dd[, "OS_time"], est_dd[, "OS_status"], as.matrix(est_dd[, -c(1, 2)]),
                                trace = TRUE, start.penalty = 500, parallel = F)

    cv.res <- cv.CoxBoost(est_dd[, "OS_time"], est_dd[, "OS_status"], as.matrix(est_dd[, -c(1, 2)]),
                          maxstepno = 500, K = fold, type = "verweij", penalty = pen$penalty, parallel = F)
    fit <- CoxBoost(est_dd[, "OS_time"], est_dd[, "OS_status"], as.matrix(est_dd[, -c(1, 2)]),
                    stepno = cv.res$optimal.step, penalty = pen$penalty)
    rid <- as.data.frame(coef(fit))
    rid$id <- rownames(rid)
    rid <- rid[which(rid$`coef(fit)` != 0), "id"]
    result <- data.frame(
      method = c(rep("CoxBoost", length(rid))),
      selected.fea = rid
    )

    selected.feature <- rbind(selected.feature, result)
    write.table(result,file=paste0(outdir,"/5.Coxboost_select_features.csv"),sep=",",row.names = F)
    return(selected.feature)
  }

  ##6.RSF
  FS.rsf<-function(est_dd,
                   outdir,
                   seed=seed){
    ##### 6.RSF ###########
    message("--- 6.RSF  ---")
    set.seed(seed)
    ns_res<-tune.nodesize(Surv(OS_time, OS_status) ~ ., est_dd)
    ###
    nodesize=ns_res$nsize.opt
    print(paste0("The optimal nodesize is ",nodesize))
    set.seed(seed)
    fit <- rfsrc(Surv(OS_time, OS_status) ~ .,
                 data = est_dd,
                 ntree = 1000, nodesize = nodesize,
                 splitrule = "logrank",
                 importance = T,
                 proximity = T,
                 forest = T,
                 seed = seed
    )
    rid <- var.select(object = fit, conservative = "high")
    rid <- rid$topvars

    result <- data.frame(
      method = c(rep("RSF", length(rid))),
      selected.fea = rid
    )

    selected.feature <- rbind(selected.feature, result)
    write.table(result,file=paste0(outdir,"/6.RSF_select_features.csv"),sep=",",row.names = F)
    return(selected.feature)
  }


  ##7.RF-RFE
  FS.rfe<-function(est_dd,
                   selected.feature,
                   fold,
                   outdir,
                   seed=seed){
    ##### 6.RFE ###########
    message("--- rf-RFE  ---")
    est_dd$OS_status <- as.factor(est_dd$OS_status)
    set.seed(seed)
    control <- rfeControl(functions = rfFuncs,
                          method = "repeatedcv",
                          number = fold,
                          repeats= 5)
    sizes <- c(1,2, 5, 10, 15,20,30,50,100)
    ###
    set.seed(seed)
    rfe_results <- rfe(x = est_dd[, -c(1:2)],
                       y = est_dd$OS_status,
                       sizes = sizes,
                       rfeControl = control)
    ###
    selected_features <- predictors(rfe_results)
    ###
    result <- data.frame(
      method = c(rep("RF-RFE", length(selected_features))),
      selected.fea = selected_features
    )

    selected.feature <- rbind(selected.feature, result)
    write.table(result,file=paste0(outdir,"/7.RF-RFE_select_features.csv"),sep=",",row.names = F)
    return(selected.feature)
  }

  ##8.stepcox
  FS.stepCox<-function(est_dd,
                       selected.feature,
                       outdir,
                       seed=seed){
    ##### 8.StepCox ###########
    message("--- 8.StepCox ---")
    all_result<-c()

    for (direction in c("both", "backward", "forward")) {
      fit <- tryCatch(
        stepAIC(coxph(Surv(OS_time, OS_status) ~ ., est_dd), direction = direction),
        error = function(e) {
          message(paste("Error in direction:", direction, "-", e$message))
          return(NULL)
        }
      )

      if (!is.null(fit)) {
        rid <- names(coef(fit))
        result <- data.frame(
          method = c(rep(paste0("StepCox", "+", direction), length(rid))),
          selected.fea = rid
        )
        all_result <- rbind(all_result, result)
      }
    }

    selected.feature <- rbind(selected.feature, all_result)
    write.table(all_result,file=paste0(outdir,"/8.stepcox_select_features.csv"),sep=",",row.names = F)
    return(selected.feature)
  }

  ##9.xgboost
  FS.xgboost<-function(est_dd,
                       selected.feature,
                       outdir,
                       seed=seed){
    ##### 9.Xgboost ###########
    message("--- 9.Xgboost  ---")
    set.seed(seed)
    train <- apply(est_dd[, -c(1)], 2, as.numeric) %>% as.data.frame()
    train_matrix <- sparse.model.matrix(OS_status ~ . - 1, data = train)
    train_label <- as.numeric(train$OS_status)
    train_fin <- list(data = train_matrix, label = train_label)
    dtrain <- xgb.DMatrix(data = train_fin$data, label = train_fin$label)
    xgb <- xgboost(
      data = dtrain, max_depth = 6, eta = 0.5,
      objective = "binary:logistic", nround = 25
    )
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
    write.table(result,file=paste0(outdir,"/9.xgboost_select_features.csv"),sep=",",row.names = F)
    return(selected.feature)
  }

  ##10. svm main function
  svmRFE.wrap <- function(test.fold, X, ...) {
      # Wrapper to run svmRFE function while omitting a given test fold
      train.data <- X[-test.fold, ]
      test.data <- X[test.fold, ]

      # Rank the features
      features.ranked <- svmRFE(train.data, ...)

      return(list(feature.ids = features.ranked, train.data.ids = row.names(train.data), test.data.ids = row.names(test.data)))
    }

  svmRFE <- function(X, k = 1, halve.above = 5000) {
      # Feature selection with Multiple SVM Recursive Feature Elimination (RFE) algorithm
      n <- ncol(X) - 1

      # Scale data up front so it doesn't have to be redone each pass
      cat("Scaling data...")
      X[, -1] <- scale(X[, -1])
      cat("Done!\n")
      flush.console()

      pb <- txtProgressBar(1, n, 1, style = 3)

      i.surviving <- 1:n
      i.ranked <- n
      ranked.list <- vector(length = n)

      # Recurse through all the features
      while (length(i.surviving) > 0) {
        if (k > 1) {
          # Subsample to obtain multiple weights vectors (i.e. mSVM-RFE)
          folds <- rep(1:k, len = nrow(X))[sample(nrow(X))]
          folds <- lapply(1:k, function(x) which(folds == x))

          # Obtain weights for each training set
          w <- lapply(folds, getWeights, X[, c(1, 1 + i.surviving)])
          w <- do.call(rbind, w)

          # Normalize each weights vector
          w <- t(apply(w, 1, function(x) x / sqrt(sum(x^2))))

          # Compute ranking criteria
          v <- w * w
          vbar <- apply(v, 2, mean)
          vsd <- apply(v, 2, sd)
          c <- vbar / vsd
        } else {
          # Only do 1 pass (i.e. regular SVM-RFE)
          w <- getWeights(NULL, X[, c(1, 1 + i.surviving)])
          c <- w * w
        }

        # Rank the features
        ranking <- sort(c, index.return = T)$ix
        if (length(i.surviving) == 1) {
          ranking <- 1
        }

        if (length(i.surviving) > halve.above) {
          # Cut features in half until less than halve.above
          nfeat <- length(i.surviving)
          ncut <- round(nfeat / 2)
          n <- nfeat - ncut

          cat("Features halved from", nfeat, "to", n, "\n")
          flush.console()

          pb <- txtProgressBar(1, n, 1, style = 3)
        } else {
          ncut <- 1
        }

        # Update feature list
        ranked.list[i.ranked:(i.ranked - ncut + 1)] <- i.surviving[ranking[1:ncut]]
        i.ranked <- i.ranked - ncut
        i.surviving <- i.surviving[-ranking[1:ncut]]

        setTxtProgressBar(pb, n - length(i.surviving))
        flush.console()
      }

      close(pb)

      return(ranked.list)
    }

  getWeights <- function(test.fold, X) {
      # Fit a linear SVM model and obtain feature weights
      train.data <- X
      if (!is.null(test.fold)) train.data <- X[-test.fold, ]

      svmModel <- svm(train.data[, -1], train.data[, 1],
                      cost = 10, cachesize = 500,
                      scale = F, type = "C-classification", kernel = "linear"
      )

      t(svmModel$coefs) %*% svmModel$SV
    }

  WriteFeatures <- function(results, input, save = T, file = "features_ranked.txt") {
      # Compile feature rankings across multiple folds
      featureID <- sort(apply(sapply(results, function(x) sort(x$feature, index.return = T)$ix), 1, mean), index = T)$ix
      avg.rank <- sort(apply(sapply(results, function(x) sort(x$feature, index.return = T)$ix), 1, mean), index = T)$x
      feature.name <- colnames(input[, -1])[featureID]
      features.ranked <- data.frame(FeatureName = feature.name, FeatureID = featureID, AvgRank = avg.rank)
      if (save == T) {
        write.table(features.ranked, file = file, quote = F, row.names = F)
      } else {
        features.ranked
      }
    }

  FeatSweep.wrap <- function(i, results, input) {
      # Wrapper to estimate generalization error across all hold-out folds, for a given number of top features
      svm.list <- lapply(results, function(x) {
        e1071::tune(svm,
                    train.x = input[x$train.data.ids, 1 + x$feature.ids[1:i]],
                    train.y = input[x$train.data.ids, 1],
                    validation.x = input[x$test.data.ids, 1 + x$feature.ids[1:i]],
                    validation.y = input[x$test.data.ids, 1],
                    # Optimize SVM hyperparamters
                    ranges = e1071::tune(svm,
                                         train.x = input[x$train.data.ids, 1 + x$feature.ids[1:i]],
                                         train.y = input[x$train.data.ids, 1],
                                         ranges  = list(gamma = 2^(-12:0), cost = 2^(-6:6))
                    )$best.par,
                    tunecontrol = tune.control(sampling = "fix")
        )$perf
      })

      error <- mean(sapply(svm.list, function(x) x$error))
      return(list(svm.list = svm.list, error = error))
    }

  ######### the main of the function ##########
  plan(multisession, workers = ncore)
  if (!is.na(seed) &
      identical(c("ID", "OS_time", "OS_status"), colnames(InputMatrix)[1:3]) &
      identical(c("ID", "OS_time", "OS_status"), common_feature[1:3]) &
      length(common_feature) > 3) {
      message("--- Data preprocessing ---")

      # Matching candidate genes to genes in each cohort
      InputMatrix <- InputMatrix[, common_feature]
      InputMatrix <- InputMatrix[!is.na(InputMatrix$OS_time) & !is.na(InputMatrix$OS_status), ]
      InputMatrix[, -c(1:3)] <- apply(InputMatrix[, -c(1:3)], 2, function(x) {
        x[is.na(x)] <- mean(x, na.rm = T)
        return(x)
      })
      est_dd <- as.data.frame(InputMatrix)[, common_feature[-1]]
      pre_var <- common_feature[-c(1:3)]

      ##### setting the pamameters ######
      seed <- seed
      iter.times <- 1000
      selected.feature <- data.frame()

      if (method=="all"){
        ### 1. Repeated Lasso  #############
        message("--- 1.Repeated lasso ---")
        selected.feature<-FS.lasso(est_dd,pre_var,iter.times,fold,selected.feature,outdir,seed=seed)

        ### 2.Enet ###########
        message("--- 2.Enet  ---")
        selected.feature<-FS.enet(est_dd,pre_var,iter.times,fold,selected.feature,outdir,seed=seed)

        ### 3. Ridge #############
        message("--- 3.Ridge ---")
        selected.feature<-FS.ridge(est_dd,pre_var,iter.times,fold,selected.feature,outdir,seed=seed)

        ### 4.Boruta ###########
        message("--- 4.Boruta  ---")
        selected.feature<-FS.boruta(est_dd,selected.feature,outdir,seed=seed)

        ### 5.coxboost ##########
        message("--- 5.coxboost  ---")
        selected.feature<-FS.coxboost(est_dd,selected.feature,fold,outdir,seed=seed)

        ### 6.rsf ##########
        message("--- 6.RSF  ---")
        selected.feature<-FS.rsf(est_dd,outdir,seed=seed)

        ##### 7.RFE ##########
        message("--- 7.RFE  ---")
        selected.feature<-FS.rfe(est_dd,selected.feature,fold,outdir,seed=seed)

        ##### 8.stepcox ##########
        message("--- 8.stepcox  ---")
        selected.feature<-FS.stepCox(est_dd,selected.feature,outdir,seed=seed)

        ##### 9.xgboost ##########
        message("--- 9.xgboost  ---")
        selected.feature<-FS.xgboost(est_dd,selected.feature,outdir,seed=seed)

        if (svm_method){
          ##### 10.SVM-REF ##########
          message("--- 10.SVM-REF  ---")
          print("This step will probably take several hours")
          input <- est_dd[, -1]
          input[,1]<-as.factor(input[,1])
          # 10CV (k-fold cross Validation)
          nfold <- fold
          nrows <- nrow(input)
          folds <- rep(1:nfold, len = nrows)[sample(nrows)]
          folds <- lapply(1:nfold, function(x) which(folds == x))
          results <- lapply(folds, svmRFE.wrap, input, k = fold, halve.above = 50)
          top.features <- WriteFeatures(results, input, save = F)
          n.features <- nrow(top.features)
          if (n.features > 300) {
            n.svm <- 300
          } else {
            n.svm <- n.features
          }
          featsweep <- lapply(1:n.features, FeatSweep.wrap, results, input)
          no.info <- min(prop.table(table(input[, 1])))
          errors <- sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))
          fea <- top.features[1:which.min(errors), "FeatureName"]
          result <- data.frame(
            method = c(rep("SVM-REF", length(fea))),
            selected.fea = fea
          )
          write.table(result,file=paste0(outdir,"/10.svm_select_features.csv"),sep=",",row.names = F)
          selected.feature <- rbind(selected.feature, result)
        }


      } else if (method=="boruta"){

        ### boruta function
        boruta_feature_select<-function(est_dd,selected.feature,outdir,seed){

          set.seed(seed)
          data<-est_dd[,-1]
          data$OS_status<-as.factor(data$OS_status)
          boruta<-Boruta(OS_status ~.,data = data,doTrace=2,maxRuns = 1000,ntree=1000)

          print("plot the importance of feature")
          boruta.variable.imp <- boruta.imp(boruta)
          head(boruta.variable.imp)
          p<-sp_boxplot(boruta.variable.imp, melted=T, xvariable = "Variable", yvariable = "Importance",
                        legend_variable = "finalDecision",x_label = "",coordinate_flip=F,
                        title="Feature importance",
                        legend_variable_order = c("shadowMax", "shadowMean", "shadowMin", "Tentative","Confirmed"),
                        xtics_angle = 90)

          print(p)
          ggsave(p,filename = paste0(outdir,"/2.boruta_feature_importance.jpg"),dpi=600,units="cm",width=10,height =8,scale = 1.5)
          write.table(boruta$finalDecision,file=paste0(outdir,"/2.boruta_finalDecision_feature.csv"),sep=",")
          result<-as.data.frame(boruta$finalDecision)
          result<-data.frame(gene=row.names(result),decision=result[,1])
          selected.feature<-result[result$decision=="Confirmed",]$gene
          return(selected.feature)
        }

        selected.feature<-boruta_feature_select(est_dd,selected.feature,
                                                outdir,seed=seed)
      } else if (method == "RSF"){
        message("--- RSF  ---")
        selected.feature<-FS.rsf(est_dd,outdir,seed=seed)
      } else if (method =="xgboost"){
        message("--- xgboost  ---")
        selected.feature<-FS.xgboost(est_dd,selected.feature,outdir,seed=seed)
      } else if (method == "rfe") {
        selected.feature<-FS.rfe(est_dd,selected.feature,fold,outdir,seed=seed)
      } else if (method == "lasso"){
        selected.feature<-FS.lasso(est_dd,pre_var,iter.times,fold,selected.feature,outdir,seed=seed)
      }

      t2<-Sys.time()
      run_time <- t2 - t1
      print(run_time)
      ##
      save(selected.feature,file=paste0(outdir,"/",sprintf("%d_selected_feature.Rdata",fold)))
      return(selected.feature)
  }

}
