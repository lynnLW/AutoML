Corefeature.Prog.Screen2 <- function(InputMatrix, ### first column ID,second OS.time, third OS, (0/1), feature list...
                                       genelist,# all, single,all_without_SVM
                                       seed = NULL,
                                       d="m") {
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
      "survivalROC",
      "future.apply",
      "plotROC",
      "survminer"
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

  ### 1000 times unicox####
  print("1000 times unicox filter")
  train<-InputMatrix[,-1]
  outTab <- NULL
  surv <- train
  for(i in 3:ncol(surv)){ # survival information (OS in this case)
    # customized function
    display.progress = function (index, totalN, breakN=20) {
      if ( index %% ceiling(totalN/breakN)  ==0  ) {
        cat(paste(round(index*100/totalN), "% ", sep=""))
      }
    }

    display.progress(index = i, totalN = ncol(surv)) # show running progression
    gene <- colnames(surv)[i]
    Mboot <- future_replicate(1000, expr = { # bootstrap for 1,000 times
      indices <- sample(rownames(surv), size = nrow(surv) * 0.8, replace = F) # extract 80% samples at each bootsratp
      data <- surv[indices,]
      fmla1 <- as.formula(Surv(data[,"OS.time"],data[,"OS"]) ~ data[,gene])
      mycox <- coxph(fmla1,data = data)
      coxResult <- summary(mycox)
      P <- coxResult$coefficients[,"Pr(>|z|)"]
    }
    )
    times <- length(Mboot[which(Mboot < 0.01)])
    outTab <- rbind(outTab,
                    cbind(gene = gene,
                          times = times))
  }
  outTab <- as.data.frame(outTab)
  outTab$times<-as.integer(outTab$times)
  filter_outTab<-outTab[outTab$times>800,]
  dir.create(paste0("feature_select/"),recursive = T)
  write.table(filter_outTab,file = paste0("feature_select/filtered_unicox_feature_1000.csv"),sep=",",row.names = F)
  features<-filter_outTab$gene

  ###
  input<-InputMatrix[,c("OS.time","OS",features)]
  selected.feature<-boruta_feature_select(input,seed=123456)
}



