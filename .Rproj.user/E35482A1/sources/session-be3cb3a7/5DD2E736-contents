  SigKMcox <- function(gene_list,
                         inputSet,
                         KM_pcutoff=0.05, # p-value
                         d="m" ##days, months, or years # KM的筛选阈值
    ) {
    print("Starting the data preprocess of KM")
    ############### data preprocess#######

    library(survival)
    library(dplyr)
    # rename column 2 3
    names(inputSet)[2]<-"OS.time"
    names(inputSet)[3]<-"OS"
    # null value replace with zero
    inputSet[is.na(inputSet)] <- 0
    # table(is.na(inputSet))
    inputSet <- inputSet %>% as.data.frame()
    inputSet$OS.time <- as.numeric(inputSet$OS.time)
    if (d=="m") {
      inputSet <- inputSet[inputSet$OS.time > 1, ] # days more than 30
    } else if (d=="y"){
      inputSet <- inputSet[inputSet$OS.time > 0.083, ]
    } else {
      inputSet <- inputSet[inputSet$OS.time > 30, ]
    }
    nsample<-nrow(inputSet)
    print(paste0(nsample," with more than 30 follow-up days for the next step"))

    # 将genelist和表达矩阵的基因名称格式统一
      gene_list <- gsub("-", ".", gene_list)
      gene_list <- gsub("_", ".", gene_list)
      colnames(inputSet)[4:ncol(inputSet)] <- gsub("-", ".", colnames(inputSet)[4:ncol(inputSet)])
      colnames(inputSet)[4:ncol(inputSet)] <- gsub("_", ".", colnames(inputSet)[4:ncol(inputSet)])

      print("Gets the intersection of genelist and expression profile")
      # 获取genelist和表达谱的交集
      comsa1 <- intersect(colnames(inputSet)[4:ncol(inputSet)], gene_list)
      # write.table(comsa1,"2.intersection_genelist_exprSet_gene.txt", row.names = F, quote = F)

      print("Processing the  input representation matrix")
      # 对输入的表达矩阵进行处理
      inputSet <- inputSet[, c("ID", "OS.time", "OS", comsa1)]

      inputSet[, c(1:2)] <- apply(inputSet[, c(1:2)], 2, as.factor)
      inputSet[, c(2:ncol(inputSet))] <- apply(inputSet[, c(2:ncol(inputSet))], 2, as.numeric)
      inputSet <- as.data.frame(inputSet)
      # rownames(inputSet) <- inputSet$ID
      print("Data preprocessing completed")
      # 自定义显示进程函数
      display.progress <- function(index, totalN, breakN = 20) {
        if (index %% ceiling(totalN / breakN) == 0) {
          cat(paste(round(index * 100 / totalN), "% ", sep = ""))
        }
      }
      ############### 单变量cox#######
      print("Stating the KM survival")


      # KM估计筛选基因
      kmoutput <- NULL

      for (i in 1:ncol(inputSet[, 4:ncol(inputSet)])) {
        display.progress(index = i, totalN = ncol(inputSet), breakN = 20)
        g <- colnames(inputSet[, 4:ncol(inputSet)])[i]
        tmp <- inputSet[, c("OS.time", "OS", g)]
        tmp$group <- ifelse(tmp[, 3] > median(tmp[, 3]), "High", "Low")
        fitd <- survdiff(Surv(OS.time, OS) ~ group, data = tmp, na.action = na.exclude)
        p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
        kmoutput <- rbind(kmoutput, data.frame(
          gene = g,
          pvalue = p.val,
          stringsAsFactors = F
        ))
      }

      print("Finished the KM selection & save result")
      skmoutput <- kmoutput[which(kmoutput$pvalue < KM_pcutoff),]
      write.table(skmoutput,paste("2.KM_",KM_pcutoff,"_result.csv",sep = ""),row.names = F, quote = F,sep=",")
      # 进行变量筛选
      selgene <- kmoutput[which(kmoutput$pvalue < KM_pcutoff), "gene"]
      return(selgene)
    }
