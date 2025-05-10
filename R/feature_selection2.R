#' Screening out the core prognostic genes with the machine learning algorithms
#'
#' A function can be used for screening out the core features from the given candidate genes with 9 machine learning algorithms.
#'
#' @param InputMatrix ### 1: ID, 2: OS_time, 3: OS_status (0/1), feature list...
#' @param genelist The candidate genes
#' @param filter_OS_time default F If keep follow up time > 30 days
#' @param meta_time When filter_OS_time=T, set follow up time months(m) or days(d) or years(y)
#' @param deg differential expression filtering
#' @param up when up=T, HR>1 in unicox analysis and log2FC>0 in differential analysis
#' @param outdir the output directory
#' @param seed  The seed
#' @param height 8
#' @param width 10
#' @return feature list
#' @export
#' @examples
#' \donttest{
#' # Requires trained data and genelist
#' data(train_data)
#' selected.feature<-feature_selection2(InputMatrix,genelist=genelist)
#' }
feature_selection2 <- function(InputMatrix,
                                genelist,
                                filter_OS_time=F,
                                meta_time="d",
                                deg=T,
                                up=F,
                                outdir='1.feature_select/',
                                seed = 123,
                                height=8,
                                width=10) {
  ### loading the packages ####
  if (T) {
    Biocductor_packages <- c(
      "survival",
      "glmnet",
      "Boruta",
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
      InputMatrix <- InputMatrix[InputMatrix$OS_time > 1, ] # days more than 30
      nsample<-nrow(InputMatrix)
      print(paste0(nsample," with more than 30 follow-up days for the next step"))
    } else if (meta_time=="y"){
      InputMatrix <- InputMatrix[InputMatrix$OS_time > 0.083, ]
      nsample<-nrow(InputMatrix)
      print(paste0(nsample," with more than 30 follow-up days for the next step"))
    } else {
      InputMatrix <- InputMatrix[InputMatrix$OS_time > 30, ]
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

    } else {
      diff_gene<-diff_gene
      print(paste0("Gets ",length(row.names(diff_gene))," differential genes with a pvalue <0.05"))
      write.csv(diff_gene,file=paste0(outdir,"/1.diff_gene_0.05.csv"))
      geneset.1<-row.names(diff_gene)
    }
  }

  ### 1000 times unicox####
  print("1000 times unicox filter")
  outTab <- NULL
  surv <- InputMatrix[,-1]
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
      fmla1 <- as.formula(Surv(data[,"OS_time"],data[,"OS_status"]) ~ data[,gene])
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
  write.table(outTab,file = paste0(outdir,"/2.unicox_feature_1000.csv"),sep=",",row.names = F)
  write.table(filter_outTab,file = paste0(outdir,"/2.filtered_unicox_feature_1000.csv"),sep=",",row.names = F)

  ###
  features<-filter_outTab$gene
  if(deg==T){
    features<-intersect(geneset.1,features)
  }
  message(paste0("---the number of candidate genes is ", length(features), " ---"))
  input<-InputMatrix[,c("OS_time","OS_status",features)]
  ###
  boruta_feature_select<-function(input,seed,outdir,width,height){
    set.seed(seed)
    data<-input[,-1]
    data$OS_status<-as.factor(data$OS_status)
    boruta<-Boruta(OS_status ~.,data = data,doTrace=2,maxRuns = 1000,ntree=1000)
    boruta
    ##
    print("plot the importance of feature")
    boruta.variable.imp <- boruta.imp(boruta)
    head(boruta.variable.imp)
    library("YSX")
    p<-sp_boxplot(boruta.variable.imp, melted=T, xvariable = "Variable", yvariable = "Importance",
                  legend_variable = "finalDecision",x_label = "",coordinate_flip=F,
                  title="Feature importance",
                  legend_variable_order = c("shadowMax", "shadowMean", "shadowMin", "Tentative","Confirmed"),
                  xtics_angle = 90)

    print(p)
    write.table(boruta$finalDecision,file=paste0(outdir,"/3.boruta_finalDecision_feature.csv"),sep=",")
    result<-as.data.frame(boruta$finalDecision)
    result<-data.frame(gene=row.names(result),decision=result[,1])
    selected.feature<-result[result$decision=="Confirmed",]$gene
    return(selected.feature)
  }
  boruta.imp <- function(x){
    library(dplyr)
    library(Boruta)
    library(caret)
    imp <- reshape2::melt(x$ImpHistory, na.rm=T)[,-1]
    colnames(imp) <- c("Variable","Importance")
    imp <- imp[is.finite(imp$Importance),]

    variableGrp <- data.frame(Variable=names(x$finalDecision),
                              finalDecision=x$finalDecision)

    showGrp <- data.frame(Variable=c("shadowMax", "shadowMean", "shadowMin"),
                          finalDecision=c("shadowMax", "shadowMean", "shadowMin"))

    variableGrp <- rbind(variableGrp, showGrp)

    boruta.variable.imp <- merge(imp, variableGrp, all.x=T)

    sortedVariable <- boruta.variable.imp %>% group_by(Variable) %>%
      summarise(median=median(Importance)) %>% arrange(median)
    sortedVariable <- as.vector(sortedVariable$Variable)


    boruta.variable.imp$Variable <- factor(boruta.variable.imp$Variable, levels=sortedVariable)

    invisible(boruta.variable.imp)
  }
  selected.feature<-boruta_feature_select(input,seed=seed,outdir,width,height)

  t2<-Sys.time()
  run_time <- t2 - t1
  print(run_time)
  save(selected.feature,file=paste0(outdir,"/selected.feature.Rdata"))
  return(selected.feature)
}



