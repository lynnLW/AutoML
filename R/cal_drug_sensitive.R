cal_drug_sensitive<-function(test.data,
                             database="CTRP2",
                             TPM=T,
                             correct="standardize" #"eb" for array,standardize  for rnaseq
                             ){
  #loading packages
  library("oncoPredict")
  library(data.table)
  library(gtools)
  library(reshape2)
  library(ggpubr)
  library(openxlsx)

  if(F){
    load("data/internal_drug_data.rdata")
  }

  ###
  if (database=="CTRP2" & TPM==T){
    training_data=internal_drug_data[['training_data']]$CTRP2_TPM_Expr
    drug_data=internal_drug_data[['drug_data']]$CTRP2
  } else if (database=="CTRP2" & TPM==F) {
    training_data=internal_drug_data[['training_data']]$CTRP2_RPKM_Expr
    drug_data=internal_drug_data[['drug_data']]$CTRP2
  } else if (database=="GDSC1") {
    training_data=internal_drug_data[['training_data']]$GDSC1_Expr
    drug_data=internal_drug_data[['drug_data']]$GDSC1
  } else if (database=="GDSC2") {
    training_data=internal_drug_data[['training_data']]$GDSC2_Expr
    drug_data=internal_drug_data[['drug_data']]$GDSC2
  }
  ###
  if (!identical(rownames(drug_data),colnames(training_data))){
    common<-intersect(rownames(drug_data),colnames(training_data))
    drug_data<-drug_data[common,]
    training_data<-training_data[,common]
    identical(rownames(drug_data),colnames(training_data))
  }
  ###
  calcPhenotype(trainingExprData = training_data,
                trainingPtype = drug_data,
                testExprData = as.matrix(test.data),
                batchCorrect = correct,  #   "eb" for array,standardize  for rnaseq
                powerTransformPhenotype = TRUE,
                removeLowVaryingGenes = 0.2,
                minNumSamples = 10,
                printOutput = TRUE,
                removeLowVaringGenesFrom = 'rawData' )

}
