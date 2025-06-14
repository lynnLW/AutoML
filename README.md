# AutoML
Construction of prognostic models using machine learning (ML) algorithms 

AutoML is used to select prognostic genes, construct predictive models, and evaluate model performance using RNA-seq and microarray data. AutoML has 11 built-in ML algorithms.

**Graphical Abstract：** ![E2F-flowchart-AutoML](https://github.com/user-attachments/assets/770a353e-4159-4453-a667-17686ad37e2c)

## Citation

**E2F targets and G2M checkpoint convergence drive prostate cancer progression - a machine learning guided prognostic framework**

Lin Wang#

## Contact

Lin Wang, PhD, [1155116558\@link.cuhk.edu.hk](1155116558@link.cuhk.edu.hk)

Institute of Trauma and Metabolism, Zhengzhou Central Hospital Affiliated to Zhengzhou University, Zhengzhou 450007, China.

Any technical question, please get in touch with Lin Wang ([1155116558\@link.cuhk.edu.hk](1155116558@link.cuhk.edu.hk)).

copyright, [LynnLab\@ZZU]

## Installation

You may install this package with:

```{r}
# options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

list.of.packages <- c("devtools", "dplyr", "stringr", "survival", "survminer", "survcomp", "aplot", "ggplot2", "ggpubr", "caret", "survivalROC", "e1071", "GSVA", "glmnet",  "msigdbr", "randomForestSRC", "plsRcox", "superpc", "gbm", "CoxBoost", "xgboost", "mboost", "oncoPredict")
#checking missing packages from the list
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
packToInst <- setdiff(list.of.packages, installed.packages())

lapply(packToInst, function(x){
  BiocManager::install(x,ask = F,update = F)
})

# You can install AutoML from GitHub:
devtools::install_github("lynnLW/AutoML")
```

## Example

### Feature selection of candidate genes:

```{r}
## load R package and internal data set
library(AutoML)
data("demo_list")
train_data<-demo_list$train_data # load example data
head(train_data[1:5,1:10])
                               ID OS.time OS       MYC   CTNNB1       JAG2    NOTCH1        DLL1      AXIN2     PSEN2
# TCGA.ZG.A9NI.01A TCGA.ZG.A9NI.01A    4.34  0 1.4807221 1.773289  0.6772469 0.3712124  0.25047264  0.2502103 0.9650362
# TCGA.ZG.A9ND.01A TCGA.ZG.A9ND.01A   13.47  0 1.1940904 1.839368  0.3133605 0.1363740  0.41909796 -0.4465003 1.1581636
# TCGA.ZG.A9N3.01A TCGA.ZG.A9N3.01A   11.47  0 1.3785982 1.729136 -0.2316246 0.5081497  0.13855600  0.2678887 0.8667854
# TCGA.ZG.A9MC.01A TCGA.ZG.A9MC.01A   14.95  0 1.5653816 1.904166  0.8083476 0.5997976  0.76862466  0.6978906 0.9582685
# TCGA.ZG.A9M4.01A TCGA.ZG.A9M4.01A   17.97  0 0.8217003 2.078277  0.2438614 0.1479537 -0.01707259  0.4449893 1.0062466

genelist<-demo_list$genelist # load gene list data
genelist
# [1] "MYC"    "CTNNB1" "JAG2"   "NOTCH1" "DLL1"   "AXIN2"  "PSEN2"  "FZD1"   "NOTCH4" "LEF1"   "AXIN1"  "NKD1"   "WNT5B" 
# [14] "CUL1"   "JAG1"   "MAML1"  "KAT2A"  "GNAI1"  "WNT6"   "PTCH1"  "NCOR2"  "DKK4"   "HDAC2"  "DKK1"   "TCF7"   "WNT1"  
# [27] "NUMB"   "ADAM17" "DVL2"   "PPARD"  "NCSTN"  "HDAC5"  "CCND2"  "FRAT1"  "CSNK1E" "RBPJ"   "FZD8"   "TP53"   "SKP2"  
# [40] "HEY2"   "HEY1"   "HDAC11"

## feature select by DEG, cox, and ML filtering
selected.feature<-feature_selection(train_data,
                                    genelist=genelist,
                                    outdir="test/feature_selection/")
head(selected.feature[1:5,]) # view selected feature
#           method selected.fea
#1           Lasso        AXIN1
#2           Lasso         JAG1
#3           Lasso        KAT2A
#4           Lasso        NCSTN
#5 Enet[alpha=0.5]        AXIN1
```
### Identified prognostic genes selected by at least 7 ML algorithms

```{r}
#example result from feature selection function
## selected.feature<-get(load(system.file("extdata", "5_selected_feature.rda", package = "AutoML"))) 
f<-top_feature_select(selected.feature = selected.feature,
                        nmethod = 7,
                        width=7.5,
                        height = 10,
                        outdir="test/feature_selection/")
# The final selected genes: 4 
```
![图片1](https://github.com/user-attachments/assets/5d593849-b46b-44b0-975e-3b02fe7b7a51)



### Construction of ML-based prognosis models

```{r}
candidate_genes<-f # The final selected genes
# modeling
model_list<-ML.survival.model(train_data,# The training data
                                candidate_genes,
                                filter_OS_time=F,
                                meta_time="m",
                                cor=F,
                                cor_threshold=0.85,
                                fold=5,
                                rep=10,
                                p=0.75,
                                deep_method = F,
                                gbm_method = F, # memory consuming
                                outdir="test/train/",
                                seed=5201314,
                                ncore=4)
```

### Model performance evaluation in the training data

``` r
# example result from ML.survival.model 
# model_list<-get(load(system.file("extdata", "10_5_model_list.rda", package = "AutoML")))
# extract cindex list
cindex_list<-lapply(model_list,function(x)x$metrics_list)
# Show C-index of each ML algorithm
cindex_rank2(cindex_list,order="valid",index="all",outdir="test/evaluation/",plot_type="boxplot")
```
![图片2](https://github.com/user-attachments/assets/9b6cf7c5-ddfd-4a2e-8a13-31b2deff5097)


```{r}
# extract model list
model.list<-lapply(model_list,function(x)x$final_model)
# Show C-index of all final ML models
cindex_rank2(model_list=model.list,outdir="test/evaluation/")
```

![图片3](https://github.com/user-attachments/assets/cf2d9c53-a331-48e7-8762-30ce7a62446c)

### Model performance evaluation in the external cohort 

```{r}
# loading final candidate genes
candidate_genes<-f 
print(candidate_genes)
#[1] "AXIN1" "JAG1"  "KAT2A" "NCSTN"

# loading external validation cohort
list_train_vali_Data<-get(load(system.file("extdata", "list_train_vali_Data", package = "AutoML"))) # The validation cohort list
# calculate c-index and time-dependent AUC values
model_auc_list<-cal_vali_index(list_train_vali_Data,candidate_genes,model.list,rep=1,outdir="test/validation/")
```

### Plot C-index in all cohorts

```{r}
# model_auc_list<-get(load(system.file("extdata", "test_index.rda", package = "AutoML"))) example result from cal_vali_index
cindex_rank(vali_auc_list = model_auc_list,index="cindex",train="Train",plot_type="barplot",outdir="test/validation/")
cindex_rank(vali_auc_list = model_auc_list,index="km_auc_1",train="Train",plot_type="barplot",outdir="test/validation/")
cindex_rank(vali_auc_list = model_auc_list,index="km_auc_2",train="Train",plot_type="barplot",outdir="test/validation/")
cindex_rank(vali_auc_list = model_auc_list,index="km_auc_3",train="Train",plot_type="barplot",outdir="test/validation/")
```

![图片4](https://github.com/user-attachments/assets/50d5875d-c4b1-4a27-a35b-9bfcb80c984d)


### Plot ROC curves in all cohorts
```{r}
roc_plot(vali_auc_list = model_auc_list,model="all",outdir="test/validation/")
```

![图片6](https://github.com/user-attachments/assets/cf120e67-c33b-400b-bc99-5e0a4f3cc79f)

### Plot KM survival curves in all cohorts
```{r}
surv_plot(vali_auc_list = model_auc_list,model="all",outdir="test/validation/")
```
![图片5](https://github.com/user-attachments/assets/51886e4f-9da9-41b7-9ccb-6a4412b02b99)


### Multicox analysis in all cohorts
```{r}
# loading meta information of each cohort
list_train_vali_meta<-get(load(system.file("extdata", "list_train_vali_meta", package = "AutoML"))) # The validation cohort list

# Extract the risk score of the best model
model_name='GBM'
rs_list=lapply(model_auc_list[[model_name]],function(x){x[[1]]$pred_df})

# Multivariate Cox analysis for each cohort
for(i in 1:length(rs_list)){
    # risk table 
    dataset_name<-names(rs_list)[i]
    print(dataset_name)
    rs<-rs_list[[i]]
    names(rs)[3]<-"rs"
    # meta table
    meta<-list_train_vali_meta[[dataset_name]]
    # multicox
    combined_meta<-merge(meta,rs,by="row.names")
    # multicox for continuous rs
    sum.cox<-generate_multicox_analysis(
      data = combined_meta,
      features = colnames(meta),
      gene = "rs",
      dataset_name = names(rs_list)[i],
      outdir = "test/multicox/GBM",
      cut_type = NULL
    )
}
```

![图片9](https://github.com/user-attachments/assets/a6c2ecf8-4a63-473c-bfbe-94b0f2cb061c)


### Comparison with published signatures
```{r}
# own signature
# model_auc_list<-get(load(system.file("extdata", "test_index.rda", package = "AutoML"))) 
own_auc_list<-list()
published_auc_list<-list()
own_auc_list[['Own']]<-model_auc_list
# published signatures modeling using ML.survival.model and cal_vali_index function
published_auc_list[['Published']]<-get(load(system.file("extdata", "published_auc_list.rda", package = "AutoML")))
indices=c("Cindex","AUC_1","AUC_2","AUC_3","AUC_5","AUC_7","AUC_10")
dir.create("test/comparison/GBM",recursive=T)
for (index in indices){
  p<-index_comp(own_auc_list =own_auc_list,
                  published_auc_list=published_auc_list,
                  model_name='GBM',
                  dataset=names(own_auc_list[[1]][[1]]),
                  index="Cindex")
  p
  ggpubr::ggexport(p,filename=paste0("test/comparison/GBM/",index,".jpg"),res=600,height =1800,width =9000)
}
```
![图片7](https://github.com/user-attachments/assets/7068733a-2230-49ec-8f49-61714b076f53)


### Functional enrichment

```{r}
# list_train_vali_Data<-get(load(system.file("extdata", "list_train_vali_Data", package = "AutoML")))  ## Full genes not only feature genes
# Predictive metabolites, metabolic and cancer hallmark pathways, immune cells scores using ssgsea
library(dplyr)
MPI_list<-list()
Immune_list<-list()
hallmark_list<-list()
metabolic_list<-list()
for (i in 1:length(list_train_vali_Data)){
  expr<-list_train_vali_Data[[i]][,-c(1:3)] %>% t() %>% as.data.frame()
  dataset_name=names(list_train_vali_Data)[i]
  gs<-geneset_cal(expr,category="H",output_dir = paste0("test/function/H/",dataset_name)) #hallmark
  hallmark_list[[i]]<-gs
  gs<-geneset_cal(expr,category="TILs",output_dir = paste0("test/function/Immune/TILs/",dataset_name)) # immune cell
  Immune_list[[i]]<-gs
  gs<-geneset_cal(expr,category="MPI",output_dir = paste0("test/function/MPI/",dataset_name)) #predictive metabolites
  MPI_list[[i]]<-gs
  gs<-geneset_cal(expr,category="Metabolic",output_dir = paste0("test/function/Metabolic/",dataset_name)) #metabolic pathways
  metabolic_list[[i]]<-gs
}
names(metabolic_list)<-names(list_train_vali_Data)
names(hallmark_list)<-names(list_train_vali_Data)
names(Immune_list)<-names(list_train_vali_Data)
names(MPI_list)<-names(list_train_vali_Data)

# Setting parallel calculations through a MulticoreParam back-end
# with workers=4 and tasks=100.
# Estimating ssGSEA scores for 11 gene sets.
# [1] "Calculating ranks..."
# [1] "Calculating absolute values from ranks..."
#   |==========================================================================================================| # 100%
# [1] "Normalizing..."
# Finished! Results saving in:test/Metabolic/Train/GSVA_Metabolic_ssgsea.csv
```

### Functional differences between high- and low-risk groups
#### Separating samples into high and low risk groups
```{r}
group_list<-list()
auc_list<-model_auc_list[['GBM']] #data from cal_vali_index function
for (i in 1:length(auc_list)){
  dataset<-names(auc_list)[i]
  rs_df<-auc_list[[i]][[1]]$pred_df
  rs_df$group<-ifelse(rs_df$pred>median(rs_df$pred),"High","Low")
  group_list[[dataset]]<-rs_df
}
save(group_list,file=paste0("test/group_list.rdata"))
```

#### Calculate the correlation between the risk score and functional pathways
```{r}
for (i in 1:length(auc_list)){
  dataset<-names(auc_list)[i]
  rs_df<-group_list[[dataset]]
  immune<-Immune_list[[dataset]]
  immune<-immune[row.names(rs_df),]
  merge_df<-cbind(rs_df$pred,immune)
  names(merge_df)[1]<-"risk_score"
  # calculate correlation
  analyze_drug_correlations(
    data = merge_df,
    risk_score_col = "risk_score",
    r_threshold = 0,
    top_n = 2,
    ncol = 2,
    outdir=paste0("test/function/immune/",dataset))
  # calculate differences
  cal_diff(immune,rs_df$group,outdir=paste0("test/function/immune/",dataset)) 
}
```
![图片8](https://github.com/user-attachments/assets/de0ec9f1-6a68-4fbd-a398-f2576288d02d)


### Calculate the drug sensitivity of each cohort

```{r}
#calculate sensitivity
library(dplyr)
library(limma)
for (i in 1:length(list_train_vali_Data)){
  dataset_name=names(list_train_vali_Data)[i]
  expr<-list_train_vali_Data[[dataset_name]][,-c(1:3)] %>% t() %>% as.data.frame() # col: samples
  cal_drug_sensitive(expr,
                    database ="CTRP2",
                    output_dir="test/drug/",
                    output_filename=paste0(dataset_name,"_drug_sensitivity"))
}
```

### Drug sensitivity differences between high- and low-risk groups

```{r}
# group_list<-get(load(system.file("extdata", "group_list.rda", package = "AutoML")))
drug_file<-list.files("test/drug/",pattern = "csv",full.names = T)
for (i in 1:length(group_list)){
  dataset<-names(group_list)[i]
  file<-drug_file[grep(dataset,drug_file)]
  result<-read.csv(file,row.names = 1)
  head(result[1:5,1:5])
  rs_df<-group_list[[dataset]]
  rs_df<-rs_df[row.names(result),]
  merge_df<-cbind(rs_df$pred,result)
  names(merge_df)[1]<-"risk_score"
  # calculate correlations
  analyze_drug_correlations(
    data = merge_df,
    risk_score_col = "risk_score",
    r_threshold = 0,
    top_n = 10,
    ncol = 2,
    outdir=paste0("test/drug/",dataset))
  # calculate differences
  cal_diff(result,rs_df$group,outdir=paste0("test/drug/",dataset))
}
```
![图片10](https://github.com/user-attachments/assets/3a3aa274-538a-4e0c-a2b2-01457333f346)

```{r}
files<-list.files(paste0("test/drug/"),pattern="cor_results.csv",full.names = T,recursive = T)
result<-data.frame()
for (i in c(1:6)){
  file=read.csv(files[i],row.names = 1)
  dataset_name<-basename(dirname(files[i]))
  file$dataset<-dataset_name
  result<-rbind(result,file)
}
names(result)[1]<-"selected.fea"
names(result)[6]<-"method"
negative<-result[result$p.value<0.05&result$r<(-0.3),]
positive<-result[result$p.value<0.05&result$r>0.3,]

top_feature_select(selected.feature = positive,nmethod = 4,width = 12,height = 12,
                   outdir =paste0("test/drug/CTRP2/up/"))
top_feature_select(selected.feature = negative,nmethod = 4,width = 12,height = 12,
                   outdir = paste0("test/drug/CTRP2/down/"))

write.csv(result,file = paste0("test/drug/CTRP2/all_cor_results.csv"))

```

![图片11](https://github.com/user-attachments/assets/234b7475-1286-4ef0-9272-e9c43926925b)

