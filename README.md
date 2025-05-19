# AutoML
Construction of prognostic models using machine learning algorithms 

AutoML is used to select prognostic genes, construct prognostic models, and evaluate model performance using RNA-seq and microarray data. AutoML has 11 built-in ML algorithms.

**Graphical Abstract：** ![E2F-flowchart-AutoML](https://github.com/user-attachments/assets/770a353e-4159-4453-a667-17686ad37e2c)

## Citation

**E2F targets and G2M checkpoint convergence drive prostate cancer progression - a machine learning guided prognostic framework**

Lin Wang#

## Contact

Lin Wang, PhD, [1155116558\@link.cuhk.edu.hk](1155116558@link.cuhk.edu.hk)

Institute of Trauma and Metabolism, Zhengzhou Central Hospital Affiliated to Zhengzhou University, Zhengzhou 450007, China.

Any technical question, please contact Lin Wang ([1155116558\@link.cuhk.edu.hk](1155116558@link.cuhk.edu.hk)).

copyright, [LynnLab\@ZZU](mailto:LynnLab@ZZU)

## Installation

You may install this package with:

```{r}
# options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

list.of.packages <- c("devtools","dplyr", "stringr", "survival", "survminer", "survcomp", "aplot", "ggplot2", "ggpubr", "caret", "survivalROC","e1071", "GSVA", "glmnet",  "msigdbr", "randomForestSRC" , "plsRcox", "superpc", "gbm", "CoxBoost", "xgboost", "mboost")
#checking missing packages from the list
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
packToInst <- setdiff(list.of.packages, installed.packages())

lapply(packToInst, function(x){
  BiocManager::install(x,ask = F,update = F)
})

# You can install AutoML from Github:
devtools::install_github("lynnLW/AutoML")
```

## Example

### Feature select of candidate genes:

```{r}
## load R package and internal data set
library(AutoML)
load("inst/extdata/train_data.Rdata") # load example data
head(InputMatrix[1:5,1:10])
                               ID OS.time OS       MYC   CTNNB1       JAG2    NOTCH1        DLL1      AXIN2     PSEN2
# TCGA.ZG.A9NI.01A TCGA.ZG.A9NI.01A    4.34  0 1.4807221 1.773289  0.6772469 0.3712124  0.25047264  0.2502103 0.9650362
# TCGA.ZG.A9ND.01A TCGA.ZG.A9ND.01A   13.47  0 1.1940904 1.839368  0.3133605 0.1363740  0.41909796 -0.4465003 1.1581636
# TCGA.ZG.A9N3.01A TCGA.ZG.A9N3.01A   11.47  0 1.3785982 1.729136 -0.2316246 0.5081497  0.13855600  0.2678887 0.8667854
# TCGA.ZG.A9MC.01A TCGA.ZG.A9MC.01A   14.95  0 1.5653816 1.904166  0.8083476 0.5997976  0.76862466  0.6978906 0.9582685
# TCGA.ZG.A9M4.01A TCGA.ZG.A9M4.01A   17.97  0 0.8217003 2.078277  0.2438614 0.1479537 -0.01707259  0.4449893 1.0062466

load("inst/extdata/genelist.Rdata") # load gene list data
genelist
# [1] "MYC"    "CTNNB1" "JAG2"   "NOTCH1" "DLL1"   "AXIN2"  "PSEN2"  "FZD1"   "NOTCH4" "LEF1"   "AXIN1"  "NKD1"   "WNT5B" 
# [14] "CUL1"   "JAG1"   "MAML1"  "KAT2A"  "GNAI1"  "WNT6"   "PTCH1"  "NCOR2"  "DKK4"   "HDAC2"  "DKK1"   "TCF7"   "WNT1"  
# [27] "NUMB"   "ADAM17" "DVL2"   "PPARD"  "NCSTN"  "HDAC5"  "CCND2"  "FRAT1"  "CSNK1E" "RBPJ"   "FZD8"   "TP53"   "SKP2"  
# [40] "HEY2"   "HEY1"   "HDAC11"

## feature select by DEG, cox, and ML filtering
selected.feature<-feature_selection(InputMatrix,
                                    genelist=genelist,
                                    outdir="1.feature_selection/")
head(selected.feature[1:5,]) # view selected feature
#           method selected.fea
#1           Lasso        AXIN1
#2           Lasso         JAG1
#3           Lasso        KAT2A
#4           Lasso        NCSTN
#5 Enet[alpha=0.5]        AXIN1

## Identified prognostic gene selected by at least 7 ML algorithms
f<-top_feature_select(selected.feature = selected.feature,
                        nmethod = 7,
                        width=7.5,
                        height = 10,
                        outdir="1.feature_selection/")
# The final selected genes: 4 
```
![1](https://github.com/user-attachments/assets/8603b5c5-1c36-44f4-b10b-e3f882b3945c)


### Construction of ML-based prognosis models

``` r
candidate_genes<-f # The final selected genes
train_data<-InputMatrix # The training data

# modeling
model.list<-ML.survival.model(train_data,
                                candidate_genes[1:2],
                                filter_OS_time=F,
                                meta_time="m",
                                cor=F,
                                cor_threshold=0.85,
                                fold=5,
                                rep=10,
                                p=0.75,
                                deep=F,
                                outdir="2.train/",
                                seed=5201314,
                                ncore=4)
```

### Model performance evaluation in the training data

``` r
load("inst/extdata/10_5_model_list.RData") # result from ML.survival.model
# extract cindex list
cindex_list<-lapply(model_list,function(x)x$metrics_list)
cindex_rank2(cindex_list,order="valid",index="all",outdir="3.figure/",plot_type="boxplot")
```

![2](https://github.com/user-attachments/assets/2da7b68d-bc1e-4165-b4f2-331f88f43339)


```{r}
# extract model list
model.list<-lapply(model_list,function(x)x$final_model)
cindex_rank2(model_list=model.list,outdir="3.figure/")
```
![image](https://github.com/user-attachments/assets/dc34422b-844c-41ee-a38e-385a18570812)



### Model performance evaluation in the external cohort 

```{r}
# loading external validation cohort
load("inst/extdata/train_features.Rdata") # The final selected genes
load("inst/extdata/list_train_vali_Data") # The validation cohort list

candidate_genes<-common_feature[4:length(common_feature)]
print(candidate_genes)
#[1] "AXIN1" "JAG1"  "KAT2A" "NCSTN"

# calculate c-index and time-dependent AUC values
model_auc_list<-cal_vali_index(list_train_vali_Data,candidate_genes,model.list,rep=1,outdir="4.test/")
```

### Plot C-index in all cohort

```{r}
cindex_rank(vali_auc_list = model_auc_list,index="cindex",train="Train",plot_type="barplot",outdir="4.test/")
cindex_rank(vali_auc_list = model_auc_list,index="km_auc_1",train="Train",plot_type="barplot",outdir="4.test/")
cindex_rank(vali_auc_list = model_auc_list,index="km_auc_2",train="Train",plot_type="barplot",outdir="4.test/")
cindex_rank(vali_auc_list = model_auc_list,index="km_auc_3",train="Train",plot_type="barplot",outdir="4.test/")
```
![3](https://github.com/user-attachments/assets/f66a09dd-9083-4cd2-b0a8-fbaf9b28a6ef)


### Plot roc curves in all cohort
```{r}
roc_plot(vali_auc_list = model_auc_list,model="all",outdir="4.test/")

```
![4](https://github.com/user-attachments/assets/ee170036-f123-4249-94e4-4ecaad3fa370)


### plot KM survival curves in all cohort
```{r}
surv_plot(vali_auc_list = model_auc_list,model="all",outdir="4.test/")
```
![5](https://github.com/user-attachments/assets/48005504-a230-4067-8182-d2436c04e44b)


### Multicox analysis in all cohorts
```{r}
# loading meta information of each cohort")
load("inst/extdata/list_train_vali_meta")

# Extract the risk score of the best model
model_name='GBM'
rs_list=lapply(model_auc_list[[model_name]],function(x){x[[1]]$pred_df})
outdir=paste0("5.multicox/GBM")
dir.create(outdir,recursive = T)

# Multicox analysis for each cohort
con_summary<-c()
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
      outdir = "5.multicox/GBM",
      cut_type = NULL
    )
}
```
![6](https://github.com/user-attachments/assets/7bb414ae-bcf3-4ab1-80e4-56e5e1e0a142)



### Comparison with published signatures
```{r}

load("inst/extdata/test_index.Rdata")
load("inst/extdata/all_model_list.rdata")
load("inst/extdata/published_auc_list.rdata")
own_auc_list[['Own']]<-model_auc_list

indexs=c("Cindex","AUC_1","AUC_2","AUC_3","AUC_5","AUC_7")
dataset<-names(own_auc_list[[1]][[1]])
outdir=paste0("../test/6.comparison/")
dir.create(paste0(outdir,"/",model_name),recursive = T)
model_name="GBM"
p<-index_comp(own_auc_list =own_auc_list,
                published_auc_list=published_auc_list,
                model_name='GBM',  
                dataset=names(own_auc_list[[1]][[1]]),
                index="Cindex")
p
ggsave(p,file=paste0(outdir,"/",model_name,"/",index,".jpg"),dpi=600,height =6,width =30,units="cm")
```
