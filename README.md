---
editor_options: 
  markdown: 
    wrap: 72
---

# ML

## Step 1.1 - Data List Preparetion

`outdir="dataset/all/mm"`

`sets=seq(1:12)`

`correct_method="none"`

`source("R/run_data_combine.R")`

`list_vali_data<-run_data_combine(list_data_surv, sets=sets,correct_method = NULL, outdir="dataset/all/mm/")`

## Step 1.2 - Meta List Preparetion

`outdir="dataset/all/mm"`
`sets=seq(1:12)`
`source("R/run_meta_combine.R")`
`list_meta_data<-run_meta_combine(list_data_surv,sets=sets,outdir="dataset/all/mm/")`

