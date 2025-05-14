#' Data preprocessing for model performance evaluation
#'
#' @param list_data_surv #column 1 Sample_ID # 2 clinical features
#' @param sets  c(1,2,3,4)
#' @param features clinical variables
#' @param outdir output path for example "Age","PSA","GS","path_T"
#' @return data list
#' @export
#'
run_meta_combine<-function(list_data_surv,
                           sets=NULL,
                           features=c("Sample_ID","Age","PSA","GS","path_T"),
                           outdir="dataset/"){
  
  ##output directory
  if(!dir.exists(outdir)){
      dir.create(outdir,recursive = T)
  }

  ##sets to combat
  if(is.null(sets)==T){
    sets=seq_along(list_data_surv)
  } else {
    sets=sets
  }
  ##list all files that need to be combined
  print(paste0(length(list_data_surv)," datasets"))

  ##sets
  list_vali_meta<-c()
  list_vali_meta<-list_data_surv[sets]
  print(paste0(length(list_vali_meta)," datasets is preparing to extract metainfo"))

  ##common meta feature
  list_vali_meta<-lapply(list_vali_meta,function(x){
    x<-x[,features]
    return(x)
  })

  for (i in 1:length(list_vali_meta)) {
    data <- list_vali_meta[[i]]

    # Set row names
    row.names(data) <- data$Sample_ID
    data <- data[, -1]  # Remove the Sample_ID column

    # Check if 'GS' column exists and contains all NA values
    if (all(is.na(data$Age))) {
      data$Age <- NULL  # Remove GS column if it has all NAs
    } else {
      # Apply conditions to create GS_group
      data$Age <- data$Age
    }

    # Check if 'GS' column exists and contains all NA values
    if (all(is.na(data$PSA))) {
      data$PSA <- NULL  # Remove GS column if it has all NAs
    } else {
      # Apply conditions to create GS_group
      data$PSA <- data$PSA
    }

    # Check if 'GS' column exists and contains all NA values
    if (all(is.na(data$GS))) {
      data$GS <- NULL  # Remove GS column if it has all NAs
    } else {
      # Apply conditions to create GS_group
      data$GS <- ifelse(data$GS >= 8, ">=8", "<8")
    }

    # Check if 'GS' column exists and contains all NA values
    if (all(is.na(data$path_T))) {
      data$path_T <- NULL  # Remove GS column if it has all NAs
    } else {
      # Apply conditions to create GS_group
      data$path_T<-gsub("T","",data$path_T)
      data$path_T <- ifelse(data$path_T >2, ">2", "<=2")
    }

    # Update the dataset in the list
    list_vali_meta[[i]] <- data
  }

  save(list_vali_meta,file=paste0(outdir,"/list_vali_meta.Rdata"))

  print("Finish the combine process")
  return(list_vali_meta)

}
