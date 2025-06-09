#' Calculate Drug Sensitivity Predictions with Custom Output Path
#'
#' Predicts drug sensitivity scores for test samples with customizable output location.
#'
#' @param test.data Gene expression matrix (genes in rows, samples in columns)
#' @param database Drug sensitivity database ("CTRP2", "GDSC1", or "GDSC2")
#' @param TPM Logical indicating whether to use TPM-normalized data (for CTRP2 only)
#' @param output_dir Directory path for output files (default: "./calcPhenotype_Output")
#' @param output_filename Name for output CSV file (without extension, default: "DrugPredictions")
#' @return Matrix of predicted drug sensitivity scores (invisibly)
#' @importFrom oncoPredict calcPhenotype
#' @importFrom utils write.csv
#' @export
cal_drug_sensitive <- function(test.data,
                               database = c("CTRP2", "GDSC1", "GDSC2","GDSC1_new","GDSC2_new"),
                               TPM = TRUE,
                               output_dir = "./calcPhenotype_Output",
                               output_filename = "DrugPredictions") {

  requireNamespace("oncoPredict")

  # validate parameter
  if (!database %in% c("CTRP2", "GDSC1", "GDSC2","GDSC1_new","GDSC2_new")) {
    stop("database must be one of: 'CTRP2', 'GDSC1', 'GDSC2','GDSC1_new','GDSC2_new'")
  }

  if (!is.matrix(test.data) && !is.data.frame(test.data)) {
    stop("test.data must be a matrix or data.frame with genes in rows and samples in columns")
  }

  # data examine
  if (nrow(test.data) < 100) {
    warning("Test data has very few genes (<100), results may be unreliable")
  }

  # output
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
  }
  output_path <- file.path(output_dir, paste0(output_filename, ".csv"))

  # loading
  data_files <- c(
    system.file("extdata", "internal_drug_data.rdata", package = "AutoML"),
    system.file("extdata", "internal_new_drug_data.rdata", package = "AutoML"),
    system.file("extdata", "internal_expr_data.rdata", package = "AutoML")
  )

  if (!all(file.exists(data_files))) {
    stop("Required data files not found in AutoML package")
  }

  # 加载数据
  if(database %in% c("CTRP2", "GDSC1", "GDSC2")){
    internal_drug_data <- get(load(data_files[1]))
    internal_expr_data <- get(load(data_files[3]))
  } else if (database %in% c("GDSC1_new", "GDSC2_new")){
    internal_drug_data <- get(load(data_files[2]))
    internal_expr_data <- get(load(data_files[3]))
  }

  # 选择训练数据
  training_data <- switch(database,
                          "CTRP2" = {
                            if (TPM) {
                              internal_drug_data$training_data$CTRP2_TPM_Expr
                            } else {
                              internal_drug_data$training_data$CTRP2_RPKM_Expr
                            }
                          },
                          "GDSC1" = internal_expr_data$GDSC1_Expr,
                          "GDSC2" = internal_expr_data$GDSC2_Expr,
                          "GDSC1_new" = internal_expr_data$GDSC1_Expr,
                          "GDSC2_new" = internal_expr_data$GDSC1_Expr
  )

  drug_data <- internal_drug_data[[database]]

  # Ensure matching samples between training and drug data
  if (!identical(rownames(drug_data), colnames(training_data))) {
    common <- intersect(rownames(drug_data), colnames(training_data))
    if (length(common) == 0) {
      stop("No matching samples found between training and drug data")
    }
    drug_data <- drug_data[common, , drop = FALSE]
    training_data <- training_data[, common, drop = FALSE]
  }

  # Calculate drug sensitivity predictions
  results <- tryCatch(
    oncoPredict::calcPhenotype(
      trainingExprData = training_data,
      trainingPtype = as.matrix(drug_data),
      testExprData = as.matrix(test.data),
      batchCorrect = "none",
      powerTransformPhenotype = T,
      removeLowVaryingGenes = 0.2,
      minNumSamples = 10,
      printOutput = TRUE,
      removeLowVaringGenesFrom = 'rawData',
      folder=F
    ),
    error = function(e) {
      stop("Drug sensitivity prediction failed: ", e$message)
    }
  )

  # 保存结果
  write.csv(results, file = output_path, row.names = TRUE)
  message("Results saved to: ", output_path)

  # return result
  invisible(results)
}
