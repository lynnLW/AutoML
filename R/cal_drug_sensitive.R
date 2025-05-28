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
                               database = c("CTRP2", "GDSC1", "GDSC2"),
                               TPM = TRUE,
                               output_dir = "./calcPhenotype_Output",
                               output_filename = "DrugPredictions") {

  # validate parameter
  if (!database %in% c("CTRP2", "GDSC1", "GDSC2")) {
    stop("database must be one of: 'CTRP2', 'GDSC1', 'GDSC2'")
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

  # loading drug data
  file_path <- system.file("extdata", "internal_drug_data.rda", package = "AutoML")
  internal_drug_data <- get(load(file_path))

  # Select appropriate training data
  training_data <- switch(database,
                          "CTRP2" = if (TPM) {
                            internal_drug_data$training_data$CTRP2_TPM_Expr
                          } else {
                            internal_drug_data$training_data$CTRP2_RPKM_Expr
                          },
                          "GDSC1" = internal_drug_data$training_data$GDSC1_Expr,
                          "GDSC2" = internal_drug_data$training_data$GDSC2_Expr
  )

  drug_data <- internal_drug_data$drug_data[[database]]

  # Ensure matching samples between training and drug data
  if (!identical(rownames(drug_data), colnames(training_data))) {
    common <- intersect(rownames(drug_data), colnames(training_data))
    if (length(common) == 0) {
      stop("No matching samples found between training and drug data")
    }
    drug_data <- drug_data[common, , drop = FALSE]
    training_data <- training_data[, common, drop = FALSE]
  }

  # revise oncopredict output
  original_dir <- getOption("oncoPredict.outputDirectory")
  options(oncoPredict.outputDirectory = output_dir)
  on.exit(options(oncoPredict.outputDirectory = original_dir))

  # Calculate drug sensitivity predictions
  results <- tryCatch(
    oncoPredict::calcPhenotype(
      trainingExprData = training_data,
      trainingPtype = drug_data,
      testExprData = as.matrix(test.data),
      batchCorrect = "none",
      powerTransformPhenotype = TRUE,
      removeLowVaryingGenes = 0.2,
      minNumSamples = 10,
      printOutput = TRUE,
      removeLowVaringGenesFrom = 'rawData'
    ),
    error = function(e) {
      stop("Drug sensitivity prediction failed: ", e$message)
    }
  )

  # output file
  original_output <- file.path(output_dir, "DrugPredictions.csv")
  if (file.exists(original_output)) {
    file.rename(original_output, output_path)
    message("Results saved to: ", output_path)
  } else {
    warning("Prediction completed but output file not found at: ", original_output)
  }

  # return result
  invisible(results)
}
