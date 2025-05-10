#' Perform Multivariable Cox Regression Analysis
#'
#' Conducts survival analysis with multiple covariates including gene expression data
#'
#' @param data Dataframe containing survival time (time), status (status), covariates and gene expression values
#' @param features Vector of covariate names to include in analysis
#' @param gene Name of target gene to analyze
#' @param dataset_name Identifier for the dataset (used in output files)
#' @param outdir Output directory for results
#' @param cut_type Discretization method for gene expression ("mean", "median", "quantile", or "best_cutoff")
#' @param plot_height Plot height in units of 400 pixels (NULL for default)
#' @param plot_width Plot width in units of 300 pixels (NULL for default)
#' @return Dataframe containing Cox regression results
#' @export
#'
#' @examples
#' \dontrun{
#' results <- generate_multicox_analysis(
#'   data = clinical_data,
#'   features = c("age", "gender", "stage"),
#'   gene = "TP53",
#'   dataset_name = "TCGA_LUAD",
#'   outdir = "./results"
#' )
#' }


generate_multicox_analysis <- function(data, features, gene, dataset_name, outdir,
                                  cut_type = NULL, plot_height = NULL,plot_width=NULL) {
  #loading package
  library(survival)
  library(survminer)
  library(dplyr)

  # log
  log_file <- file.path(outdir, "analysis_errors.log")
  if (!file.exists(log_file)) file.create(log_file)

  # main function ------------------------------------------------------------

  #' check the input
  validate_input <- function(data, features, gene) {
    required_cols <- c("time", "status", features, gene)
    missing_cols <- setdiff(required_cols, colnames(data))
    if (length(missing_cols) > 0) {
      stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
    }
  }

  #' preprocessing of continual or discrete
  preprocess_data <- function(data, gene, cut_type) {
    if (!is.null(cut_type)) {
      ##best cutpoint
      cutpoint <- surv_cutpoint(data,time = "time",event = "status",variables = gene)
      ##other
      cutoff <- switch(
        cut_type,
        "mean" = mean(data[[gene]], na.rm = TRUE),
        "median" = median(data[[gene]], na.rm = TRUE),
        "quantile" = quantile(data[[gene]], probs = 0.75, na.rm = TRUE),
        "best_cutoff"=as.numeric(cutpoint$cutpoint$cutpoint)
      )
      data[[gene]] <- factor(
        ifelse(data[[gene]] > cutoff, "High", "Low"),
        levels = c("Low", "High")
      )
    }
    return(data)
  }

  #' checking separation
  check_separation <- function(data, gene, features) {
    # categorical variables
    if (!is.numeric(data[[gene]])) {
      tab <- table(data[[gene]], data$status)
      if (any(rowSums(tab) == 0)) {
        stop(paste("Complete separation detected in", gene))
      }
    }

    if (!is.numeric(data[[gene]])) {
      tab <- table(data[[gene]], data$status)
      if (any(rowSums(tab) == 0)) {
        stop(paste("Complete separation detected in", gene))
      }
    }

    # continuous
    for (var in features) {
      if (!is.numeric(data[[var]])) {
        tab <- table(data[[var]], data$status)
        if (any(rowSums(tab) == 0)) {
          stop(paste("Complete separation detected in", gene))
        }
      } else {
        if (var(data[[var]], na.rm = TRUE) < 1e-6) {
          stop(paste("Near-zero variance in", var))
        }
      }
    }
  }

  #' modeling
  fit_robust_cox <- function(data, features, gene, log_file) {
    formula_str <- paste(
      "Surv(time, status) ~",
      paste(c(features,gene), collapse = " + ")
    )

    withCallingHandlers(
      {
        model <- coxph(as.formula(formula_str), data = data)
        if (any(is.infinite(coef(model)))) {
          stop("Infinite coefficients detected")
        }
        model
      },
      warning = function(w) {
        if (grepl("Loglik converged before variable", conditionMessage(w))) {
          log_msg <- paste(Sys.time(), "[Warn]", "Cox convergence issue:", conditionMessage(w))
          write(log_msg, log_file, append = TRUE)
          invokeRestart("muffleWarning")
        }
      }
    )
  }

  #' extract result
  extract_cox_results <- function(model, dataset_name) {
    summary(model)$coefficients %>%
      as.data.frame() %>%
      tibble::rownames_to_column("variable") %>%
      dplyr::mutate(dataset = dataset_name)
  }

  #' saving results
  save_cox_results <- function(result, dataset_name, outdir, cut_type) {
    suffix <- ifelse(is.null(cut_type), "continuous", paste0("discrete_", cut_type))
    filename <- file.path(outdir, paste0(dataset_name, "_", suffix, "_cox.csv"))
    write.csv(result, filename, row.names = FALSE)
  }

  #' generate forest plot
  generate_cox_plot <- function(model,data,dataset_name,width,
                                height,cut_type,outdir,log_file) {
    tryCatch({
      #source("R/ggforest2.R")
      p <- ggforest2(
        model,
        data = data,
        cpositions = c(0.02, 0.15, 0.35),
        main = dataset_name
      ) +
        ggplot2::theme(
          legend.position = "none",
          panel.grid = element_blank(),
          title = element_text(colour = "black", size = 12),
          axis.title = element_text(colour = "black", size = 12),
          axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 10)
        )
      print(p)
      suffix <- ifelse(is.null(cut_type), "", paste0("_", cut_type))
      output_file <- file.path(outdir, paste0(dataset_name, suffix, "_cox_plot.jpg"))
      dynamic_height <- 1 * length(features)
      height <- ifelse(is.null(height),dynamic_height, height)
      width <- ifelse(is.null(width),12,width)
      ggpubr::ggexport(
        p,
        filename = output_file,
        res = 600,
        width = width *300,
        height = height * 400
      )
    }, error = function(e) {
      log_msg <- paste(Sys.time(), "[Plot Error]", dataset_name, ":", conditionMessage(e))
      write(log_msg, log_file, append = TRUE)
      message(log_msg)
    })
  }

  #
  tryCatch({
    # parameter
    validate_input(data, features, gene)

    # preprocessing
    processed_data <- preprocess_data(data, gene, cut_type)

    # checking separation----
    check_separation(processed_data, gene, features)

    # modeling----
    cox_model <- fit_robust_cox(processed_data, features, gene, log_file)
    if (is.null(cox_model)) return(NULL)

    # result
    result <- extract_cox_results(cox_model, dataset_name)

    # output
    save_cox_results(result, dataset_name, outdir, cut_type)

    # view plot
    # forest plot----
    generate_cox_plot(cox_model, processed_data, dataset_name,plot_width,
                      plot_height,cut_type,outdir,log_file)
    return(result)
  }, error = function(e) {
    log_msg <- paste(Sys.time(), "[Error]", dataset_name, ":", conditionMessage(e))
    write(log_msg, file = log_file, append = TRUE)
    message(log_msg)
    return(NULL)
  })
}
