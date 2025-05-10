#' Calculate risk Score for New Dataset
#'
#' Computes risk score using trained survival analysis models
#'
#' @param newdata Dataframe containing time-to-event data with 'time' and 'status' columns
#' @param model Trained survival analysis model object
#' @param model_name Model specification. Supported models:
#'   "Lasso", "Ridge", "Enet", "RFRSF", "GBM", "CoxBoost", "plsRcox",
#'   "XGBoost", "BlackBoost", "DeepHit", "DeepSurv", "SurvivalSVM"
#' @return risk score dataframe
#' @export
#' @examples
#' \donttest{
#' # Requires pre-trained model
#' library(survival)
#' data(example_data)
#' cox_model <- coxph(Surv(time, status) ~ ., data = example_data)
#' pred_df <- cal_pred(
#'   newdata = example_data,
#'   model = cox_model,
#'   model_name = "CoxBoost"
#' )
#' }
cal_pred<-function(newdata,model,model_name){
  # Parameter validation ----------------------------------------------------
  required_cols <- c("time", "status")
  if (!all(required_cols %in% colnames(newdata))) {
    stop("Input data must contain 'time' and 'status' columns")
  }

  valid_models <- c("Lasso", "Ridge", "Enet", "RFRSF", "GBM", "CoxBoost",
                    "plsRcox", "XGBoost", "BlackBoost", "DeepHit",
                    "DeepSurv", "SurvivalSVM","GLMBoost","SuperPC")
  model_name <- match.arg(model_name, choices = valid_models)

  if (is.null(model)) stop("Must provide valid trained model")

  # loading package---------------------------------------------------------
  if (!requireNamespace(c("survcomp","survival"), quietly = TRUE)) {
    stop("Package survcomp required: install.packages('survcomp')")
  }

  pred_coxb<-c()
  pred_df<-data.frame()
  ##
  if (model_name %in% c("Lasso","Enet","Ridge")){
    pred_coxb =as.numeric(predict(model,type='link',
                                  newx=as.matrix(newdata[,3:ncol(newdata)]),
                                  s = model$lambda.min))
  } else if (model_name %in% c("RFRSF","SurvivalSVM")){
    pred_coxb=as.numeric(predict(model,newdata)$predicted)

  } else if (model_name %in% c("GBM","GLMBoost","Survreg")){
    pred_coxb=as.numeric(predict(model,newdata))

  } else if (model_name %in% c("CoxBoost","plsRcox")){
    pred_coxb <- as.numeric(predict(model,newdata=newdata[,-c(1,2)],
                                    newtime=newdata[,'time'],
                                    newstatus=newdata[,'status'],
                                    type="lp"))

  } else if (model_name %in% c("SuperPC")){
    test <- list(x=t(newdata[,-c(1,2)]),
                 y=newdata$time,
                 censoring.status=newdata$status,
                 featurenames=colnames(newdata)[-c(1,2)])

    pred_coxb <-superpc::superpc.predict(model$fit,
                                 model$train_data,
                                 test,
                                 threshold = model$cv.fit$thresholds[which.max(model$cv.fit[["scor"]][1,])],
                                 n.components = 1)$v.pred
    pred_coxb <- as.numeric(pred_coxb)
  } else if (model_name %in% c("XGBoost")){
    pred_coxb = as.numeric(predict(model,data.matrix(newdata[,-c(1:2)])))

  } else if (model_name %in% c("BlackBoost")){
    pred_coxb <- predict(model,
                         newdata = newdata,
                         type = "link",
                         lambda = 1)
    pred_coxb<-as.numeric(pred_coxb)
  } else if (model_name %in% c("DeepHit")){
    train_fold<-data.frame(model$y,model$x)
    params=model$best_param[[1]]
    set_seed(123)
    re_model <- survivalmodels::deephit(formula = Surv(time, status) ~ ., data = train_fold,
                        frac = 0.2, activation = "relu",
                        num_nodes = c(params$nd1, params$nd2, params$nd3, params$nd4),
                        early_stopping = TRUE, batch_size = params$batch,
                        epochs = params$epochs, learning_rate = params$lr)
    pred_coxb=as.numeric(predict(re_model,newdata,type = "risk"))
  } else if (model_name %in% c("DeepSurv")){
    train_fold<-data.frame(model$y,model$x)
    params=model$best_param[[1]]
    set_seed(123)
    re_model <- survivalmodels::deepsurv(formula = Surv(time, status) ~ ., data = train_fold,
                         frac = 0.2, activation = "relu",
                         num_nodes = c(params$nd1, params$nd2, params$nd3, params$nd4),
                         early_stopping = TRUE, batch_size = params$batch,
                         epochs = params$epochs, learning_rate = params$lr)
    pred_coxb=as.numeric(predict(re_model,newdata,type = "risk"))
  }

  pred_df=data.frame(row.names = row.names(newdata),
                     "time"=newdata$time,
                     "status"=newdata$status,
                     "pred"=pred_coxb)
  return(pred_df)
}
