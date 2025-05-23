% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cal_pred.R
\name{cal_pred}
\alias{cal_pred}
\title{Calculate risk Score for New Dataset}
\usage{
cal_pred(newdata, model, model_name)
}
\arguments{
\item{newdata}{Dataframe containing time-to-event data with 'time' and 'status' columns}

\item{model}{Trained survival analysis model object}

\item{model_name}{Model specification. Supported models:
"Lasso", "Ridge", "Enet", "RFRSF", "GBM", "CoxBoost", "plsRcox",
"XGBoost", "BlackBoost", "DeepHit", "DeepSurv", "SurvivalSVM"}
}
\value{
risk score dataframe
}
\description{
Computes risk score using trained survival analysis models
}
