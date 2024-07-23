ML.boruta<-function(est_dd,
                    seed=123456
    ){
    ##### 3.Boruta ###########
    set.seed(seed)

    message("--- 3.Boruta  ---")
    boruta <- Boruta(
      x = as.matrix(est_dd[,-1]), y = as.factor(est_dd[, c(2)]), pValue = 0.01, mcAdj = T,
      maxRuns = 11
    )
    boruta.imp <- function(x) {
      imp <- reshape2::melt(x$ImpHistory, na.rm = T)[, -1]
      colnames(imp) <- c("Variable", "Importance")
      imp <- imp[is.finite(imp$Importance), ]

      variableGrp <- data.frame(
        Variable = names(x$finalDecision),
        finalDecision = x$finalDecision
      )

      showGrp <- data.frame(
        Variable = c("shadowMax", "shadowMean", "shadowMin"),
        finalDecision = c("shadowMax", "shadowMean", "shadowMin")
      )

      variableGrp <- rbind(variableGrp, showGrp)

      boruta.variable.imp <- merge(imp, variableGrp, all.x = T)

      sortedVariable <- boruta.variable.imp %>%
        group_by(Variable) %>%
        summarise(median = median(Importance)) %>%
        arrange(median)
      sortedVariable <- as.vector(sortedVariable$Variable)


      boruta.variable.imp$Variable <- factor(boruta.variable.imp$Variable, levels = sortedVariable)

      invisible(boruta.variable.imp)
    }

    boruta.variable.imp <- boruta.imp(boruta)
    # head(boruta.variable.imp)
    boruta.finalVars <- data.frame(Item = getSelectedAttributes(boruta, withTentative = T), Type = "Boruta")

    result <- data.frame(
      method = c(rep("Boruta", length(boruta.finalVars$Item))),
      selected.fea = boruta.finalVars$Item
    )

    selected.feature <- rbind(selected.feature, result)
    write.table(result,file="3.boruta_select_features.csv",sep=",",row.names = F)
    return(selected.feature)
}
