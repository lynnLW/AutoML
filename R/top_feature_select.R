top_feature_select <- function(object, # output of ML.Corefeature.Prog.Screen
                              col = NULL, # color value for segment and point
                              top = NULL, # top number gene to show
                              top_select=NULL,
                              nmethod=NULL
                              ) {
  library(ggplot2)

  if (is.null(col) == T) {
    col <- c("#B09C85", "#E18727")
  } else {
    col <- col
  }

  if (is.null(nmethod) == T) {
    nmethod <- 5
  } else {
    nmethod <- nmethod
  }

  if (is.null(top) == T) {
    top <- 20
  } else {
    top <- top
  }

  tmp <- as.data.frame(table(object$selected.fea))
  tmp$Var1 <- gsub("\\.", "-", tmp$Var1)
  tmp <- tmp[order(tmp$Freq, decreasing = T), ]
  ##
  tmp$Var1 <- factor(tmp$Var1, levels = tmp$Var1)
  colnames(tmp)[2] <- "Frequence"
  jpeg(filename="3.top_feature_selection.jpg",res=600,width = 12,height = 10,units = "cm")
  p1 <- ggplot(tmp[(nrow(tmp) - top + 1):nrow(tmp), ], aes(x = Var1, y = Frequence)) +
    geom_segment(aes(x = Var1, xend = Var1, y = 0, yend = Frequence), color = col[1]) +
    geom_point(aes(size = Frequence), color = col[2], alpha = 0.5) +
    DOSE::theme_dose(10) +
    # scale_color_manual(values=dataset_col[t],name="Cohort")+
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
      # legend.position = "",
      plot.title = element_text(hjust = 0.5),
      panel.background = element_rect(fill = "white")
    ) +
    # scale_y_continuous(position = "right")+
    labs(y = "Frequence of screening", x = "", title = paste("Top", top, "selected genes")) +
    coord_flip()
  print(p1)
  dev.off()
  ##select feature
  if(top_select==T){
    final.feature<-as.character(tmp[1:top,]$Var1)
    print(paste0("Finally, ",length(final.feature)," features were selected for model training"))
  } else {
    final.feature<-as.character(tmp[tmp$Frequence>6,]$Var1)
    print(paste0("Finally, ",length(final.feature)," features were selected for model training"))
  }
  return(final.feature)
  save(final.feature,file="3.final_feature.Rdata")
}
