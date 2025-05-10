generat_rfrsf_importance<-function(fit,prefix,outdir){
  importance=as.data.frame(fit$importance)
  names(importance)[1]<-"Importance"
  importance$gene<-row.names(importance)

  library(ggplot2)

  # order by importance
  importance <- importance[order(importance$Importance), ]

  # as factor
  importance$gene <- factor(importance$gene, levels = importance$gene)

  # bar plot, fill = "steelblue"
  p <- ggplot(importance, aes(x = Importance, y = gene)) +
      geom_bar(stat = "identity",fill = "steelblue") +
      geom_text(aes(label = sprintf("%.3f", Importance)),
              position = position_dodge(width = 0.9),
              hjust = -0.1,
              size = 3,
              color = "black") +
      theme_bw(base_size = 12) +
      theme(
        axis.text.y = element_text(size = 10,color="black"),
        axis.text.x = element_text(size = 10,color="black"),
        axis.title.x = element_text(size = 12,color="black"),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        panel.grid = element_blank()
      ) +
      labs(
        title = "",
        x = "Importance",
        y = ""
      )

    # print plot
    print(p)

    # saving result
    ggsave(paste0(outdir,"/",prefix,"_importance_plot.jpg"), plot = p, width =5, height = 5, dpi = 600)
    importance <- importance[order(importance$Importance,decreasing = TRUE), ]
    importance$gene <- factor(importance$gene, levels = importance$gene)
    write.csv(importance,file=paste0(outdir,"/",prefix,"_importance.csv"))
    return(importance)
}
