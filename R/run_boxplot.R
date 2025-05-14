#' boxplot of differential expression between two groups
#'
#' @param expr_matrix expr dataframe col features rows accession
#' @param genelist a gene list
#' @param grouping grouping datafarme
#' @param orders order for group
#' @param dataset_name  title
#' @param method "t.test" or "wilcox"
#' @param width plot width
#' @param height plot height
#' @param outdir the output directory
#' @import ggpubr
#' @import dplyr
#' @importFrom utils combn
generate_two_boxplot<-function(expr_matrix, ### expr dataframe col features rows accession
                         genelist, ### gene or geneset
                         grouping, # grouping information col1:accession col2:group
                         orders=NULL, # order for group
                         dataset_name=NULL,
                         method="t.test",
                         width=5.5,
                         height=7,
                         outdir="gene_expression/boxplot/"){

  ##
  if (is.null(dataset_name)){
    dataset_name=gene
  } else {
    dataset_name=dataset_name
  }
  ##
  if (is.null(orders)){
    orders=NULL
  } else {
    orders=orders
  }

  ## output path
  if (! dir.exists(outdir)){
      dir.create(outdir,recursive = T)
  }

  ##matrix row.names: genesymbol; colnames: samples
  ##grouping info column 1: accession; column 2: group;
  data<-expr_matrix
  colnames(grouping)[1]<-"accession"
  data<- data[grouping$accession, , drop = FALSE]

  ## group infomation
  group_info<-data.frame(row.names = grouping$accession,group=grouping$group)
  group_info<-na.omit(group_info)
  ## combine with group_info
  expr <- data[, genelist, drop = FALSE]
  expr<-expr[row.names(group_info), , drop=FALSE]
  expr1<-cbind(expr,group_info)
  ## colors
  colors = c("#d6d69b","#95CADA","#193E8F","#E26844","#9fc3d5",
             "#9497C5","#2a347a","#FFCC77",
             "#95A54A","#CC6E0E")
  ## groups

  if (is.null(orders)) {
    grouplist<-c(unique(expr1$group))
    print("The dataset include group")
    print(unlist(grouplist))
    pairwise_comparisons <- combn(unique(grouplist), 2, simplify = FALSE)
    ## boxplot
    for(gene in genelist){
      p<-ggviolin(data=expr1,
                  x="group",
                  y=gene,
                  xlab = "",
                  fill="group",
                  color="white",
                  width = 0.8,
                  alpha = 1,
                  ylab=paste0(gene),
                  legend="none",
                  palette ="Paired",
                  title=dataset_name)
      p<-p+
        geom_boxplot(aes(x = .data$group, y = expr1[[gene]]),
                     fill = "white",
                     color="black",
                     width = 0.07,
                     size = 0.2,
                     outlier.shape = NA) +
        theme_bw() +
        theme(
          legend.position = "none",
          panel.grid = element_blank(),
          title = element_text(colour = "black", size = 10),
          axis.title.x = element_text(colour = "black", size = 10),
          axis.title.y = element_text(colour = "black", size = 10),
          axis.text.x = element_text(colour = "black", size = 10),
          axis.text.y = element_text(colour = "black", size = 8)
        )
      p<-p+stat_compare_means(comparisons = pairwise_comparisons,size = 3,method=method)
      p
      ggsave(p,filename=paste0(outdir,"/",gene,"_",dataset_name,".jpg"),width = width,height = height,
             units ="cm",dpi=600)
    }

  } else {
      grouplist<-c(orders)
      print("The dataset include group")
      print(unlist(grouplist))
      pairwise_comparisons <- combn(unique(grouplist), 2, simplify = FALSE)
      expr2<-expr1[expr1$group %in% grouplist,]
      ## boxplot
      for(gene in genelist){
        p<-ggviolin(data=expr2,
                    x="group",
                    y=gene,
                    xlab = "",
                    fill="group",
                    color="white",
                    order=orders,
                    width = 0.8,
                    alpha = 1,
                    ylab=paste0(gene),
                    legend="none",
                    palette ="Paired",
                    title=dataset_name)
        p<-p+
          geom_boxplot(aes(x = .data$group, y = expr2[[gene]]),
                       fill = "white",
                       color="black",
                       width = 0.07,
                       size = 0.2,
                       outlier.shape = NA) +
          theme_bw() +
          theme(
            legend.position = "none",
            panel.grid = element_blank(),
            title = element_text(colour = "black", size = 10),
            axis.title.x = element_text(colour = "black", size = 10),
            axis.title.y = element_text(colour = "black", size = 10),
            axis.text.x = element_text(colour = "black", size = 10),
            axis.text.y = element_text(colour = "black", size = 8)
          )
        p<-p+stat_compare_means(comparisons = pairwise_comparisons,size = 3,method=method)
        p
        ggsave(p,filename=paste0(outdir,"/",gene,"_",dataset_name,".jpg"),width = width,height = height,
               units ="cm",dpi=600)
    }
    }

}
