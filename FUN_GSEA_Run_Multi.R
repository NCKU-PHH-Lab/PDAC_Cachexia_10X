
GSEA_Run_Multi <- function(CCMarker_SPA.lt,
                           GeneSets = Pathway.all.MM, TopNum = 10,
                           Save.Path = Subfolder.Path, FileName = "/PBMC_GSEA_EnrichPlot.pdf"
                           ) {

  GSEA_Large.lt <- list()
  GSEA_Large.df <- as.data.frame(matrix(nrow=0,ncol=10))
  colnames(GSEA_Large.df) <- c("GeneType","PhenoType","pathway","pval","padj","log2err","ES", "NES" ,"size","leadingEdge")
  GSEA_Large.df.TOP <- GSEA_Large.df


  pdf(file = paste0(Save.Path, FileName),width = 15, height = 7 )

  for(i in 1:length(CellType.list)){

    gseaDat <- CCMarker_SPA.lt[[paste0(CellType.list[i])]][["CCMarker.All"]]
    gseaDat <- data.frame(row.names(gseaDat),gseaDat)
    colnames(gseaDat)[[1]] <- c("Gene")
    ranks <- gseaDat$avg_log2FC
    names(ranks) <- gseaDat$Gene
    # head(ranks)
    # barplot(sort(ranks, decreasing = T))

    GSEA_Large.Output <- GSEA_Run_LargeGeneSet(ranks,GeneSets,TopNum)

    fgseaRes <- GSEA_Large.Output[["fgseaRes"]]
    # head(fgseaRes[order(padj, -abs(NES)), ], n=10)

    pathwaysH <- GSEA_Large.Output[["Pathway.all.list"]]

    # plot.new()
    # plotEnrichment(pathwaysH[[as.character(fgseaRes$pathway[1])]], ranks)

    topPathways <- GSEA_Large.Output[["topPathways"]]

    library(ggplot2)
    plot.new()
    plotGseaTable(pathwaysH[topPathways$pathway],
                  ranks,
                  fgseaRes,
                  gseaParam = 0.5) + title( paste0("PBMC.",CellType.list[i]), adj = 0, line =3)

    plotEnrichment_Pos1 <- plotEnrichment(pathwaysH[[as.character(topPathways[1,1])]], ranks)+ labs(title= paste0("PBMC.",CellType.list[i],": ",as.character(topPathways[1,1])))
    #plotEnrichment_Pos1
    plotEnrichment_Neg1 <- plotEnrichment(pathwaysH[[as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])]], ranks)+ labs(title= paste0("PBMC.",CellType.list[i],": ",as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])))
    #plotEnrichment_Neg1

    Sum <- list(gseaDat,ranks,pathwaysH,fgseaRes,plotEnrichment_Pos1,plotEnrichment_Neg1)
    names(Sum) <- c("gseaDat","ranks","pathwaysH","fgseaRes","plotEnrichment_Pos1","plotEnrichment_Neg1")
    GSEA_Large.lt[[i]] <- Sum
    names(GSEA_Large.lt)[[i]] <- paste0(CellType.list[i])

    fgseaRes2 <- data.frame(paste0(CellType.list[i]),fgseaRes)
    colnames(fgseaRes2)[[1]] <- c("PhenoType")
    GSEA_Large.df <- rbind(GSEA_Large.df,fgseaRes2 )

    topPathways2 <- data.frame(paste0(CellType.list[i]),topPathways)
    colnames(topPathways2)[[1]] <- c("PhenoType")
    GSEA_Large.df.TOP <- rbind(GSEA_Large.df.TOP, topPathways2)

    rm(gseaDat,ranks,pathwaysH,fgseaRes,fgseaRes2,plotEnrichmen,Sum,topPathways,topPathways2)

  }

  dev.off()

  #####------------------------ Output ------------------------ #####
  OUTPUT <- list(GSEA_Large.lt, GSEA_Large.df, GSEA_Large.df.TOP)
  names(OUTPUT) <- c("GSEA_Large.lt","GSEA_Large.df","GSEA_Large.df.TOP")

return(OUTPUT)
}


