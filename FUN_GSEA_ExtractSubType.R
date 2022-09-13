GSEA_ExtractSubType = function(GSEA_Large.Sum.TOP.S,
                               KeyWordSet.lt = list(Mode = "KWSet", KW = c("CD4+T","CD8+T","T")), ## Mode =c("KWSet","Grep")
                               GSEA_Color = list(high = "#ef476f",mid = "white",low = "#0077b6",
                               Save.Path = "")
){

  ## Extract Cell Type
  if(KeyWordSet.lt[["Mode"]]=="KWSet"){
    GSEA_Sub.df <- GSEA_Large.Sum.TOP.S[GSEA_Large.Sum.TOP.S$PhenoType %in% KeyWordSet.lt[["KW"]],]

  }else{
    GSEA_Sub.df <- GSEA_Large.Sum.TOP.S[grep(KeyWordSet.lt[["KW"]],GSEA_Large.Sum.TOP.S$PhenoType),]

  }

  ## Bubble plot
  BBPlot_Sub <- ggplot(GSEA_Sub.df,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) +
              geom_point() +scale_size_area(max_size = 7)+
              scale_colour_gradient2(low = GSEA_Color.lt[["low"]], mid = GSEA_Color.lt[["mid"]], high = GSEA_Color.lt[["high"]],
                                     guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ theme_bw()+
              theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

  BBPlot_SubB <- BBPlot_Sub %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,OL_Thick = 1.5,
                                               XtextSize=15 ,  YtextSize=10,AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)
  # BBPlot_SubB <- BBPlot_SubB +theme(axis.title.y=element_blank(),
  #                  axis.text.y=element_blank(),
  #                  axis.ticks.y=element_blank())

  BBPlot_SubB

  ## Sort the plot
  BBPlot_SubB_Sort <- BBPlot_SubB %>%
                insert_left(GSEA_ggplot_SPA.lt[["Y_Order"]],width = 0.2)
  BBPlot_SubB_Sort

  ## Export pdf
  pdf(file = paste0(Save.Path,"/PBMC_GSEA_Bubble_SPA_SubType_T.pdf"),width = 17, height = 7 )
    BBPlot_SubB
    BBPlot_SubB_Sort
  dev.off()

  ## Output
  OUTPUT <- list(GSEA_Sub.df, BBPlot_SubB, BBPlot_SubB_Sort)
  names(OUTPUT) <- c("GSEA_Sub.df", "BBPlot_SubB", "BBPlot_SubB_Sort")

  return(OUTPUT)

}
