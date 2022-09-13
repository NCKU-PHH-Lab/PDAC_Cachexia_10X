GSEA_ExtractSubType = function(GSEA_Large.Sum.TOP.S,
                               KeyWordSet.lt = list(Mode = "KWSet", KW = c("CD4+T","CD8+T","T")), ## Mode =c("KWSet","Grep")
                               GSEA_Color = list(high = "#ef476f",mid = "white",low = "#0077b6")
){

  ## Extract Cell Type
  if(KeyWordSet.lt[["Mode"]]=="KWSet"){
    GSEA_Sub.df <- GSEA_Large.Sum.TOP.S[GSEA_Large.Sum.TOP.S$PhenoType %in% KeyWordSet.lt[["KW"]],]

  }else{
    GSEA_Sub.df <- GSEA_Large.Sum.TOP.S[grep(KeyWordSet.lt[["KW"]],GSEA_Large.Sum.TOP.S$PhenoType),]

  }


  # ## T Cell
  # # GSEA_Sub.df <- GSEA_Large.Sum.TOP.S[grep("T",GSEA_Large.Sum.TOP.S$PhenoType),]
  # GSEA_Sub.df <- GSEA_Large.Sum.TOP.S[GSEA_Large.Sum.TOP.S$PhenoType %in% c("CD4+T","CD8+T","T"),]

  BBPlot_T <- ggplot(GSEA_Sub.df,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) +
    geom_point() +scale_size_area(max_size = 7)+
    scale_colour_gradient2(low = GSEA_Color.lt[["low"]], mid = GSEA_Color.lt[["mid"]], high = GSEA_Color.lt[["high"]],
                           guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

  BBPlot_T

  BBPlot_TB <- BBPlot_T %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,OL_Thick = 1.5,
                                           XtextSize=15 ,  YtextSize=10,AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)
  # BBPlot_TB <- BBPlot_TB +theme(axis.title.y=element_blank(),
  #                  axis.text.y=element_blank(),
  #                  axis.ticks.y=element_blank())
  BBPlot_TB

  BBPlot_TB1 <- BBPlot_TB %>%
    insert_left(GSEA_ggplot_SPA.lt[["Y_Order"]],width = 0.2)
  BBPlot_TB1


  pdf(file = paste0(Subfolder.Path,"/PBMC_GSEA/PBMC_GSEA_Bubble_SPA_SubType_T.pdf"),width = 17, height = 7 )
  BBPlot_TB
  BBPlot_TB1
  dev.off()

  return(GSEA_BBPlot.lt)
}
