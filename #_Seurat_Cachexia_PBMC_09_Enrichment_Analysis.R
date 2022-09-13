##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)

##### ***Load data if necessary*** #####
# load("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-09-09_PBMC_Main/08_2_Find_CCmarker_in_different_Cell_type_and_VolcanoPlot(SPA).RData")
load("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-09-09_PBMC_Main/08_2_Find_CCmarker_in_different_Cell_type_and_VolcanoPlot(SPA).RData")

##### ***Current path and new folder setting if necessary*** #####
Save.Path = paste0(Sys.Date(),"_","PBMC_GSEA")
dir.create(Save.Path)
SampleType = "PBMC"

## Create folder
Subfolder.Path = paste0(Save.Path,"/",SampleType,"_GSEA")
if (!dir.exists(Subfolder.Path)){
  dir.create(Subfolder.Path)
}


##### 09_0 GSEA Analysis (Geneset Prepare) #####
  # https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html
  # install # https://bioconductor.org/packages/release/bioc/html/GSEABase.html
  library(tidyverse)
  library(fgsea)
  source("FUN_GSEA_LargeGeneSet.R")
  source("FUN_HSsymbol2MMsymbol.R")
  source("FUN_GSEA_ggplot.R")

  #### Load Geneset ####
  ## Geneset from GSEA
  # Pathway.all <- read.delim(paste0(getwd(),"/Pathway.all.v7.4.symbols.gmt"),header = F)
  Pathway.all.MM <- read.delim2(paste0(getwd(),"/GSEA_Geneset/Customized_GSEAGenesets_Pathway3D_Mm_ComB.gmt"),
                                col.names = 1:max(count.fields(paste0(getwd(),"/GSEA_Geneset/Customized_GSEAGenesets_Pathway3D_Mm_ComB.gmt"))),
                                header = F,sep = "\t")

  #### Genename convert ####


  #### Update gene name ####

# ##### save.image #####
# save.image(paste0(Subfolder.Path,"/09_0_GSEA_Analysis_(Geneset Prepare).RData"))

##### 09_1 GSEA Analysis (SPA) #####
## Create folder
dir.create(paste0(Subfolder.Path))


## GSEA analysis
GSEA_Large <- list()
GSEA_Large.df <- as.data.frame(matrix(nrow=0,ncol=10))
colnames(GSEA_Large.df) <- c("GeneType","PhenoType","pathway","pval","padj","log2err","ES", "NES" ,"size","leadingEdge")
GSEA_Large.df.TOP <- GSEA_Large.df

pdf(file = paste0(Subfolder.Path, "/PBMC_GSEA_SPA.pdf"),width = 15, height = 7 )

for(i in 1:length(CellType.list)){

  gseaDat <- CCMarker_SPA.lt[[paste0(CellType.list[i])]][["CCMarker.All"]]
  gseaDat <- data.frame(row.names(gseaDat),gseaDat)
  colnames(gseaDat)[[1]] <- c("Gene")
  ranks <- gseaDat$avg_log2FC
  names(ranks) <- gseaDat$Gene
  # head(ranks)
  # barplot(sort(ranks, decreasing = T))

  GSEA_Large.Output <- FUN_GSEA_LargeGeneSet(ranks,Pathway.all.MM,10)

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
  GSEA_Large[[i]] <- Sum
  names(GSEA_Large)[[i]] <- paste0(CellType.list[i])

  fgseaRes2 <- data.frame(paste0(CellType.list[i]),fgseaRes)
  colnames(fgseaRes2)[[1]] <- c("PhenoType")
  GSEA_Large.df <- rbind(GSEA_Large.df,fgseaRes2 )

  topPathways2 <- data.frame(paste0(CellType.list[i]),topPathways)
  colnames(topPathways2)[[1]] <- c("PhenoType")
  GSEA_Large.df.TOP <- rbind(GSEA_Large.df.TOP, topPathways2)

  rm(gseaDat,ranks,pathwaysH,fgseaRes,fgseaRes2,plotEnrichmen,Sum,topPathways,topPathways2)

}

dev.off()

## GSEA_Large.Sum.TOP ##
GSEA_Large.Sum.TOP <- rbind(GSEA_Large.df.TOP)
GSEA_Large.Sum.TOP <- GSEA_Large.Sum.TOP[,!colnames(GSEA_Large.Sum.TOP) %in% c("leadingEdge")]
write.table(GSEA_Large.Sum.TOP, file=paste0(Subfolder.Path,"/PBMC_GSEA_Pathway_LargeTOP_SPA.txt"),sep="\t",
            row.names=F, quote = FALSE)

##### Bubble plot #####
library(ggplot2)
library(scales)
GSEA_Color.lt = list(high = "#ef476f",mid = "white",low = "#0077b6")

GSEA_Large.Sum.TOP$PhenoType <- factor(GSEA_Large.Sum.TOP$PhenoType,
                                       levels = Cell_Type_Order.set)

GSEA_ggplot_SPA.lt <- GSEA_ggplot(GSEA_Large.Sum.TOP, NES_Th = 1.5, padj_Th = 0.01)
GSEA_Large.Sum.TOP.S <- GSEA_ggplot_SPA.lt[["GSEA_TOP.df"]]

# GSEA_Large.Sum.TOP.S <- GSEA_Large.Sum.TOP[abs(GSEA_Large.Sum.TOP$NES) > 1,]
# GSEA_Large.Sum.TOP.S <- GSEA_Large.Sum.TOP.S[abs(GSEA_Large.Sum.TOP.S$padj) < 0.05,]

# GSEA_Large.Sum.TOP.S <- GSEA_Large.Sum.TOP[abs(GSEA_Large.Sum.TOP$padj) < 0.25,]
# GSEA_Large.Sum.TOP.S <- GSEA_Large.Sum.TOP.S[abs(GSEA_Large.Sum.TOP.S$pval) < 0.05,]

pdf(file = paste0(Subfolder.Path,"/PBMC_GSEA_Bubble_SPA.pdf"),width = 17, height = 12 )
  GSEA_ggplot_SPA.lt[["BBPlot_Ori"]]
  GSEA_ggplot_SPA.lt[["BBPlot"]]
  GSEA_ggplot_SPA.lt[["BBPlot2"]]
  GSEA_ggplot_SPA.lt[["BBPlotB1"]]
  GSEA_ggplot_SPA.lt[["BBPlotB1"]]
dev.off()


##### Extract SubType #####
source("FUN_GSEA_ExtractSubType.R")

## T Cell
GSEA_T.lt <- GSEA_ExtractSubType(GSEA_Large.Sum.TOP.S,
                                 KeyWordSet.lt = list(Mode = "KWSet", KW = c("CD4+T","CD8+T","T")),
                                 OrderSet = GSEA_ggplot_SPA.lt[["Y_Order"]],
                                 GSEA_Color = GSEA_Color.lt,
                                 Save.Path = paste0(Subfolder.Path),
                                 FileName = "/PBMC_GSEA_Bubble_SPA_SubType_T.pdf")

## Mac
GSEA_Mac.lt <- GSEA_ExtractSubType(GSEA_Large.Sum.TOP.S,
                                 KeyWordSet.lt = list(Mode = "Grep", KW = c("Mac")),
                                 OrderSet = GSEA_ggplot_SPA.lt[["Y_Order"]],
                                 GSEA_Color = GSEA_Color.lt,
                                 Save.Path = paste0(Subfolder.Path),
                                 FileName = "/PBMC_GSEA_Bubble_SPA_SubType_Mac.pdf")


rm(p2,p3,BBPlotB1,BBPlotB2,BBPlotB,BBPlot_Cluster,df1.1.clust.Pheno,df1.1.clust.Pathway,
   df1.1,df1,BBPlot,BBPlot_Mac,BBPlot_MacB,BBPlot_T,BBPlot_TB)

##### save.image #####
save.image(paste0(Subfolder.Path,"/09_1_GSEA_Analysis_(SPA).RData"))

##### 09_2 GSEA Analysis (SSA_MAle) #####
GSEA_Large_Male <- list()
GSEA_Large_Male.df <- as.data.frame(matrix(nrow=0,ncol=10))
colnames(GSEA_Large_Male.df) <- c("GeneType","PhenoType","pathway","pval","padj","log2err","ES", "NES" ,"size","leadingEdge")
GSEA_Large_Male.df.TOP <- GSEA_Large_Male.df


pdf(file = paste0(Subfolder.Path, "/PBMC_GSEA_SSA_Male.pdf"),width = 15, height = 7 )

for(i in 1:length(CellType.list)){

  gseaDat <- CCMarker_Male.lt[[paste0(CellType.list[i])]][["CCMarker.All"]]
  gseaDat <- data.frame(row.names(gseaDat),gseaDat)
  colnames(gseaDat)[[1]] <- c("Gene")
  ranks <- gseaDat$avg_log2FC
  names(ranks) <- gseaDat$Gene
  # head(ranks)
  # barplot(sort(ranks, decreasing = T))


  #GSEA_Large_Male.Output <- FUN_GSEA_Large_MaleGeneSet(ranks,Pathway.all,10)
  GSEA_Large_Male.Output <- FUN_GSEA_LargeGeneSet(ranks,Pathway.all.MM,10)

  fgseaRes <- GSEA_Large_Male.Output[["fgseaRes"]]
  # head(fgseaRes[order(padj, -abs(NES)), ], n=10)

  pathwaysH <- GSEA_Large_Male.Output[["Pathway.all.list"]]

  # plot.new()
  # plotEnrichment(pathwaysH[[as.character(fgseaRes$pathway[1])]], ranks)

  topPathways <- GSEA_Large_Male.Output[["topPathways"]]

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
  GSEA_Large_Male[[i]] <- Sum
  names(GSEA_Large_Male)[[i]] <- paste0(CellType.list[i])

  fgseaRes2 <- data.frame(paste0(CellType.list[i]),fgseaRes)
  colnames(fgseaRes2)[[1]] <- c("PhenoType")
  GSEA_Large_Male.df <- rbind(GSEA_Large_Male.df,fgseaRes2 )

  topPathways2 <- data.frame(paste0(CellType.list[i]),topPathways)
  colnames(topPathways2)[[1]] <- c("PhenoType")
  GSEA_Large_Male.df.TOP <- rbind(GSEA_Large_Male.df.TOP,topPathways2 )

  rm(gseaDat,ranks,pathwaysH,fgseaRes,fgseaRes2,plotEnrichmen,Sum,topPathways,topPathways2)

}

dev.off()

## GSEA_Large_Male.Sum.TOP ##
GSEA_Large_Male.Sum.TOP <- rbind(GSEA_Large_Male.df.TOP)
GSEA_Large_Male.Sum.TOP <- GSEA_Large_Male.Sum.TOP[,!colnames(GSEA_Large_Male.Sum.TOP) %in% c("leadingEdge")]
write.table(GSEA_Large_Male.Sum.TOP, file=paste0(Subfolder.Path,"/PBMC_GSEA_Pathway_LargeTOP_SSA_Male.txt"),sep="\t",
            row.names=F, quote = FALSE)



##### Bubble plot #####
library(ggplot2)
library(scales)

GSEA_Large_Male.Sum.TOP$PhenoType <- factor(GSEA_Large_Male.Sum.TOP$PhenoType,
                                            levels = Cell_Type_Order.set)

GSEA_ggplot_SSA_Male.lt <- GSEA_ggplot(GSEA_Large_Male.Sum.TOP,NES_Th = 1.5, padj_Th = 0.01)
GSEA_Large_Male.Sum.TOP.S <- GSEA_ggplot_SSA_Male.lt[["GSEA_TOP.df"]]

pdf(file = paste0(Subfolder.Path,"/PBMC_GSEA_Bubble_SSA_Male.pdf"),width = 17, height = 12 )
  GSEA_ggplot_SSA_Male.lt[["BBPlot_Ori"]]
  GSEA_ggplot_SSA_Male.lt[["BBPlot"]]
  GSEA_ggplot_SSA_Male.lt[["BBPlot2"]]
  GSEA_ggplot_SSA_Male.lt[["BBPlotB1"]]
  GSEA_ggplot_SSA_Male.lt[["BBPlotB1"]]
dev.off()


##### Extract SubType #####
source("FUN_GSEA_ExtractSubType.R")

## T Cell
GSEA_T_Male.lt <- GSEA_ExtractSubType(GSEA_Large_Male.Sum.TOP.S,
                                 KeyWordSet.lt = list(Mode = "KWSet", KW = c("CD4+T","CD8+T","T")),
                                 OrderSet = GSEA_ggplot_SSA_Male.lt[["Y_Order"]],
                                 GSEA_Color = GSEA_Color.lt,
                                 Save.Path = paste0(Subfolder.Path),
                                 FileName = "/PBMC_GSEA_Bubble_SSA_Male_SubType_T.pdf")

## Mac
GSEA_Mac_Male.lt <- GSEA_ExtractSubType(GSEA_Large_Male.Sum.TOP.S,
                                   KeyWordSet.lt = list(Mode = "Grep", KW = c("Mac")),
                                   OrderSet = GSEA_ggplot_SSA_Male.lt[["Y_Order"]],
                                   GSEA_Color = GSEA_Color.lt,
                                   Save.Path = paste0(Subfolder.Path),
                                   FileName = "/PBMC_GSEA_Bubble_SSA_Male_SubType_Mac.pdf")

rm(p2,p3,BBPlotB1,BBPlotB2,BBPlotB,BBPlot_Cluster,df1.1.clust.Pheno,df1.1.clust.Pathway,
   df1.1,df1,BBPlot,BBPlot_Mac,BBPlot_MacB,BBPlot_T,BBPlot_TB)

##### save.image #####
save.image(paste0(Subfolder.Path,"/09_2_GSEA_Analysis_(SSA_Male).RData"))

##### 09_3 GSEA Analysis (SSA_Female) #####
GSEA_Large_Female <- list()
GSEA_Large_Female.df <- as.data.frame(matrix(nrow=0,ncol=10))
colnames(GSEA_Large_Female.df) <- c("GeneType","PhenoType","pathway","pval","padj","log2err","ES", "NES" ,"size","leadingEdge")
GSEA_Large_Female.df.TOP <- GSEA_Large_Female.df


pdf(file = paste0(Subfolder.Path, "/PBMC_GSEA_SSA_Female.pdf"),width = 15, height = 7 )

for(i in 1:length(CellType.list)){

  gseaDat <- CCMarker_Female.lt[[paste0(CellType.list[i])]][["CCMarker.All"]]
  gseaDat <- data.frame(row.names(gseaDat),gseaDat)
  colnames(gseaDat)[[1]] <- c("Gene")
  ranks <- gseaDat$avg_log2FC
  names(ranks) <- gseaDat$Gene
  # head(ranks)
  # barplot(sort(ranks, decreasing = T))


  #GSEA_Large_Female.Output <- FUN_GSEA_Large_FemaleGeneSet(ranks,Pathway.all,10)
  GSEA_Large_Female.Output <- FUN_GSEA_LargeGeneSet(ranks,Pathway.all.MM,10)

  fgseaRes <- GSEA_Large_Female.Output[["fgseaRes"]]
  # head(fgseaRes[order(padj, -abs(NES)), ], n=10)

  pathwaysH <- GSEA_Large_Female.Output[["Pathway.all.list"]]

  # plot.new()
  # plotEnrichment(pathwaysH[[as.character(fgseaRes$pathway[1])]], ranks)

  topPathways <- GSEA_Large_Female.Output[["topPathways"]]

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
  GSEA_Large_Female[[i]] <- Sum
  names(GSEA_Large_Female)[[i]] <- paste0(CellType.list[i])

  fgseaRes2 <- data.frame(paste0(CellType.list[i]),fgseaRes)
  colnames(fgseaRes2)[[1]] <- c("PhenoType")
  GSEA_Large_Female.df <- rbind(GSEA_Large_Female.df,fgseaRes2 )

  topPathways2 <- data.frame(paste0(CellType.list[i]),topPathways)
  colnames(topPathways2)[[1]] <- c("PhenoType")
  GSEA_Large_Female.df.TOP <- rbind(GSEA_Large_Female.df.TOP,topPathways2 )

  rm(gseaDat,ranks,pathwaysH,fgseaRes,fgseaRes2,plotEnrichmen,Sum,topPathways,topPathways2)

}

dev.off()

## GSEA_Large_Female.Sum.TOP ##
GSEA_Large_Female.Sum.TOP <- rbind(GSEA_Large_Female.df.TOP)
GSEA_Large_Female.Sum.TOP <- GSEA_Large_Female.Sum.TOP[,!colnames(GSEA_Large_Female.Sum.TOP) %in% c("leadingEdge")]
write.table(GSEA_Large_Female.Sum.TOP, file=paste0(Subfolder.Path,"/PBMC_GSEA_Pathway_LargeTOP_SSA_Female.txt"),sep="\t",
            row.names=F, quote = FALSE)

##### Bubble plot #####
library(ggplot2)
library(scales)

GSEA_Large_Female.Sum.TOP$PhenoType <- factor(GSEA_Large_Female.Sum.TOP$PhenoType,
                                              levels = Cell_Type_Order.set)

GSEA_ggplot_SSA_Female.lt <- GSEA_ggplot(GSEA_Large_Female.Sum.TOP,NES_Th = 1.5, padj_Th = 0.01)
GSEA_Large_Female.Sum.TOP.S <- GSEA_ggplot_SSA_Female.lt[["GSEA_TOP.df"]]

pdf(file = paste0(Subfolder.Path,"/PBMC_GSEA_Bubble_SSA_Female.pdf"),width = 17, height = 12 )
GSEA_ggplot_SSA_Female.lt[["BBPlot_Ori"]]
GSEA_ggplot_SSA_Female.lt[["BBPlot"]]
GSEA_ggplot_SSA_Female.lt[["BBPlot2"]]
GSEA_ggplot_SSA_Female.lt[["BBPlotB1"]]
GSEA_ggplot_SSA_Female.lt[["BBPlotB1"]]
dev.off()


##### Extract SubType #####
source("FUN_GSEA_ExtractSubType.R")

## T Cell
GSEA_T_Female.lt <- GSEA_ExtractSubType(GSEA_Large_Female.Sum.TOP.S,
                                      KeyWordSet.lt = list(Mode = "KWSet", KW = c("CD4+T","CD8+T","T")),
                                      OrderSet = GSEA_ggplot_SSA_Female.lt[["Y_Order"]],
                                      GSEA_Color = GSEA_Color.lt,
                                      Save.Path = paste0(Subfolder.Path),
                                      FileName = "/PBMC_GSEA_Bubble_SSA_Female_SubType_T.pdf")

## Mac
GSEA_Mac_Female.lt <- GSEA_ExtractSubType(GSEA_Large_Female.Sum.TOP.S,
                                        KeyWordSet.lt = list(Mode = "Grep", KW = c("Mac")),
                                        OrderSet = GSEA_ggplot_SSA_Female.lt[["Y_Order"]],
                                        GSEA_Color = GSEA_Color.lt,
                                        Save.Path = paste0(Subfolder.Path),
                                        FileName = "/PBMC_GSEA_Bubble_SSA_Female_SubType_Mac.pdf")

rm(p2,p3,BBPlotB1,BBPlotB2,BBPlotB,BBPlot_Cluster,df1.1.clust.Pheno,df1.1.clust.Pathway,
   df1.1,df1,BBPlot,BBPlot_Mac,BBPlot_MacB,BBPlot_T,BBPlot_TB)


##### save.image #####
save.image(paste0(Subfolder.Path,"/09_3_GSEA_Analysis_(SSA_Female).RData"))

##### 09_4 GSEA Analysis (SSA) #####
GSEA_Large_Male.df.TOP2 <- GSEA_Large_Male.df.TOP
GSEA_Large_Female.df.TOP2 <- GSEA_Large_Female.df.TOP
GSEA_Large.df.TOP2 <- GSEA_Large.df.TOP

GSEA_Large_Male.df.TOP2$PhenoType <- paste0("M_", GSEA_Large_Male.df.TOP2$PhenoType)
GSEA_Large_Female.df.TOP2$PhenoType <- paste0("F_", GSEA_Large_Female.df.TOP2$PhenoType)
GSEA_Large.df.TOP2$PhenoType <- paste0("SPA_", GSEA_Large.df.TOP2$PhenoType)


GSEA_Large_SumTOP_Sex.df <- rbind(GSEA_Large.df.TOP2,GSEA_Large_Male.df.TOP2,GSEA_Large_Female.df.TOP2)
GSEA_Large_SumTOP_Sex.df <- GSEA_Large_SumTOP_Sex.df[,!colnames(GSEA_Large_SumTOP_Sex.df) %in% c("leadingEdge")]

GSEA_ggplot_SSA.lt <- GSEA_ggplot(GSEA_Large_SumTOP_Sex.df,NES_Th = 1.5, padj_Th = 0.01)
GSEA_Large_SumTOP_Sex.df.S <- GSEA_ggplot_SSA.lt[["GSEA_TOP.df"]]

write.table(GSEA_Large_SumTOP_Sex.df, file=paste0(Subfolder.Path,"/PBMC_GSEA_Pathway_3Dataset_SSA.txt"),sep="\t",
            row.names=F, quote = FALSE)

##### Extract SubType (PBMC) #####
## SPA
GSEA_Large_SumTOP_SPA.df.S <- GSEA_Large_SumTOP_Sex.df.S[grep("SPA",GSEA_Large_SumTOP_Sex.df.S$PhenoType),]

BBPlot_SPA <- ggplot(GSEA_Large_SumTOP_SPA.df.S,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) +
  geom_point() +
  scale_size_area(max_size = 7)+
  scale_colour_gradient2(low = GSEA_Color.lt[["low"]], mid = GSEA_Color.lt[["mid"]], high = GSEA_Color.lt[["high"]],
                         guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

BBPlot_SPA

BBPlot_SPAB <- BBPlot_SPA %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                                             XtextSize=15 ,  YtextSize=10,AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)
BBPlot_SPAB

BBPlot_SPAB1 <- BBPlot_SPAB %>%
  insert_left( GSEA_ggplot_SSA.lt[["Y_Order"]],width = 0.2)
BBPlot_SPAB1

## Female
GSEA_Large_SumTOP_F.df.S <- GSEA_Large_SumTOP_Sex.df.S[grep("F_",GSEA_Large_SumTOP_Sex.df.S$PhenoType),]

BBPlot_F <- ggplot(GSEA_Large_SumTOP_F.df.S,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) +
  geom_point() +
  scale_size_area(max_size = 7)+
  scale_colour_gradient2(low = GSEA_Color.lt[["low"]], mid = GSEA_Color.lt[["mid"]], high = GSEA_Color.lt[["high"]],
                         guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

BBPlot_F

BBPlot_F_B <- BBPlot_F %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                                          XtextSize=15 ,  YtextSize=10,AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)
BBPlot_F_B

BBPlot_F_B1 <- BBPlot_F_B %>%
  insert_left( GSEA_ggplot_SSA.lt[["Y_Order"]],width = 0.2)
BBPlot_F_B1


## Male
GSEA_Large_SumTOP_M.df.S <- GSEA_Large_SumTOP_Sex.df.S[grep("M_",GSEA_Large_SumTOP_Sex.df.S$PhenoType),]

BBPlot_M <- ggplot(GSEA_Large_SumTOP_M.df.S,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) +
  geom_point() +
  scale_size_area(max_size = 7)+
  scale_colour_gradient2(low = GSEA_Color.lt[["low"]], mid = GSEA_Color.lt[["mid"]], high = GSEA_Color.lt[["high"]],
                         guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+ theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

BBPlot_M

BBPlot_M_B <- BBPlot_M %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                                          XtextSize=15 ,  YtextSize=10,AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)
BBPlot_M_B

BBPlot_M_B1 <- BBPlot_M_B %>%
  insert_left( GSEA_ggplot_SSA.lt[["Y_Order"]],width = 0.2)
BBPlot_M_B1


source("FUN_GSEA_ExtractSubType.R")

## T Cell
GSEA_T_Sum.lt <- GSEA_ExtractSubType(GSEA_Large_SumTOP_Sex.df.S,
                                     KeyWordSet.lt = list(Mode = "KWSet", KW = c("F_CD4+T","F_CD8+T","F_T",
                                                                                 "M_CD4+T","M_CD8+T","M_T",
                                                                                 "SPA_CD4+T","SPA_CD8+T","SPA_T")),
                                     OrderSet = GSEA_ggplot_SSA.lt[["Y_Order"]],
                                     GSEA_Color = GSEA_Color.lt,
                                     Save.Path = paste0(Subfolder.Path),
                                     FileName = "/PBMC_GSEA_Bubble_Sum_SubType_T.pdf")

BBPlot_T_B <- GSEA_T_Sum.lt[["BBPlot_SubB"]]
BBPlot_T_B1 <- GSEA_T_Sum.lt[["BBPlot_SubB_Sort"]]

## Mac
GSEA_Mac_Sum.lt <- GSEA_ExtractSubType(GSEA_Large_SumTOP_Sex.df.S,
                                       KeyWordSet.lt = list(Mode = "Grep", KW = c("Mac")),
                                       OrderSet = GSEA_ggplot_SSA.lt[["Y_Order"]],
                                       GSEA_Color = GSEA_Color.lt,
                                       Save.Path = paste0(Subfolder.Path),
                                       FileName = "/PBMC_GSEA_Bubble_Sum_SubType_Mac.pdf")
BBPlot_Mac_B <- GSEA_Mac_Sum.lt[["BBPlot_SubB"]]
BBPlot_Mac_B1 <- GSEA_Mac_Sum.lt[["BBPlot_SubB_Sort"]]

## Mast
GSEA_Mast_Sum.lt <- GSEA_ExtractSubType(GSEA_Large_SumTOP_Sex.df.S,
                                        KeyWordSet.lt = list(Mode = "Grep", KW = c("Mast")),
                                        OrderSet = GSEA_ggplot_SSA.lt[["Y_Order"]],
                                        GSEA_Color = GSEA_Color.lt,
                                        Save.Path = paste0(Subfolder.Path),
                                        FileName = "/PBMC_GSEA_Bubble_Sum_SubType_Mast.pdf")

BBPlot_Mast_B <- GSEA_Mast_Sum.lt[["BBPlot_SubB"]]
BBPlot_Mast_B1  <- GSEA_Mast_Sum.lt[["BBPlot_SubB_Sort"]]

## Neu
GSEA_Neu_Sum.lt <- GSEA_ExtractSubType(GSEA_Large_SumTOP_Sex.df.S,
                                        KeyWordSet.lt = list(Mode = "Grep", KW = c("Neu")),
                                        OrderSet = GSEA_ggplot_SSA.lt[["Y_Order"]],
                                        GSEA_Color = GSEA_Color.lt,
                                        Save.Path = paste0(Subfolder.Path),
                                        FileName = "/PBMC_GSEA_Bubble_Sum_SubType_Neu.pdf")

BBPlot_Neu_B <- GSEA_Neu_Sum.lt[["BBPlot_SubB"]]
BBPlot_Neu_B1  <- GSEA_Neu_Sum.lt[["BBPlot_SubB_Sort"]]


##### Export Bubble plot #####

pdf(file = paste0(Subfolder.Path,"/PBMC_GSEA_Bubble_Sum.pdf"),width = 35, height = 17 )
GSEA_ggplot_SSA.lt[["BBPlotB1"]]
  BBPlotB1
  BBPlot_SPAB
  BBPlot_SPAB1
  BBPlot_F_B
  BBPlot_F_B1
  BBPlot_M_B
  BBPlot_M_B1
  BBPlot_T_B
  BBPlot_T_B1
  BBPlot_Mac_B
  BBPlot_Mac_B1
  BBPlot_Mast_B
  BBPlot_Mast_B1
  BBPlot_Neu_B
  BBPlot_Neu_B1
dev.off()



##### Bubble plot #####
GSEA_Large_SumTOP_Sex.df.S <- GSEA_Large_SumTOP_Sex.df[abs(GSEA_Large_SumTOP_Sex.df$NES) > 1,]
GSEA_Large_SumTOP_Sex.df.S <- GSEA_Large_SumTOP_Sex.df.S[abs(GSEA_Large_SumTOP_Sex.df.S$padj) < 0.05,]
# GSEA_Large_SumTOP_Sex.df.S <- GSEA_Large_SumTOP_Sex.df[abs(GSEA_Large_SumTOP_Sex.df$padj) < 0.25,]
# GSEA_Large_SumTOP_Sex.df.S <- GSEA_Large_SumTOP_Sex.df.S[abs(GSEA_Large_SumTOP_Sex.df.S$pval) < 0.05,]
library(ggplot2)
library(scales)

pdf(file = paste0(Subfolder.Path,"/PBMC_GSEA_Bubble_SSA.pdf"),width = 17, height = 20 )


ggplot(GSEA_Large_SumTOP_Sex.df.S,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) +
  geom_point() + scale_colour_gradient2(low = "#0077b6", mid = "white", high = "#ef476f",
                                        guide = "colourbar",midpoint = 0)
ggplot(GSEA_Large_SumTOP_Sex.df.S,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) +
  geom_point() + scale_colour_gradient2(low = "#04873f", mid = "white", high = "#e3672d",
                                        guide = "colourbar",midpoint = 0)


BBPlot <- ggplot(GSEA_Large_SumTOP_Sex.df.S,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) +
  geom_point() + scale_colour_gradient2(low = "#ef476f", mid = "white", high = "#0077b6",
                                        guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")

BBPlot %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                          XtextSize=10,  YtextSize=5,AxisTitleSize=1, AspRat=1, XaThick=0.8, YaThick=0.8 ,OL_Thick = 1.5)

dev.off()


rm(p2,p3,BBPlotB1,BBPlotB2,BBPlotB,BBPlot_Cluster,df1.1.clust.Pheno,df1.1.clust.Pathway,
   df1.1,df1,BBPlot,BBPlot_Mac,BBPlot_MacB,BBPlot_T,BBPlot_TB)

##### save.image #####
save.image(paste0(Subfolder.Path,"/09_4_GSEA_Analysis_(SSA).RData"))

