##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)

##### ***Load data if necessary*** #####
# load("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-09-09_SC_Main/08_2_Find_CCmarker_in_different_Cell_type_and_VolcanoPlot(SPA).RData")
load("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-09-09_SC_Main/08_2_Find_CCmarker_in_different_Cell_type_and_VolcanoPlot(SPA).RData")

##### ***Current path and new folder setting if necessary*** #####
Save.Path = paste0(Sys.Date(),"_","SC_GSEA")
dir.create(Save.Path)
SampleType = "SC"

## Create folder
Subfolder.Path = paste0(Save.Path,"/",SampleType,"_GSEA")
if (!dir.exists(Subfolder.Path)){
  dir.create(Subfolder.Path)
}

##### Load packages ####
# https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html
# install # https://bioconductor.org/packages/release/bioc/html/GSEABase.html
library(tidyverse)
library(fgsea)
library(aplot)

##### Function setting  #####
## Call function
source("FUN_GSEA_Run_LargeGeneSet.R")
source("FUN_HSsymbol2MMsymbol.R")
source("FUN_GSEA_ggplot.R")
source("FUN_GSEA_Run_Multi.R")
source("FUN_GSEA_ExtractSubType.R")

##### Parameter setting  #####
GSEATopNum = 10
NES_Th_Set = 1.5

##### 09_0 GSEA Analysis (Geneset Prepare) #####
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

##*******************************************************************************************************************##
##### 09_1 GSEA Analysis (SPA) #####
## Create folder
dir.create(paste0(Subfolder.Path))

## GSEA analysis
GSEA_SPA.lt <- GSEA_Run_Multi(CCMarker_SPA.lt,
                              GeneSets = Pathway.all.MM, TopNum = GSEATopNum,
                              Save.Path = Subfolder.Path, FileName = "/SC_GSEA_SPA_EnrichPlot.pdf")

GSEA_Large.df.TOP <- GSEA_SPA.lt[["GSEA_Large.df.TOP"]]
GSEA_Large.df.TOP <- GSEA_Large.df.TOP[,!colnames(GSEA_Large.df.TOP) %in% c("leadingEdge")]
write.table(GSEA_Large.df.TOP, file=paste0(Subfolder.Path,"/SC_GSEA_Pathway_LargeTOP_SPA.txt"),sep="\t",
            row.names=F, quote = FALSE)

##### Bubble plot #####
library(ggplot2)
library(scales)
GSEA_Color.lt = list(high = "#ef476f",mid = "white",low = "#0077b6")

GSEA_Large.df.TOP$PhenoType <- factor(GSEA_Large.df.TOP$PhenoType,
                                       levels = Cell_Type_Order.set)

GSEA_ggplot_SPA.lt <- GSEA_ggplot(GSEA_Large.df.TOP, NES_Th = NES_Th_Set, padj_Th = 0.01)
GSEA_Large.df.TOP.S <- GSEA_ggplot_SPA.lt[["GSEA_TOP.df"]]

# GSEA_Large.df.TOP.S <- GSEA_Large.df.TOP[abs(GSEA_Large.df.TOP$NES) > 1,]
# GSEA_Large.df.TOP.S <- GSEA_Large.df.TOP.S[abs(GSEA_Large.df.TOP.S$padj) < 0.05,]

# GSEA_Large.df.TOP.S <- GSEA_Large.df.TOP[abs(GSEA_Large.df.TOP$padj) < 0.25,]
# GSEA_Large.df.TOP.S <- GSEA_Large.df.TOP.S[abs(GSEA_Large.df.TOP.S$pval) < 0.05,]

pdf(file = paste0(Subfolder.Path,"/SC_GSEA_Bubble_SPA.pdf"),width = 17, height = 12 )
  GSEA_ggplot_SPA.lt[["BBPlot_Ori"]]
  GSEA_ggplot_SPA.lt[["BBPlot"]]
  GSEA_ggplot_SPA.lt[["BBPlot2"]]
  GSEA_ggplot_SPA.lt[["BBPlotB1"]]
dev.off()


##### Extract SubType #####
source("FUN_GSEA_ExtractSubType.R")

# ## T Cell
# GSEA_T.lt <- GSEA_ExtractSubType(GSEA_Large.df.TOP,
#                                  KeyWordSet.lt = list(Mode = "KWSet", KW = c("CD4+T","CD8+T","T")),
#                                  OrderSet = GSEA_ggplot_SPA.lt[["Y_Order"]],
#                                  GSEA_Color = GSEA_Color.lt,
#                                  Save.Path = paste0(Subfolder.Path),
#                                  FileName = "/SC_GSEA_Bubble_SPA_SubType_T.pdf")

## Mac
GSEA_Mac.lt <- GSEA_ExtractSubType(GSEA_Large.df.TOP.S,
                                 KeyWordSet.lt = list(Mode = "Grep", KW = c("Mac")),
                                 OrderSet = GSEA_ggplot_SPA.lt[["Y_Order"]],
                                 GSEA_Color = GSEA_Color.lt,
                                 Save.Path = paste0(Subfolder.Path),
                                 FileName = "/SC_GSEA_Bubble_SPA_SubType_Mac.pdf")

## Duc
GSEA_Duc.lt <- GSEA_ExtractSubType(GSEA_Large.df.TOP.S,
                                   KeyWordSet.lt = list(Mode = "Grep", KW = c("Duc")),
                                   OrderSet = GSEA_ggplot_SPA.lt[["Y_Order"]],
                                   GSEA_Color = GSEA_Color.lt,
                                   Save.Path = paste0(Subfolder.Path),
                                   FileName = "/SC_GSEA_Bubble_SPA_SubType_Duc.pdf")

## Fib
GSEA_Fib.lt <- GSEA_ExtractSubType(GSEA_Large.df.TOP.S,
                                   KeyWordSet.lt = list(Mode = "Grep", KW = c("Fib")),
                                   OrderSet = GSEA_ggplot_SPA.lt[["Y_Order"]],
                                   GSEA_Color = GSEA_Color.lt,
                                   Save.Path = paste0(Subfolder.Path),
                                   FileName = "/SC_GSEA_Bubble_SPA_SubType_Fib.pdf")


rm(p2,p3,BBPlotB1,BBPlotB2,BBPlotB,BBPlot_Cluster,df1.1.clust.Pheno,df1.1.clust.Pathway,
   df1.1,df1,BBPlot,BBPlot_Mac,BBPlot_MacB,BBPlot_T,BBPlot_TB)

##### save.image #####
save.image(paste0(Subfolder.Path,"/09_1_GSEA_Analysis_(SPA).RData"))

##*******************************************************************************************************************##
##### 09_2 GSEA Analysis (SSA_Male) #####
## GSEA analysis
GSEA_SSA_Male.lt <- GSEA_Run_Multi(CCMarker_Male.lt,
                              GeneSets = Pathway.all.MM, TopNum = GSEATopNum,
                              Save.Path = Subfolder.Path, FileName = "/SC_GSEA_SSA_Male_EnrichPlot.pdf")

GSEA_Large_Male.df.TOP <- GSEA_SSA_Male.lt[["GSEA_Large.df.TOP"]]
GSEA_Large_Male.df.TOP <- GSEA_Large_Male.df.TOP[,!colnames(GSEA_Large_Male.df.TOP) %in% c("leadingEdge")]
write.table(GSEA_Large_Male.df.TOP, file=paste0(Subfolder.Path,"/SC_GSEA_Pathway_LargeTOP_SSA_Male.txt"),sep="\t",
            row.names=F, quote = FALSE)

##### Bubble plot #####
library(ggplot2)
library(scales)

GSEA_Large_Male.df.TOP$PhenoType <- factor(GSEA_Large_Male.df.TOP$PhenoType,
                                            levels = Cell_Type_Order.set)

GSEA_ggplot_SSA_Male.lt <- GSEA_ggplot(GSEA_Large_Male.df.TOP,NES_Th = NES_Th_Set, padj_Th = 0.01)
GSEA_Large_Male.df.TOP.S <- GSEA_ggplot_SSA_Male.lt[["GSEA_TOP.df"]]

pdf(file = paste0(Subfolder.Path,"/SC_GSEA_Bubble_SSA_Male.pdf"),width = 17, height = 12 )
  GSEA_ggplot_SSA_Male.lt[["BBPlot_Ori"]]
  GSEA_ggplot_SSA_Male.lt[["BBPlot"]]
  GSEA_ggplot_SSA_Male.lt[["BBPlot2"]]
  GSEA_ggplot_SSA_Male.lt[["BBPlotB1"]]
dev.off()

##### Extract SubType #####
source("FUN_GSEA_ExtractSubType.R")

# ## T Cell
# GSEA_T_Male.lt <- GSEA_ExtractSubType(GSEA_Large_Male.df.TOP.S,
#                                  KeyWordSet.lt = list(Mode = "KWSet", KW = c("CD4+T","CD8+T","T")),
#                                  OrderSet = GSEA_ggplot_SSA_Male.lt[["Y_Order"]],
#                                  GSEA_Color = GSEA_Color.lt,
#                                  Save.Path = paste0(Subfolder.Path),
#                                  FileName = "/SC_GSEA_Bubble_SSA_Male_SubType_T.pdf")

## Mac
GSEA_Mac_Male.lt <- GSEA_ExtractSubType(GSEA_Large_Male.df.TOP.S,
                                   KeyWordSet.lt = list(Mode = "Grep", KW = c("Mac")),
                                   OrderSet = GSEA_ggplot_SSA_Male.lt[["Y_Order"]],
                                   GSEA_Color = GSEA_Color.lt,
                                   Save.Path = paste0(Subfolder.Path),
                                   FileName = "/SC_GSEA_Bubble_SSA_Male_SubType_Mac.pdf")

## Duc
GSEA_Duc_Male.lt <- GSEA_ExtractSubType(GSEA_Large_Male.df.TOP.S,
                                   KeyWordSet.lt = list(Mode = "Grep", KW = c("Duc")),
                                   OrderSet = GSEA_ggplot_SPA.lt[["Y_Order"]],
                                   GSEA_Color = GSEA_Color.lt,
                                   Save.Path = paste0(Subfolder.Path),
                                   FileName = "/SC_GSEA_Bubble_SSA_Male_SubType_Duc.pdf")

## Fib
GSEA_Fib_Male.lt <- GSEA_ExtractSubType(GSEA_Large_Male.df.TOP.S,
                                   KeyWordSet.lt = list(Mode = "Grep", KW = c("Fib")),
                                   OrderSet = GSEA_ggplot_SPA.lt[["Y_Order"]],
                                   GSEA_Color = GSEA_Color.lt,
                                   Save.Path = paste0(Subfolder.Path),
                                   FileName = "/SC_GSEA_Bubble_SSA_Male_SubType_Fib.pdf")

rm(p2,p3,BBPlotB1,BBPlotB2,BBPlotB,BBPlot_Cluster,df1.1.clust.Pheno,df1.1.clust.Pathway,
   df1.1,df1,BBPlot,BBPlot_Mac,BBPlot_MacB,BBPlot_T,BBPlot_TB)

##### save.image #####
save.image(paste0(Subfolder.Path,"/09_2_GSEA_Analysis_(SSA_Male).RData"))

##*******************************************************************************************************************##
##### 09_3 GSEA Analysis (SSA_Female) #####
## GSEA analysis
GSEA_SSA_Female.lt <- GSEA_Run_Multi(CCMarker_Female.lt,
                                   GeneSets = Pathway.all.MM, TopNum = GSEATopNum,
                                   Save.Path = Subfolder.Path, FileName = "/SC_GSEA_SSA_Female_EnrichPlot.pdf")

GSEA_Large_Female.df.TOP <- GSEA_SSA_Female.lt[["GSEA_Large.df.TOP"]]
GSEA_Large_Female.df.TOP <- GSEA_Large_Female.df.TOP[,!colnames(GSEA_Large_Female.df.TOP) %in% c("leadingEdge")]
write.table(GSEA_Large_Female.df.TOP, file=paste0(Subfolder.Path,"/SC_GSEA_Pathway_LargeTOP_SSA_Female.txt"),sep="\t",
            row.names=F, quote = FALSE)

##### Bubble plot #####
library(ggplot2)
library(scales)

GSEA_Large_Female.df.TOP$PhenoType <- factor(GSEA_Large_Female.df.TOP$PhenoType,
                                              levels = Cell_Type_Order.set)

GSEA_ggplot_SSA_Female.lt <- GSEA_ggplot(GSEA_Large_Female.df.TOP,NES_Th = NES_Th_Set, padj_Th = 0.01)
GSEA_Large_Female.df.TOP.S <- GSEA_ggplot_SSA_Female.lt[["GSEA_TOP.df"]]

pdf(file = paste0(Subfolder.Path,"/SC_GSEA_Bubble_SSA_Female.pdf"),width = 17, height = 12 )
  GSEA_ggplot_SSA_Female.lt[["BBPlot_Ori"]]
  GSEA_ggplot_SSA_Female.lt[["BBPlot"]]
  GSEA_ggplot_SSA_Female.lt[["BBPlot2"]]
  GSEA_ggplot_SSA_Female.lt[["BBPlotB1"]]
dev.off()


##### Extract SubType #####
source("FUN_GSEA_ExtractSubType.R")

# ## T Cell
# GSEA_T_Female.lt <- GSEA_ExtractSubType(GSEA_Large_Female.df.TOP.S,
#                                       KeyWordSet.lt = list(Mode = "KWSet", KW = c("CD4+T","CD8+T","T")),
#                                       OrderSet = GSEA_ggplot_SSA_Female.lt[["Y_Order"]],
#                                       GSEA_Color = GSEA_Color.lt,
#                                       Save.Path = paste0(Subfolder.Path),
#                                       FileName = "/SC_GSEA_Bubble_SSA_Female_SubType_T.pdf")

## Mac
GSEA_Mac_Female.lt <- GSEA_ExtractSubType(GSEA_Large_Female.df.TOP.S,
                                        KeyWordSet.lt = list(Mode = "Grep", KW = c("Mac")),
                                        OrderSet = GSEA_ggplot_SSA_Female.lt[["Y_Order"]],
                                        GSEA_Color = GSEA_Color.lt,
                                        Save.Path = paste0(Subfolder.Path),
                                        FileName = "/SC_GSEA_Bubble_SSA_Female_SubType_Mac.pdf")

## Duc
GSEA_Duc_Female.lt <- GSEA_ExtractSubType(GSEA_Large_Female.df.TOP.S,
                                        KeyWordSet.lt = list(Mode = "Grep", KW = c("Duc")),
                                        OrderSet = GSEA_ggplot_SPA.lt[["Y_Order"]],
                                        GSEA_Color = GSEA_Color.lt,
                                        Save.Path = paste0(Subfolder.Path),
                                        FileName = "/SC_GSEA_Bubble_SSA_Female_SubType_Duc.pdf")

## Fib
GSEA_Fib_Female.lt <- GSEA_ExtractSubType(GSEA_Large_Female.df.TOP.S,
                                        KeyWordSet.lt = list(Mode = "Grep", KW = c("Fib")),
                                        OrderSet = GSEA_ggplot_SPA.lt[["Y_Order"]],
                                        GSEA_Color = GSEA_Color.lt,
                                        Save.Path = paste0(Subfolder.Path),
                                        FileName = "/SC_GSEA_Bubble_SSA_Female_SubType_Fib.pdf")


rm(p2,p3,BBPlotB1,BBPlotB2,BBPlotB,BBPlot_Cluster,df1.1.clust.Pheno,df1.1.clust.Pathway,
   df1.1,df1,BBPlot,BBPlot_Mac,BBPlot_MacB,BBPlot_T,BBPlot_TB)

##### save.image #####
save.image(paste0(Subfolder.Path,"/09_3_GSEA_Analysis_(SSA_Female).RData"))

##*******************************************************************************************************************##
##### 09_4 GSEA Analysis (SSA) #####
GSEA_Large_Male.df.TOP2 <- GSEA_Large_Male.df.TOP
GSEA_Large_Female.df.TOP2 <- GSEA_Large_Female.df.TOP
GSEA_Large.df.TOP2 <- GSEA_Large.df.TOP

GSEA_Large_Male.df.TOP2$PhenoType <- paste0("M_", GSEA_Large_Male.df.TOP2$PhenoType)
GSEA_Large_Female.df.TOP2$PhenoType <- paste0("F_", GSEA_Large_Female.df.TOP2$PhenoType)
GSEA_Large.df.TOP2$PhenoType <- paste0("SPA_", GSEA_Large.df.TOP2$PhenoType)


GSEA_Large_SumTOP_Sex.df <- rbind(GSEA_Large.df.TOP2,GSEA_Large_Male.df.TOP2,GSEA_Large_Female.df.TOP2)
GSEA_Large_SumTOP_Sex.df <- GSEA_Large_SumTOP_Sex.df[,!colnames(GSEA_Large_SumTOP_Sex.df) %in% c("leadingEdge")]

GSEA_ggplot_SSA.lt <- GSEA_ggplot(GSEA_Large_SumTOP_Sex.df,NES_Th = NES_Th_Set, padj_Th = 0.01)
GSEA_Large_SumTOP_Sex.df.S <- GSEA_ggplot_SSA.lt[["GSEA_TOP.df"]]

write.table(GSEA_Large_SumTOP_Sex.df, file=paste0(Subfolder.Path,"/SC_GSEA_Pathway_3Dataset_All.txt"),sep="\t",
            row.names=F, quote = FALSE)


##### Extract SubType #####
source("FUN_GSEA_ExtractSubType.R")

## SPA
GSEA_SPA_Sum.lt <- GSEA_ExtractSubType(GSEA_Large_SumTOP_Sex.df.S,
                                        KeyWordSet.lt = list(Mode = "Grep", KW = c("SPA")),
                                        OrderSet = GSEA_ggplot_SSA.lt[["Y_Order"]],
                                        GSEA_Color = GSEA_Color.lt,
                                        Save.Path = paste0(Subfolder.Path),
                                        FileName = "/SC_GSEA_Bubble_Sum_SubType_SPA.pdf",
                                        PDFwidth = 35, PDFheight = 17)

BBPlot_SPA_B <- GSEA_SPA_Sum.lt[["BBPlot_SubB"]]
BBPlot_SPA_B1 <- GSEA_SPA_Sum.lt[["BBPlot_SubB_Sort"]]

## Male
GSEA_Male_Sum.lt <- GSEA_ExtractSubType(GSEA_Large_SumTOP_Sex.df.S,
                                        KeyWordSet.lt = list(Mode = "Grep", KW = c("M_")),
                                        OrderSet = GSEA_ggplot_SSA.lt[["Y_Order"]],
                                        GSEA_Color = GSEA_Color.lt,
                                        Save.Path = paste0(Subfolder.Path),
                                        FileName = "/SC_GSEA_Bubble_Sum_SubType_Male.pdf",
                                        PDFwidth = 35, PDFheight = 17)

BBPlot_M_B <- GSEA_Male_Sum.lt[["BBPlot_SubB"]]
BBPlot_M_B1 <- GSEA_Male_Sum.lt[["BBPlot_SubB_Sort"]]

## Female
GSEA_Female_Sum.lt <- GSEA_ExtractSubType(GSEA_Large_SumTOP_Sex.df.S,
                                        KeyWordSet.lt = list(Mode = "Grep", KW = c("F_")),
                                        OrderSet = GSEA_ggplot_SSA.lt[["Y_Order"]],
                                        GSEA_Color = GSEA_Color.lt,
                                        Save.Path = paste0(Subfolder.Path),
                                        FileName = "/SC_GSEA_Bubble_Sum_SubType_Female.pdf",
                                        PDFwidth = 35, PDFheight = 17)

BBPlot_F_B <- GSEA_Female_Sum.lt[["BBPlot_SubB"]]
BBPlot_F_B1 <- GSEA_Female_Sum.lt[["BBPlot_SubB_Sort"]]


# ## T Cell
# GSEA_T_Sum.lt <- GSEA_ExtractSubType(GSEA_Large_SumTOP_Sex.df.S,
#                                      KeyWordSet.lt = list(Mode = "KWSet", KW = c("F_CD4+T","F_CD8+T","F_T",
#                                                                                  "M_CD4+T","M_CD8+T","M_T",
#                                                                                  "SPA_CD4+T","SPA_CD8+T","SPA_T")),
#                                      OrderSet = GSEA_ggplot_SSA.lt[["Y_Order"]],
#                                      GSEA_Color = GSEA_Color.lt,
#                                      Save.Path = paste0(Subfolder.Path),
#                                      FileName = "/SC_GSEA_Bubble_Sum_SubType_T.pdf",
#                                      PDFwidth = 35, PDFheight = 17 )

BBPlot_T_B <- GSEA_T_Sum.lt[["BBPlot_SubB"]]
BBPlot_T_B1 <- GSEA_T_Sum.lt[["BBPlot_SubB_Sort"]]

## Mac
GSEA_Mac_Sum.lt <- GSEA_ExtractSubType(GSEA_Large_SumTOP_Sex.df.S,
                                       KeyWordSet.lt = list(Mode = "Grep", KW = c("Mac")),
                                       OrderSet = GSEA_ggplot_SSA.lt[["Y_Order"]],
                                       GSEA_Color = GSEA_Color.lt,
                                       Save.Path = paste0(Subfolder.Path),
                                       FileName = "/SC_GSEA_Bubble_Sum_SubType_Mac.pdf",
                                       PDFwidth = 35, PDFheight = 17 )
BBPlot_Mac_B <- GSEA_Mac_Sum.lt[["BBPlot_SubB"]]
BBPlot_Mac_B1 <- GSEA_Mac_Sum.lt[["BBPlot_SubB_Sort"]]

## Duc
GSEA_Duc_Sum.lt <- GSEA_ExtractSubType(GSEA_Large_SumTOP_Sex.df.S,
                                        KeyWordSet.lt = list(Mode = "Grep", KW = c("Duc")),
                                        OrderSet = GSEA_ggplot_SSA.lt[["Y_Order"]],
                                        GSEA_Color = GSEA_Color.lt,
                                        Save.Path = paste0(Subfolder.Path),
                                        FileName = "/SC_GSEA_Bubble_Sum_SubType_Duc.pdf",
                                        PDFwidth = 35, PDFheight = 17 )

BBPlot_Duc_B <- GSEA_Duc_Sum.lt[["BBPlot_SubB"]]
BBPlot_Duc_B1  <- GSEA_Duc_Sum.lt[["BBPlot_SubB_Sort"]]

## Fib
GSEA_Fib_Sum.lt <- GSEA_ExtractSubType(GSEA_Large_SumTOP_Sex.df.S,
                                       KeyWordSet.lt = list(Mode = "Grep", KW = c("Fib")),
                                       OrderSet = GSEA_ggplot_SSA.lt[["Y_Order"]],
                                       GSEA_Color = GSEA_Color.lt,
                                       Save.Path = paste0(Subfolder.Path),
                                       FileName = "/SC_GSEA_Bubble_Sum_SubType_Fib.pdf",
                                       PDFwidth = 35, PDFheight = 17 )

BBPlot_Fib_B <- GSEA_Fib_Sum.lt[["BBPlot_SubB"]]
BBPlot_Fib_B1  <- GSEA_Fib_Sum.lt[["BBPlot_SubB_Sort"]]


##### Export Bubble plot #####
BBPlotB2 <- GSEA_ggplot_SSA.lt[["BBPlot"]]
BBPlotB2 <- BBPlotB2 + theme(axis.title.y=element_blank(),
                            axis.text.y=element_blank(),
                            axis.ticks.y=element_blank())
## Sort the plot
BBPlotB2 <- BBPlotB2 %>% insert_left(GSEA_ggplot_SSA.lt[["Y_Order"]], width = 0.2)
BBPlotB2

pdf(file = paste0(Subfolder.Path,"/SC_GSEA_Bubble_Sum.pdf"),width = 35, height = 17 )
  GSEA_ggplot_SSA.lt[["BBPlot"]]
  GSEA_ggplot_SSA.lt[["BBPlotB1"]]
  BBPlotB2

  BBPlot_SPA_B
  BBPlot_SPA_B1

  BBPlot_F_B
  BBPlot_F_B1
  BBPlot_M_B
  BBPlot_M_B1
  # BBPlot_T_B
  # BBPlot_T_B1
  BBPlot_Mac_B
  BBPlot_Mac_B1
  BBPlot_Fib_B
  BBPlot_Fib_B1
  BBPlot_Duc_B
  BBPlot_Duc_B1
dev.off()


##### Bubble plot #####
GSEA_Large_SumTOP_Sex.df.S <- GSEA_Large_SumTOP_Sex.df[abs(GSEA_Large_SumTOP_Sex.df$NES) > 1,]
GSEA_Large_SumTOP_Sex.df.S <- GSEA_Large_SumTOP_Sex.df.S[abs(GSEA_Large_SumTOP_Sex.df.S$padj) < 0.05,]
# GSEA_Large_SumTOP_Sex.df.S <- GSEA_Large_SumTOP_Sex.df[abs(GSEA_Large_SumTOP_Sex.df$padj) < 0.25,]
# GSEA_Large_SumTOP_Sex.df.S <- GSEA_Large_SumTOP_Sex.df.S[abs(GSEA_Large_SumTOP_Sex.df.S$pval) < 0.05,]
library(ggplot2)
library(scales)

pdf(file = paste0(Subfolder.Path,"/SC_GSEA_Bubble_SSA.pdf"),width = 17, height = 20 )


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
