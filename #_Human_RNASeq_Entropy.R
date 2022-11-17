##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)

##### Load Package ######
library(tidyverse)
library(compiler)
library(ComplexHeatmap)
library(circlize)

##### Load data ######
GE.df <- read.delim("HumanRNA_12RNASeqCutoffIS1.gct.txt",header=T, stringsAsFactors = FALSE)
TarGene.df <- read.delim("HumanRNA_TargetGene.txt",header=T, stringsAsFactors = FALSE)
HumanSampleAnno.df <- read.delim("HumanRNA_SampleAnno.txt",header=T, stringsAsFactors = FALSE)

##### Order setting #####
# Order_SampleID <- HumanSampleAnno.df$SampleID %>% sort()
# Order_SampleID <- c("B1", "B3", "B4", "B6", "B7", "B8", "P1", "P2", "P3", "P4", "S5", "S6")
Order_SampleID <- c("B6", "B3", "B4", "B8", "B1", "B7", "S5", "S6", "P1", "P2", "P3", "P4")

## Order GE col
GE.df <- GE.df[,c(colnames(GE.df)[1:2],Order_SampleID)]
## Order HumanSampleAnno.df
HumanSampleAnno.df <- left_join(data.frame(SampleID=Order_SampleID),HumanSampleAnno.df)

##### Filter data #####
## Filter by TarGene
GE_TarGene.df <- GE.df[GE.df$Description %in% TarGene.df$Gene,]
row.names(GE_TarGene.df) <- GE_TarGene.df$NAME
colnames(GE_TarGene.df)


## Filter by Expression
GE_Flt.set <- ""
for (i in 1:nrow(GE_TarGene.df)){
  GE_Flt.set[i]<- any(GE_TarGene.df[i,-1:-2] > 1)
  }


GE_TarGene.df <- GE_TarGene.df[which(GE_Flt.set =="TRUE" ),]


##### Entropy base heatmap #####
## PMID: 23664764
## Ref: https://pubmed.ncbi.nlm.nih.gov/23664764/

# GE_TarGene.df$Sum <- rowSums(GE_TarGene.df[,-1:-2])
GE_MTX.df <-  GE_TarGene.df[,-1:-2]
GE_MTX_RowSum <- rowSums(GE_MTX.df)


GE_MTX_DSum.df <- as.data.frame(matrix(data = NA,nrow = nrow(GE_MTX.df),ncol = ncol(GE_MTX.df)))
colnames(GE_MTX_DSum.df) <- colnames(GE_MTX.df)
row.names(GE_MTX_DSum.df) <- row.names(GE_MTX.df)

GE_MTX_DSum.df <- GE_MTX.df

for (i in 1:nrow(GE_MTX.df)) {

  for (j in 1:ncol(GE_MTX.df)) {

    GE_MTX_DSum.df[i,j] <- GE_MTX.df[i,j]/GE_MTX_RowSum[i]

  }

}


GE_MTX_DSum_Log2.df <- log2(GE_MTX_DSum.df+1)

HFactor <- rowSums(GE_MTX_DSum.df*GE_MTX_DSum_Log2.df)




##***************************************************************************##
### HF filter and sort

### HF filter
GE_TarGene_Sum.df <- GE_TarGene.df
GE_TarGene_Sum.df$HF <- HFactor

GE_TarGene_Sum_HFFlt.df <- GE_TarGene_Sum.df[GE_TarGene_Sum.df$HF > 0.5,]


## Sort
# GE_TarGene_Sum_HFFltTTT.df <- GE_TarGene_Sum_HFFlt.df[order(GE_TarGene_Sum_HFFlt.df$S5),]
# GE_TarGene_Sum_HFFltTTT.df <- GE_TarGene_Sum_HFFlt.df[order(GE_TarGene_Sum_HFFlt.df[,3]),]

# # GE_TarGeneTTT.df <- GE_TarGene.df[order(GE_TarGene.df[,3],GE_TarGene.df[,4],decreasing = T),]
##Bug# GE_TarGeneTTT.df  <- apply(X = GE_TarGene.df, MARGIN = 2, FUN = sort)
##Bug# GE_TarGene.df <- apply(X = GE_TarGene.df, MARGIN = 2, FUN = sort)

GE_TarGene.df <- GE_TarGene.df[order(GE_TarGene.df[,3],GE_TarGene.df[,4],GE_TarGene.df[,5],
                                     GE_TarGene.df[,6],GE_TarGene.df[,7],GE_TarGene.df[,8],
                                     decreasing = T),]

##****************************************************************************##
## Create MTX
GE_TarGene.mx <- GE_TarGene.df[,-1:-2] %>% t()
GE_TarGene_Log.mx <- log10(GE_TarGene.mx+1e-5 %>% as.numeric())
quantile(GE_TarGene.mx)

MTX <- GE_TarGene_Log.mx


# # str <- ""
# # for (i in 3:ncol(GE_TarGene.df)) {
# #   str <- paste0(str,"GE_TarGene.df[,",i,"],")
# # }
# #
# # rm(i)
# #
# # formals(order)[names(str)] <- str
# #
# # GE_TarGeneTTT.df <- do.call("arrange", list(as.formula(str), .data=GE_TarGene.df[,-1:-2]))
# # GE_TarGeneTTT.df <- GE_TarGene.df[do.call("order", as.formula(str)),]
#
# # f <- paste0(TarGene[i]," ~ Cachexia") # f <- "Vwf ~ Cachexia"
# # SummaryTable_Temp.df <- do.call("compare_means", list(as.formula(f), data=Anno_Temp.df, group.by = "celltype", method = "wilcox.test"))
# # rm(f)



##***************************************************************************##

##### ComplexHeatmap #####
## https://jokergoo.github.io/ComplexHeatmap-reference/book/
library(ComplexHeatmap)
library(circlize)

## Try Heatmap
# Heatmap(GE_TarGene.df[,-1:-2])
# Heatmap(GE_TarGene.df[,-1:-2] %>% t())


#### Set row annotation ####
## Color setting
colType <- c("#4a8f8f","#a1338d","#e86020")
names(colType) <- c("Serum","Blood","Pancreatic juice")

ha_row = rowAnnotation(
  Type = HumanSampleAnno.df$Type,
  col = list(Type = colType),
  show_legend = T,
  show_annotation_name = T
)

#### Set col annotation ####
# GE_TarGene.df$Description %>% unique() %>% length()
GeneList <- GE_TarGene.df$Description %>% unique()

## Color setting
library(ggsci)
library(ggplot2)
# vignette( "ggsci") #Check the color setting
col3 = 	pal_d3("category20", alpha = 0.7)(length(GeneList))
col3 = 	c(pal_d3("category20", alpha = 0.7)(20),
          pal_d3("category20b", alpha = 0.7)(20),
          pal_d3("category20c", alpha = 0.7)(20),
          pal_d3("category10", alpha = 0.7)(10))
colCT <- col3[1:length(GeneList)]
names(colCT) <- GeneList

colHF <- c("#ffffff","#076e69","#0c5451")


ha_column_T = HeatmapAnnotation(
  HF = HFactor,
  GeneName = GE_TarGene.df$Description,  # anno_colum.df$sample_type,
  col = list(HF = colorRamp2(c(0,0.8,1),colHF),
             GeneName = colCT),
  show_legend = T,
  show_annotation_name = F
)

draw(ha_column_T)

# ## Color setting
# sample = c("#2267a4", "#3d85c6", "#d5a6bd", "#c27ba0")
# names(sample) <- Anno_Cell.df$sample %>% unique()
#
# library(ggsci)
# library(ggplot2)
# # vignette( "ggsci") #Check the color setting
# col3 = 	pal_d3("category20", alpha = 0.7)(length(CellType.Order))
# colCT <- col3[1:length(CellType.Order)]
# names(colCT) <- CellType.Order
#
# ha_column_T = HeatmapAnnotation(
#   Sample = Anno_Cell.df[,"sample"],  # anno_colum.df$sample_type,
#   Cachexia = Anno_Cell.df[,"Cachexia"], # anno_colum.df$gender,
#   Celltype = Anno_Cell.df[,"celltype"],
#   col = list(Sample = sample,
#              Celltype = colCT ,#pal_npg(), #colCT ,
#              Cachexia = c("EOCX"="#5b517d", "PreCX"="#a095c7")), #,"Medium"="#b57545"
#   show_legend = T,
#   show_annotation_name = F
# )




### Gene Expression Heatmap
# col_HMap2 <- c("#d0e0e3","#45818e")
# col_HMap2 <- c("#5a6ce0", "#EEEEEE","#d63c8e")
# col_HMap2 <- c("#5a6ce0", "#ffffff","#d63c8e")
# col_HMap2 <- c("#646466","#5a6ce0","#d63c8e")
col_HMap2 <- c("#ffffff","#d63c8e","#6e073d") #,"#9e135d"
col_HMap2 <- c("#ffffff","#ffffff","#ffffff","#666363","#1f1e1e") #,"#9e135d"


## Order by HF
Heatmap(
  MTX  ,
  cluster_rows = F,
  cluster_columns = F,
  # column_order = order(Anno_Cell.df$celltype,Anno_Cell.df$Cachexia),
  column_order = order(HFactor),
  row_order = Order_SampleID,
  show_column_names = F,
  show_row_names = T,
  name = "Log10GeneExp",
  # ## set color
  # col = col_HMap2,
  col = colorRamp2(c(min(MTX),(min(MTX)/2), 0, (max(MTX)/2), max(MTX)), col_HMap2),
  # col = colorRamp2(c(-max(GE_TarGene_Log.mx), 0, max(GE_TarGene_Log.mx)), col_HMap2),

  show_heatmap_legend = T,
  use_raster = F,
  top_annotation = ha_column_T,
  right_annotation = ha_row,
  # width = ncol( TarGeneAnno.mtx)*unit(5, "mm"), # length(CellType.Order)*unit(6, "mm")
  # height = nrow( TarGeneAnno.mtx)*unit(5, "mm") # length(CellType.Order)*unit(15, "mm"),
) -> P.Heatmap_GeneEXP

P.Heatmap_GeneEXP %>% print


## Order by cell type
Heatmap(
  MTX  ,
  cluster_rows = F,
  cluster_columns = F,
  # column_order = order(Anno_Cell.df$celltype,Anno_Cell.df$Cachexia),,
  show_column_names = F,
  show_row_names = T,
  name = "Log10GeneExp",
  # ## set color
  # col = col_HMap2,
  col = colorRamp2(c(min(MTX),(min(MTX)/2), 0, (max(MTX)/2), max(MTX)), col_HMap2),
  # col = colorRamp2(c(-max(GE_TarGene_Log.mx), 0, max(GE_TarGene_Log.mx)), col_HMap2),

  show_heatmap_legend = T,
  use_raster = F,
  top_annotation = ha_column_T,
  right_annotation = ha_row,
  # width = ncol( TarGeneAnno.mtx)*unit(5, "mm"), # length(CellType.Order)*unit(6, "mm")
  # height = nrow( TarGeneAnno.mtx)*unit(5, "mm") # length(CellType.Order)*unit(15, "mm"),
) -> P.Heatmap_GeneEXP2

P.Heatmap_GeneEXP2 %>% print

## Clustering
Heatmap(
  MTX  ,
  cluster_rows = T,
  cluster_columns = T,
  # column_order = order(Anno_Cell.df$celltype,Anno_Cell.df$Cachexia),
  show_column_names = F,
  show_row_names = T,
  name = "Log10GeneExp",
  # ## set color
  # col = col_HMap2,
  col = colorRamp2(c(min(MTX),(min(MTX)/2), 0, (max(MTX)/2), max(MTX)),  col_HMap2),
  # col = colorRamp2(c(-max(GE_TarGene_Log.mx), 0, max(GE_TarGene_Log.mx)), col_HMap2),

  show_heatmap_legend = T,
  use_raster = F,
  top_annotation = ha_column_T,
  right_annotation = ha_row,
  # width = ncol( TarGeneAnno.mtx)*unit(5, "mm"), # length(CellType.Order)*unit(6, "mm")
  # height = nrow( TarGeneAnno.mtx)*unit(5, "mm") # length(CellType.Order)*unit(15, "mm"),
) -> P.Heatmap_GeneEXP3

P.Heatmap_GeneEXP3 %>% print


## Check
TarGene.df$Gene %>% unique() %>% sort()
GeneList %>% sort()



##***************************************************************************##
##### Current path and new folder setting*  #####
Version = paste0(Sys.Date(),"_Cachexia_12HumanRNAseq_Entropy")
SaveCC.Path = Version
dir.create(SaveCC.Path)


##### Export PDF #####
pdf(file = paste0(SaveCC.Path,"/",Sys.Date(),"_", "Cachexia_12HumanRNAseq.pdf"),width = 15, height = 12 )
P.Heatmap_GeneEXP
P.Heatmap_GeneEXP2
P.Heatmap_GeneEXP3
dev.off()

##### Export TSV #####
write.table( GE_TarGene.df ,
             file = paste0(SaveCC.Path,"/",Sys.Date(),"_", "Cachexia_12HumanRNAseq_HeatmapMatrix.tsv"),
             sep = "\t",
             quote = F,
             row.names = F
)

