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


##### Filter data #####
GE_TarGene.df <- GE.df[GE.df$Description %in% TarGene.df$Gene,]
row.names(GE_TarGene.df) <- GE_TarGene.df$NAME


##### ComplexHeatmap #####
## https://jokergoo.github.io/ComplexHeatmap-reference/book/
library(ComplexHeatmap)
library(circlize)

Heatmap(GE_TarGene.df[,-1:-2])
Heatmap(GE_TarGene.df[,-1:-2] %>% t())


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
GE_TarGene.df$Description %>% unique() %>% length()
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


ha_column_T = HeatmapAnnotation(
  GeneName = GE_TarGene.df$Description,  # anno_colum.df$sample_type,
  col = list(GeneName = colCT),
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
col_HMap2 <- c("#ffffff","#a19fa0","#1f1e1e") #,"#9e135d"

GE_TarGene.mx <- GE_TarGene.df[,-1:-2] %>% t()
GE_TarGene_Log.mx <- log10(GE_TarGene.mx+1e-5 %>% as.numeric())
quantile(GE_TarGene.mx)

MTX <- GE_TarGene_Log.mx

Heatmap(
  MTX  ,
  cluster_rows = T,
  cluster_columns = F,
  # column_order = order(Anno_Cell.df$celltype,Anno_Cell.df$Cachexia),
  show_column_names = F,
  show_row_names = T,
  name = "Log10GeneExp",
  # ## set color
  # col = col_HMap2,
  col = colorRamp2(c(min(MTX), 0, max(MTX)), col_HMap2),
  # col = colorRamp2(c(-max(GE_TarGene_Log.mx), 0, max(GE_TarGene_Log.mx)), col_HMap2),

  show_heatmap_legend = T,
  use_raster = F,
  top_annotation = ha_column_T,
  right_annotation = ha_row,
  # width = ncol( TarGeneAnno.mtx)*unit(5, "mm"), # length(CellType.Order)*unit(6, "mm")
  # height = nrow( TarGeneAnno.mtx)*unit(5, "mm") # length(CellType.Order)*unit(15, "mm"),
) -> P.Heatmap_GeneEXP

P.Heatmap_GeneEXP %>% print

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
  col = colorRamp2(c(min(MTX), 0, max(MTX)), col_HMap2),
  # col = colorRamp2(c(-max(GE_TarGene_Log.mx), 0, max(GE_TarGene_Log.mx)), col_HMap2),

  show_heatmap_legend = T,
  use_raster = F,
  top_annotation = ha_column_T,
  right_annotation = ha_row,
  # width = ncol( TarGeneAnno.mtx)*unit(5, "mm"), # length(CellType.Order)*unit(6, "mm")
  # height = nrow( TarGeneAnno.mtx)*unit(5, "mm") # length(CellType.Order)*unit(15, "mm"),
) -> P.Heatmap_GeneEXP2

P.Heatmap_GeneEXP2 %>% print


## Check
TarGene.df$Gene %>% unique() %>% sort()
GeneList %>% sort()


##### Export PDF #####
pdf(file = paste0("Cachexia_12HumanRNAseq.pdf"),width = 20, height = 12 )
P.Heatmap_GeneEXP
P.Heatmap_GeneEXP2
dev.off()

##### Export TSV #####
write.table( GE_TarGene.df ,
             file = paste0("Cachexia_12HumanRNAseq_HeatmapMatrix.tsv"),
             sep = "\t",
             quote = F,
             row.names = F
)

