##### data loading #####  
raw.df <- read.csv(
#  "D:/KLC/### Datasets/Cachexia_differential-gene-heatmap/20210923/202109018_PBMC/PBMC_Cachexia.Marker_Sum_LogFC_M.txt",
  "D:/Dropbox/Chang Charlene/##_GitHub/0-R/##_PHH_Lab/2021_Cachexia/20211121_PBMC/PBMC_Cachexia.Marker_Sum_LogFC_M.txt",
  
  sep = "\t"
)
raw.df[2,3] <- "CD4T"
raw.df[2,1] <- "M.CD4T"
raw.df[3,3] <- "CD8T"
raw.df[3,1] <- "M.CD8T"
raw.df[12,3] <- "CD4T"
raw.df[12,1] <- "F.CD4T"
raw.df[13,3] <- "CD8T"
raw.df[13,1] <- "F.CD8T"



Sig.df <- read.csv(
#  "D:/KLC/### Datasets/Cachexia_differential-gene-heatmap/20210923/202109018_PBMC/Cachexia_PBMC_GSEAforCheck.txt",
  "D:/Dropbox/Chang Charlene/##_GitHub/0-R/##_PHH_Lab/2021_Cachexia/20211121_PBMC/Cachexia_PBMC_GSEAforCheck2.txt",
  
  sep = "\t"
)

Sig.df[2,2] <- "M.CD4T"
Sig.df[3,2] <- "M.CD8T"
Sig.df[14,2] <- "M.CD4T"
Sig.df[15,2] <- "M.CD8T"
Sig.df[26,2] <- "F.CD4T"
Sig.df[27,2] <- "F.CD8T"
Sig.df[38,2] <- "F.CD4T"
Sig.df[39,2] <- "F.CD8T"


# arrange_celltype <- c("MacM0_F","MacM0_M","MacM1_F","MacM1_M","MacM2_F","MacM2_M","Mon_F","Mon_M",
#                       "Neu_F","Neu_M","B_F","B_M","NK_F","NK_M","B_F","B_M",
#                       "CD4T_F","CD4T_M","CD8T_F","CD8T_M","NK_F","NK_M","Ery_F","Ery_M","Thr_F","Thr_M")

arrange_celltype <- c("MacM0_F","MacM0_M","MacM2_F","MacM2_M","Mon_F","Mon_M",
                      "Neu_F","Neu_M","B_F","B_M", "CD4T_F","CD4T_M","CD8T_F",
                      "CD8T_M", "NK_F","NK_M","Ery_F","Ery_M","Thr_F","Thr_M")


##### data arrangement #####
library(magrittr)
library(dplyr)

# arrange raw.df sample
row.names(raw.df) <- paste(raw.df$CellType,raw.df$Sex,sep = "_")
raw.df <- raw.df[arrange_celltype,]
row.names(raw.df) <- c(1:nrow(raw.df)) # recover the serial number of row



# generate arrange.df
arrange.df <- data.frame(
  "gene"=colnames(raw.df)[c(-1:-3)],
  "count"=rep(NA,ncol(raw.df)-3), # non-zero value count
  "order"=rep(NA,ncol(raw.df)-3), # sum of serial number of row
  "sum"=rep(NA,ncol(raw.df)-3), # sum of expression level
  "co_exp"=rep(NA,ncol(raw.df)-3) # count the gene co-express in F and M
)

# generate count and order value
for(i in 1:nrow(arrange.df)){
  
  ## count
  arrange.df$count[i] <- c(raw.df[,arrange.df$gene[i]]!=0) %>% 
    as.numeric() %>% 
    sum()
  
  ## order
  set <- row.names(raw.df)[raw.df[,arrange.df$gene[i]]!=0] %>%
    as.numeric() 
  
  arrange.df$order[i] <- c(nrow(raw.df)+1-set) %>% 
    exp() %>% 
    sum()
  rm(set)
  
  ## sum
  arrange.df$sum[i] <- c(raw.df[,arrange.df$gene[i]]) %>% 
    sum()
  
  ## co_exp
  
  if(arrange.df$count[i]==2){
    odd <- seq(from =1, to = nrow(raw.df), by = 2)
    even <- odd+1
    arrange.df$co_exp[i] <- c(c(raw.df[odd,arrange.df$gene[i]]!=0) & c(raw.df[even,arrange.df$gene[i]]!=0)) %>% 
      as.numeric() %>% 
      sum()
    rm(odd,even)
  }else{
    arrange.df$co_exp[i] <- 0
  }
  
}
rm(i)


# arrangement and exclude sample have zero count
arrange.df <- arrange(
  arrange.df,
  -co_exp,
  count,
  -order,
  -sum
)
arrange.df <- c(arrange.df$count>0) %>% 
  arrange.df[.,]

# generate arranged df
df <- raw.df[c("Sex","CellType",arrange.df$gene)]
df <- df[-1] %>% 
  mutate("CellType"=paste(df$CellType,df$Sex,sep = "_"))
df <- data.frame(df[-1],row.names = df[[1]])

# arrange sample
df <- df[arrange_celltype,]
mt <- as.matrix(df)

# mt.p.df <- data.frame(matrix(NA,4,ncol(mt)))
# colnames(mt.p.df) <- colnames(mt)
# 
# row.names(mt.p.df) <- c("MacM1_F","MacM1_M" ,"Th_F","Th_M")
# mt <- rbind(mt,mt.p.df)
# 
# arrange_celltype2 <- c("MacM0_F","MacM0_M","MacM1_F","MacM1_M","MacM2_F","MacM2_M","Mon_F","Mon_M",
#                       "Neu_F","Neu_M","B_F","B_M", "CD4T_F","CD4T_M","CD8T_F",
#                       "CD8T_M", "Th_F","Th_M","NK_F","NK_M","Ery_F","Ery_M","Thr_F","Thr_M")
# mt <- mt[match(arrange_celltype2, rownames(mt)),]
##### generate side annotation plot #####
library(tidyr)
library(magrittr)
library(dplyr)

Anno.df <- Sig.df[c(-3,-6,-7,-9)]
Anno.df <- separate(Anno.df, PhenoType,c("gender","celltype")) %>% 
  data.frame(.,"pheno_type" = paste(.$celltype, .$gender, sep = "_")) 
#Anno.df$celltype <- factor(Anno.df$celltype , levels = c("Mac","Mon","Neu","B","T","NK","Ery","Thr"))
Anno.df$celltype <- factor(Anno.df$celltype , levels = c("MacM0","MacM1","MacM2","Mon","Neu","B","CD4T","CD8T","Th","NK","Ery","Thr"))

Anno.df$gender <- factor(Anno.df$gender , levels = c("F","M"))
Anno.df$GeneType <- factor(Anno.df$GeneType , levels = c("Pos","Neg"))

# NA_df <- Anno.df[1:2,]
# NA_df$GeneType <- c("Neg","Neg")
# NA_df$gender <- rep("F",times=2)
# NA_df$celltype <- c("T","B")
# NA_df$pheno_type <- c("T_F","B_F")
# NA_df$pval <- c(NA,NA)
# NA_df$padj <- c(NA,NA)
# NA_df$NES <- c(NA,NA)

NA_df <- Anno.df[1,]
NA_df$GeneType <- c("Neg")
NA_df$gender <- rep("F",times=1)
NA_df$celltype <- c("B")
NA_df$pheno_type <- c("B_F")
NA_df$pval <- c(NA)
NA_df$padj <- c(NA)
NA_df$NES <- c(NA)

Anno.df <- bind_rows(Anno.df,NA_df)

Anno.df <- arrange(
  Anno.df,
  GeneType,
  celltype,
  gender
)
rm(Sig.df)


Pos_ES <- Anno.df[Anno.df$GeneType=="Pos",] 
Pos_ES <- Pos_ES[match(arrange_celltype, Pos_ES$pheno_type),]
Neg_ES <- Anno.df[Anno.df$GeneType=="Neg",] 
Neg_ES <- Neg_ES[match(arrange_celltype, Neg_ES$pheno_type),]
Pos_ES.mt <- Pos_ES$NES %>% as.matrix()
Neg_ES.mt <- Neg_ES$NES %>% as.matrix()

# ## Anno.df <- Anno.df[match(arrange_celltype, Anno.df$pheno_type),]
# Pos_ES.mt <- Anno.df$NES[Anno.df$GeneType=="Pos"] %>% 
#   as.matrix()
# 
# Neg_ES.mt <- Anno.df$NES[Anno.df$GeneType=="Neg"] %>% 
#   as.matrix()

##### making heatmap #####
library(ComplexHeatmap)
library(circlize)
library(ggplot2)

lgd_ES <- Legend(
  title = "NES", 
  col_fun = colorRamp2(c(-2, -1, 0, 1, 2), c("#5CADAD","#B3D9D9", "#E0E0E0","#DAB1D5", "#AE57A4")), 
  at = c(-2, -1, 0, 1, 2), 
  legend_height = unit(2, "cm"),
  direction = "horizontal"
)

lgd_hp <- Legend(
  title = "log2FC",
  col_fun = colorRamp2(
    c(-4 ,-2, 0, 2, 4),
    c("#0000C6","#9393FF","#FFFFFF","#FF2D2D","#AE0000")
  ),
  at = c(-4 ,-2, 0, 2, 4),
  legend_height = unit(2, "cm"),
  title_position = "topleft",
  direction = "horizontal"
)

fig1 <- rowAnnotation(
  Pos = Pos_ES.mt,
  Neg = Neg_ES.mt,
  col = list(
    Pos=colorRamp2(c(-2, -1, 0, 1, 2), c("#5CADAD","#B3D9D9", "#E0E0E0","#DAB1D5", "#AE57A4")),
    Neg=colorRamp2(c(-2, -1, 0, 1, 2), c("#5CADAD","#B3D9D9", "#E0E0E0","#DAB1D5", "#AE57A4"))
  ),
  show_legend = F,
  na_col = "white"
) + 
  Heatmap(
    mt,
    column_title = "PBMC_cachexia marker",
    column_title_side = "top",
    cluster_rows = F,
    cluster_columns = F,
    show_column_names = F,
    name = "mt",
    # set color
    col = colorRamp2(
      c(-5 ,-3,-1, 0, 1, 3, 5),
      c("#0000C6","#9393FF","#FBFBFF","#FFFFFF","#FFECEC","#FF2D2D","#AE0000")
    ),
    show_heatmap_legend = F
  ) 
draw(fig1, annotation_legend_side = "bottom", annotation_legend_list = list(lgd_hp,lgd_ES))  

##### Generate output file #####
output.df <- df[arrange.df$co_exp>0]
write.table(
  data.frame("Phenotype"=row.names(output.df),df[arrange.df$co_exp>0]) ,
  file = paste0(getwd(),"/20211121_PBMC/PBMC_PairedMarker.tsv"),
  sep = "\t",
  quote = F,
  row.names = F
)

