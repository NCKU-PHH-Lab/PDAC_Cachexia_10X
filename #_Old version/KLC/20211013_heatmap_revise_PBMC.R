##### data loading #####  
raw.df <- read.csv(
  "D:/KLC/### Datasets/Cachexia_differential-gene-heatmap/20210923/202109018_PBMC/PBMC_Cachexia.Marker_Sum_LogFC_M.txt",
  sep = "\t"
)
Sig.df <- read.csv(
  "D:/KLC/### Datasets/Cachexia_differential-gene-heatmap/20210923/202109018_PBMC/Cachexia_PBMC_GSEAforCheck.txt",
  sep = "\t"
)
arrange_celltype <- c("Mac_F","Mac_M","Mon_F","Mon_M","Neu_F","Neu_M","B_F","B_M","T_F","T_M","NK_F","NK_M","Ery_F","Ery_M","Thr_F","Thr_M")



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

##### generate side annotation plot #####
library(tidyr)
library(magrittr)
library(dplyr)

Anno.df <- Sig.df[c(-3,-6,-7,-9)]
Anno.df <- separate(Anno.df, PhenoType,c("gender","celltype")) %>% 
  data.frame(.,"pheno_type" = paste(.$celltype, .$gender, sep = "_")) 
Anno.df$celltype <- factor(Anno.df$celltype , levels = c("Mac","Mon","Neu","B","T","NK","Ery","Thr"))
Anno.df$gender <- factor(Anno.df$gender , levels = c("F","M"))
Anno.df$GeneType <- factor(Anno.df$GeneType , levels = c("Pos","Neg"))

NA_df <- Anno.df[1:2,]
NA_df$GeneType <- c("Neg","Neg")
NA_df$gender <- rep("F",times=2)
NA_df$celltype <- c("T","B")
NA_df$pheno_type <- c("T_F","B_F")
NA_df$pval <- c(NA,NA)
NA_df$padj <- c(NA,NA)
NA_df$NES <- c(NA,NA)

Anno.df <- bind_rows(Anno.df,NA_df)

Anno.df <- arrange(
  Anno.df,
  GeneType,
  celltype,
  gender
)
rm(Sig.df)

Pos_ES.mt <- Anno.df$NES[Anno.df$GeneType=="Pos"] %>% 
  as.matrix()

Neg_ES.mt <- Anno.df$NES[Anno.df$GeneType=="Neg"] %>% 
  as.matrix()

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
  col = colorRamp2(
    c(-4 ,-2, 0, 2, 4),
    c("#0000C6","#9393FF","#FFFFFF","#FF2D2D","#AE0000")
  ),
  at = c(-4 ,-2, 0, 2, 4),
  title = "log2FC",
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
  file = paste0(getwd(),"/output/PBMC_PairedMarker.tsv"),
  sep = "\t",
  quote = F,
  row.names = F
)

