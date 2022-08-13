##### RNA-seq analysis in R #####
# https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html
# install # https://bioconductor.org/packages/release/bioc/html/GSEABase.html
library(fgsea)
source("FUN_GSEA_LargeGeneSet.R")
source("FUN_HSsymbol2MMsymbol.R")


# Geneset from GSEA
# Pathway.all <- read.delim(paste0(PathName,"/Pathway.all.v7.4.symbols.gmt"),header = F)
Pathway.all <- read.delim2(paste0(PathName,"/GSEA_Geneset_Pathway_AllIndex.txt"),
                           col.names = 1:max(count.fields(paste0(PathName,"/GSEA_Geneset_Pathway_AllIndex.txt"))),
                           header = F,sep = "\t")

Pathway.all.MM = as.data.frame(matrix(nrow=nrow(Pathway.all),ncol=ncol(Pathway.all)*1.5))
for (i in 1:nrow(Pathway.all)) {
  #Pathway.all[,i] <- data.frame(colnames(Pathway.all)[i]=Pathway.all[,i]) %>% HSsymbol2MMsymbol(.,colnames(Pathway.all)[i])
  PathwayN <- data.frame(Pathway.all[i,3:ncol(Pathway.all)]) %>% t() 
  colnames(PathwayN)="Test"
  PathwayN <- HSsymbol2MMsymbol(PathwayN,"Test")
  Pathway.all.MM[i,1:length(unique(PathwayN$MM.symbol))] <- unique(PathwayN$MM.symbol)
  
}

Pathway.all.MM <- data.frame(Pathway.all[,1:2],Pathway.all.MM)
colnames(Pathway.all.MM) <- seq(1:ncol(Pathway.all.MM))

# assign(paste0("marrow_sub_DucT2_TOP2ACenter_T", i),marrow_sub_DucT2_TOP2ACenter_Tn)
# assign(colnames(Pathway.all)[i],Pathway.all[,i])

##### GSEA.Large_Female #####
GSEA.Large_Female <- list()
GSEA.Large_Female.df <- as.data.frame(matrix(nrow=0,ncol=10))
colnames(GSEA.Large_Female.df) <- c("GeneType","PhenoType","pathway","pval","padj","log2err","ES", "NES" ,"size","leadingEdge")
GSEA.Large_Female.df.TOP <- GSEA.Large_Female.df

for(i in 1:length(CellType.list)){
  
  gseaDat <- Cachexia.Marker.Female[[paste0(CellType.list[i])]][["Cachexia.Marker.All"]]
  gseaDat <- data.frame(row.names(gseaDat),gseaDat)
  colnames(gseaDat)[[1]] <- c("Gene")
  ranks <- gseaDat$avg_log2FC
  names(ranks) <- gseaDat$Gene
  # head(ranks)
  # barplot(sort(ranks, decreasing = T))
  
  
  #GSEA.Large.Output <- FUN_GSEA_LargeGeneSet(ranks,Pathway.all,10)
  GSEA.Large.Output <- FUN_GSEA_LargeGeneSet(ranks,Pathway.all.MM,10)
  
  fgseaRes <- GSEA.Large.Output[["fgseaRes"]]
  # head(fgseaRes[order(padj, -abs(NES)), ], n=10)
  
  pathwaysH <- GSEA.Large.Output[["Pathway.all.list"]] 
  
  # plot.new()
  # plotEnrichment(pathwaysH[[as.character(fgseaRes$pathway[1])]], ranks)
  
  topPathways <- GSEA.Large.Output[["topPathways"]]
  
  library(ggplot2)
  plot.new()
  plotGseaTable(pathwaysH[topPathways$pathway], 
                ranks, 
                fgseaRes, 
                gseaParam = 0.5) + title( paste0("SC.Female.",CellType.list[i]), adj = 0, line =3)
  
  plotEnrichment_Pos1 <- plotEnrichment(pathwaysH[[as.character(topPathways[1,1])]], ranks)+ labs(title= paste0("SC.Female.",CellType.list[i],": ",as.character(topPathways[1,1])))
  #plotEnrichment_Pos1
  plotEnrichment_Neg1 <- plotEnrichment(pathwaysH[[as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])]], ranks)+ labs(title= paste0("SC.Female.",CellType.list[i],": ",as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])))
  #plotEnrichment_Neg1
  
  Sum <- list(gseaDat,ranks,pathwaysH,fgseaRes,plotEnrichment_Pos1,plotEnrichment_Neg1)
  names(Sum) <- c("gseaDat","ranks","pathwaysH","fgseaRes","plotEnrichment_Pos1","plotEnrichment_Neg1")
  GSEA.Large_Female[[i]] <- Sum
  names(GSEA.Large_Female)[[i]] <- paste0("Female.",CellType.list[i])
  
  fgseaRes2 <- data.frame(paste0("F.",CellType.list[i]),fgseaRes)
  colnames(fgseaRes2)[[1]] <- c("PhenoType")
  GSEA.Large_Female.df <- rbind(GSEA.Large_Female.df,fgseaRes2 )
  
  topPathways2 <- data.frame(paste0("F.",CellType.list[i]),topPathways)
  colnames(topPathways2)[[1]] <- c("PhenoType")
  GSEA.Large_Female.df.TOP <- rbind(GSEA.Large_Female.df.TOP,topPathways2 )
  
  rm(gseaDat,ranks,pathwaysH,fgseaRes,fgseaRes2,plotEnrichmen,Sum,topPathways,topPathways2)
  
}

##### GSEA.Large_Male #####
library(fgsea)
source("FUN_GSEA_LargeGeneSet.R")
source("FUN_HSsymbol2MMsymbol.R")

GSEA.Large_Male <- list()
GSEA.Large_Male.df <- as.data.frame(matrix(nrow=0,ncol=10))
colnames(GSEA.Large_Male.df) <- c("GeneType","PhenoType","pathway","pval","padj","log2err","ES", "NES" ,"size","leadingEdge")
GSEA.Large_Male.df.TOP <- GSEA.Large_Male.df

for(i in 1:length(CellType.list)){
  
  gseaDat <- Cachexia.Marker.Male[[paste0(CellType.list[i])]][["Cachexia.Marker.All"]]
  gseaDat <- data.frame(row.names(gseaDat),gseaDat)
  colnames(gseaDat)[[1]] <- c("Gene")
  ranks <- gseaDat$avg_log2FC
  names(ranks) <- gseaDat$Gene
  # head(ranks)
  # barplot(sort(ranks, decreasing = T))
  
  
  #GSEA.Large.Output <- FUN_GSEA_LargeGeneSet(ranks,Pathway.all,10)
  GSEA.Large.Output <- FUN_GSEA_LargeGeneSet(ranks,Pathway.all.MM,10)
  
  fgseaRes <- GSEA.Large.Output[["fgseaRes"]]
  # head(fgseaRes[order(padj, -abs(NES)), ], n=10)
  
  pathwaysH <- GSEA.Large.Output[["Pathway.all.list"]] 
  
  # plot.new()
  # plotEnrichment(pathwaysH[[as.character(fgseaRes$pathway[1])]], ranks)
  
  topPathways <- GSEA.Large.Output[["topPathways"]]
  
  library(ggplot2)
  plot.new()
  plotGseaTable(pathwaysH[topPathways$pathway], 
                ranks, 
                fgseaRes, 
                gseaParam = 0.5) + title( paste0("SC.Male.",CellType.list[i]), adj = 0, line =3)
  
  plotEnrichment_Pos1 <- plotEnrichment(pathwaysH[[as.character(topPathways[1,1])]], ranks)+ labs(title= paste0("SC.Male.",CellType.list[i],": ",as.character(topPathways[1,1])))
  #plotEnrichment_Pos1
  plotEnrichment_Neg1 <- plotEnrichment(pathwaysH[[as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])]], ranks)+ labs(title= paste0("SC.Male.",CellType.list[i],": ",as.character(topPathways[length(as.data.frame(topPathways)[,1]),1])))
  #plotEnrichment_Neg1
  
  Sum <- list(gseaDat,ranks,pathwaysH,fgseaRes,plotEnrichment_Pos1,plotEnrichment_Neg1)
  names(Sum) <- c("gseaDat","ranks","pathwaysH","fgseaRes","plotEnrichment_Pos1","plotEnrichment_Neg1")
  GSEA.Large_Male[[i]] <- Sum
  names(GSEA.Large_Male)[[i]] <- paste0("Male.",CellType.list[i])
  
  fgseaRes2 <- data.frame(paste0("M.",CellType.list[i]),fgseaRes)
  colnames(fgseaRes2)[[1]] <- c("PhenoType")
  GSEA.Large_Male.df <- rbind(GSEA.Large_Male.df,fgseaRes2 )
  
  topPathways2 <- data.frame(paste0("M.",CellType.list[i]),topPathways)
  colnames(topPathways2)[[1]] <- c("PhenoType")
  GSEA.Large_Male.df.TOP <- rbind(GSEA.Large_Male.df.TOP,topPathways2 )
  
  rm(gseaDat,ranks,pathwaysH,fgseaRes,fgseaRes2,plotEnrichmen,Sum,topPathways,topPathways2)
  
}



##### GSEA.Large.Sum.TOP #####
GSEA.Large.Sum.TOP <- rbind(GSEA.Large_Male.df.TOP,GSEA.Large_Female.df.TOP)
GSEA.Large.Sum.TOP <- GSEA.Large.Sum.TOP[,!colnames(GSEA.Large.Sum.TOP) %in% c("leadingEdge")]
write.table(GSEA.Large.Sum.TOP, file=paste0(PathName,"/",RVersion,"/Cachexia_SC_GSEA_Pathway_LargeTOP.txt"),sep="\t",
            row.names=F, quote = FALSE)


##### Bubble plot #####
GSEA.Large.Sum.TOP.S <- GSEA.Large.Sum.TOP[abs(GSEA.Large.Sum.TOP$NES) > 1,]
GSEA.Large.Sum.TOP.S <- GSEA.Large.Sum.TOP.S[abs(GSEA.Large.Sum.TOP.S$padj) < 0.05,]
# GSEA.Large.Sum.TOP.S <- GSEA.Large.Sum.TOP[abs(GSEA.Large.Sum.TOP$padj) < 0.25,]
# GSEA.Large.Sum.TOP.S <- GSEA.Large.Sum.TOP.S[abs(GSEA.Large.Sum.TOP.S$pval) < 0.05,]
library(ggplot2)
library(scales)
ggplot(GSEA.Large.Sum.TOP.S,aes(x=PhenoType, y = pathway, color = padj, size = abs(NES))) + 
  geom_point() +
  scale_color_distiller()

ggplot(GSEA.Large.Sum.TOP.S,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) + 
  geom_point() +
  scale_colour_gradient2(low = "#2d76e3", mid = "white", high = "#e32dc2", 
                         guide = "colourbar",midpoint = 0)

BBPlot <- ggplot(GSEA.Large.Sum.TOP.S,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) + 
  geom_point() +
  scale_colour_gradient2(low = "#04873f", mid = "white", high = "#e3672d", 
                         guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")
BBPlot %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                          XtextSize=10,  YtextSize=5,AxisTitleSize=1, AspRat=2, XaThick=0.8, YaThick=0.8)

