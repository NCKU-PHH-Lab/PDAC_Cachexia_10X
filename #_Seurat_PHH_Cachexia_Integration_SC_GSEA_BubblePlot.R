# https://zhuanlan.zhihu.com/p/352143854



#####  Current path and new folder setting ##### 
PathName = setwd(getwd())
RVersion = "20211201_SC"
# dir.create(paste0(PathName,"/",RVersion))

GSEA.df <- read.csv(paste0(PathName,"/",RVersion,"/Cachexia_SC_GSEA_Pathway_LargeTOP_Ch_FDR005.txt"),
                    sep = "\t")

library(ggplot2)
library("magrittr")
source("FUN_Beautify_ggplot.R")

##### Original Bubble Plot #####
BBPlot <- ggplot(GSEA.df,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) + 
  geom_point() +
  scale_colour_gradient2(low = "#2d76e3", mid = "white", high = "#e32dc2", 
                         guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")+
  scale_size_area(max_size = 3)

BBPlot

BBPlot %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                          XtextSize=10,  YtextSize=5,AxisTitleSize=1, AspRat=2, XaThick=0.8, YaThick=0.8)

##
BBPlot2 <- ggplot(GSEA.df,aes(x=pathway , y = PhenoType, color = NES, size = -log10(padj))) + 
  geom_point() +
  scale_colour_gradient2(low = "#2d76e3", mid = "white", high = "#e32dc2", 
                         guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")

BBPlot2

BBPlot2 %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                          XtextSize=10,  YtextSize=5,AxisTitleSize=1, AspRat=2, XaThick=0.8, YaThick=0.8)



##### Clustering #####
df1 <- reshape2::dcast(GSEA.df,PhenoType~pathway,value.var = "NES")
rownames(df1)<-df1$id

df1.1<-df1[,2:ncol(df1)]

df1.1.clust.Pheno<-hclust(dist(df1.1))
df1.1.clust.Pathway<-hclust(dist(t(df1.1)))


library(ggtree)
library(dplyr)
PhenoType_Order <- data.frame(No=row.names(df1),PhenoType=df1[,1]) 
PhenoType_Order$No <- as.numeric(PhenoType_Order$No)
PhenoType_Treeclust_Order <- data.frame(No=df1.1.clust.Pheno[["order"]])
PhenoType_Treeclust_Order <- left_join(PhenoType_Treeclust_Order,PhenoType_Order) 
df1.1.clust.Pheno[["order"]] <- PhenoType_Treeclust_Order$PhenoType
df1.1.clust.Pheno[["labels"]] <- PhenoType_Treeclust_Order$PhenoType

p2 <- ggtree(df1.1.clust.Pheno)
p2+
  geom_tiplab()+
  xlim(NA,7)
p2.2<-p2+
  geom_tiplab()+
  xlim(NA,7)+
  #geom_tiplab(angle=90)+
  #theme_tree2()+
  layout_dendrogram()
p2.2

p3<-ggtree(df1.1.clust.Pathway)
p3+
  geom_tiplab()+
  xlim(NA,7)
## Plot
library(aplot)
BBPlot_Cluster<- BBPlot %>%
                 insert_left(p3,width = 0.2)


BBPlot_Cluster

BBPlotB <- BBPlot %>% 
  BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                   XtextSize=10,  YtextSize=5,AxisTitleSize=1, AspRat=2, XaThick=0.6, YaThick=0.6)

BBPlotB %>%
  insert_left(p3,width = 0.2)


BBPlotB %>%
   insert_left(p3,width = 0.2)%>%
   insert_top(p2+layout_dendrogram(),height = 0.2)

pdf(file = paste0(PathName,"/",RVersion,"/SC_GSEA_Bubble.pdf"),width = 17, height = 7 )
BBPlotB %>%
  insert_left(p3,width = 0.2)

BBPlotB %>%
  insert_left(p3,width = 0.2)%>%
  insert_top(p2+layout_dendrogram(),height = 0.2)
dev.off()


##### Extract SubType #####

## Fib
GSEA_Fib.df <- GSEA.df[grep("Fib",GSEA.df$PhenoType),]

BBPlot_Fib <- ggplot(GSEA_Fib.df,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) + 
  geom_point() +
  scale_size_area(max_size = 7)+
  scale_colour_gradient2(low = "#2d76e3", mid = "white", high = "#e32dc2", 
                         guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")

BBPlot_Fib

BBPlot_FibB <- BBPlot_Fib %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                          XtextSize=15 ,  YtextSize=12,AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)

pdf(file = paste0(PathName,"/",RVersion,"/SC_GSEA_Bubble_SubType_Fib.pdf"),width = 17, height = 7 )
BBPlot_FibB %>%
  insert_left(p3,width = 0.2)
dev.off()

## Duc
GSEA_Duc.df <- GSEA.df[grep("Duc",GSEA.df$PhenoType),]

BBPlot_Duc <- ggplot(GSEA_Duc.df,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) + 
  geom_point() +
  scale_size_area(max_size = 5)+
  scale_colour_gradient2(low = "#2d76e3", mid = "white", high = "#e32dc2", 
                         guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")

BBPlot_Duc

BBPlot_DucB <- BBPlot_Duc %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                                             XtextSize=15,  YtextSize=12, AxisTitleSize=1, AspRat=2, XaThick=0.8, YaThick=0.8)

pdf(file = paste0(PathName,"/",RVersion,"/SC_GSEA_Bubble_SubType_Duc.pdf"),width = 17, height = 22 )
BBPlot_DucB %>%
  insert_left(p3,width = 0.2)
dev.off()

## Mac
GSEA_Mac.df <- GSEA.df[grep("Mac",GSEA.df$PhenoType),]

BBPlot_Mac <- ggplot(GSEA_Mac.df,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) + 
  geom_point() +
  scale_size_area(max_size = 5)+
  scale_colour_gradient2(low = "#2d76e3", mid = "white", high = "#e32dc2", 
                         guide = "colourbar",midpoint = 0)+ theme(legend.position = "bottom")

BBPlot_Mac

BBPlot_MacB <- BBPlot_Mac %>% BeautifyggPlot(LegPos  = "bottom",LegBox = "horizontal",LegDir="horizontal", xangle =90,
                                             XtextSize=15,  YtextSize=10, AxisTitleSize=1, AspRat=4, XaThick=0.8, YaThick=0.8)

pdf(file = paste0(PathName,"/",RVersion,"/SC_GSEA_Bubble_SubType_Mac.pdf"),width = 17, height = 20 )
BBPlot_MacB %>%
  insert_left(p3,width = 0.2)
dev.off()