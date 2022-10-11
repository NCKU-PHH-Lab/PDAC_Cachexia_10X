SampleType = "PBMC"

## INTCHG: Interchangeable
## SubType Setting
if(SampleType == "PBMC"){
  # For PBMC
  scRNA.SeuObj <- PBMC.combined

}else{
  # For SC
  scRNA.SeuObj <- SC.combined

}


##### Beautify the Figures #####

library(ggplot2)

##### UMAP sample umap #####
p1 <- DimPlot(scRNA.SeuObj, reduction = "umap", group.by = "sample")
p1 +
  # scale_colour_gradient2(low = "black", mid = "#3528c7", high = "#ff85e9",
  #                        guide = "colourbar",midpoint = 0, labs(fill ="Exp")) +
  #            ggtitle(Main)+
  # ggtitle(paste0(Main,"(",Sub_Name,")"))+
  theme(axis.text.x = element_text(face="bold",  size=17),
        axis.text.y = element_text(face="bold",size=17),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(0.9),face="bold"),
        plot.title = element_text(color="black", size=20,
                                  face="bold.italic",hjust = 0.1,vjust =-10), # margin = margin(t = 0.5, b = -7),
        #     plot.background = element_rect(fill = 'chartreuse'),
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"),
        legend.background = element_rect(fill = alpha("white", 0.5)),
        #      legend.position = c(0.1, 0.18),
        #     plot.text = element_text(size = 20),
        aspect.ratio=1) + #square plot
  theme(axis.line.x = element_line(colour = "black", size = 1.2),
        axis.line.y = element_line(colour = "black", size = 1.2))+
  theme(legend.position = c(0.75, 0.15))

##### UMAP sample tsne #####
p1_2 <- DimPlot(scRNA.SeuObj, reduction = "tsne", group.by = "sample")
p1_2 +
  # scale_colour_gradient2(low = "black", mid = "#3528c7", high = "#ff85e9",
  #                        guide = "colourbar",midpoint = 0, labs(fill ="Exp")) +
  #            ggtitle(Main)+
  # ggtitle(paste0(Main,"(",Sub_Name,")"))+
  theme(axis.text.x = element_text(face="bold",  size=17),
        axis.text.y = element_text(face="bold",size=17),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(0.9),face="bold"),
        plot.title = element_text(color="black", size=20,
                                  face="bold.italic",hjust = 0.1,vjust =-10), # margin = margin(t = 0.5, b = -7),
        #     plot.background = element_rect(fill = 'chartreuse'),
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"),
        legend.background = element_rect(fill = alpha("white", 0.5)),
        #      legend.position = c(0.1, 0.18),
        #     plot.text = element_text(size = 20),
        aspect.ratio=1) + #square plot
  theme(axis.line.x = element_line(colour = "black", size = 1.2),
        axis.line.y = element_line(colour = "black", size = 1.2))+
  theme(legend.position = c(0.75, 0.1))


##### UMAP Cluster #####
p2 <- DimPlot(scRNA.SeuObj, label = TRUE, repel = TRUE,label.size = 4,ncol = 3)
p2 +
  # scale_colour_gradient2(low = "black", mid = "#3528c7", high = "#ff85e9",
  #                        guide = "colourbar",midpoint = 0, labs(fill ="Exp")) +
  #            ggtitle(Main)+
  # ggtitle(paste0(Main,"(",Sub_Name,")"))+
  theme(axis.text.x = element_text(face="bold",  size=17),
        axis.text.y = element_text(face="bold",size=17),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(0.9),face="bold"),
        plot.title = element_text(color="black", size=20,
                                  face="bold.italic",hjust = 0.1,vjust =-10), # margin = margin(t = 0.5, b = -7),
        #     plot.background = element_rect(fill = 'chartreuse'),
        #legend.direction = "horizontal",
        #legend.box = "horizontal",
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=10,face="bold"),
        legend.background = element_rect(fill = alpha("white", 0.5)),
        #      legend.position = c(0.1, 0.18),
        #     plot.text = element_text(size = 20),
        aspect.ratio=1) + #square plot
  theme(axis.line.x = element_line(colour = "black", size = 1.2),
        axis.line.y = element_line(colour = "black", size = 1.2))+
  theme(legend.position = c(1, 0.50)) # legends position
  # guides(col = guide_legend(ncol = 3)) + # multiple row/col legends
  # theme(legend.key.size = unit(1, 'cm'))
  # guides(color = guide_legend(override.aes = list(size = 3) ) ) + # https://aosmith.rbind.io/2020/07/09/ggplot2-override-aes/

##### UMAP sample #####
p3 <- DimPlot(scRNA.SeuObj, label = TRUE, reduction = "umap", ncol = 2,split.by = "sample",label.size = 4)
p3 +
  # scale_colour_gradient2(low = "black", mid = "#3528c7", high = "#ff85e9",
  #                        guide = "colourbar",midpoint = 0, labs(fill ="Exp")) +
  #            ggtitle(Main)+
  # ggtitle(paste0(Main,"(",Sub_Name,")"))+
  theme(axis.text.x = element_text(face="bold",  size=14),
        axis.text.y = element_text(face="bold",size=14),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(1),face="bold"),
        strip.text = element_text(size=15), ### https://github.com/satijalab/seurat/issues/2471
        plot.title = element_text(color="black", size=20,
                                  face="bold.italic",hjust = 0.1,vjust =-10), # margin = margin(t = 0.5, b = -7),
        #     plot.background = element_rect(fill = 'chartreuse'),
        legend.direction = "vertical",  # vertical, horizontal
        legend.box = "vertical", # vertical, horizontal
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=10,face="bold"),
        legend.background = element_rect(fill = alpha("white", 0.5)),
        #      legend.position = c(0.1, 0.18),
        #     plot.text = element_text(size = 20),
        aspect.ratio=1) + #square plot
  theme(axis.line.x = element_line(colour = "black", size = 1.2),
        axis.line.y = element_line(colour = "black", size = 1.2))+
  theme(legend.position = c(1, 0.5)) # legends position

##### UMAP Cachexia #####
p4 <- DimPlot(scRNA.SeuObj, label = TRUE, reduction = "umap", ncol = 2,split.by = "Cachexia",label.size = 5)
p4 +
  # scale_colour_gradient2(low = "black", mid = "#3528c7", high = "#ff85e9",
  #                        guide = "colourbar",midpoint = 0, labs(fill ="Exp")) +
  #            ggtitle(Main)+
  # ggtitle(paste0(Main,"(",Sub_Name,")"))+
  theme(axis.text.x = element_text(face="bold",  size=14),
        axis.text.y = element_text(face="bold",size=14),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(0.9),face="bold"),
        plot.title = element_text(color="black", size=20,
                                  face="bold.italic",hjust = 0.1,vjust =-10), # margin = margin(t = 0.5, b = -7),
        #     plot.background = element_rect(fill = 'chartreuse'),
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.title = element_text(size=15, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.background = element_rect(fill = alpha("white", 0.5)),
        #      legend.position = c(0.1, 0.18),
        #     plot.text = element_text(size = 20),
        aspect.ratio=1) + #square plot
  theme(axis.line.x = element_line(colour = "black", size = 1.0),
        axis.line.y = element_line(colour = "black", size = 1.0))+
  theme(legend.position = c(0.35, 0.9)) # legends position

##### UMAP Sex #####
p5 <- DimPlot(scRNA.SeuObj, label = TRUE, reduction = "umap", ncol = 2,split.by = "Sex",label.size = 5)
p5 +
  # scale_colour_gradient2(low = "black", mid = "#3528c7", high = "#ff85e9",
  #                        guide = "colourbar",midpoint = 0, labs(fill ="Exp")) +
  #            ggtitle(Main)+
  # ggtitle(paste0(Main,"(",Sub_Name,")"))+
  theme(axis.text.x = element_text(face="bold",  size=14),
        axis.text.y = element_text(face="bold",size=14),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(0.9),face="bold"),
        plot.title = element_text(color="black", size=20,
                                  face="bold.italic",hjust = 0.1,vjust =-10), # margin = margin(t = 0.5, b = -7),
        #     plot.background = element_rect(fill = 'chartreuse'),
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.title = element_text(size=15, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12, face="bold"),
        legend.background = element_rect(fill = alpha("white", 0.5)),
        #      legend.position = c(0.1, 0.18),
        #     plot.text = element_text(size = 20),
        aspect.ratio=1) + #square plot
  theme(axis.line.x = element_line(colour = "black", size = 1.0),
        axis.line.y = element_line(colour = "black", size = 1.0))+
  theme(legend.position = c(0.35, 0.9)) # legends position

##### UMAP seurat_clusters #####
p6 <- DimPlot(scRNA.SeuObj, label = TRUE ,group.by = "seurat_clusters", label.size = 5)
p6 +
  # scale_colour_gradient2(low = "black", mid = "#3528c7", high = "#ff85e9",
  #                        guide = "colourbar",midpoint = 0, labs(fill ="Exp")) +
  #            ggtitle(Main)+
  # ggtitle(paste0(Main,"(",Sub_Name,")"))+
  theme(axis.text.x = element_text(face="bold",  size=17),
        axis.text.y = element_text(face="bold",size=17),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(0.9),face="bold"),
        plot.title = element_text(color="black", size=20,
                                  face="bold.italic",hjust = 0.1,vjust =-10), # margin = margin(t = 0.5, b = -7),
        #     plot.background = element_rect(fill = 'chartreuse'),
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=10,face="bold"),
        legend.background = element_rect(fill = alpha("white", 0.5)),
        #      legend.position = c(0.1, 0.18),
        #     plot.text = element_text(size = 20),
        aspect.ratio=1) + #square plot
  theme(axis.line.x = element_line(colour = "black", size = 1.2),
        axis.line.y = element_line(colour = "black", size = 1.2))+
  theme(legend.position = c(0.6, 0.10)) # legends position

##### UMAP CellType #####
p7 <- DimPlot(scRNA.SeuObj, label = TRUE, reduction = "umap",group.by = "celltype", label.size = 5)
p7 +
  # scale_colour_gradient2(low = "black", mid = "#3528c7", high = "#ff85e9",
  #                        guide = "colourbar",midpoint = 0, labs(fill ="Exp")) +
  #            ggtitle(Main)+
  # ggtitle(paste0(Main,"(",Sub_Name,")"))+
  theme(axis.text.x = element_text(face="bold",  size=17),
        axis.text.y = element_text(face="bold",size=17),
        axis.line = element_line(colour = "darkblue", size = 2, linetype = "solid"),
        axis.title = element_text(size = rel(0.9),face="bold"),
        plot.title = element_text(color="black", size=20,
                                  face="bold.italic",hjust = 0.1,vjust =-10), # margin = margin(t = 0.5, b = -7),
        #     plot.background = element_rect(fill = 'chartreuse'),
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"),
        legend.background = element_rect(fill = alpha("white", 0.5)),
        #      legend.position = c(0.1, 0.18),
        #     plot.text = element_text(size = 20),
        aspect.ratio=1) + #square plot
  theme(axis.line.x = element_line(colour = "black", size = 1.2),
        axis.line.y = element_line(colour = "black", size = 1.2))+
  theme(legend.position = c(0.9, 0.8))


##### Bubble TypeID #####
Bubble1 <- DotPlot(scRNA.SeuObj, features = markers.to.plot, cols = c("#1f2041","#4b3f72","#119da4", "#19647e"), dot.scale = 8, split.by = "TypeID") +
  RotatedAxis()

Bubble1 +
  # scale_colour_gradient2(low = "black", mid = "#3528c7", high = "#ff85e9",
  #                        guide = "colourbar",midpoint = 0, labs(fill ="Exp")) +
  #            ggtitle(Main)+
  # ggtitle(paste0(Main,"(",Sub_Name,")"))+
  theme(axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",size=10),
        axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.title = element_text(size = rel(0.9),face="bold"),
        plot.title = element_text(color="black", size=20,
                                  face="bold.italic",hjust = 0.1,vjust =-10), # margin = margin(t = 0.5, b = -7),
        #     plot.background = element_rect(fill = 'chartreuse'),
        legend.title = element_text(size=11, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"),
        legend.background = element_rect(fill = alpha("white", 0.5)) )+
        theme(legend.position="bottom",legend.box = "horizontal")# +#,
        #      legend.position = c(0.1, 0.18),
        #     plot.text = element_text(size = 20),
        # aspect.ratio=1) + #square plot
        # theme(axis.line.x = element_line(colour = "black", size = 1.2),
        #       axis.line.y = element_line(colour = "black", size = 1.2))+
        # theme(legend.position = c(0.9, 0.8))

##### Bubble Cachexia #####
Bubble2 <- DotPlot(scRNA.SeuObj, features = markers.to.plot, cols = c("#cb1b16","#003f88"), dot.scale = 8, split.by = "Cachexia") +
    RotatedAxis()

Bubble2 +
    theme(axis.text.x = element_text(face="bold",  size=10),
          axis.text.y = element_text(face="bold",size=10),
          axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
          axis.title = element_text(size = rel(0.9),face="bold"),
          plot.title = element_text(color="black", size=20,
                                    face="bold.italic",hjust = 0.1,vjust =-10), # margin = margin(t = 0.5, b = -7),
          #     plot.background = element_rect(fill = 'chartreuse'),
          legend.title = element_text(size=11, color = "black", face="bold"),
          legend.text = element_text(colour="black", size=12,face="bold"),
          legend.background = element_rect(fill = alpha("white", 0.5)) )+
          theme(legend.position="bottom",legend.box = "horizontal")# +#,

##### Bubble Sex #####
Bubble3 <- DotPlot(scRNA.SeuObj, features = markers.to.plot, cols = c( "#5b8e7d","#7b2cbf"), dot.scale = 8, split.by = "Sex") +
  RotatedAxis()

Bubble3 +
  theme(axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",size=10),
        axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.title = element_text(size = rel(0.9),face="bold"),
        plot.title = element_text(color="black", size=20,
                                  face="bold.italic",hjust = 0.1,vjust =-10), # margin = margin(t = 0.5, b = -7),
        #     plot.background = element_rect(fill = 'chartreuse'),
        legend.title = element_text(size=11, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"),
        legend.background = element_rect(fill = alpha("white", 0.5)) )+
        theme(legend.position="bottom",legend.box = "horizontal")# +#,


##### Bubble plot #####
ggplot(GSEA.Large.Sum.TOP.S,aes(x=PhenoType, y = pathway, color = NES, size = -log10(padj))) +
  geom_point() +
  scale_colour_gradient2(low = "#04873f", mid = "white", high = "#e3672d",
                         guide = "colourbar",midpoint = 0)+
  theme(axis.text.x = element_text(face="bold",  size=13, angle = 60,vjust =0.5),
        axis.text.y = element_text(face="bold",size=10),
        axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.title = element_text(size = rel(1.5),face="bold"),
        plot.title = element_text(color="black", size=20,
                                  face="bold.italic",hjust = 0.1,vjust =-10), # margin = margin(t = 0.5, b = -7),
        #     plot.background = element_rect(fill = 'chartreuse'),
        legend.title = element_text(size=11, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"),
        legend.background = element_rect(fill = alpha("white", 0.5)) )+ scale_size(range = c(2,10))+
  theme(panel.background = element_rect(fill = "#edebeb", colour = "black",
                                        size = 1))+ theme(legend.position = "bottom")

##### Bubble Candidate genes #####
Bubble4 <- DotPlot(scRNA.SeuObj, features = Cachexia_Marker_HM, cols = c("#e864a4", "#c90c5e","#4169e1","#092e9c"), dot.scale = 8, split.by = "TypeID") +
  RotatedAxis()
Bubble4+
  theme(axis.text.x = element_text(face="bold",  size=10),
        axis.text.y = element_text(face="bold",size=10),
        axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.title = element_text(size = rel(0.9),face="bold"),
        plot.title = element_text(color="black", size=20,
                                  face="bold.italic",hjust = 0.1,vjust =-10), # margin = margin(t = 0.5, b = -7),
        #     plot.background = element_rect(fill = 'chartreuse'),
        legend.title = element_text(size=11, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"),
        legend.background = element_rect(fill = alpha("white", 0.5)) )+
  theme(legend.position="bottom",legend.box = "horizontal")# +#,

