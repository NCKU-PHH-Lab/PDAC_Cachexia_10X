##### Beautify the Figures #####

library(ggplot2)
p1 <- DimPlot(PBMC.combined, reduction = "umap", group.by = "TypeID")
p1 + 
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
        legend.title = element_text(size=12, color = "black", face="bold"),
        legend.text = element_text(colour="black", size=12,face="bold"),
        legend.background = element_rect(fill = alpha("white", 0.5)),
        #      legend.position = c(0.1, 0.18),
        #     plot.text = element_text(size = 20),
        aspect.ratio=1) + #square plot
  theme(axis.line.x = element_line(colour = "black", size = 0.8),
        axis.line.y = element_line(colour = "black", size = 0.8))+
  theme(legend.position = c(0.75, 0.15))
