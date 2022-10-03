# ## INTCHG: Interchangeable
#   # For PBMC
#   scRNA.SeuObj <- PBMC.combined
#   SampleType = "PBMC"
#
#   ## For SC
#   # scRNA.SeuObj <- SC.combined
#   # SampleType = "SC"

source("FUN_HSsymbol2MMsymbol.R")

##### HALLMARK_TNFA_SIGNALING_VIA_NFKB #####
  HALLMARK_TNFA_SIGNALING_VIA_NFKB <- read.delim(paste0(getwd(),"/GSEA_Geneset/HALLMARK_TNFA_SIGNALING_VIA_NFKB.txt"),header=T, stringsAsFactors = FALSE)
  HALLMARK_TNFA_SIGNALING_VIA_NFKB <- data.frame(HALLMARK_TNFA_SIGNALING_VIA_NFKB[-1,])
  colnames(HALLMARK_TNFA_SIGNALING_VIA_NFKB) <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB")

  HALLMARK_TNFA_SIGNALING_VIA_NFKB_Mouse <- HSsymbol2MMsymbol(HALLMARK_TNFA_SIGNALING_VIA_NFKB,"HALLMARK_TNFA_SIGNALING_VIA_NFKB")
  HALLMARK_TNFA_SIGNALING_VIA_NFKB_Mouse2 <- unique(HALLMARK_TNFA_SIGNALING_VIA_NFKB_Mouse[,2])
  # FeaturePlot(scRNA.SeuObj, features = HALLMARK_TNFA_SIGNALING_VIA_NFKB_Mouse2[1:24], min.cutoff = "q9")
  # FeaturePlot(scRNA.SeuObj, features = HALLMARK_TNFA_SIGNALING_VIA_NFKB_Mouse2[25:50], min.cutoff = "q9")
  # FeaturePlot(scRNA.SeuObj, features =c("Crip1","Arg1","Mt2","Tspo"), min.cutoff = "q9")

  scRNA.SeuObj <- AddModuleScore(object = scRNA.SeuObj, features = as.data.frame(HALLMARK_TNFA_SIGNALING_VIA_NFKB_Mouse2),
                                  ctrl = 5, name = 'HALLMARK_TNFA_SIGNALING_VIA_NFKB')
  FeaturePlot(object = scRNA.SeuObj, features = 'HALLMARK_TNFA_SIGNALING_VIA_NFKB1', split.by = "Cachexia")+
    scale_colour_gradient2(low ="#0077b6", mid = "white", high = "#ef476f",
                           guide = "colourbar",midpoint = 0, labs(fill ="Exp"))

  FeaturePlot(object = scRNA.SeuObj, features = 'HALLMARK_TNFA_SIGNALING_VIA_NFKB1', split.by = "sample",
              cols = c("#0077b6", "white",  "#ef476f"))
  plots <- VlnPlot(scRNA.SeuObj, features = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB1"), split.by = "sample", group.by = "celltype",
                   pt.size = 0, combine = FALSE,cols = c("#ef476f", "#e87993",  "#0077b6",  "#3c9ccf"))
  wrap_plots(plots = plots, ncol = 1)

  plots <- VlnPlot(scRNA.SeuObj, features = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB1"), split.by = "Cachexia", group.by = "celltype",
                   pt.size = 0, combine = FALSE,cols = c("#ef476f",  "#0077b6"))
  p1 <- wrap_plots(plots = plots, ncol = 1) %>%  BeautifyggPlot(.,AspRat=1,LegPos = c(-0.08, -0.08),AxisTitleSize=1.7,
                                                        XtextSize=18,  YtextSize=,18, xangle = 90,
                                                        LegTextSize = 15)  +
                   theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

  p1

  ###
##### REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION #####
  REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION <- read.delim(paste0(getwd(),"/GSEA_Geneset/REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION.txt"),header=T, stringsAsFactors = FALSE)
  REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION <- data.frame(REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION[-1,])
  colnames(REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION) <- c("REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION")

  REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION_Mouse <- HSsymbol2MMsymbol(REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION,"REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION")
  REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION_Mouse2 <- unique(REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION_Mouse[,2])
  # FeaturePlot(scRNA.SeuObj, features = REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION_Mouse2[1:24], min.cutoff = "q9")
  # FeaturePlot(scRNA.SeuObj, features = REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION_Mouse2[25:50], min.cutoff = "q9")
  # FeaturePlot(scRNA.SeuObj, features =c("Crip1","Arg1","Mt2","Tspo"), min.cutoff = "q9")

  scRNA.SeuObj <- AddModuleScore(object = scRNA.SeuObj, features = as.data.frame(REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION_Mouse2),
                                  ctrl = 5, name = 'REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION')
  FeaturePlot(object = scRNA.SeuObj, features = 'REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION1', split.by = "Cachexia")+
    scale_colour_gradient2(low ="#0077b6", mid = "white", high = "#ef476f",
                           guide = "colourbar",midpoint = 0, labs(fill ="Exp"))

  FeaturePlot(object = scRNA.SeuObj, features = 'REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION1', split.by = "sample",
              cols = c("#0077b6", "white",  "#ef476f"))
  plots2 <- VlnPlot(scRNA.SeuObj, features = c("REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION1"), split.by = "sample", group.by = "celltype",
                   pt.size = 0, combine = FALSE,cols = c("#ef476f", "#e87993",  "#0077b6",  "#3c9ccf"))
  wrap_plots(plots2 = plots2, ncol = 1)

  plots2 <- VlnPlot(scRNA.SeuObj, features = c("REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION1"), split.by = "Cachexia", group.by = "celltype",
                   pt.size = 0, combine = FALSE,cols = c("#ef476f",  "#0077b6"))
  p2 <- wrap_plots(plots2 = plots2, ncol = 1) %>%  BeautifyggPlot(.,AspRat=1,LegPos = c(-0.08, -0.08),AxisTitleSize=1.7,
                                                                XtextSize=18,  YtextSize=,18, xangle = 90,
                                                                LegTextSize = 15)  +
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

  p2


##### REACTOME_ANTIMICROBIAL_PEPTIDES #####
  REACTOME_ANTIMICROBIAL_PEPTIDES <- read.delim(paste0(getwd(),"/GSEA_Geneset/REACTOME_ANTIMICROBIAL_PEPTIDES.txt"),header=T, stringsAsFactors = FALSE)
  REACTOME_ANTIMICROBIAL_PEPTIDES <- data.frame(REACTOME_ANTIMICROBIAL_PEPTIDES[-1,])
  colnames(REACTOME_ANTIMICROBIAL_PEPTIDES) <- c("REACTOME_ANTIMICROBIAL_PEPTIDES")

  REACTOME_ANTIMICROBIAL_PEPTIDES_Mouse <- HSsymbol2MMsymbol(REACTOME_ANTIMICROBIAL_PEPTIDES,"REACTOME_ANTIMICROBIAL_PEPTIDES")
  REACTOME_ANTIMICROBIAL_PEPTIDES_Mouse2 <- unique(REACTOME_ANTIMICROBIAL_PEPTIDES_Mouse[,2])
  # FeaturePlot(scRNA.SeuObj, features = REACTOME_ANTIMICROBIAL_PEPTIDES_Mouse2[1:24], min.cutoff = "q9")
  # FeaturePlot(scRNA.SeuObj, features = REACTOME_ANTIMICROBIAL_PEPTIDES_Mouse2[25:50], min.cutoff = "q9")
  # FeaturePlot(scRNA.SeuObj, features =c("Crip1","Arg1","Mt2","Tspo"), min.cutoff = "q9")

  scRNA.SeuObj <- AddModuleScore(object = scRNA.SeuObj, features = as.data.frame(REACTOME_ANTIMICROBIAL_PEPTIDES_Mouse2),
                                  ctrl = 5, name = 'REACTOME_ANTIMICROBIAL_PEPTIDES')
  FeaturePlot(object = scRNA.SeuObj, features = 'REACTOME_ANTIMICROBIAL_PEPTIDES1', split.by = "Cachexia")+
    scale_colour_gradient2(low ="#0077b6", mid = "white", high = "#ef476f",
                           guide = "colourbar",midpoint = 0, labs(fill ="Exp"))

  FeaturePlot(object = scRNA.SeuObj, features = 'REACTOME_ANTIMICROBIAL_PEPTIDES1', split.by = "sample",
              cols = c("#0077b6", "white",  "#ef476f"))
  plots3 <- VlnPlot(scRNA.SeuObj, features = c("REACTOME_ANTIMICROBIAL_PEPTIDES1"), split.by = "sample", group.by = "celltype",
                    pt.size = 0, combine = FALSE,cols = c("#ef476f", "#e87993",  "#0077b6",  "#3c9ccf"))
  wrap_plots(plots3 = plots3, ncol = 1)

  plots3 <- VlnPlot(scRNA.SeuObj, features = c("REACTOME_ANTIMICROBIAL_PEPTIDES1"), split.by = "Cachexia", group.by = "celltype",
                    pt.size = 0, combine = FALSE,cols = c("#ef476f",  "#0077b6"))
  p3 <- wrap_plots(plots3 = plots3, ncol = 1) %>%  BeautifyggPlot(.,AspRat=1,LegPos = c(-0.08, -0.08),AxisTitleSize=1.7,
                                                                   XtextSize=18,  YtextSize=,18, xangle = 90,
                                                                   LegTextSize = 15)  +
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

  p3


##### REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION #####
  REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION <- read.delim(paste0(getwd(),"/GSEA_Geneset/REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION.txt"),header=T, stringsAsFactors = FALSE)
  REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION <- data.frame(REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION[-1,])
  colnames(REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION) <- c("REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION")

  REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION_Mouse <- HSsymbol2MMsymbol(REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION,"REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION")
  REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION_Mouse2 <- unique(REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION_Mouse[,2])
  # FeaturePlot(scRNA.SeuObj, features = REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION_Mouse2[1:24], min.cutoff = "q9")
  # FeaturePlot(scRNA.SeuObj, features = REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION_Mouse2[25:50], min.cutoff = "q9")
  # FeaturePlot(scRNA.SeuObj, features =c("Crip1","Arg1","Mt2","Tspo"), min.cutoff = "q9")

  scRNA.SeuObj <- AddModuleScore(object = scRNA.SeuObj, features = as.data.frame(REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION_Mouse2),
                                  ctrl = 5, name = 'REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION')
  FeaturePlot(object = scRNA.SeuObj, features = 'REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION1', split.by = "Cachexia")+
    scale_colour_gradient2(low ="#0077b6", mid = "white", high = "#ef476f",
                           guide = "colourbar",midpoint = 0, labs(fill ="Exp"))

  FeaturePlot(object = scRNA.SeuObj, features = 'REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION1', split.by = "sample",
              cols = c("#0077b6", "white",  "#ef476f"))
  plots4 <- VlnPlot(scRNA.SeuObj, features = c("REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION1"), split.by = "sample", group.by = "celltype",
                    pt.size = 0, combine = FALSE,cols = c("#ef476f", "#e87993",  "#0077b6",  "#3c9ccf"))
  wrap_plots(plots4 = plots4, ncol = 1)

  plots4 <- VlnPlot(scRNA.SeuObj, features = c("REACTOME_PLATELET_ACTIVATION_SIGNALING_AND_AGGREGATION1"), split.by = "Cachexia", group.by = "celltype",
                    pt.size = 0, combine = FALSE,cols = c("#ef476f",  "#0077b6"))
  p4 <- wrap_plots(plots4 = plots4, ncol = 1) %>%  BeautifyggPlot(.,AspRat=1,LegPos = c(-0.08, -0.08),AxisTitleSize=1.7,
                                                                   XtextSize=18,  YtextSize=,18, xangle = 90,
                                                                   LegTextSize = 15)  +
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

  p4

##### REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2 #####
  REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2 <- read.delim(paste0(getwd(),"/GSEA_Geneset/REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2.txt"),header=T, stringsAsFactors = FALSE)
  REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2 <- data.frame(REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2[-1,])
  colnames(REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2) <- c("REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2")

  REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2_Mouse <- HSsymbol2MMsymbol(REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2,"REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2")
  REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2_Mouse2 <- unique(REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2_Mouse[,2])
  # FeaturePlot(scRNA.SeuObj, features = REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2_Mouse2[1:24], min.cutoff = "q9")
  # FeaturePlot(scRNA.SeuObj, features = REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2_Mouse2[25:50], min.cutoff = "q9")
  # FeaturePlot(scRNA.SeuObj, features =c("Crip1","Arg1","Mt2","Tspo"), min.cutoff = "q9")

  scRNA.SeuObj <- AddModuleScore(object = scRNA.SeuObj, features = as.data.frame(REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2_Mouse2),
                                  ctrl = 5, name = 'REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2')
  FeaturePlot(object = scRNA.SeuObj, features = 'REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA21', split.by = "Cachexia")+
    scale_colour_gradient2(low ="#0077b6", mid = "white", high = "#ef476f",
                           guide = "colourbar",midpoint = 0, labs(fill ="Exp"))

  FeaturePlot(object = scRNA.SeuObj, features = 'REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA21', split.by = "sample",
              cols = c("#0077b6", "white",  "#ef476f"))
  plots5 <- VlnPlot(scRNA.SeuObj, features = c("REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA21"), split.by = "sample", group.by = "celltype",
                    pt.size = 0, combine = FALSE,cols = c("#ef476f", "#e87993",  "#0077b6",  "#3c9ccf"))
  wrap_plots(plots5 = plots5, ncol = 1)

  plots5 <- VlnPlot(scRNA.SeuObj, features = c("REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA21"), split.by = "Cachexia", group.by = "celltype",
                    pt.size = 0, combine = FALSE,cols = c("#ef476f",  "#0077b6"))
  p5 <- wrap_plots(plots5 = plots5, ncol = 1) %>%  BeautifyggPlot(.,AspRat=1,LegPos = c(-0.08, -0.08),AxisTitleSize=1.7,
                                                                   XtextSize=18,  YtextSize=,18, xangle = 90,
                                                                   LegTextSize = 15)  +
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

  p5


##### REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS #####
  REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS <- read.delim(paste0(getwd(),"/GSEA_Geneset/REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS.txt"),header=T, stringsAsFactors = FALSE)
  REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS <- data.frame(REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS[-1,])
  colnames(REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS) <- c("REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS")

  REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS_Mouse <- HSsymbol2MMsymbol(REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS,"REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS")
  REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS_Mouse2 <- unique(REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS_Mouse[,2])
  # FeaturePlot(scRNA.SeuObj, features = REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS_Mouse2[1:24], min.cutoff = "q9")
  # FeaturePlot(scRNA.SeuObj, features = REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS_Mouse2[25:50], min.cutoff = "q9")
  # FeaturePlot(scRNA.SeuObj, features =c("Crip1","Arg1","Mt2","Tspo"), min.cutoff = "q9")

  scRNA.SeuObj <- AddModuleScore(object = scRNA.SeuObj, features = as.data.frame(REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS_Mouse2),
                                  ctrl = 5, name = 'REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS')
  FeaturePlot(object = scRNA.SeuObj, features = 'REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS1', split.by = "Cachexia")+
    scale_colour_gradient2(low ="#0077b6", mid = "white", high = "#ef476f",
                           guide = "colourbar",midpoint = 0, labs(fill ="Exp"))

  FeaturePlot(object = scRNA.SeuObj, features = 'REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS1', split.by = "sample",
              cols = c("#0077b6", "white",  "#ef476f"))
  plots6 <- VlnPlot(scRNA.SeuObj, features = c("REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS1"), split.by = "sample", group.by = "celltype",
                    pt.size = 0, combine = FALSE,cols = c("#ef476f", "#e87993",  "#0077b6",  "#3c9ccf"))
  wrap_plots(plots6 = plots6, ncol = 1)

  plots6 <- VlnPlot(scRNA.SeuObj, features = c("REACTOME_CREATION_OF_C4_AND_C2_ACTIVATORS1"), split.by = "Cachexia", group.by = "celltype",
                    pt.size = 0, combine = FALSE,cols = c("#ef476f",  "#0077b6"))
  p6 <- wrap_plots(plots6 = plots6, ncol = 1) %>%  BeautifyggPlot(.,AspRat=1,LegPos = c(-0.08, -0.08),AxisTitleSize=1.7,
                                                                  XtextSize=18,  YtextSize=,18, xangle = 90,
                                                                  LegTextSize = 15)  +
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

  p6

##### KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION #####
  KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION <- read.delim(paste0(getwd(),"/GSEA_Geneset/KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION.txt"),header=T, stringsAsFactors = FALSE)
  KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION <- data.frame(KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION[-1,])
  colnames(KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION) <- c("KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION")

  KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION_Mouse <- HSsymbol2MMsymbol(KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION,"KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION")
  KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION_Mouse2 <- unique(KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION_Mouse[,2])
  # FeaturePlot(scRNA.SeuObj, features = KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION_Mouse2[1:24], min.cutoff = "q9")
  # FeaturePlot(scRNA.SeuObj, features = KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION_Mouse2[25:50], min.cutoff = "q9")
  # FeaturePlot(scRNA.SeuObj, features =c("Crip1","Arg1","Mt2","Tspo"), min.cutoff = "q9")

  scRNA.SeuObj <- AddModuleScore(object = scRNA.SeuObj, features = as.data.frame(KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION_Mouse2),
                                  ctrl = 5, name = 'KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION')
  FeaturePlot(object = scRNA.SeuObj, features = 'KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION1', split.by = "Cachexia")+
    scale_colour_gradient2(low ="#0077b6", mid = "white", high = "#ef476f",
                           guide = "colourbar",midpoint = 0, labs(fill ="Exp"))

  FeaturePlot(object = scRNA.SeuObj, features = 'KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION1', split.by = "sample",
              cols = c("#0077b6", "white",  "#ef476f"))
  plots7 <- VlnPlot(scRNA.SeuObj, features = c("KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION1"), split.by = "sample", group.by = "celltype",
                    pt.size = 0, combine = FALSE,cols = c("#ef476f", "#e87993",  "#0077b6",  "#3c9ccf"))
  wrap_plots(plots7 = plots7, ncol = 1)

  plots7 <- VlnPlot(scRNA.SeuObj, features = c("KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION1"), split.by = "Cachexia", group.by = "celltype",
                    pt.size = 0, combine = FALSE,cols = c("#ef476f",  "#0077b6"))
  p7 <- wrap_plots(plots7 = plots7, ncol = 1) %>%  BeautifyggPlot(.,AspRat=1,LegPos = c(-0.08, -0.08),AxisTitleSize=1.7,
                                                                  XtextSize=18,  YtextSize=,18, xangle = 90,
                                                                  LegTextSize = 15)  +
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

  p7


##### REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX #####
  REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX <- read.delim(paste0(getwd(),"/GSEA_Geneset/REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX.txt"),header=T, stringsAsFactors = FALSE)
  REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX <- data.frame(REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX[-1,])
  colnames(REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX) <- c("REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX")

  REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX_Mouse <- HSsymbol2MMsymbol(REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX,"REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX")
  REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX_Mouse2 <- unique(REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX_Mouse[,2])
  # FeaturePlot(scRNA.SeuObj, features = REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX_Mouse2[1:24], min.cutoff = "q9")
  # FeaturePlot(scRNA.SeuObj, features = REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX_Mouse2[25:50], min.cutoff = "q9")
  # FeaturePlot(scRNA.SeuObj, features =c("Crip1","Arg1","Mt2","Tspo"), min.cutoff = "q9")

  scRNA.SeuObj <- AddModuleScore(object = scRNA.SeuObj, features = as.data.frame(REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX_Mouse2),
                                  ctrl = 5, name = 'REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX')
  FeaturePlot(object = scRNA.SeuObj, features = 'REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX1', split.by = "Cachexia")+
    scale_colour_gradient2(low ="#0077b6", mid = "white", high = "#ef476f",
                           guide = "colourbar",midpoint = 0, labs(fill ="Exp"))

  FeaturePlot(object = scRNA.SeuObj, features = 'REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX1', split.by = "sample",
              cols = c("#0077b6", "white",  "#ef476f"))
  plots8 <- VlnPlot(scRNA.SeuObj, features = c("REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX1"), split.by = "sample", group.by = "celltype",
                    pt.size = 0, combine = FALSE,cols = c("#ef476f", "#e87993",  "#0077b6",  "#3c9ccf"))
  wrap_plots(plots8 = plots8, ncol = 1)

  plots8 <- VlnPlot(scRNA.SeuObj, features = c("REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX1"), split.by = "Cachexia", group.by = "celltype",
                    pt.size = 0, combine = FALSE,cols = c("#ef476f",  "#0077b6"))
  p8 <- wrap_plots(plots8 = plots8, ncol = 1) %>%  BeautifyggPlot(.,AspRat=1,LegPos = c(-0.08, -0.08),AxisTitleSize=1.7,
                                                                  XtextSize=18,  YtextSize=,18, xangle = 90,
                                                                  LegTextSize = 15)  +
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

  p8



##### Export PDF file #####
  pdf(file = paste0(Save.Path,"/PBMC_GSEA_Violin.pdf"),
      width = 7, height = 7 )
    p1
    p2
    p3
    p4
    p5
    p6
    p7
    p8
  dev.off() # graphics.off()





