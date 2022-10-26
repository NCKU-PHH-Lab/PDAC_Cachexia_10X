## Tutorial:Gene Set Enrichment Analysis (fgsea)
## Ref: https://github.com/hamidghaedi/Enrichment-Analysis
## Ref: https://www.biostars.org/p/467197/
## RNA-seq analysis in R
## Ref: https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

  Save.Path <- c("D:/Dropbox/##_GitHub/##_PHH_Lab/PDAC_Cachexia_10X/2022-10-17_PBMC_Main")
  # SampleType = "PBMC"

##### Load Packages #####
  # if(!require("tidyverse")) install.packages("tidyverse")
  # library(tidyverse)

  #### Basic installation ####
  ## Check whether the installation of those packages is required from basic
  Package.set <- c("tidyverse","ggplot2","Seurat","SeuratData","patchwork","plyr","eoffice","DT")
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      install.packages(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)


  #### BiocManager installation ####
  ## Set the desired organism
  # organism = "org.Hs.eg.db" ## c("org.Hs.eg.db","org.Mm.eg.db")   ##  c("org.Dm.eg.db")

  ## Check whether the installation of those packages is required from BiocManager
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  Package.set <- c("clusterProfiler","enrichplot","pathview") # c(organism,"fgsea","clusterProfiler","enrichplot","pathview")
  for (i in 1:length(Package.set)) {
    if (!requireNamespace(Package.set[i], quietly = TRUE)){
      BiocManager::install(Package.set[i])
    }
  }
  ## Load Packages
  lapply(Package.set, library, character.only = TRUE)
  rm(Package.set,i)

  #### GitHub installation ####
  if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
  devtools::install_github("lionel-/ggstance")
  library(ggstance)

##### Function setting #####
  ## Call function
  source("FUN_Beautify_ggplot.R")

##### Load RData* #####
  load(paste0(Save.Path,"/09_1_GSEA_Analysis_(SPA).RData"))

  ## INTCHG: Interchangeable
  ## SubType Setting
    if(SampleType == "PBMC"){
      # For PBMC
      scRNA.SeuObj <- PBMC.combined

    }else if(SampleType == "SC"){
      # For SC
      scRNA.SeuObj <- SC.combined

    }


  #### Clean up Object ####
  rm(list=setdiff(ls(), c("scRNA.SeuObj","SampleType","Save.Path",
                          str_subset(objects(), pattern = "GSEA"))))
  ## Save Ori
  scRNA_Ori.SeuObj <- scRNA.SeuObj

  ## Clean up data (Delete other type)
  scRNA.SeuObj <- scRNA.SeuObj[,!grepl("Other", scRNA.SeuObj@meta.data[["celltype"]] )]

  ## Clean up data (Delete other type)

#####***************************************************************************#####
#####*  Plot by previous results *#####
##### Set condition #####
  SubType = "Mac3"
  Set_GSEA_FDR <- 0.01
  Set_GSEA_NES <- 1

  Set_GE_FDR <- 0.01
  Set_GE_LogFC <- 1

  ## Barplot setting
  NumGenesetsPlt=15

##### Extract df #####
  ## Extract by SubType
  GSEA_Sub.df <- GSEA_Large.df.TOP[GSEA_Large.df.TOP$PhenoType %in% SubType,]

  ## Filter by FDR & NES
  GSEA_Sub.df <- GSEA_Sub.df[GSEA_Sub.df$padj <= Set_GSEA_FDR & abs(GSEA_Sub.df$NES) > 1,]

##### Plot #####
  Barplot <- ggplot(GSEA_Sub.df, aes(NES, fct_reorder(pathway, NES), fill = padj), showCategory=(NumGenesetsPlt*2)) +
    geom_barh(stat='identity') +
    scale_fill_continuous(low = "#d45772", high = "#3b74bf", guide = guide_colorbar(reverse=TRUE)) +
    theme_minimal() + ylab(NULL)

  Barplot <- Barplot %>% BeautifyggPlot(LegPos = c(0.9, 0.15), AxisTitleSize=1.7, YtextSize=11,OL_Thick = 1.5)
  Barplot

#####***************************************************************************#####
#####*  Rerun GSEA *#####
##### Load Package #####
  library(DESeq2)
  library(org.Hs.eg.db)
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(fgsea)
  library(ggplot2)
  library(reshape2)
  library(ComplexHeatmap)
  library(circlize)

##### Function setting #####
  ## Call function
  source("FUN_Beautify_ggplot.R")
  source("FUN_Find_Markers.R")
  source("FUN_VolcanoPlot.R")
  source("FUN_ggPlot_vline.R")
  source("FUN_GSEA_ANAL.R")
  source("FUN_DistrPlot.R")

##### Import setting and Import #####
  ## Import GSEA gene sets
  # InputGSEA <- "GSEA_Geneset_Pathway_3Database_WithoutFilter.txt"
  # InputGSEA <- "m5_go_bp_v0_3_symbols.gmt"  # InputGSEA <- "m2.all.v0.3.symbols.gmt"

  InputGSEA <- "m2.all.v0.3.symbols.gmt"
  InFOLName_GSEA <- "Input_Genesets"
  # Pathway.all <- read.delim2(paste0(getwd(),"/",InFOLName_GSEA,"/",InputGSEA),
  #                            col.names = 1:max(count.fields(paste0(getwd(),"/",InFOLName_GSEA,"/",InputGSEA))),
  #                            header = F,sep = "\t")
  Pathway.all <- gmtPathways(paste0(getwd(),"/",InFOLName_GSEA,"/",InputGSEA))

##### Conditions setting* #####
  Group_Mode <- "GoupByPheno"   # c("GoupByPheno","GoupByGeneExp")
  TarGene_name <- "Chil3"
  PhenoGrp_name1 <- "Cachexia"
  PhenoGrp_name2 <- c("EOCX","PreCX")

  GeneExpSet.lt <- list(GeneExpMode = "Mean", # c("Mean","Mean1SD","Mean2SD","Mean3SD","Median","Quartile","Customize"))
                        UpCutoff = 1, LowerCutoff = 1)

  if(Group_Mode == "GoupByGeneExp"){
    ## Group by GeneExp
    AnnoSet.lt <- list(GroupType = TarGene_name, GroupCompare = c("High","Low") )   ## DEG by GeneExp group
  }else{
    ## Group by Pheno
    AnnoSet.lt <- list(GroupType = PhenoGrp_name1, GroupCompare = PhenoGrp_name2 )
  }

  Thr.lt <- list(LogFC = c("logFC",1), pVal = c("PValue",0.05) )

##### Current path and new folder setting* #####
  ProjectName = "CC10X"

  if(Group_Mode == "GoupByGeneExp")
  {
    # ExportAnno = "Chil3Mean_PathM2"
    ExportAnno = paste0(TarGene_name,"_",GeneExpSet.lt$GeneExpMode,"_",SubType)
  }else
  {
    # ExportAnno = "Recur2Prim"
    ExportAnno = paste0(Group_Mode,"_",paste0(PhenoGrp_name2[1],PhenoGrp_name2[2],"_",SubType))
  }

  ExportName = paste0(ProjectName,"_",SampleType,"_",ExportAnno)

  Version = paste0(Sys.Date(),"_",ProjectName,"_",SampleType,"_", ExportAnno)
  SaveSub.Path = paste0(getwd(),"/",Version)
  ## Create new folder
  if (!dir.exists(SaveSub.Path)){
    dir.create(SaveSub.Path)
  }

##### Update the genename ####
  ## Update the genename ##* Take very long time
  UpdateGene <- "No"  # UpdateGene <- c("Yes","No")
  if(UpdateGene == "Yes"){
    row.names(GeneExp.df) <- UpdateSymbolList(row.names(GeneExp.df))
  }

#************************************************************************************************************************#
##### Data preprocess setting #####
  ## Extract data from scRNA.SeuObj
  ## Old version (Without normalizaiton) ## GeneExp.df <- scRNA.SeuObj@assays[["RNA"]]@counts %>% as.data.frame()
  GeneExp.df <- GetAssayData(scRNA.SeuObj, assay = "RNA", slot = "data") %>% as.data.frame() # normalized data matrix

  Anno.df <- scRNA.SeuObj@meta.data
  Anno.df <- data.frame(ID=row.names(Anno.df), Anno.df)
  Anno.df <- left_join(data.frame("ID"=colnames(GeneExp.df)),
                       Anno.df)
  row.names(Anno.df) <- Anno.df[,1]

  ## Select Pheno column
  Anno_Ori.df <- Anno.df
  colnames(Anno.df)

  # PhenoColKeep.set <- c("X_INTEGRATION","X_PATIENT","histological_type","sample_type","gender")
  # Anno.df <- Anno.df[,c(PhenoColKeep.set)]
  # colnames(Anno.df)
  #
  # head(Anno.df)

  # ## Select Pheno row
  # PhenoRowKeep.set <- list(col="Cachexia" ,row=c("EOCX"))
  # Anno.df <- Anno.df[Anno.df[,PhenoRowKeep.set[["col"]]] %in% PhenoRowKeep.set[["row"]], ]
  #
  # GeneExp.df <- GeneExp.df[,colnames(GeneExp.df) %in% Anno.df[,1] ]
  # rm(PhenoRowKeep.set)

  ## Select Pheno row2
  PhenoRowKeep.set <- list(col="celltype" ,row=c(SubType))
  Anno.df <- Anno.df[Anno.df[,PhenoRowKeep.set[["col"]]] %in% PhenoRowKeep.set[["row"]], ]

  GeneExp.df <- GeneExp.df[,colnames(GeneExp.df) %in% Anno.df[,1] ]
  rm(PhenoRowKeep.set)

  # ## Delete specific cell type
  # ## Clean up data
  # Anno.df <- Anno.df[!grepl("Other", Anno.df$celltype),]
  # GeneExp.df <- GeneExp.df[,colnames(GeneExp.df) %in% Anno.df$ID]



#************************************************************************************************************************#
##### Visualization #####
  source("FUN_DistrPlot.R")
  ##### Group by gene expression 1: CutOff by total  #####
  Plot.DistrPlot <- FUN_DistrPlot(GeneExp.df,
                                  TarGeneName = TarGene_name, GroupSet = GeneExpSet.lt,
                                  Save.Path = SaveSub.Path, ExportName = ExportName)
  Plot.DistrPlot_SD_Q <- Plot.DistrPlot[["TGeneDen_SD_Q.p"]]
  Plot.DistrPlot_SD_Q


#************************************************************************************************************************#
##### Grouping #####
  if(Group_Mode == "GoupByGeneExp")
  {
    source("FUN_Group_GE.R")
    ##### Group by gene expression 1: CutOff by total  #####
    GeneExp_group.set <- FUN_Group_GE(GeneExp.df, Anno.df,
                                      TarGeneName = TarGene_name, GroupSet = GeneExpSet.lt,
                                      Save.Path = SaveSub.Path, ExportName = ExportName)
    Anno.df <- GeneExp_group.set[["AnnoNew.df"]]
    GeneExp_high.set <- GeneExp_group.set[["GeneExp_high.set"]]
    GeneExp_low.set <- GeneExp_group.set[["GeneExp_low.set"]]

    ##### Group by gene expression 2: CutOff by Comparison #####
    ## FUN Comparison (Visualization and value)


  }

  ##### Group by phenotype #####


#************************************************************************************************************************#
##### Run Enrichment analysis in R #####
  #### Run DEG ####
  source("FUN_DEG_Analysis.R")
  DEG_ANAL.lt <- FUN_DEG_Analysis(GeneExp.df, Anno.df,
                              GroupType = AnnoSet.lt[["GroupType"]], GroupCompare = AnnoSet.lt[["GroupCompare"]],
                              ThrSet = Thr.lt,
                              TarGeneName = TarGene_name, GroupMode = GeneExpSet.lt, SampleID = "ID",
                              Save.Path = SaveSub.Path, ExportName = ExportName, AnnoName = "AvB")
  DE_Extract.df <- DEG_ANAL.lt[["DE_Extract.df"]]


  #### Old version ####
  # #### Run GSEA ####
  # source("FUN_GSEA_ANAL.R")
  #
  # GSEA_Result.lt <- FUN_GSEA_ANAL(DE_Extract.df, CMGeneSet = Pathway.all,
  #                                 NumGenesetsPlt=15,
  #                                 TarGeneName = TarGene_name,
  #                                 ThrSet = Thr.lt, Species = "Homo sapiens", # Speices type can check by msigdbr_species()
  #                                 Save.Path = SaveSub.Path, ExportName = ExportName, AnnoName = "Path")


  #************************************************************************************************************************#

  ## Tutorial:Gene Set Enrichment Analysis (fgsea)
  ## Ref: https://github.com/hamidghaedi/Enrichment-Analysis
  ## Ref: https://www.biostars.org/p/467197/
  summary(DE_Extract.df)
  res <- DE_Extract.df
  res$SYMBOL <- res$Gene

  ranks <- DE_Extract.df$logFC
  names(ranks) <- DE_Extract.df$Gene

  fgseaRes <- fgsea(Pathway.all, ranks, minSize=15, maxSize = 500)

  # Tidy the results:
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES)) # order by normalized enrichment score (NES)

  # To see what genes are in each of these pathways:
  gene.in.pathway <- Pathway.all %>%
    enframe("pathway", "SYMBOL") %>%
    unnest(cols = c(SYMBOL)) %>%
    inner_join(res, by="SYMBOL")

  #______________________VISUALIZATION______________________________#

  #__________bar plot _______________#
  # Plot the normalized enrichment scores.
  #Color the bar indicating whether or not the pathway was significant:
  fgseaResTidy$adjPvalue <- ifelse(fgseaResTidy$padj <= 0.05, "significant", "non-significant")
  cols <- c("non-significant" = "grey", "significant" = "red")
  ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
    geom_col() +
    scale_fill_manual(values = cols) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="Hallmark pathways Enrichment Score from GSEA")

  ##### Plot #####
  NumGenesetsPlt=10
  fgseaResTidy_Ori <- fgseaResTidy
  fgseaResTidy <- fgseaResTidy[fgseaResTidy$padj <= Set_GSEA_FDR & abs(fgseaResTidy$NES) >  Set_GSEA_NES,]



  Barplot <- ggplot(fgseaResTidy, aes(NES, fct_reorder(pathway, NES), fill = padj), showCategory=(NumGenesetsPlt*2)) +
    geom_barh(stat='identity') +
    scale_fill_continuous(low = "#d45772", high = "#3b74bf", guide = guide_colorbar(reverse=TRUE)) +
    theme_minimal() + ylab(NULL)

  Barplot <- Barplot %>% BeautifyggPlot(LegPos = c(0.9, 0.15), AxisTitleSize=1.7, YtextSize=11,OL_Thick = 1.5)
  Barplot


  #__________Enrichment  Plot_______#
  # Enrichment plot for E2F target gene set
  plotEnrichment(pathway = Pathway.all[["JARDIM_PERASSI_TRIPLE_NEGATIVE_BREAST_CANCER_MOUSE_XENOGRAFT_MELATONIN_UP"]], ranks)
  graphics.off()

  plotGseaTable(Pathway.all[fgseaRes$pathway[fgseaRes$padj < 0.05]], ranks, fgseaRes,
                gseaParam=0.5)
  graphics.off()


  #________ Heatmap Plot_____________#
  # pathways with significant enrichment score
  sig.path <- fgseaResTidy$pathway[fgseaResTidy$adjPvalue == "significant"]
  sig.gen <- unique(na.omit(gene.in.pathway$SYMBOL[gene.in.pathway$pathway %in% sig.path]))

  ### create a new data-frame that has '1' for when a gene is part of a term, and '0' when not
  h.dat <- dcast(gene.in.pathway[, c(1,2)], SYMBOL~pathway)
  rownames(h.dat) <- h.dat$SYMBOL
  h.dat <- h.dat[, -1]

  h.dat <- h.dat[rownames(h.dat) %in% sig.gen, ]
  h.dat <- h.dat[, colnames(h.dat) %in% sig.path]

  h.dat_ori <- h.dat
  for (i in 1:nrow( h.dat)) {
    for (j in 1:ncol( h.dat)) {
      if(is.na( h.dat[i,j])){
         h.dat[i,j] = 0
      }else{
         h.dat[i,j] = 1
      }
    }
  }

  h.dat_rownames <- rownames(h.dat)
  h.dat <- data.frame(apply(h.dat, 2, function(x) as.numeric(as.character(x))))
  rownames(h.dat) <- h.dat_rownames

  # keep those genes with 3  or more occurnes
  table(data.frame(rowSums(h.dat)))

  h.dat <- h.dat[data.frame(rowSums(h.dat)) >= 3, ]

  #
  topTable <- res[res$SYMBOL %in% rownames(h.dat), ]
  rownames(topTable) <- topTable$SYMBOL
  # match the order of rownames in toptable with that of h.dat
  topTableAligned <- topTable[which(rownames(topTable) %in% rownames(h.dat)),]
  topTableAligned <- topTableAligned[match(rownames(h.dat), rownames(topTableAligned)),]
  all(rownames(topTableAligned) == rownames(h.dat))

  topTableAligned <- topTableAligned[topTableAligned$FDR <= Set_GE_FDR &
                                       abs(topTableAligned$logFC) >= Set_GE_LogFC,]
  h.dat <- h.dat[rownames(h.dat) %in% topTableAligned$Gene,]

  # colour bar for -log10(adjusted p-value) for sig.genes
  dfMinusLog10FDRGenes <- data.frame(-log10(
    topTableAligned[which(rownames(topTableAligned) %in% rownames(h.dat)), 'FDR']))
  dfMinusLog10FDRGenes[dfMinusLog10FDRGenes == 'Inf'] <- 0

  # colour bar for fold changes for sigGenes
  dfFoldChangeGenes <- data.frame(
    topTableAligned[which(rownames(topTableAligned) %in% rownames(h.dat)), 'logFC'])

  # merge both
  dfGeneAnno <- data.frame(dfMinusLog10FDRGenes, dfFoldChangeGenes)
  colnames(dfGeneAnno) <- c('Gene score', 'Log2FC')
  dfGeneAnno[,2] <- ifelse(dfGeneAnno$Log2FC > 0, 'Up-regulated',
                           ifelse(dfGeneAnno$Log2FC < 0, 'Down-regulated', 'Unchanged'))
  colours <- list(
    'Log2FC' = c('Up-regulated' = 'royalblue', 'Down-regulated' = 'yellow', 'Unchanged' ='black'))
  haGenes <- rowAnnotation(
    df = dfGeneAnno,
    col = colours,
    width = unit(1,'cm'),
    annotation_name_side = 'top')

  # Now a separate color bar for the GSEA enrichment padj. This will
  # also contain the enriched term names via annot_text()

  # colour bar for enrichment score from fgsea results
  dfEnrichment <- fgseaRes[, c("pathway", "NES")]
  dfEnrichment <- dfEnrichment[dfEnrichment$pathway %in% colnames(h.dat)]
  dd <- dfEnrichment$pathway
  dfEnrichment <- dfEnrichment[, -1]
  rownames(dfEnrichment) <- dd
  colnames(dfEnrichment) <- 'Normalized\n Enrichment score'
  haTerms <- HeatmapAnnotation(
    df = dfEnrichment,
    Term = anno_text(
      colnames(h.dat),
      rot = 45,
      just = 'right',
      gp = gpar(fontsize = 12)),
    annotation_height = unit.c(unit(1, 'cm'), unit(8, 'cm')),
    annotation_name_side = 'left')
  # now generate the heatmap
  hmapGSEA <- Heatmap(h.dat,
                      name = 'GSEA hallmark pathways enrichment',
                      split = dfGeneAnno[,2],
                      col = c('0' = 'white', '1' = 'forestgreen'),
                      rect_gp = gpar(col = 'grey85'),
                      cluster_rows = TRUE,
                      show_row_dend = TRUE,
                      row_title = 'Top Genes',
                      row_title_side = 'left',
                      row_title_gp = gpar(fontsize = 11, fontface = 'bold'),
                      row_title_rot = 90,
                      show_row_names = TRUE,
                      row_names_gp = gpar(fontsize = 11, fontface = 'bold'),
                      row_names_side = 'left',
                      row_dend_width = unit(35, 'mm'),
                      cluster_columns = TRUE,
                      show_column_dend = TRUE,
                      column_title = 'Enriched terms',
                      column_title_side = 'top',
                      column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
                      column_title_rot = 0,
                      show_column_names = FALSE,
                      show_heatmap_legend = FALSE,
                      clustering_distance_columns = 'euclidean',
                      clustering_method_columns = 'ward.D2',
                      clustering_distance_rows = 'euclidean',
                      clustering_method_rows = 'ward.D2',
                      bottom_annotation = haTerms)


  draw(hmapGSEA + haGenes,
       heatmap_legend_side = 'right',
       annotation_legend_side = 'right')


  tiff("GSEA_enrichment_2.tiff", units="in", width=15, height=22, res=400)
  draw(hmapGSEA + haGenes,
       heatmap_legend_side = 'right',
       annotation_legend_side = 'right')
  dev.off()





  #### Run ORA ####
  ## FUN ORA

#************************************************************************************************************************#
# ##### Build files for GSEA official input #####
#   source("FUN_GSEA_ForOFFL.R")
#
#   FUN_GSEA_ForOFFL(GeneExp.df, Group1 = GeneExp_high.set, Group2 = GeneExp_low.set,
#                    GroupMode = Group_Mode,
#                    TarGeneName = TarGene_name, GeneExpSet = GeneExpSet.lt,
#                    Save.Path = SaveSub.Path, ExportName = ExportName,
#                    AnnoName="Recur2Prim")
#
# ##### Build files for Metascape official input #####
#
#
#
# #### Save RData ####
#   save.image(paste0(Save.Path,"/GseaGo_",ExportName,".RData"))




