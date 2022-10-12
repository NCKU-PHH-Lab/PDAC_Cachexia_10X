##### Load Packages  #####
  library(infercnv)

##### Input Ref #####
  Input_Ref.Path <- paste0(getwd(),"/Input_inferCNV_Ref/mm10_genomic_mapinfo_one.tsv")

##### Current path and new folder setting #####
  SaveinferCNV.Path <- paste0(Save.Path,"/", SampleType, "inferCNV")

  Ref_group.Set = c("Fib1_PreCX","Fib2_PreCX","Fib3_PreCX")

# ##### Example #####
#
#   infercnv_obj = CreateInfercnvObject(raw_counts_matrix=system.file("extdata", "oligodendroglioma_expression_downsampled.counts.matrix.gz", package = "infercnv"),
#                                       annotations_file=system.file("extdata", "oligodendroglioma_annotations_downsampled.txt", package = "infercnv"),
#                                       delim="\t",
#                                       gene_order_file=system.file("extdata", "gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt", package = "infercnv"),
#                                       ref_group_names=c("Microglia/Macrophage","Oligodendrocytes (non-malignant)"))
#
#   infercnv_obj = infercnv::run(infercnv_obj,
#                                cutoff=1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
#                                out_dir=tempfile(),
#                                cluster_by_groups=TRUE,
#                                denoise=TRUE,
#                                HMM=TRUE)

#######################
  ## Exp
  SC_Exp.df <- as.data.frame(SC.combined@assays[["RNA"]]@counts)
  SC_Exp.mx <- as.matrix(SC_Exp.df)

  ## Ano
  SC_Ano_CT.df <- data.frame(celltype = SC.combined@meta.data[["celltype.Cachexia"]])
  row.names(SC_Ano_CT.df) <- SC.combined@assays[["RNA"]]@counts@Dimnames[[2]]
  colnames(SC_Ano_CT.df)<- ""
  SC_Ano_CT.mx <- as.matrix(SC_Ano_CT.df)
  #

  # SC_Ano_ID.df <- data.frame(ID = SC.combined@assays[["RNA"]]@counts@Dimnames[[2]])
  #
  # SC_Ano.df <- cbind(SC_Ano_ID.df,SC_Ano_CT.df)
  # row.names(SC_Ano.df)<- SC.combined@assays[["RNA"]]@counts@Dimnames[[2]]
  # #colnames(SC_Ano.df)<- ""
  # # SC_Ano.df$celltype <- as.character(SC_Ano.df$celltype)
  # # SC_Ano.df$ID <- as.character(SC_Ano.df$ID)
  # SC_Ano.mx <- as.matrix(SC_Ano.df)
  # # colnames(SC_Ano.mx)[1] <- 1
  # # colnames(SC_Ano.mx)[2] <- 2
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix = SC_Exp.mx,
                                      annotations_file = SC_Ano_CT.mx,
                                      delim="\t",
                                      gene_order_file = Input_Ref.Path,
                                      ref_group_names = Ref_group.Set)

  infercnv_obj = infercnv::run(infercnv_obj,
                                cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                out_dir=SaveinferCNV.Path, # out_dir=tempfile(),
                                cluster_by_groups=TRUE,
                                denoise=TRUE,
                                HMM=TRUE,
                                hclust_method='ward.D2')
