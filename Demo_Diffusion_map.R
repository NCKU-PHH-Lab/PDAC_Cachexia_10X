## Ref: https://broadinstitute.github.io/2019_scWorkshop/functional-pseudotime-analysis.html
## Ref: https://rdrr.io/cran/diffusionMap/

##### Load Packages ####
library(tidyverse)
library(Seurat)

GeneExp.df <- GetAssayData(scRNA.SeuObj, assay = "RNA", slot = "data") %>% as.data.frame() # normalized data matrix
GeneExp.mtx <- GeneExp.df %>% as.matrix()

# install.packages("diffusionMap")
# library("diffusionMap")
# dm <- diffuse(t(GeneExp.mtx[1:100,1:100]))


# install.packages("gam")
# library(gam)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("slingshot")
# library(slingshot)
# library(Seurat)
# install.packages("corrplot")
# library(corrplot)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("destiny")
library(destiny)

dm <- DiffusionMap(t(GeneExp.mtx[,1:1000]))
dm <- DiffusionMap(t(GeneExp.mtx))


# Plot diffusion component 1 vs diffusion component 2 (DC1 vs DC2).
tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                  DC2 = eigenvectors(dm)[, 2])
ggplot(tmp, aes(x = DC1, y = DC2)) +
  geom_point() +# scale_color_tableau() +
  xlab("Diffusion component 1") +
  ylab("Diffusion component 2") +
  theme_classic()

