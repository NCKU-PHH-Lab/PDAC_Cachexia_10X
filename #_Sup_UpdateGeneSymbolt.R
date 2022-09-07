## UpdateSymbolList: Get updated synonyms for gene symbols
## Ref: https://rdrr.io/cran/Seurat/man/UpdateSymbolList.html
library(Seurat)

## Not run:
GeneSybmolThesarus(symbols = c("FAM64A"))

## End(Not run)

## Not run:
UpdateSymbolList(symbols = cc.genes$s.genes)
UpdateSymbolList(symbols = "SEPT1")
UpdateSymbolList(symbols = "1-Sep")

## End(Not run)
PBMC.combined@assays[["RNA"]]@counts@Dimnames[[1]] <- UpdateSymbolList(PBMC.combined@assays[["RNA"]]@counts@Dimnames[[1]])

A <- PBMC.combined@assays[["RNA"]]@counts@Dimnames[[1]]
B <- UpdateSymbolList(PBMC.combined@assays[["RNA"]]@counts@Dimnames[[1]])

summary(A==B)

