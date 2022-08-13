  

  GO.OBO <- read.delim2(paste0(getwd(),"/go.obo"))
  GO.OBO.S <- data.frame(GO.OBO[1:1000,])
  GO.OBO.S_Try_id <- GO.OBO.S[grep("id:", GO.OBO.S$GO.OBO.1.1000...), ]
  GO.OBO.S_Try_name <- GO.OBO.S[grep("name:", GO.OBO.S$GO.OBO.1.1000...), ]