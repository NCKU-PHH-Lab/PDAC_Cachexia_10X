##### Script info #####
  # Project: cachexia scRNA analysis
  # Data from Charlene
  # Aim: to generate heatmap of EOCX and PreCX marker pattern
  # Notice the setting when you see the symbol "*"

##### Packages #####
  library(tidyverse)
  library(ComplexHeatmap)
  library(ggplot2)
  library(circlize)

##### Input setting* #####
  # The directory path of raw dataset for plot
  file.dir <- "2022-09-09_Results_1stSubmission/datasets_for_Heatmap"

  # The file format
  file.form <- ".csv" # default: csv

##### Data import #####
  file.set <- list.files(file.dir) %>% str_detect(file.form) %>% list.files(file.dir)[.]
  raw_data.lt <- list()
  if(file.form == ".tsv"){
    sep = "\t"
  }else{
    sep = ","
  }
  for (f in file.set) {
    raw_data.lt[[f]] <- read.csv(
      paste0(file.dir,"/",f),
      sep =sep
    )
    rm(f)
  }
  rm(sep)

  # rename the data frame
  names(raw_data.lt) <- seq(1,length(raw_data.lt))
  for (i in 1:length(raw_data.lt)) {
    colnames(raw_data.lt[[i]])[1:4] <- c("Cell_type", "State", "Strategy", "ID")
    rm(i)
  }

##### Arrange setting* #####

  # exclude unwanted cell_type*
  unwanted.set <- c("Other1", "Other2")
  raw_data.lt[[1]] <- raw_data.lt[[1]][!raw_data.lt[[1]]$Cell_type %in% unwanted.set,]
  rm(unwanted.set)

  # generate the unique cell_type
  unarrange.lt <- list()
  for (i in 1:length(raw_data.lt)) {
    unarrange.lt[[i]] <- raw_data.lt[[i]]["Cell_type"] %>% unique()
    rm(i)
  }

  # check the unique cell_type
  for(i in 1:length(unarrange.lt)){
    print(unarrange.lt[[i]])
    rm(i)
  } # It's helpful for arranging the cell type

  # generate the arrange list*
  arrange.lt <- list(
    "Cell_type_1" = c("All", "T", "CD4+T", "CD8+T", "B", "Mac1", "Mac2", "Mac3",
                      "Mast", "NK", "Neu", "Ery"),
    "Cell_type_2" = c("All", "Fib1", "Fib2", "Fib3", "Duc1", "Duc2", "Duc3", "Duc4",
                      "Duc5", "Duc6", "Mac1", "Mac2", "Mac3", "Mac4", "Mac5"),
    "State" = c("EOCX", "PreCX"),
    "Strategy" = c("SPA", "I", "F", "M")
  )
  ## check
  for(i in 1:length(unarrange.lt)){
    print(sort(unarrange.lt[[i]][[1]])==sort(arrange.lt[[i]]))
    rm(i)
  } # The printed result should be all "TRUE"
  rm(unarrange.lt)
##### Data preprocess #####

  # factor transform of the features

  for (i in 1:length(raw_data.lt)) {
    raw_data.lt[[i]]$Cell_type <- factor(raw_data.lt[[i]]$Cell_type,
                                         levels = arrange.lt[paste0("Cell_type_",i)][[1]])
    raw_data.lt[[i]]$State <- factor(raw_data.lt[[i]]$State,
                                         levels = arrange.lt["State"][[1]])
    raw_data.lt[[i]]$Strategy <- factor(raw_data.lt[[i]]$Strategy,
                                         levels = arrange.lt["Strategy"][[1]])
    rm(i)
  }

  # generate dataframe for heatmap making
  hp_data.lt <- list()

  ## Total
    hp_data.lt[["Total"]] <- list()
    for (i in 1:length(raw_data.lt)) {
      temp.df <- raw_data.lt[[i]][raw_data.lt[[i]]$Cell_type=="All",] %>%
        arrange(Strategy, State)
      row.names(temp.df) <- temp.df$ID
      hp_data.lt[["Total"]][[i]] <- as.matrix(temp.df[-(1:4)]) %>% t()
      rm(i, temp.df)
    }

  ## Different strategies
    for (STR in arrange.lt$Strategy) {
      hp_data.lt[[STR]] <- list()

      for (i in 1:length(raw_data.lt)) {
        temp.df <- raw_data.lt[[i]][raw_data.lt[[i]]$Strategy==STR & raw_data.lt[[i]]$Cell_type!="All",] %>%
          arrange(Cell_type, State)
        row.names(temp.df) <- temp.df$ID
        hp_data.lt[[STR]][[i]] <- as.matrix(temp.df[-(1:4)]) %>% t()
        rm(i, temp.df)
      }
      rm(STR)
    }

  ## Ranking the gene set
    ranked_gene.lt <- list()
    for (i in 1:length(hp_data.lt[["Total"]])) {
      ranked_gene.lt[[i]] <- hp_data.lt[["Total"]][[i]] %>% as.data.frame() %>% abs() %>%
        arrange(across(starts_with("All"),desc)) %>% row.names()
      rm(i)
    }

    for (STR in names(hp_data.lt)) {
      for (i in 1:length(hp_data.lt[[STR]])) {
        hp_data.lt[[STR]][[i]] <- hp_data.lt[[STR]][[i]][ranked_gene.lt[[i]],]
        rm(i)
      }
      rm(STR)
    }

##### Making Heatmap* #####
    # set color*
    col_fun_T = colorRamp2(c(-10, 0, 10), c("#0077b6", "#FFFFFF", "#ef476f"))
    col_fun_S = colorRamp2(c(-1, 0, 1), c("#0077b6", "#FFFFFF", "#ef476f"))
    is_gray_fence = TRUE

    # making plot
    ht_opt$message = FALSE
    for (i in 1:length(file.set)) {
      ## set file name
      file_name <- gsub(file.form, ".pdf", file.set[[i]]) %>%
        paste0("plot/",Sys.Date(), "_",.)
      ## generate pdf
      pdf(file_name)
      for(j in 1:length(hp_data.lt)){
        ### color check
        if(j == 1){
          col_fun <- col_fun_T
        }else{
          col_fun <- col_fun_S
        }
        ### gray fence check

        if(is_gray_fence & j!=1){
          mt <- hp_data.lt[[j]][[i]]
          mt[,c(seq(3,ncol(mt),4),seq(4,ncol(mt),4))][mt[,c(seq(3,ncol(mt),4),seq(4,ncol(mt),4))]==0] <- NA
        }else{
          mt <- hp_data.lt[[j]][[i]]
        }


        ### plotting
        hp <- Heatmap(
          mt,
          na_col = "#F0F0F0",
          cluster_rows = F,
          cluster_columns = F,
          show_column_names = F,
          show_row_names = F,
          name = "count",
          col =col_fun,
          show_heatmap_legend = F
        )
        draw(hp)
        hp <- Heatmap(
          mt,
          na_col = "#F0F0F0",
          column_title = names(hp_data.lt)[[j]],
          column_title_side = "top",
          cluster_rows = F,
          cluster_columns = F,
          show_column_names = T,
          show_row_names = F,
          name = "count",
          col =col_fun,
          show_heatmap_legend = T
        )
        draw(hp)
        rm(j,mt,hp)
      }
      graphics.off()
      rm(i)
    }
