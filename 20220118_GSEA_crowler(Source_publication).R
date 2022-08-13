##### Version info #####
  # _                           
  # platform       x86_64-w64-mingw32          
  # arch           x86_64                      
  # os             mingw32                     
  # system         x86_64, mingw32             
  # status                                     
  # major          4                           
  # minor          0.5                         
  # year           2021                        
  # month          03                          
  # day            31                          
  # svn rev        80133                       
  # language       R                           
  # version.string R version 4.0.5 (2021-03-31)
  # nickname       Shake and Throw  

##### Crawler function #####
GSEA_crowler <- function(GeneSet.lt, core){
  
  require(parallel)
  require(magrittr)
  require(utils)

  cl <- makeCluster(core)
  
  Crowler <- function(GeneSet){
    require(rvest)
    require(xml2)
    # print(names(GeneSet)[1])
    
    GeneSet.path	 <- '/html/body/div[3]/div[4]/div/table/tbody/tr/td[3]/table/tbody/tr[1]/td'
    Exact_source.path	 <- '/html/body/div[3]/div[4]/div/table/tbody/tr/td[3]/table/tbody/tr[7]/td'
   
    
    output.df <- data.frame(
      "Gene_Set" = rep(NA,length(GeneSet)),
      "Exact_source" = rep(NA,length(GeneSet))
    )
    
    for(j in 1:length(GeneSet)){
     
      url <- paste0("https://www.gsea-msigdb.org/gsea/msigdb/cards/",GeneSet[j],".html")
      html_content <- read_html(url)
      
      GeneSet.term <- html_content %>% html_nodes(xpath = GeneSet.path) %>% html_text2()
      Exact_source.term <- html_content %>% html_nodes(xpath = Exact_source.path) %>% html_text2()
      
      # ## Ch
      # require(rvest)
      # require(V8)
      # # https://datascienceplus.com/scraping-javascript-rendered-web-content-using-r/
      # # Create a new v8 context
      # ct <- v8()
      # Exact_source.term <- read_html(ct$eval(gsub('document.write',
      #                                             '',Exact_source.term))) %>% 
      #                                             html_text2()
      
      
      if(length(GeneSet.term)==0){
        GeneSet.term <- NA
      }
      if(length(Exact_source.term)==0){
        Exact_source.term <- NA
      }
      output.df[j,1] <- GeneSet.term[1]
      output.df[j,2] <- Exact_source.term[1]
    }
    return(output.df) 
  }
  
  output.df <- parLapply(
    cl, GeneSet.lt, Crowler
  )
  output.df <- do.call(rbind, lapply(output.df, data.frame))
  
  stopCluster(cl)
  return(output.df)
}

##### Test #####
  test.lt <- read.delim2(
    "c5.go.v7.5.symbols/c5.go.v7.5.symbols.gmt",
    col.names = 1:max(count.fields("c5.go.v7.5.symbols/c5.go.v7.5.symbols.gmt")),
    header = F, sep = "\t"
  )[[1]] %>% as.list()
  
  test.df <- GSEA_crowler(test.lt[1:10], core = 4)
