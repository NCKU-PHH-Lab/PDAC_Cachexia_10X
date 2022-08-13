# https://www.itread01.com/content/1549283773.html

##### Presetting ######
  rm(list = ls()) # Clean variable
  memory.limit(150000)

##### Load Packages  ##### 
  library(XML)
  library(RCurl)

##### #####  
  myheader <- c(
    "User-Agent"="Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/48.0.2564.97 Safari/537.36",
    "Accept"="text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8",
    "Accept-Language"="gzip, deflate, sdch",
    "Connection"="keep-alive",
    "Accept-Charset"="GB2312,utf-8;q=0.7,*;q=0.7"
  )
  
  # url='https://www.gsea-msigdb.org/gsea/msigdb/cards/GOBP_ACID_SECRETION.html'
  url='https://portal.gdc.cancer.gov/legacy-archive/files/76bf1f8c-bc30-4530-9767-508e34c04fa8'

  
  d <- debugGatherer()
  web <- getURL(url, httpheader = myheader, debugfunction = d$update, verbose = T,ssl.verifyhost=FALSE,ssl.verifypeer=FALSE)
doc <- htmlParse(web,encoding = "UTF-8")
project_title <- sapply(getNodeSet(doc,"//h2[@class='title']//a"),xmlValue)

web2 <- getURL(url, debugfunction = d$update, verbose = T,ssl.verifyhost=FALSE,ssl.verifypeer=FALSE)
doc2 <- htmlParse(web2,encoding = "UTF-8")




doc.df <- readHTMLTable(web, header = TRUE)
doc.df <- htmlTreeParse(web)

## KLC
library(parallel)
library(magrittr)
library(utils)
library(rvest)
library(xml2)

GeneSet.path	 <- '/html/body/div[3]/div[4]/div/table/tbody/tr/td[3]/table/tbody/tr[1]/td'
Exact_source.path	 <- '/html/body/div[3]/div[4]/div/table/tbody/tr/td[3]/table/tbody/tr[7]/td'


html_content <- read_html(url)
GeneSet.term <- html_content %>% html_nodes(xpath = GeneSet.path) %>% html_text2()
Exact_source.term <- html_content %>% html_nodes(xpath = Exact_source.path) %>% html_text2()


html_content2 <- read_html(url, httpheader = myheader)
html_content3 <- read_html(web)

GeneSet.term <- html_content3 %>% html_nodes(xpath = GeneSet.path) %>% html_text2()
Exact_source.term <- html_content3 %>% html_nodes(xpath = Exact_source.path) %>% html_text2()
