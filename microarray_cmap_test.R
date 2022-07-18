library(GEOquery)
library(limma)

objects <- read.csv("GSE51068_objects\\rdata_objects.csv", header = F)
anno <- read.csv('GSE51068_probes.txt', sep = '\t')
rownames(anno) <-anno$ID
for (object in objects$V1)
{
  load(paste0('GSE51068_objects\\', object))
  dm <- dim(fit1$coefficients)
  tbl <- topTable(fit1, coef = 2, n=dm[1], p.value=0.05)
  if(dim(tbl)[1] == 0)
  {
    next
  }
  
  tbl$gene_id<- anno[rownames(tbl),]$Entrez.Gene
  tbl$gene_name<- anno[rownames(tbl),]$Gene.Symbol
  write.table(tbl, paste0("GSE51068_objects_results\\",object,".csv"), quote=F,
              row.names=F,
              col.names = T,
              sep = '\t')
}

################### GSE60408 ###################

objects <- read.csv("GSE60408_objects\\rdata_objects.csv", header = F)
anno <- read.csv('GSE60408_probes.txt', sep = '\t')
rownames(anno) <-anno$ID
for (object in objects$V1)
{
  load(paste0('GSE60408_objects\\', object))
  dm <- dim(fit$coefficients)
  tbl <- topTable(fit, coef = 3, n=dm[1], p.value=0.05)
  if(dim(tbl)[1] == 0)
  {
    next
  }
  
  tbl$gene_id<- anno[rownames(tbl),]$Entrez.Gene
  tbl$gene_name<- anno[rownames(tbl),]$Gene.Symbol
  write.table(tbl, paste0("GSE60408_objects_results\\",object,".csv"), quote=F,
              row.names=F,
              col.names = T,
              sep = '\t')
}

################### GSE102006 ###################

objects <- read.csv("GSE102006_objects\\rdata_objects.csv", header = F)
anno <- read.csv('GSE102006_probes.txt', sep = '\t')
rownames(anno) <-anno$ID
for (object in objects$V1)
{
  load(paste0('GSE102006_objects\\', object))
  dm <- dim(fit2$coefficients)
  tbl <- topTable(fit2, coef = 2, n=dm[1], p.value=0.05)
  if(dim(tbl)[1] == 0)
  {
    next
  }
  
  tbl$gene_id<- as.character(anno[rownames(tbl),]$GENE)
  tbl$gene_name <- as.character(anno[rownames(tbl),]$GENE_SYMBOL)
  write.table(tbl, paste0("GSE102006_objects_results\\",object,".csv"), quote=F,
              row.names=F,
              col.names = T,
              sep = '\t')
}


###################
library(knitr)
library(tinytex)
library(httr)
library(jsonlite)
library(htmltools)

## Signature as a list of up and down genes
```{r up and down genes}
```{r connected signatures by uploading signature file}

setwd("GSE60408_ilincs")
files <- read.csv("files", header = F)$V1

for(sigFile in files)
{
  apiUrl<-"http://www.ilincs.org/api/SignatureMeta/upload"
  req <- POST(apiUrl, body=list(file=httr::upload_file(sigFile)))
  signatureFile <- httr::content(req)$status$fileName[[1]]
  apiUrl <- paste("http://www.ilincs.org/api/ilincsR/signatureEnrichment?sigFile=",signatureFile,"&library=LIB_5&metadata=TRUE",sep="")
  req <- GET(apiUrl)
  json <- httr::content(req, as = "text")
  iLincsConnectedCompoundPerturbations <- fromJSON(json)$enrichment
  write.table(iLincsConnectedCompoundPerturbations, paste0("results_signatureEnrichment\\",sigFile, "_ilincs_results.csv"),
              quote=F,
              row.names=F,
              col.names = T,
              sep = '\t')
}


files <- read.csv("files", header = F)$V1

for(sigFile in tail(files, length(files)-51))
{
  apiUrl="http://www.ilincs.org/api/ilincsR/findConcordancesSC"
  
  df<-read.csv(sigFile)
  up_genes <- df[df$Value_LogDiffExp == 1,]$Name
  down_genes <- df[df$Value_LogDiffExp == -1,]$Name
  if(length(up_genes)<2)
  {
    next
  }
  
  if(length(down_genes)<2)
  {
    next
  }
  topUpRegulatedGenes <- list(genesUp=up_genes)
  topDownregulatedGenes <- list(genesDown=down_genes)
  req <- POST("http://www.ilincs.org/api/ilincsR/findConcordancesSC", body = list(mode="UpDn",metadata=TRUE,signatureProfile = c(topUpRegulatedGenes, topDownregulatedGenes)),encode = "json")
  ilincsUpDnConnectedSignatures <- data.table::rbindlist(httr::content(req)$concordanceTable, use.names = TRUE, fill = TRUE)

  write.table(ilincsUpDnConnectedSignatures, paste0("results_findConcordancesSC\\",sigFile, "_ilincs_results.csv"),
              quote=F,
              row.names=F,
              col.names = T,
              sep = '\t')
}
