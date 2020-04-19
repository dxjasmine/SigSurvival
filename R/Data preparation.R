#package required
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")

library(GEOquery)

# datasets preparation
constructDB <- function(){
  
  mydb <- read.csv("my datasets",header = T,sep = "\t",stringsAsFactors = F)
  #remove sample size <30
  keep <- as.numeric(dataset_matrix[,3])>30
  dataset_matrix <- dataset_matrix[keep,]
  
  data_db <- list()
  for (i in seq_along(mydb$GEO.id)){
    
    print(i)
    gse = getGEO(mydb$GEO.id[i], GSEMatrix =TRUE)
    x <- exprs(gse[[1]])
    
    # remove Affymetrix control probes
    #x <- x[-grep('^AFFX', rownames(x)),]
    
    # transform the expression data to Z scores
    x <- t(scale(t(x)))
    
    # extract information of interest from the phenotype data (pdata)
    idx <- grep("*:ch1",colnames(pData(gse[[1]])))
    
    
    metadata <- data.frame(pData(gse[[1]])[,idx],
                           row.names = rownames(pData(gse[[1]])))
    
    
    # filter the Z-scores expression data to match the samples in our pdata
    x <- x[,which(colnames(x) %in% rownames(metadata))]
    
    # check that sample names match exactly between pdata and Z-scores 
    all((colnames(x) == rownames(metadata)) == TRUE)
    ## [1] TRUE
    
    # create a merged pdata and Z-scores object
    coxdata <- data.frame(metadata, t(x))
    data_db[[i]] <- coxdata
  }
  names(data_db) <- mydb$GEO.id
}

#change some column name of datasets to make them easier to process in later query
cleaning_myDB <- function(data_db){
  clean_data_db <- list()
  gseid <- "GSE41271"
  coxdata <- data_db[[gseid]]
  #calculate survival time
  sel <- which(is.na(coxdata$last.follow.up.survival.ch1))
  coxdata <- coxdata[-sel,]
  coxdata$last.follow.up.survival.ch1
  coxdata$survival.time.in.year.ch1 <- (as.Date(as.character(coxdata$last.follow.up.survival.ch1), format="%Y-%m-%d") - 
                                          as.Date(as.character(coxdata$date.of.surgery.ch1), format="%Y-%m-%d"))/365
  #censor event
  coxdata$status.ch1 <- rep(1,length(coxdata$last.follow.up.survival.ch1))
  clean_data_db[[gseid]] <- coxdata
  #===================================
  gseid <- "GSE50081"
  coxdata <- data_db[[gseid]]
  clean_data_db[[gseid]] <- coxdata
  #=============
  gseid<- "GSE14814"
  coxdata <- data_db[[gseid]]
  coxdata <- coxdata[which(coxdata$Post.Surgical.Treatment.ch1=="OBS"),]
  discard <- c("dss.survival.time.ch1", "DSS.time.ch1","DSS.status.ch1" )
  coxdata <- coxdata[,!names(coxdata) %in% discard]
  clean_data_db[[gseid]] <- coxdata
  
}



#=================

library(hgu133a.db)
library(annotate)
library(illuminaHumanv3.db)
list_of_gene_symbol<- c("ANXA1", "APCS", "C6", "COL14A1", "CYB5R3", "FAM82A1", "GBE1", "GNAI1", "HLA-DQA1", "IL16", "LTBP2", "MMS19", "MXRA5", "SEC13")

geneToAffy<- function(list_of_gene_symbol){
  x <- hgu133aSYMBOL
  # Get the probe identifiers - gene symbol mappings
  mapped_probes <- mappedkeys(x)
  # Convert to a dataframe
  genesym.affy_probeid <- as.data.frame(x[mapped_probes])
  #head(genesym.probeid)
  results <- genesym.affy_probeid[which(genesym.affy_probeid$symbol %in% list_of_gene_symbol),]
  
  return(results)
}

AffyToGene<- function(list_of_affy){
  x <- hgu133aSYMBOL
  # Get the probe identifiers - gene symbol mappings
  mapped_probes <- mappedkeys(x)
  # Convert to a dataframe
  genesym.affy_probeid <- as.data.frame(x[mapped_probes])
  #head(genesym.probeid)
  results <- genesym.affy_probeid[which(genesym.affy_probeid$probe_id %in% list_of_affy),]
  
  return(results)
}



#Convert from illumina id to  gene symbol

illumToGene<- function(list_of_illum_id){
  
  
  x <- illuminaHumanv3ALIAS2PROBE
  # Get the probe identifiers that are mapped to an ACCNUM
  mapped_probes <- mappedkeys(x)
  # Convert to a list
  genesym.illumina_probeid <- as.data.frame(x[mapped_probes])
  results <- genesym.illumina_probeid[which(genesym.illumina_probeid$symbol %in% list_of_illum_id),]
  return(results)
  
}

#Convert from gene symbol to illumina id
GeneToillum<- function(list_of_gene){
  
  
  
  x <- illuminaHumanv3ALIAS2PROBE
  # Get the probe identifiers that are mapped to an ACCNUM
  mapped_probes <- mappedkeys(x)
  # Convert to a list
  genesym.illumina_probeid <- as.data.frame(x[mapped_probes])
  results <- genesym.illumina_probeid[which(genesym.illumina_probeid$alias_symbol %in% list_of_gene),]
  return(results)
  
}