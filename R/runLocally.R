list_of_gene_symbol <- genesym.illumina_probeid$alias_symbol
gene_to_affy_result <- geneToAffy(list_of_gene_symbol = list_of_gene_symbol)

sigset <- read.csv("50081sig",header = F,stringsAsFactors = F)
sigset2 <- read.csv("sig compare",header = T,stringsAsFactors = F)

names(sigset) <- "sig"

l <- sigset$sig
r1 <- AffyToGene(l)
r2 <- GeneToillum(r1$symbol)
sig2 <- data.frame(sig= r1$probe_id,geneid = r1$symbol,stringsAsFactors = F)

sigtocompare <- data.frame(sig= sigset2$probe_id,geneid = sigset2$symbol,stringsAsFactors = F)
gseid<- "GSE14814" 
cleaning_myDB(data_db)
result <- getCoxEfficient(gseid = gseid)
coxdata<- result[[2]]
res <- result[[1]]
#change signature
sur_data_result <- getSurvData(coxdata ,res = res,sig_set = sig2)
sur_data <- sur_data_result[[1]]
sig_set <-sur_data_result[[2]]
surv_plot <-plotSurv(sur_data)

coxph_plot <- plotCoxph(sig_set,sur_data = sur_data)
#======
l <- sigset2$symbol
r2 <- GeneToillum(l)
sig2 <- data.frame(sig= r2$probe_id,geneid = r2$alias_symbol,stringsAsFactors = F)
sigtocompare <- data.frame(sig= sigset2$probe_id,geneid = sigset2$symbol,stringsAsFactors = F)
gseid <- "GSE50081" 
cleaning_myDB(data_db)
result <- getCoxEfficient(gseid = gseid)
coxdata<- result[[2]]
res <- result[[1]]
#change signature
sur_data_result <- getSurvData(coxdata ,res = res,sig_set = sigtocompare)
sur_data <- sur_data_result[[1]]
sig_set <-sur_data_result[[2]]
surv_plot <-plotSurv(sur_data)

coxph_plot <- plotCoxph(sig_set,sur_data = sur_data)

