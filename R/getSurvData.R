# tidy column type - survival status
recordstatus <- function(x){
  if(grepl("^[dDy]",x )) {rs = 1}
  if(grepl("^[aAn]",x )){rs = 0}
  return (rs)
}

NYstatus <- function(x){
  if(grepl("N",x )) {rs = 1}
  if(grepl("Y",x )){rs = 0}
  return (rs)
}

getCoxEfficient <- function(gseid){
  coxdata <- clean_data_db[[gseid]]
  
  #remove unknown status 
  status_name <- names(coxdata)[grep("status",names(coxdata))][1]
  
  
  if (gseid == "GSE37745" | gseid == "gse37745") {
    status_name <- "dead.ch1"
  }
  if(length(unique(coxdata[,status_name]))>2){
    discard <- which(coxdata[,status_name]=="Not available" | coxdata[,status_name]=="NTL" )
    coxdata<- coxdata[-discard,]
  }
  
  
  #record status for censored event
  if(typeof(coxdata[,status_name]) =="double"){
    print("here")
    coxdata$record_status.ch1 <- coxdata[,status_name]
  }else{
    for (i in 1:length(rownames(coxdata))) {
      coxdata$record_status.ch1[i] <- recordstatus(coxdata[,status_name][i]) 
    }
  }
  print("123")
  
  
  
  
  # numberic type conversion
  surv_time_name <- names(coxdata)[grep("time",names(coxdata))][1]
  coxdata[,surv_time_name] <- as.numeric(coxdata[,surv_time_name])
  
  
  
  #2, test each gene independently via Cox regression
  f <- paste('Surv(',surv_time_name,',','record_status.ch1',')' ,'~ [*]')
  print(f)
  
  fname <- paste(gseid,surv_time_name,sep = " ",collapse = "")
  fcheckname <- paste("./",fname,sep="",collapse = "")
  if(!file.exists(fcheckname)){
    res <- RegParallel(
      data = coxdata,
      formula = f,
      FUN = function(formula, data)
        coxph(formula = formula,
              data = data,
              ties = 'breslow',
              singular.ok = TRUE),
      FUNtype = 'coxph',
      variables = colnames(coxdata)[(length(idx)+1):ncol(coxdata)],
      blocksize = 2000,
      cores = 2,
      nestedParallel = FALSE,
      conflevel = 95)
    
    write.csv(res,fname)
  }else{
    print("456")
    res <- read.csv(fname,header = T,stringsAsFactors = F)
  }
  
  r<- list()
  r[[1]]<-res
  r[[2]]<-coxdata
  return(r)
}


getSurvData <- function(coxdata,res,sig_set,type = "survival"){
  #get signature risk score
  sig2 <- sig_set
  #sig2 <- data.frame(sig= sigcompare$affyid,gene= sigcompare$symbol)
  for (i in 1:length(sig2$sig)) {
    print(i)
    b_i <- grep(sig2$sig[i],res$Term)
    sig2$coeff[i] <- res[b_i,]$Beta
    sig2$logrank[i]<- res[b_i,]$LogRank
    sig2$pvalue[i]<- res[b_i,]$P
    sig2$hr[i]<- res[b_i,]$HR
    
  }
  sig3 <- data.frame(geneid = unique(sig2$geneid),stringsAsFactors = F)
  
  for (i in seq_along(sig3$geneid)) {
    cur <- sig2[which(sig2$geneid == sig3$geneid[i]),]
    sig3$sig[i] <- cur$sig[1]
    sig3$coeff[i] <- max(cur$coeff)
    sig3$logrank[i]<- max(cur$logrank)
    sig3$pvalue[i]<- max(cur$pvalue)
    sig3$hr[i]<- max(cur$hr)
  }
  sig2 <- sig3
  
  
  
  #prepare survial data
  col_sel <- grep("ch1",colnames(coxdata))
  
  sur_data2 <- coxdata[,col_sel]
  for (i in 1:length(sig2$sig)) {
    sig_index <- grep(sig2$sig[i],names(coxdata))
    probe_name <- as.character(sig2$geneid[i])
    sur_data2[,probe_name] <- coxdata[,sig_index]
    
  }
  
  
  
  
  #calculate risk score
  sur_data3 <- sur_data2
  rs_sum <- numeric(length = length(rownames(sur_data3)))
  for (i in 1:length(sig2$coeff)) {
    cur_gene <- sur_data3[,sig2$geneid[i]]
    rs_sum <- rs_sum +  cur_gene*sig2$coeff[i]
  }
  
  for (i in 1: length(rownames(sur_data3))) {
    risk_score <- rs_sum[i]
    sur_data3$rs.ch1[i] <- risk_score
  }
  
  
  highExpr <- median(sur_data3$rs)
  
  for (i in 1: length(rownames(sur_data3))) {
    sur_data3$risk_status.ch1[i] <- ifelse(sur_data3$rs[i] >= highExpr, 'High','low')
  }
  r<- list()
  r[[1]]<-sur_data3
  r[[2]]<-sig2
  return(r)
}