plotSurv<- function(sur_data,type="survival"){
  surv_time_name <- names(sur_data)[grep("time",names(sur_data))][1]
  
  surv_object <- Surv(time = sur_data[,surv_time_name],
                      event = sur_data[,"record_status.ch1"])
  
  fit1 <- survfit(surv_object~risk_status.ch1, data = sur_data)
  p <- ggsurvplot(fit1, 
                  data = sur_data, 
                  pval = TRUE,
                  risk.table = T)
  return(p)
  
  
}

plotCoxph<- function(sig_set,sur_data){
  gene_s <- paste0(sig_set$geneid,sep = "+",collapse = "")
  gene_s <- gsub('.{1}$', '', gene_s)
  
  surv_object <- Surv(time = sur_data[,surv_time_name],
                      event = sur_data[,"record_status.ch1"])
  
  model <- coxph( formula = as.formula(paste('surv_object ~',gene_s)) ,data = sur_data)
  result <- exp(cbind(coef(model),confint(model)))
  
  p <-ggforest(model,data = sur_data)
  exp(cbind(coef(model),confint(model)))
  #univ_formulas <- sapply(gene_s,function(x)as.formula(paste('surv_object ~',x)))
  
  #univ_models <- lapply(univ_formulas, function(x){coxph(x,data=sur_data3)})
  #univ_results <- lapply(univ_models,function(x){return(exp(cbind(coef(x),confint(x))))})
  return(p)
  
}
