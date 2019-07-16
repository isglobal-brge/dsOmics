lmFeatureDS <- function(cpgIndex, model, molecular.data, pheno.data){
  
  # ASSIGN FUNCTION
  
  sel <- molecular.data$feature[cpgIndex]
  
  mt <- as.formula(model)
  vars <- all.vars(mt)
  
  data.sel <- data.frame(as.numeric(t(molecular.data[which(molecular.data$feature==sel),-c(which(names(molecular.data)=='feature'))])), pheno.data[, vars])
  colnames(data.sel) <- c(sel, vars)
  rownames(data.sel) <- c()

  return(data.sel)
  
}
