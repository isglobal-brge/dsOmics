lmFeatureDS <- function(cpgIndex, model, molecular.data, pheno.data){
  
  # ASSIGN FUNCTION
  
  #Extracting name of feature (CpG site)
  sel <- molecular.data$feature[cpgIndex]
  
  mt <- as.formula(model)
  vars <- all.vars(mt)
  
  #Extracting expression values from molecular.data for feature given by "sel"
  exprValues = as.numeric(t(molecular.data[which(molecular.data$feature==sel), 
                                           -c(which(names(molecular.data)=='feature'))]))
  
  #Adding sample names to corresponding expression values
  m.data <- as.data.frame(cbind(colnames(molecular.data)[-1], exprValues))
  colnames(m.data) <- c("sample", sel)
  
  #Converting expression values from factor back to numeric
  m.data[, sel] = as.numeric(as.character(m.data[, sel]))
  
  #Extracting values from pheno.data for specified phenotypic covariates
  p.data <- pheno.data[, c("sample", vars)]
  colnames(p.data) <- c("sample", vars)
  
  #Joining methylation and pheno data by "sample" column
  data.sel = merge(m.data, p.data, by="sample")
  
  rownames(data.sel) <- c()

  return(data.sel)
  
}
