ds.lmFeature <- function(i, model, molecular.data, pheno.data, datasources){

  mt <- as.formula(model)
  vars <- all.vars(mt)
  
  ans <- matrix(0, ncol=3, nrow=length(i))
  colnames(ans) <- c("beta", "s.e.", "p.value")
  rowNames <- c()
  
  for (k in 1:length(i)){
  
  cally <- paste0("lmFeatureDS(i=", i[k], ", ", "model=", model, ", " ,"molecular.data=", molecular.data, ", ", "pheno.data=", pheno.data,")")
  datashield.assign(datasources, 'data.sel', as.symbol(cally))
  
  dep <- paste(vars, collapse="+")
  colnames <- ds.colnames('data.sel')
  feature <- colnames$study1[1]
  mm <- as.formula(paste(feature, "~", dep))
  
  mod <- ds.glm.o(mm, family='gaussian', data='data.sel')
  
  ans[k,] <- mod$coefficients[2, c(1,2,4)]
  rowNames[k] <- feature
  
  }
  
  rownames(ans) <- rowNames
  
  return(ans)
  
}