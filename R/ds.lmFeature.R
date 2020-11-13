ds.lmFeature <- function(i, model, molecular.data, pheno.data, datasources){

  mt <- as.formula(model)
  vars <- all.vars(mt)
  
  cally <- paste0("lmFeatureDS(i=", i, ", ", "model=", model, ", " ,"molecular.data=", molecular.data, ", ", "pheno.data=", pheno.data,")")
  datashield.assign(datasources, 'data.sel', as.symbol(cally))
  
  dep <- paste(vars, collapse="+")
  colnames <- ds.colnames('data.sel')
  feature <- colnames$study1[1]
  mm <- as.formula(paste(feature, "~", dep))
  
  mod <- ds.glm.o(mm, family='gaussian', data='data.sel')
  
  ans <- matrix(mod$coefficients[1, c(1,2,4)], nrow=1)
  colnames(ans) <- c("beta", "s.e.", "p.value")
  rownames(ans) <- feature
  return(ans)
  
}