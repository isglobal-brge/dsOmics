lmFeature <- function(feature, vars, data, cellEstim){
  cally <- paste("selFeature(", data, ", feature=", deparse(feature), ", vars=", deparse(vars), ")")
  datashield.assign(conns, 'dat', as.symbol(cally))
  
  if (!ds.isNA(cellEstim)[[1]]){
      ds.cbind(c("dat", "cell.counts"), newobj="dat")
  }
    
  mm <- as.formula(paste(feature, "~ ."))
  mod <- ds.glm(mm, family='gaussian', data='dat', viewIter = FALSE)
  metrics <- as.data.frame(mod$coefficients[2, c(1,2,4)])
  names(metrics) <- feature
  return(metrics)
}

glmSNP <- function(snp, vars, data, family="binomial"){
  y <- vars[1]
  dep <- paste(c(snp, vars[-1]), collapse="+")
  mm <- as.formula(paste(y, "~", dep))
  mod <- ds.glm(mm, family=family, data=data)
  metrics <- as.data.frame(mod$coefficients[2, c(1,2,4)])
  names(metrics) <- snp
  return(metrics)
}
