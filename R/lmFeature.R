lmFeatureDS <- function(feature, vars, data, cellEstim, connections){
  cally <- paste("selFeatureDS(", data, ", feature=", deparse(feature), ", vars=", deparse(vars), ")")
  datashield.assign(connections, 'dat', as.symbol(cally))
  
  if (!ds.isNA(cellEstim)[[1]]){
      ds.cbind(c("dat", "cell.counts"), newobj="dat")
  }
    
  mm <- as.formula(paste(feature, "~ ", paste(ds.colnames('dat')[[1]][-1], collapse="+")))
  mod <- ds.glm(mm, family='gaussian', data='dat', viewIter = FALSE)
  metrics <- as.data.frame(mod$coefficients[2, c(1,2,4)])
  names(metrics) <- feature
  return(metrics)
}

glmSNPDS <- function(snp, snps, gds, covars, vars, connections, family="binomial"){
  i <- which(snp==snps)
  cally <- paste0("selSNP(", gds, ", i=", i, ", covars=", covars,
                  ", vars=", deparse(vars), ")")
  datashield.assign(connections, 'dat', as.symbol(cally))
  y <- vars[1]
  dep <- paste(c("snp", vars[-1]), collapse="+")
  mm <- as.formula(paste(y, "~", dep))
  mod <- try(ds.glm(mm, family=family, data='dat', 
                datasources=connections, viewIter = FALSE), TRUE)
  if (inherits(mod, "try-error")) {
    metrics <- data.frame(a=rep(NA,3))
    row.names(metrics) <- c("Estimate", "Std. Error", "p-value")  
  }
  else if (mod$errorMessage!="No errors"){
    metrics <- data.frame(a=rep(NA,3))
    row.names(metrics) <- c("Estimate", "Std. Error", "p-value")  
  }
  else {
    metrics <- as.data.frame(mod$coefficients[2, c(1,2,4)])
  }
  colnames(metrics) <- snp
  return(metrics)
}


limmaDS <- function(model, eSet, sva){
  design <- stats::model.matrix(stats::as.formula(model), 
                                data=Biobase::pData(eSet))
  ee <- Biobase::exprs(eSet)
  fit <- limma::lmFit(eSet, design)
  fit <- limma::eBayes(fit)
  ans <- limma::topTable(fit, n=Inf)
  return(ans)
}

