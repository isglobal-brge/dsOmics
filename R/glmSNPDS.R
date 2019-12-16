#'
#' @title Fits a GLM to assess association between the outcome (binary variable) and a given SNP
#' @description To be supplied
#' @param snp
#' @param snps
#' @param gds
#' @param covars
#' @param vars
#' @param connections
#' @param family
#' @return a vector with effect estimates, standard error and associated p-value
#' @author Gonzalez, JR.
#'
#' @export 

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