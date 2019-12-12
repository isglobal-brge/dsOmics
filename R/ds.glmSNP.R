##' Performing a logistic regression analysis on pooled data from multiple studies for every SNP
##' 
##' Function fits a generalized linear model to gentic data for each SNP
##' in the data sets considered, using user specified outcome and phenotypic variabes as covariates
##' Outputs a matrix containing a beta value, standard error and p-value for each SNP
##' @title Logistic regression analysis of pooled data for each SNP site in study
##' @param snps an optional parameter input as a vector of integer values which indicates the names of
##' SNPs (rs numbers) sites that should be analysed. If missing all SNPs are analysed
##' @param model list of phenotypic variables to use as covariates in the regression analysis in the form: 
##' "outcome ~ covar1 + covar2 +  ... + covarN"
##' @param data name of the DataSHIELD object to which the genotype (snpMatrix) and phenotypic data (data.frame) has been assigned
##' @param datasources Opal object or list of opal objects denoting the opal server(s) information
##' @param mc.cores optional parameter that allows the user to specify the number of CPU cores to use during 
##' parallel processing. Argument can only be > 1 when the function is run on a linux machine
##' @export
##' @examples
##' 

ds.glmSNP <- function(snps=NULL, model, eSet, 
                      type.p.adj='fdr', mc.cores = 1){
  
  
  mt <- as.formula(model)
  vars <- all.vars(mt)
  
  # Setting the SNPs to loop over as the total number of 
  # SNPs in the studies if no SNPs are specified
  if(is.null(features)){
    ff <- datashield.aggregate(conns, "featureNamesDS(ES)")
    features <- Reduce(intersect, ff)
  }
  
  
  
  #Setting the number of SNPs to loop over as the total number of 
  #features in the studies if no rs are specified
  
  
  if(is.null(snps)){
    ff <- nFeatures[[1]]["Features"]
    features <- ds.colnames(eSet.data)[[1]][1:ff]
  }

  cally <- paste0("ds.snpMatrix(ES, snp='", snps, "')")
  datashield.assign(conns, 'dd', as.symbol(cally))
  
  
  ans <- t(as.data.frame(parallel::mclapply(snps, glmSNP, vars=vars,
                                            data='dd', mc.cores = mc.cores))) 
  
  #Ordering matrix by ascending p.value
  ans <- ans[order(ans[, 'p-value']), ]
  
  #Calculating adjusted p-values
  ans <- cbind(ans, p.adj = p.adjust(ans[,'p.value'],  method=type.p.adj))
  
  class(ans) <- c("dsGlmSNP", class(ans))
  
  return(ans)
  
}
