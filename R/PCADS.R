#' Title
#'
#' @param genoData 
#' @param rs 
#' @param means 
#' @param sd_hw 
#'
#' @return
#' @export
#'
#' @examples
PCADS <- function(genoData, pca_rs, pca_means, pca_sd_hw){
  
  # TODO passar el snpBlock com un argument!!!!!!
  # TODO what do we do with the missing values of the genotype??? if set to NA, sva(x) fails!!!
  # TODO format error messages of the disclosure check!
  
  if(!is.null(pca_rs)){
    dt <- tibble(rs = strsplit(pca_rs, split=",")[[1]], 
                 means = as.numeric(strsplit(pca_means, split=",")[[1]]), 
                 sd_hw = as.numeric(strsplit(pca_sd_hw, split=",")[[1]]))
    stand <- T
  } else {stand <- F}
  
  results <- Reduce(cbind, lapply(genoData, function(i){
    if(inherits(i, "GenotypeData")){
      geno <- GWASTools::GenotypeBlockIterator(i, snpBlock=1000L)
    } else if (inherits(i, "GdsGenotypeReader")) {
      gds_gd <- GWASTools::GenotypeData(i)
      geno <- GWASTools::GenotypeBlockIterator(gds_gd, snpBlock=1000L)
    } else {
      stop()
    }
    gfile <- openfn.gds(geno@data@filename, allow.duplicate = T, allow.fork = T)
    n <- index.gdsn(gfile, "genotype")
    iterations <- length(geno@snpFilter)
    rs <- GWASTools::getVariable(geno, "snp.rs.id")
    
    if(stand){
      svd_partial <- Reduce(cbind, lapply(1:iterations, function(x){
        snp.sel <- geno@snpFilter[[x]]
        genoData_temp <- readex.gdsn(n, list(NULL, snp.sel))#, .value = 3, .substitute = NA)
        if(!disclosureCheck(genoData_temp)){
          dt_temp <- dt[which(dt$rs %in% rs[snp.sel]),]
          genoData_temp <- t((t(genoData_temp) - dt_temp$means) / dt_temp$sd_hw)
          svd_partial <- svdPartial(genoData_temp)
        }else{
          stop()
        }
      }))
    } else {
      svd_partial <- Reduce(cbind, lapply(1:iterations, function(x){
        snp.sel <- geno@snpFilter[[x]]
        genoData_temp <- readex.gdsn(n, list(NULL, snp.sel))#, .value = 3, .substitute = NA)
        if(!disclosureCheck(genoData_temp)){
          svd_partial <- svdPartial(genoData_temp)
        }else{
          stop()
        }
      }))
    }
    
    return(svd_partial)
    
  }))
  
  results <- svdPartial(results)
  
  return(results)
  
}

#' Title
#'
#' @param x 
#'
#' @return
#'
#' @examples
svdPartial <- function(x){
  
  ss <- svd(x)
  ans <- sweep(ss$u, 2, FUN="*", ss$d)
  return(ans)
  
}

#' Title
#'
#' @param dataframe 
#'
#' @return
#' @export
#'
#' @examples
disclosureCheck <- function(dataframe){
  errorMessage <- FALSE
  
  # The PCA is potentially disclosive when the calculation of the covariate matrix can be disclosive.
  # For that reason, the disclosive check of the funcion dsBase::covDS is copied, if it passes the check
  # a PCA is calculated using an off-the-shelve PCA function
  
  #############################################################
  #MODULE 1: CAPTURE THE nfilter SETTINGS
  thr <- dsBase::listDisclosureSettingsDS()
  nfilter.tab <- as.numeric(thr$nfilter.tab)
  nfilter.glm <- as.numeric(thr$nfilter.glm)
  #nfilter.subset <- as.numeric(thr$nfilter.subset)
  #nfilter.string <- as.numeric(thr$nfilter.string)
  #############################################################
  
  # names of the variables
  cls <- colnames(dataframe)
  
  # number of the input variables
  N.vars <- ncol(dataframe)
  
  ######################
  # DISCLOSURE CONTROLS
  ######################
  
  ##############################################################
  # FIRST TYPE OF DISCLOSURE TRAP - TEST FOR OVERSATURATION
  # TEST AGAINST nfilter.glm								
  ##############################################################
  
  varcov.saturation.invalid <- 0
  if(N.vars > (nfilter.glm * nrow(dataframe))){
    
    varcov.saturation.invalid <- 1
    
    errorMessage <- TRUE
  }
  
  # CHECK X MATRIX VALIDITY
  # Check no dichotomous X vectors with between 1 and nfilter.tab value
  # observations at either level
  
  X.mat <- as.matrix(dataframe)
  
  Xpar.invalid <- rep(0, N.vars)
  
  for(pj in 1:N.vars){
    unique.values.noNA <- unique((X.mat[,pj])[stats::complete.cases(X.mat[,pj])])
    if(length(unique.values.noNA)==2){
      tabvar <- table(X.mat[,pj])[table(X.mat[,pj])>=1] #tabvar COUNTS N IN ALL CATEGORIES WITH AT LEAST ONE OBSERVATION
      min.category <- min(tabvar)
      if(min.category < nfilter.tab){
        Xpar.invalid[pj] <- 1
      }
    }
  }
  
  # if any of the vectors in X matrix is invalid then the function returns all the
  # outputs by replacing their values with NAs
  
  if(is.element('1', Xpar.invalid)==TRUE & varcov.saturation.invalid==0){
    
    errorMessage <- TRUE
    
  }
  
  return(errorMessage)
  
}
