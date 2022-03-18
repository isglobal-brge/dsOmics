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
PCADS <- function(genoData, pca_rs, pca_means, pca_sd_hw, snpBlock){
  
  # TODO passar el snpBlock com un argument!!!!!!
  # TODO what do we do with the missing values of the genotype??? if set to NA, sva(x) fails!!!
  # TODO format error messages of the disclosure check!
  
  if(!is.null(pca_rs)){
    dt <- tibble(rs = strsplit(pca_rs, split=",")[[1]], 
                 means = as.numeric(strsplit(pca_means, split=",")[[1]]), 
                 sd_hw = as.numeric(strsplit(pca_sd_hw, split=",")[[1]]))
    stand <- T
  } else {stand <- F}
  
  # results <- Reduce(cbind, lapply(genoData, function(i){
  #   if(inherits(i, "GenotypeData")){
  #     geno <- GWASTools::GenotypeBlockIterator(i, snpBlock=snpBlock)
  #   } else if (inherits(i, "GdsGenotypeReader")) {
  #     gds_gd <- GWASTools::GenotypeData(i)
  #     geno <- GWASTools::GenotypeBlockIterator(gds_gd, snpBlock=snpBlock)
  #   } else {
  #     stop()
  #   }
  #   gfile <- openfn.gds(geno@data@filename, allow.duplicate = T, allow.fork = T)
  #   n <- index.gdsn(gfile, "genotype")
  #   iterations <- length(geno@snpFilter)
  #   rs <- GWASTools::getVariable(geno, "snp.rs.id")
  #   
  #   if(stand){
  #     # svd_partial <- Reduce(cbind, lapply(1:iterations, function(x){
  #     #   snp.sel <- geno@snpFilter[[x]]
  #     #   genoData_temp <- readex.gdsn(n, list(NULL, snp.sel))#, .value = 3, .substitute = NA)
  #     #   dt_temp <- dt[which(dt$rs %in% rs[snp.sel]),]
  #     #   genoData_temp <- t((t(genoData_temp) - dt_temp$means) / dt_temp$sd_hw)
  #     #   svd_partial <- svdPartial(genoData_temp)
  #     # }))
  #   } else {
  #     browser()
  #     svd_partial <- Reduce(cbind, lapply(1:iterations, function(x){
  #       snp.sel <- geno@snpFilter[[x]]
  #       genoData_temp <- readex.gdsn(n, list(NULL, snp.sel))#, .value = 3, .substitute = NA)
  #       svd_partial <- svdPartial(genoData_temp)
  #     }))
  #   }
  #   return(svd_partial)
  # }))
  # browser()
  results <- do.call(cbind, lapply(genoData, function(x){
    GWASTools::getVariable(x, "genotype")
  }))

  results <- svdPartial(t(results))
  
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

#' @title Add PCA results to geno
#' 
#' @description Adds the results of the block method to the geno, creating a new dt
#'
#' @param object \code{geno} to add the PCA results
#' @param pca \code{raw} Serialized PCA object
#'
#' @return \code{dt}
#' @export

gneo_pca_pooled_addPCDS <- function(geno, pca, ncomp){

  # Extract genotype data
  genoSNPS <- do.call(cbind, lapply(geno, function(x){
    GWASTools::getVariable(x, "genotype")
  }))
  
  # Unserialize PCA object
  pca <- matrix(pca, ncol = ncomp)
  browser()
  # Calculate individual princ comp values
  pca_coord <- genoSNPS %*% pca
  colnames(pca_coord) <- paste0("Dim.", 1:ncol(pca))
  ### ADDD ROWNAMES TOO!!!!!
  return(data.frame(pca_coord))
  
}
