#' @title Subset genotype
#' 
#' @description Subset genotype by SNPs linked with ethcity differences from 
#'
#' @param genoData \code{GenotypeData} or \code{GdsGenotypeReader} object
#' @param snp_list \code{character} Use \code{"ethnic_snps"}
#'
#' @return \code{GenotypeData} or \code{GdsGenotypeReader} object
#' @export

subsetGenoDS <- function(genoData, snp_list){
  
  genoData_og <- get(genoData, envir = parent.frame())

  if(snp_list == "ethnic_snps"){
    data("ethnic_snps", package = "dsOmics")
    snp_subset <- ethnic_snps
  } else {
    stop("Not implemented a differnt argument for `snp_list` other than: 'ethnic_snps'")
  }

  original_rs <- GWASTools::getVariable(genoData_og, "snp.rs.id")
  original_ids <- GWASTools::getVariable(genoData_og, "snp.id")
  ids_of_interest <- original_ids[which(original_rs %in% snp_subset)]
  
  #############################################################
  # CAPTURE THE nfilter SETTINGS
  thr <- dsBase::listDisclosureSettingsDS()
  nfilter.subset <- as.numeric(thr$nfilter.subset)
  #############################################################
  
  if(length(ids_of_interest) <= nfilter.subset){
    stop("The resulting subset has length <= ", 
         nfilter.subset, "(DataSHIELD nfilter.subset)")
  }
  
  new_f <- tempfile()
  
  if(inherits(genoData_og, "GenotypeData")){
    gdsSubset2(genoData_og@data@filename, new_f,
              sample.include=NULL, snp.include=ids_of_interest,
              sub.storage=NULL,
              compress="LZMA_RA",
              verbose=TRUE,
              allow.fork=TRUE)
  } else if(inherits(genoData_og, "GdsGenotypeReader")){
    gdsSubset2(genoData_og@handler$filename, new_f,
              sample.include=NULL, snp.include=ids_of_interest,
              sub.storage=NULL,
              compress="LZMA_RA",
              verbose=TRUE,
              allow.fork=TRUE)
  } else {
    stop('Provided genotype object is not of class "GenotypeData" or "GdsGenotypeReader"')
  }

  assign(genoData, genoData_og, envir = parent.frame())
  
  new_gds <- GWASTools::GdsGenotypeReader(new_f, allow.fork = TRUE)
  
  # Add pheno information if GenotypeData object
  if(inherits(genoData_og, "GenotypeData")){
    new_gds <- GWASTools::GenotypeData(new_gds, scanAnnot = genoData_og@scanAnnot)
  }
  
  return(new_gds)
  
}