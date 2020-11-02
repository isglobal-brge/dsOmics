#' @title Allelic frequency
#' 
#' @description Calculates the frequency of the A allele
#'
#' @param genoData \code{GenotypeData} \code{\link{GenotypeData}} object
#' @param sexcol \code{character} Name of the sex column on the covariates file used to create the 
#' \code{\link{GenotypeData}} object
#' @param male \code{character} Name of the male indicator of the sex column on the covariates file used to create the 
#' \code{\link{GenotypeData}} object. (Note that it is case sensitive so it's not the same \code{male} than \code{Male})
#' @param female \code{character} Name of the female indicator of the sex column on the covariates file used to create the 
#' \code{\link{GenotypeData}} object. (Note that it is case sensitive so it's not the same \code{female} than \code{Female})
#'
#' @return A matrix with a row for each SNP. Columns "M" for males, "F" for females, and "all" for all scans give 
#' frequencies of the A allele. Sample size for males, females, and all is returned as "n.M", "n.F", and "n", respectively.
#'  "MAF" is the minor allele frequency over all scans.
#' @export

alleleFrequencyDS <- function(genoData, sexcol, male, female){
  
  if(inherits(genoData, "GenotypeData")){
    genoData@scanAnnot@sexCol <- sexcol
    genoData@scanAnnot@data[,sexcol][genoData@scanAnnot@data[,sexcol] == male] <- "M"
    genoData@scanAnnot@data[,sexcol][genoData@scanAnnot@data[,sexcol] == female] <- "F"
    
    return(GWASTools::alleleFrequency(genoData, verbose = FALSE))
  }
  else(stop(paste0("Object of incorrect type [", class(genoData), "] alleleFrequency requires object of type GenotypeData")))
}