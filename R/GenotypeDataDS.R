#' @title Class GenotypeData. 
#' 
#' @description Container for storing genotype data from a GWAS together with the metadata associated with the subjects (e.g. phenotypes and covariates) and SNPs involved in the study. This is a wrapper of \code{GenotypeData} function from GWASTools package.
#' @param x ..
#' @param covars ...
#' @param columnId ...
#' @return ...
#' @author Gonzalez, JR.
#'
#' @export

GenotypeDataDS <- function(x, covars, columnId){
  names(covars)[columnId] <- "scanID"
  scanAnnot <- GWASTools::ScanAnnotationDataFrame(data.frame(covars))
  geno <- GWASTools::GenotypeData(x, scanAnnot = scanAnnot)
  geno
}