#' @title Class GenotypeData. 
#' 
#' @description Container for storing genotype data from a GWAS together with the metadata associated with the subjects (e.g. phenotypes and covariates) and SNPs involved in the study. This is a wrapper of \code{GenotypeData} function from GWASTools package.
#' @param x an eSet-class object having beta values in the experimental data.
#' @param feature a character vector indicanting the selected feature.
#' @param vars a character vector with the variables to be selected.
#' @return a data frame with samples as rows and selected features and variables as columns.
#' @author Gonzalez, JR.
#'
#' @export

GenotypeDataDS <- function(x, covars, columId, ...){
  g <- GWASTools::GdsGenotypeReader(x)
  names(covars)[columnId] <- "scanID"
  scanAnnot <- GWASTools::ScanAnnotationDataFrame(covars)
  geno <- GenotypeData(g, scanAnnot = scanAnnot)
  close(g)
  geno
}
