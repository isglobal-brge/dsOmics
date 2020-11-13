#' @title Subset the information of a given SNP in a GDS object. 
#' 
#' @description This function creates a data frame with a given SNP and variables from a GDS and data.frame objects.
#' @param gds the GDS object
#' @param i the position of the SNP (e.g. column) in the SNP-matrix
#' @return a data frame with samples as rows and selected SNP and variables as columns.
#' @author Gonzalez, JR.
#'
#' @export

selSNPDS <- function(genoData, i) {
  g <- as.numeric(GWASTools::getGenotypeSelection(genoData, i))
  vv <- getScanVariableNames(genoData)
  ans <- data.frame(snp=g, getScanVariable(genoData, vv))
  ans
}
