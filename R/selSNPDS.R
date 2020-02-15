#' @title Subset the information of a given SNP in a GDS object. 
#' 
#' @description This function creates a data frame with a given SNP and variables from a GDS and data.frame objects.
#' @param gds the GDS object
#' @param i the position of the SNP (e.g. column) in the SNP-matrix
#' @param covars the data.frame with the covariates
#' @return a data frame with samples as rows and selected SNP and variables as columns.
#' @author Gonzalez, JR.
#'
#' @export

selSNPDS <- function(gds, i, covars) {
  g <- as.numeric(SNPRelate::snpgdsGetGeno(gds, snp.id = i))
  ans <- data.frame(snp=g, covars)
  ans
}
