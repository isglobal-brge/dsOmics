#' @title Create a list of candidate SNPs based on specified criteria
#' @description This function is similar to function \code{snpgdsSelectSNPDS} in the 'SNPRelate' package.
#' @param x an object of class \code{SNPGDSFileClass}, a SNP GDS file; or characters to specify the file name of SNP GDS
#' @param sample.id
#' @param snp.id
#' @param ... other parameters of  \code{snpgdsSelectSNPDS} function
#' @return The function returns an integer matrix with values 0, 1, 2 or NA representing the number of reference allele when with.id=FALSE; or list(genotype, sample.id, snp.id) when with.id=TRUE. The orders of sample and SNP IDs in the genotype matrix are actually consistent with sample.id and snp.id in the GDS file, which may not be as the same as the arguments sampel.id and snp.id specified by users.
#' @author Gonzalez, JR.
#'
#' @export

snpgdsSelectSNPDS <- function(x, sample.id=NULL, snp.id=NULL, ...){
  if (!is.null(sample.id))
    sample.id <- unlist(strsplit(sample.id, split=","))
  if (!is.null(snp.id))
    snp.id <- unlist(strsplit(snp.id, split=","))
  SNPRelate::snpgdsSelectSNP(x, sample.id=sample.id, snp.id=snp.id, ...)
} 