#' @title Hardy-Weinberg Equilibrium testing
#' 
#' @description This function performs exact Hardy-Weinberg Equilibrium testing (using Fisher's Test) 
#' over a selection of SNPs. It also counts genotype, calculates allele frequencies, 
#' and calculates inbreeding coefficients.
#'
#' @param genoData 
#' @param sexcol \code{character} Name of the sex column on the covariates file used to create the 
#' \code{\link{GenotypeData}} object
#' @param male \code{character} Name of the male indicator of the sex column on the covariates file used to create the 
#' \code{\link{GenotypeData}} object. (Note that it is case sensitive so it's not the same \code{male} than \code{Male})
#' @param female \code{character} Name of the female indicator of the sex column on the covariates file used to create the 
#' \code{\link{GenotypeData}} object. (Note that it is case sensitive so it's not the same \code{female} than \code{Female})
#' @param geno.counts \code{bool} if \code{TRUE}, genotype counts are returned in the output data.frame
#' @param snpStart \code{numeric} index of the first SNP to analyze, defaults to first SNP (\code{NULL})
#' @param snpEnd \code{numeric} index of the last SNP to analyze, defaults to last SNP (\code{NULL})
#' @param block.size \code{numeric}  number of SNPs to read in at once
#' @param permute \code{bool} logical indicator for whether to permute alleles before calculations
#'
#' @return
#' @export
#'
#' @import dplyr
#'
#' @examples
exactHWEDS <- function(genoData, sexcol, male, female, geno.counts, snpStart, snpEnd, block.size, permute){
  
  if(inherits(genoData, "GenotypeData")){
    genoData@scanAnnot@sexCol <- sexcol
    genoData@scanAnnot@data[,sexcol][genoData@scanAnnot@data[,sexcol] == male] <- "M"
    genoData@scanAnnot@data[,sexcol][genoData@scanAnnot@data[,sexcol] == female] <- "F"
    ans <- GWASTools::exactHWE(genoData = genoData, 
                    geno.counts = geno.counts,
                    snpStart = snpStart,
                    snpEnd = snpEnd,
                    block.size = block.size, 
                    verbose = FALSE,
                    permute = permute)
    return(as_tibble(ans))
  }
  else{
    stop(paste0("Object of incorrect type [", class(genoData), "] exactHWE requires object of type GenotypeData"))
  }
}