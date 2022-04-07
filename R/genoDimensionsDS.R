#' @title Get main dimensions of Genotype data
#' 
#' @description Get the number of SNPs, number of scans and number of chromosomes on the genotype file
#'
#' @param x \code{GenotypeData} or \code{GdsGenotypeReader} object
#'
#' @return \code{list} with results
#' @export

genoDimensionsDS <- function(x){
  if(inherits(x, c("GenotypeData", "GdsGenotypeReader"))){
    snp_number <- GWASTools::nsnp(x)
    scan_number <- GWASTools::nscan(x)
    chromosomes <- unique(GWASTools::getChromosome(x))
    return(list(snp_number = snp_number, scan_number = scan_number,
                chromosomes = chromosomes))
  } else {
    stop('[', x, '] Is not a GenotypeData or GdsGenotypeReader')
  }
}