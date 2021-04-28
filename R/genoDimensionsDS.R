#' Title
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
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