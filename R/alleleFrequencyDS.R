#' @title Allelic frequency
#' 
#' @description Calculates the frequency of the A allele
#'
#' @param genoData \code{GenotypeData} \code{\link{GenotypeData}} object
#'
#' @return A matrix with a row for each SNP. Columns "M" for males, "F" for females, and "all" for all scans give 
#' frequencies of the A allele. Sample size for males, females, and all is returned as "n.M", "n.F", and "n", respectively.
#'  "MAF" is the minor allele frequency over all scans.
#' @export
#' 
#' @import dplyr

alleleFrequencyDS <- function(genoData){
  if(inherits(genoData, "GenotypeData")){
    ans <- GWASTools::alleleFrequency(genoData, verbose = FALSE)
    rs <- GWASTools::getVariable(genoData, "snp.rs.id")
    return(tibble::as_tibble(ans) %>% tibble::add_column(rs=rs) %>% dplyr::relocate(rs))
  }
  else(stop(paste0("Object of incorrect type [", class(genoData), "] alleleFrequency requires object of type GenotypeData")))
}