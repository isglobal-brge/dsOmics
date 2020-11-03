#' @title Get names of chromosomes
#' 
#' @description Get the names of the chromosomes of a GenotypeData object
#'
#' @param genoData \code{GenotypeData} object
#'
#' @return A \code{list} with: \cr
#' -autosomes: integer codes for the autosomes \cr
#' -Xchr: integer code for the X chromosome \cr
#' -pseudoautosomalRegionXY: integer code for the pseudoautosomal region \cr
#' -Ychr: integer code for the Y chromosome \cr
#' -mitochondrial: integer code for mitochondrial SNPs \cr
#' 
#' @export

getChromosomeNamesDS <- function(genoData){
  
  autosomeCode <- autosomeCode(genoData) # Returns the integer codes for the autosomes.
  XchromCode <- XchromCode(genoData) # Returns the integer code for the X chromosome.
  XYchromCode <- XYchromCode(genoData) # Returns the integer code for the pseudoautosomal region.
  YchromCode <- YchromCode(genoData) # Returns the integer code for the Y chromosome.
  MchromCode <- MchromCode(genoData) # Returns the integer code for mitochondrial SNPs.
  
  return(list(autosomes = autosomeCode, Xchr = XchromCode, pseudoautosomalRegionXY = XYchromCode, 
              Ychr = YchromCode, mitochondrial = MchromCode))
  
}