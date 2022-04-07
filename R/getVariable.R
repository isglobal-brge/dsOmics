#' @title Get slot of GDS object
#' 
#' @description Available slots are limited to: \code{snp.id}, \code{snp.rs.id}, \code{snp.position}, 
#' \code{snp.chromosome} and \code{snp.allele}
#'
#' @param x \code{GenotypeData} or \code{GdsGenotypeReader} object
#' @param slot \code{character} slot to retrieve
#'
#' @return Vector of characters
#' @export

getVariable <- function(x, slot){
  if(slot %in% c("snp.id", "snp.rs.id", "snp.position", "snp.chromosome", "snp.allele")){
    if(inherits(x, "list")){
      results <- do.call(c,lapply(x, function(i){
        GWASTools::getVariable(i, slot)
      }))
    } else {
      results <- GWASTools::getVariable(x, slot)
    }
    return(results)
  } else {stop('[', slot, '] Not valid slot')}
}