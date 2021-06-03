#' Title
#'
#' @param x 
#' @param slot 
#'
#' @return
#' @export
#'
#' @examples
getVariable <- function(x, slot){
  if(slot %in% c("snp.id", "snp.rs.id", "snp.position", "snp.chromosome", "snp.allele")){
    return(GWASTools::getVariable(x, slot))
  } else {stop('[', slot, '] Not valid slot')}
}