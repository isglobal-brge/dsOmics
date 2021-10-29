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
    results <- do.call(c,lapply(x, function(i){
      GWASTools::getVariable(i, slot)
    }))
    return(results)
  } else {stop('[', slot, '] Not valid slot')}
}