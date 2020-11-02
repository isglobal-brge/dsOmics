#' Title
#'
#' @param genoData 
#'
#' @return
#' @export
#'
#' @examples
alleleFrequencyDS <- function(genoData, sexcol, male, female){
  
  
  if(inherits(genoData, "GenotypeData")){
    return(list(class(genoData@scanAnnot@data), genoData@scanAnnot@data))
    genoData@scanAnnot@sexCol <- sexcol
    genoData@scanAnnot@data[,sexcol][genoData@scanAnnot@data[,sexcol] == male] <- "M"
    genoData@scanAnnot@data[,sexcol][genoData@scanAnnot@data[,sexcol] == female] <- "F"
    
    return(GWASTools::alleleFrequency(genoData, verbose = FALSE))
  }
  else(stop(paste0("Object of incorrect type [", class(genoData), "] alleleFrequency requires object of type GenotypeData")))
  
}