#' Title
#'
#' @param genoData 
#' @param geno.counts 
#' @param snpStart 
#' @param snpEnd 
#' @param block.size 
#' @param permute 
#'
#' @return
#' @export
#'
#' @examples
exactHWEDS <- function(genoData, sexcol, male, female, geno.counts, snpStart, snpEnd, block.size, permute){
  
  if(inherits(genoData, "GenotypeData")){
    genoData@scanAnnot@sexCol <- sexcol
    genoData@scanAnnot@data[,sexcol][genoData@scanAnnot@data[,sexcol] == male] <- "M"
    genoData@scanAnnot@data[,sexcol][genoData@scanAnnot@data[,sexcol] == female] <- "F"
    
    ans <- exactHWE(genoData = genoData, 
                    geno.counts = geno.counts,
                    snpStart = snpStart,
                    snpEnd = snpEnd,
                    block.size = block.size, 
                    verbose = FALSE,
                    permute = permute)
    return(ans)
  }
  else{
    stop(paste0("Object of incorrect type [", class(genoData), "] exactHWE requires object of type GenotypeData"))
  }
}