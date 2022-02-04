#' Title
#'
#' @param gds 
#'
#' @return
#' @export
#'
#' @examples
extractPhenoFromGDSDS <- function(gds){
  if(GWASTools::hasScanAnnotation(gds)){
    pheno <- GWASTools::getScanAnnotation(gds)@data
    return(pheno)
  } else {
    stop('The provided gds object does not have an ScanAnnotation slot (no phenotypes)')
  }
}