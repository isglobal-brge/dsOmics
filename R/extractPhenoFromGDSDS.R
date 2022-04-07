#' @title Extract phenotype table from GenotypeData object
#'
#' @param gds \code{GenotypeData} object
#'
#' @return \code{data.table} of the phenotypes
#' @export

extractPhenoFromGDSDS <- function(gds){
  if(GWASTools::hasScanAnnotation(gds)){
    pheno <- GWASTools::getScanAnnotation(gds)@data
    return(pheno)
  } else {
    stop('The provided gds object does not have an ScanAnnotation slot (no phenotypes)')
  }
}