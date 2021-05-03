#' Title
#'
#' @param x 
#' @param group 
#' @param geno 
#'
#' @return
#' @export
#'
#' @examples
plotPCASNPSDS <- function(x, geno, group){
  
  if(!is.null(geno)){
    if(!inherits(geno, 'GenotypeData')){
      stop('Geno object is not of class "GenotypeData"')
    }
  }
  
  #############################################################
  #MODULE 1: CAPTURE THE nfilter SETTINGS                     #
  thr <- dsBase:::listDisclosureSettingsDS()                  #
  nfilter.levels <- as.numeric(thr$nfilter.levels)            #
  #############################################################
  
  feno <- GWASTools::getScanAnnotation(geno)@data
  if(is.null(feno[[group]])){
    stop('Grouping variable [', group, '] is not inside of the genotype file. Check the column
         names of the covariates file used to create the GenotypeData file')
  }
  factor.levels.present.in.source <- levels(factor(feno[[group]]))
  num.levels<-length(factor.levels.present.in.source)
  max.allowed.levels<-length(feno[[group]])*nfilter.levels
  
  if(num.levels>max.allowed.levels)
  {
    error.message<-
      paste0("FAILED: this variable has too many levels and may be disclosive. The ds.asFactor function allows no more than ",
             max.allowed.levels," levels in this particular study. This variable has ",num.levels)
    return(list(error.message=error.message))
  } else{
    if(!all(unlist(feno$scanID) %in% GWASTools::getScanID(x))){
      stop('The covariates file does not contain all the individuals on the GDS')
    } else {
      return(feno[match(GWASTools::getScanID(x), unlist(feno$scanID)),][[group]])
    }
  }
}