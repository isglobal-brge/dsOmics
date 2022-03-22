#' Title
#'
#' @param object 
#' @param group 
#'
#' @return
#' @export
#'
#' @examples
plotPCADS <- function(object, geno, group){

  if(inherits(geno, "GdsGenotypeReader") & !missing(group)){
    stop("Grouping is not available for the dataset used to perform the PCA, make sure to perform the PCA with an object created by the function ds.GenotypeData() if grouping is desired on the plot.")
  }
  
  #############################################################
  #MODULE 1: CAPTURE THE nfilter SETTINGS                     #
  thr <- dsBase:::listDisclosureSettingsDS()                  #
  nfilter.levels.density <- as.numeric(thr$nfilter.levels) #
  #############################################################
  
  feno <- GWASTools::getScanAnnotation(geno)@data
  if(is.null(feno[[group]])){
    stop('Grouping variable [', group, '] is not inside of the genotype file. Check the column names of the covariates file used to create the GenotypeData file')
  }
  
  factor.levels.present.in.source <- levels(factor(feno[[group]]))
  num.levels<-length(factor.levels.present.in.source)
  max.levels.by.density<-nfilter.levels.density*length(feno[[group]])
  
  
  if(num.levels>max.levels.by.density)
  {
    error.message<-
      paste0("FAILED: this variable has too many levels and may be disclosive. The ds.asFactor function allows no more than ",
             max.allowed.levels," levels in this particular study. This variable has ",num.levels)
    return(list(error.message=error.message))
  } 
  else{
    if(!all(unlist(feno$scanID) %in% GWASTools::getScanID(geno))){
      stop('The covariates file does not contain all the individuals on the GDS')
    } else {
      return(feno[match(GWASTools::getScanID(geno), unlist(feno$scanID)),][[group]])
    }
  }
}