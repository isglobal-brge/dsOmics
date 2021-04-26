#' @title Subset the information of a given SNP in a GDS object. 
#' 
#' @description This function creates a data frame with a given SNP and variables from a GDS and data.frame objects.
#' @param genoData the GDS object
#' @param i the position of the SNP (e.g. column) in the SNP-matrix
#' @return a data frame with samples as rows and selected SNP and variables as columns.
#' @author Gonzalez, JR.
#'
#' @export

selSNPDS <- function(genoData, i, strat_variable=NULL, level_number=NULL){ #output_var_factor = NULL) {
  g <- as.numeric(GWASTools::getGenotypeSelection(genoData, i))
  vv <- GWASTools::getScanVariableNames(genoData)
  ans <- data.frame(snp=g, GWASTools::getScanVariable(genoData, vv))
  if(!is.null(strat_variable)){
    #############################################################
    #MODULE 1: CAPTURE THE nfilter SETTINGS                     #
    thr <- listDisclosureSettingsDS()                           #
    nfilter.levels <- as.numeric(thr$nfilter.levels)            #
    #############################################################
    
    factor.levels.present.in.source <- levels(factor(ans[[strat_variable]]))
    num.levels<-length(factor.levels.present.in.source)
    max.allowed.levels<-length(ans[[strat_variable]])*nfilter.levels
    
    if(num.levels>max.allowed.levels)
    {
      error.message<-
        paste0("FAILED: this variable has too many levels and may be disclosive. The ds.asFactor function allows no more than ",
               max.allowed.levels," levels in this particular study. This variable has ",num.levels)
      return(list(error.message=error.message))
    } else{
      ans <- ans[ans[[strat_variable]] == level_number,]
    }
  }
  ans
}
