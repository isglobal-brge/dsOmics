#' @title Retrieve feature names from eSets.
#' @description This function is similar to the Biobase function \code{featureNames}.
#' @param x Object, possibly derived from eSet-class or AnnotatedDataFrame.
#' @return a caracter vector uniquely identifying each feature
#' @author Gonzalez, JR.
#'
#' @export

varLabelsDS <- function(x){
  if(inherits(x, "ExpressionSet"))
    return(Biobase::varLabels(x))
  else if (inherits(x, c("SummarizedExperiment","RangedSummarizedExperiment")))
    return(colnames(SummarizedExperiment::colData(x)))
  else if(inherits(x, "GenotypeData")){
    return(Biobase::varLabels(GWASTools::getScanAnnotation(x)))
  }
  else{
    stop("implements the proper method")
  }
    
} 