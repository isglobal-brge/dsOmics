#' @title Retrieve feature names from eSets.
#' @description This function is similar to the Biobase function \code{featureNames}.
#' @param x Object, possibly derived from eSet-class or AnnotatedDataFrame.
#' @return a caracter vector uniquely identifying each feature
#' @author Gonzalez, JR.
#'
#' @export

featureNamesDS <- function(x){
  if(inherits(x, "ExpresionSet"))
    Biobase::featureNames(x)
  else if (inherits(x, "RangedSummarizedExperiment"))
    SummarizedExperiment::assayNames(x)
  else
    stop("implements the proper method")
} 