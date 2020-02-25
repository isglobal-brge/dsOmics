#' @title Retrieve names of annotation from eSets.
#' @description This function is similar to the Biobase function \code{featureNames}.
#' @param x Object, possibly derived from eSet-class or AnnotatedDataFrame.
#' @return a caracter vector uniquely identifying each feature
#' @author Gonzalez, JR.
#'
#' @export

fvarLabelsDS <- function(x){
  if(inherits(x, "ExpressionSet"))
    return(Biobase::fvarLabels(x))
  else if (inherits(x, "RangedSummarizedExperiment"))
    stop("implements the proper method")
  else
    stop("implements the proper method")
} 