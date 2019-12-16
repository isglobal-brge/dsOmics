#'
#' @title Retrieve feature and sample names from eSets.
#' @description This function is similar to the Biobase function 'featureNames'.
#' @param x Object, possibly derived from eSet-class or AnnotatedDataFrame.
#' @return a caracter vector uniquely identifying each feature
#' @author Gonzalez, JR.
#'
#' @export

featureNamesDS <- function(x){
  Biobase::featureNames(x)
} 