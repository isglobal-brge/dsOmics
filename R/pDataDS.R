#'
#' @title Retrieve information on experimental phenotypes recorded in eSet and ExpressionSet-derived classes.
#' @description This function is similar to the Biobase function 'pData'.
#' @param x Object, possibly derived from eSet-class or AnnotatedDataFrame.
#' @return a data frame with samples as rows, variables as columns
#' @author Gonzalez, JR.
#'
#' @export
 
pDataDS <- function(x) {
  Biobase::pData(x)
}