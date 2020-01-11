#' @title Retrieve information on features recorded in eSet-derived objects.
#' 
#' @description This function is similar to the Biobase function \code{featureData}.
#' @param object Object, possibly derived from \code{eSet-class}.
#' @return an object containing information on both variable values and variable meta-data.
#' @author Gonzalez, JR.
#'
#' @export
#' 
featureDataDS <- function(x){
  Biobase::featureData(x)
} 
