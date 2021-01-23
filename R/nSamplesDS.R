#' @title Retrieve number of samples from eSets and RangedSummarizedExperiments.
#' @description This function is similar to the Biobase function \code{dim}.
#' @param x Object, possibly derived from eSet-class or RangedSummarizedExperiment.
#' @return a numeric value
#' @author Gonzalez, JR.
#'
#' @export

nSamplesDS <- function(x){
  if(inherits(x, c("ExpressionSet", "SummarizedExperiment","RangedSummarizedExperiment")))
    return(dim(x)[2])
  else
    stop("implements the proper method")
} 