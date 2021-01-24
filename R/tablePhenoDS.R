#' @title Retrieve the number of counts at each factor levels for metadata covariates of eSets and RangedSummarizedExperiments.
#' @description This function is similar to the base function \code{table}.
#' @param x Object, possibly derived from eSet-class or RangedSummarizedExperiment.
#' @return a numeric value for each factor level
#' @author Gonzalez, JR.
#'
#' @export

tablePhenoDS <- function(x){
  if(inherits(x, c("ExpressionSet", "SummarizedExperiment","RangedSummarizedExperiment")))
    return(table(x))
  else
    stop("implements the proper method")
} 