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
  else if (inherits(x, "RangedSummarizedExperiment")){
    df <- as.data.frame(SummarizedExperiment::rowRanges(x), 
                        stringsAsFactors = FALSE)
    colnames(df)[1] <- "chromosome"
    df[, 1] <- as.character(df[, 1])
    return(colnames(df))
  } else if(inherits(x, "SummarizedExperiment")){
    df <- as.data.frame(SummarizedExperiment::rowData(x), 
                        stringsAsFactors = FALSE)
    return(colnames(df))
  }
  else
    stop("implements the proper method")
} 