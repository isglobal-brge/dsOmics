#' @title Filter potential CpG outliers
#'
#' @param x \code{ExpressionSet} or \code{SummarizedExperiment / RangedSummarizedExperiment} object 
#' to which perform the filtering
#' @param pct \code{numeric} Tail and head quantiles that will be considered as 
#' outliers; for example: a \code{pct=0.125} will use the function 
#' \code{matrixStats::rowQuantiles(probs = c(0.125, 1-0.125))} to detect outliers.
#'
#' @return \code{ExpressionSet} with the outliers marked as NAs
#' @export

removeOutliersDS <- function(x, pct){
  if(inherits(x, "ExpressionSet")){
    probes <- Biobase::exprs(x)
  }
  else if (inherits(x, c("SummarizedExperiment","RangedSummarizedExperiment"))){
    probes <- SummarizedExperiment::assay(x)
  }
  
  rows <- matrixStats::rowQuantiles(probes, probs = c(pct, 1-pct), na.rm = T)
  maskL<- probes < rows[,1]
  maskU<- probes > rows[,2]
  probes[maskL] <- NA
  probes[maskU] <- NA

  data.frame(outliers.lower, outliers.upper, n)
  
  if(inherits(x, "ExpressionSet")){
    Biobase::exprs(x) <- probes
  }
  else if (inherits(x, c("SummarizedExperiment","RangedSummarizedExperiment"))){
    SummarizedExperiment::assay(x) <- probes
  }
  return(x)
}