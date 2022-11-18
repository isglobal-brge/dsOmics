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
  
  quantiles <- matrixStats::rowQuantiles(probes, probs=c(pct,1-pct), na.rm=T)
  low <- quantiles[,1]
  upper <- quantiles[,2]
  
  outliers.lower <- rowSums(probes < low, na.rm=T)
  outliers.upper <- rowSums(probes > upper, na.rm=T)
  
  idx <- which(probes < low, arr.ind=T)
  probes[idx] <- low[idx[,1]]
  
  idx <- which(probes > upper, arr.ind=T)
  probes[idx] <- upper[idx[,1]]
  
  if(inherits(x, "ExpressionSet")){
    Biobase::exprs(x) <- probes
  }
  else if (inherits(x, c("SummarizedExperiment","RangedSummarizedExperiment"))){
    SummarizedExperiment::assay(x) <- probes
  }
  return(x)
}