#' @title Filter potential CpG outliers
#'
#' @param x \code{ExpressionSet} or \code{SummarizedExperiment / RangedSummarizedExperiment} object 
#' to which perform the filtering
#'
#' @return
#' @export
#'
#' @examples
removeOutliersDS <- function(x){
  if(inherits(x, "ExpressionSet")){
    probes <- Biobase::exprs(x)
  }
  else if (inherits(x, c("SummarizedExperiment","RangedSummarizedExperiment"))){
    probes <- SummarizedExperiment::assay(x)
  }
  rowIQR<- matrixStats::rowIQRs(probes, na.rm = TRUE)
  row2575 <- matrixStats::rowQuantiles(probes, probs = c(0.25, 0.75), na.rm = T)
  maskL<- probes < row2575[,1] - 3 * rowIQR
  maskU<- probes > row2575[,2] + 3 * rowIQR
  probes[maskL] <- NA
  probes[maskU] <- NA
  if(inherits(x, "ExpressionSet")){
    Biobase::exprs(x) <- probes
  }
  else if (inherits(x, c("SummarizedExperiment","RangedSummarizedExperiment"))){
    SummarizedExperiment::assay(x) <- probes
  }
  return(x)
}