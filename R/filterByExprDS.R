#' @title Filter Genes By Expression Level
#' 
#' @description Determine which genes have sufficiently large counts to be retained in a statistical analysis. This is a similar function to \code{edgeR::filterByExpr}
#'
#'
#' @param object \code{eSet}, \code{RangedSummarizedExperiment} on the server
#'
#' @return A \code{RangedSummarizedExperiment} object with filtered genes
#' @export

filterByExprDS <- function(object){
  keep <- edgeR::filterByExpr(object)
  object.filt <- object[keep,]
  return(object.filt)
}