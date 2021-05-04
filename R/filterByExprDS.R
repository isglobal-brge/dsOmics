#' @title Filter Genes By Expression Level
#' 
#' @description Determine which genes have sufficiently large counts to be retained in a statistical analysis. This is a similar function to \code{edgeR::filterByExpr}
#'
#'
#' @param object \code{eSet}, \code{RangedSummarizedExperiment} on the server
#' @param group vector or factor giving group membership for a oneway layout, if appropriate.
#' 
#' @return A \code{RangedSummarizedExperiment} object with filtered genes
#' @export

filterByExprDS <- function(object, group){
  if(!is.null(group)){
    if(group %in% varLabelsDS(object)){
      keep <- edgeR::filterByExpr(object, 
                                  group=object[[group]])
    } else{
      stop('Group [', group, '] is not present on the object. Check "ds.varLabels()"')
    }
    
  } else{
    keep <- edgeR::filterByExpr(object)
  }
  object.filt <- object[keep,]
  return(object.filt)
}