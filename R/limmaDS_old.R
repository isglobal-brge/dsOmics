#'
#' @title Differential expression analysis using limma
#' @description To be supplied
#' @param model
#' @param Set
#' @param sva
#' @return a matrix with genes ordered by p-value
#' @author Gonzalez, JR.
#'
#' @export 
limmaDS <- function(model, Set, sva){
  
  if(inherits(Set, "ExpressionSet"))
    exprsData <- Biobase::pData(Set)
  else if (inherits(Set, "RangedSummarizedExperiment")) {
    exprsData <- SummarizedExperiment::colData(Set)
    Set <- t(SummarizedExperiment::assay(Set, 1))
  }
  else
    stop("Set must be a 'eSet' or 'rse' object")
  
  design <- stats::model.matrix(stats::as.formula(model), 
                                data=exprsData)
  if (isTRUE(sva)){
    ds.cbind(c('design', 'sva'), newobj='design')
  }
  
  fit <- limma::lmFit(Set, design)
  fit <- limma::eBayes(fit)
  ans <- limma::topTable(fit, n=Inf)
  return(ans)
}
