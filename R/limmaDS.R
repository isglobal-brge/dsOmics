#'
#' @title Differential expression analysis using limma
#' @description To be supplied
#' @param model
#' @param eSet
#' @param sva
#' @return a matrix with genes ordered by p-value
#' @author Gonzalez, JR.
#'
#' @export 
limmaDS <- function(model, eSet, sva){
  design <- stats::model.matrix(stats::as.formula(model), 
                                data=Biobase::pData(eSet))
  ee <- Biobase::exprs(eSet)
  fit <- limma::lmFit(eSet, design)
  fit <- limma::eBayes(fit)
  ans <- limma::topTable(fit, n=Inf)
  return(ans)
}
