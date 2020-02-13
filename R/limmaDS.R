#'
#' @title Differential expression analysis using limma
#' @description To be supplied
#' @param Set
#' @param variable_names
#' @param covariable_names
#' @param sva
#' @param annotCols
#' @return a matrix with genes ordered by p-value
#' @author Gonzalez, JR.
#'
#' @export 
#' 
limmaDS <- function(Set, variable_names, covariable_names, type, sva, fNames=NULL){
  
  if (type==2){
    if(inherits(Set, "ExpressionSet")){
      Set.counts <- Biobase::exprs(Set)
      pheno <- Biobase::pData(Set)
    }
    else if (inherits(Set, "RangedSummarizedExperiment")){
      Set.counts <- SummarizedExperiment::assay(Set)
      pheno <- SummarizedExperiment::colData(Set)
    }
    ff <- paste("~", paste(c(variable_names, covariable_names), collapse="+")) 
    design <- model.matrix(formula(ff), data=pheno)
    v <- limma::voom(Set.counts, design = design)$E
    if(inherits(Set, "ExpressionSet"))
      Biobase::exprs(Set) <- v
    else if (inherits(Set, "RangedSummarizedExperiment"))
      SummarizedExperiment::assay(Set) <- v
  }
    
  res <- MEAL::runPipeline(set = Set, 
                           variable_names = variable_names,
                           covariable_names = covariable_names,
                           sva=sva)
  ans <- MEAL::getProbeResults(res, fNames=fNames)
  return(ans)
}
