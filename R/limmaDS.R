#'
#' @title Differential expression analysis using limma
#' @description Performs differential expression analysis using LIMMA
#' @param Set either a \code{ExpressionSet} or a \code{RangedSummarizedExperiment}
#' @param variable_names grouping variable used to perform differential expression analysis
#' @param covariable_names name of variables used in the adjusted models
#' @param sva should differential expression analysis be adjusted by SVA?
#' @param annotCols variables from the annotation data used in the output
#' @return a matrix with genes ordered by p-value
#' @author Gonzalez, JR.
#'
#' @export 
#' 
limmaDS <- function(Set, variable_names, covariable_names, type, sva, annotCols=NULL){
  
  if (type==2){
    if(inherits(Set, "ExpressionSet")){
      Set.counts <- Biobase::exprs(Set)
      pheno <- Biobase::pData(Set)
    }
    else if (inherits(Set, "RangedSummarizedExperiment")){
      Set.counts <- SummarizedExperiment::assay(Set)
      pheno <- SummarizedExperiment::colData(Set)
    }
    covariable_names <- unlist(strsplit(covariable_names, split=","))
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
  ans <- dplyr::as_tibble(MEAL::getProbeResults(res, fNames=annotCols))
  return(ans)
}
