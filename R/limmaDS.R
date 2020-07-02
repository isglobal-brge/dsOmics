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
#' @import dplyr
#' @export 
#' 
limmaDS <- function(Set, variable_names, covariable_names, type, contrasts, levels, coef, sva, annotCols=NULL){
  
   Set<-eval(parse(text=Set))
  
  if (!is.null(covariable_names))
    covariable_names <- unlist(strsplit(covariable_names, split=","))
  if (!is.null(annotCols))
    annotCols <- unlist(strsplit(annotCols, split=","))
 
  
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
  
  if (!is.null(annotCols)){
    annotCols <- unlist(strsplit(annotCols, split=","))
  }
  
  if(!is.null(contrasts) & levels != "design")
  { 
    levels <- unlist(strsplit(levels, split=","))
    contrasts <- unlist(strsplit(contrasts, split=","))
    colnames(design)<-levels
    contrasts<-limma::makeContrasts(contrasts = contrasts,levels = levels)
  }
    
  res <- MEAL::runPipeline(set = Set, 
                           variable_names = variable_names,
                           covariable_names = covariable_names,
                           sva=sva)
  temp <- MEAL::getProbeResults(res, fNames=annotCols, coef = coef, contrast = contrasts)
  ans <- as_tibble(temp) %>% tibble::add_column(.before=1, id=rownames(temp)) %>%
    select(id, tail(names(.), length(annotCols)), everything())
  return(ans)
}
