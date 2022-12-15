#'
#' @title Differential expression analysis using limma on the server-side
#' @description Performs differential expression analysis using LIMMA
#' @param Set either a \code{ExpressionSet} or a \code{RangedSummarizedExperiment}
#' @param variable_names grouping variable used to perform differential expression analysis
#' @param covariable_names name of variables used in the adjusted models
#' @param type ...
#' @param contrasts ...
#' @param levels ...
#' @param coef ...
#' @param sva should differential expression analysis be adjusted by SVA?
#' @param annotCols variables from the annotation data used in the output
#' @param method String indicating the method used in the regression (e.g. lmFit function of limma: "ls" or 
#' "robust". (Default: "ls") 
#' @param robust Logical value indicating whether robust method is applied in the eBayes function of limma. Default is FALSE.
#' @param normalization String indicating the normalize method used when using voom for RNAseq data (see normalized.method argument in limma::vomm)
#' @param voomQualityWeights Logical value indicating whether limma::voomWithQualityWeights should be used instead of
#' limma::voom. 
#' @param big Logical value indicating whether SmartSVA should be used instead of SVA 
#' (TRUE recommended for methylation or when having large number of samples). 
#' 
#' 
#' @return a matrix with genes ordered by p-value
#' @author Gonzalez, JR.
#' 
#' @import dplyr
#' @export 
#' 
limma2DS <- function(Set, variable_names, covariable_names, type, contrasts, 
                    levels, coef, sva, annotCols, method, robust, normalization,
                    voomQualityWeights, big, sort.by){
  
  res <- dsOmics::limmaDS(Set, variable_names, covariable_names, type, contrasts, 
                   levels, coef, sva, annotCols, method, robust, normalization,
                   voomQualityWeights, big, sort.by)
  return(res)
}
