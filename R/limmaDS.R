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
#' 
limmaDS <- function(Set, variable_names, covariable_names, sva){
  res <- MEAL::runPipeline(set = Set, 
                           variable_names = variable_names,
                           covariable_names = covariable_names
                           sva=sva)
  ans <- MEAL::getProbeResults(res)
  return(ans)
}
