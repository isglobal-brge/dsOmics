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
limmaDS <- function(model, Set, sva){
  res <- runDiffMeanAnalysis(set = Set, model = model, sva=sva)
  ans <- getProbeResults(res)
  return(ans)
}
