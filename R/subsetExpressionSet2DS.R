#' @title Subset ExpressionSet
#' 
#' @description 
#'
#' @return Subseted \code{ExpressionSet}
#' @export

subsetExpressionSet2DS <- function(eSet, n_ids){
  
  neweSet <- eSet[,1:n_ids]
  
  return(neweSet)
  
}
