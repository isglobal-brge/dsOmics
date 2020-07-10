#'
#' @title Differential expression analysis on the server-side
#' @description Performs the DGE based on
#'  the negative binomial distribution as a model for count variability,
#'  including empirical Bayes methods, exact tests, 
#'  and generalized linear models.
#' @details This function perform a DGE analysis 
#' for \code{SummarizedExperiment} input 
#' 
#' @author L. Abarrategui for DataSHILED development team
#' @export
#'

edgeRDS<-function(vars, set)
{ 
  set<-eval(parse(text=set))
  
  #design formula
  vars <- unlist(strsplit(vars, split=","))
  ff <- paste("~", paste(c(vars), collapse="+"))
  
  
  return(res)
}

#AGGREGATE FUNCTION
#edgeRDS