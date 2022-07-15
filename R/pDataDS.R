#' @title Retrieve information on experimental phenotypes recorded in eSet and ExpressionSet-derived classes
#'
#' @param x Object (on the study server), possibly derived from \link{eSet-class} or \link{AnnotatedDataFrame}
#'
#' @return  Phenotypes \code{data.frame}
#' @export

pDataDS <- function(x){
  return(Biobase::pData(x))
}