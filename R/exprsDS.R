#'
#' @title Retrieve information on molecular data recorded in eSet and ExpressionSet-derived classes.
#' @description This function is similar to the Biobase function 'exprs'.
#' @param x Object, possibly derived from eSet-class.
#' @return a matrix features as rows, individuals/samples as columns
#' @author Gonzalez, JR.
#'
#' @export 
 
exprsDS <- function(x) {
  Biobase::exprs(x)
}
