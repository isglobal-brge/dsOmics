#' Title
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
make_valid_column_namesDS <- function(x){
  colnames(x) <- make.names(colnames(x))
  return(x)
}