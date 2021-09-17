#' @title Rename column names
#' @description Passes the column names through \code{make.names} and rewrites them
#'
#' @param x \code{data.frame} Table to change the column names
#'
#' @return \code{data.frame} with the column names updated
#' @export
#'
#' @examples
make_valid_column_namesDS <- function(x){
  colnames(x) <- make.names(colnames(x))
  return(x)
}