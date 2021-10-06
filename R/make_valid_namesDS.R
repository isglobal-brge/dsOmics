#' @title Rename column names and character variables
#'
#' @description Passes the column names and the character variables through \code{make.names}, this creates 
#' names without any special characters so there are no problems for accessing any columns or any category 
#' when using the DataSHIELD parser.
#'
#' @param x \code{data.frame} Entry table to modify
#'
#' @return \code{data.frame} with the column names and character columns updated
#' @export
#'
#' @examples
make_valid_namesDS <- function(x){
  x <- eval(parse(text=x), envir = parent.frame())
  colnames(x) <- make.names(colnames(x))
  x <- x %>% rowwise() %>% mutate(across(where(is.character), ~ if(!is.na(.x)){make.names(.x)}else{NA}))
  return(x)
}