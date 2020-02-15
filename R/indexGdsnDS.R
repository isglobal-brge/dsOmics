#' @title Get the specific GDS node
#' @description This function is similar to function \code{index.gdsn} in the 'gdsfmt' package.
#' @param x an object of class \code{gdsn.class} (a GDS node), or \code{gds.class} (a GDS file)
#' @param ... other parameters of  \code{index.gdsn} function
#' 
#' @return An object of class \code{gdsn.class} for the specified node
#' @author Gonzalez, JR.
#'
#' @export

indexGdsnDS <- function(x, ...){
  gdsfmt::index.gdsn(x, ...)
} 