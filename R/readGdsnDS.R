#' @title Read data file of a GDS node.
#' @description This function is similar to function \code{read.gdsn} in the 'gdsfmt' package.
#' @param x an object of class \code{gdsn.class}, a GDS node
#' @param ... other parameters of  \code{read.gdsn} function
#' 
#' @return An array, list or data.frame
#' @author Gonzalez, JR.
#'
#' @export

readGdsnDS <- function(x, ...){
  gdsfmt::read.gdsn(x, ...)
} 