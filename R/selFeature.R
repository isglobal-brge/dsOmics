#'
#' @title Subset the information of a given feature in a eSet object. 
#' @description This function creates a data frame with a given feature and variables of a eSet object.
#' @param x an eSet-class object having beta values in the experimental data.
#' @param feature a character vector indicanting the selected feature.
#' @param vars a character vector with the variables to be selected.
#' @return a data frame with samples as rows and selected features and variables as columns.
#' @author Gonzalez, JR.
#'

selFeature <- function(x, feature, vars){
  ans <- data.frame(x[feature, ])[, c(feature, vars)]
  ans
}
