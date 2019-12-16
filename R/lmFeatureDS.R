#'
#' @title Fits a LM to assess association between the features (outcome) and grouping variable (e.g case/control, condition, ...)
#' @description To be supplied
#' @param feature
#' @param vars
#' @param data
#' @param cellEstim
#' @param connections
#' @return a vector with effect estimates, standard error and associated p-value
#' @author Gonzalez, JR.
#'
#' @export 

lmFeatureDS <- function(feature, vars, data, cellEstim, connections){
  cally <- paste("selFeatureDS(", data, ", feature=", deparse(feature), ", vars=", deparse(vars), ")")
  datashield.assign(connections, 'dat', as.symbol(cally))
  
  if (!ds.isNA(cellEstim)[[1]]){
      ds.cbind(c("dat", "cell.counts"), newobj="dat")
  }
    
  mm <- as.formula(paste(feature, "~ ", paste(ds.colnames('dat')[[1]][-1], collapse="+")))
  mod <- ds.glm(mm, family='gaussian', data='dat', viewIter = FALSE)
  metrics <- as.data.frame(mod$coefficients[2, c(1,2,4)])
  names(metrics) <- feature
  return(metrics)
}