#'
#' @title Fits a LM to assess association between the features (outcome) and grouping variable (e.g case/control, condition, ...)
#' @description To be supplied
#' @param feature
#' @param vars
#' @param data
#' @param cellEstim
#' @return a vector with effect estimates, standard error and associated p-value
#' @author Gonzalez, JR.
#'
#' @export 

lmFeatureDS <- function(feature, vars, eSet, cellEstim){
  
  if (!ds.isNA(cellEstim)){
      ds.cbind(c("dat", "cell.counts"), newobj="dat")
  }
    
  dat <- data.frame(eSet[feature, ])[, c(feature, vars)]
  mm <- as.formula(paste(feature, "~ ", 
                         paste(colnames(dat)[-1], collapse="+")))
  mod <- ds.glm(mm, family='gaussian', data='dat', viewIter = FALSE)
  metrics <- as.data.frame(mod$coefficients[2, c(1,2,4)])
  names(metrics) <- feature
  return(metrics)
}
