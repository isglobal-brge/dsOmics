#'
#' @title Surrogate variable estimates 
#' @description This function a data frame with surrogate variables used to correct for unwanted variability
#' @param x an eSet-class object having beta values in the experimental data.
#' @param feature a character vector indicanting the selected feature.
#' @param vars a character vector with the variables to be selected.
#' @return a data frame with samples as rows and selected features and variables as columns.
#' @author Gonzalez, JR.
#'
#' @export

svaDS <- function(vars, eSet){
  edata <- Biobase::exprs(eSet)
  pheno <- Biobase::pData(eSet)
  mm <- paste(" ~ ", paste(vars, collapse=" + "))
  
  # think about missings
  mod <- stats::model.matrix( as.formula(mm), data=pheno)
  if(nrow(mod)!=nrow(pheno))
    stop("There are missing values. Select complete cases or impute them.")
  mod0 <- stats::model.matrix( ~ 1, data=pheno)
  
  
  n.sv <- try(sva::num.sv(edata, mod))
  if (!inherits(n.sv, "try-error"))
    return(NA)
  else {
    svobj <- try(sva::sva(edata, mod, mod0, n.sv=n.sv), TRUE)
     if (!inherits(svobj, "try-error"))
      return(svobj$sv)
     else 
    return(NA)
  }
}
