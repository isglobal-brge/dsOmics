#'
#' @title Estimate cell counts from methylation data 
#' @description This function uses the function 'meffil.estimate.cell.counts.from.betas' from pacakge 'meffil'
#' @param x an eSet-class object having beta values in the experimental data.
#' @param cellTypeRef Character string name of the cell type reference to use for estimating cell counts. 
#' See meffil.list.cell.type.references() for a list of available references in the package 'meffil'.#' 
#' @return a data frame with samples as rows and variables as columns having the cell count estimates
#' @author Gonzalez, JR.
#' @export
#'

cellCounts <- function(x, cellTypeRef="blood gse35069 complete"){
  xx <- Biobase::exprs(x)
  ans <- try( meffil::meffil.estimate.cell.counts.from.betas(xx,
                                                      cell.type.reference = cellTypeRef, verbose = FALSE), TRUE)
  if (!inherits(ans, "try-error"))
    return(ans)
  else 
    return(NA)
}



