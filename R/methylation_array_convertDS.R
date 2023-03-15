#' @title Convert a Methylation ExpressionSet from EPIC to 450k Array
#'
#' @description This function subsets a methylation ExpressionSet from an EPIC array (850k probes) to a 450k array or vice versa.
#'
#' @param eSet A methylation ExpressionSet object containing either EPIC or 450k array data.
#' @param objective_array A character string specifying the target array type. Can be either "450k" or "epic".
#'
#' @return An ExpressionSet object containing the subset of probes based on the target array type.
#' @export
#'
#' @examples
#' # Load example data as 'eSet'
#' # Convert the eSet from EPIC array to 450k array
#' eSet_450k <- methylation_array_convertDS(eSet, objective_array = "450k")
#' # Convert the eSet from 450k array to EPIC array
#' eSet_epic <- methylation_array_convertDS(eSet, objective_array = "epic")

methylation_array_convertDS <- function(eSet, objective_array){

  if(objective_array == "450k"){
    cpgs <- rownames(IlluminaHumanMethylation450kanno.ilmn12.hg19::SNPs.Illumina)
  } else if (objective_array == "epic"){
    cpgs <- rownames(IlluminaHumanMethylationEPICanno.ilm10b2.hg19::SNPs.Illumina)
  } else {
    stop()
  }

  eSet_converted <- eSet[rownames(eSet) %in% cpgs, ]
  return(eSet_converted)
}