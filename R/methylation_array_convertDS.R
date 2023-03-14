#' @title Convert methylation arrays
#'
#' @param eSet `ExpressionSet` to convert
#' @param objective_array `character` Array to which conver (`"450k"` or `"epic"`)
#'
#' @return `ExpressionSet` converted
#' @export
#'
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