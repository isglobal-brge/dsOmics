#' @title Filter Genes By Expression Level
#' 
#' @description Determine which genes have sufficiently large counts to be retained in a statistical analysis. This is a similar function to \code{edgeR::filterByExpr}
#'
#'
#' @param object \code{eSet}, \code{RangedSummarizedExperiment} on the server
#' @param group vector or factor giving group membership for a oneway layout, if appropriate.
#' 
#' @return A \code{RangedSummarizedExperiment} object with filtered genes
#' @export

RNAseqPreprocDS <- function(object, group){
  
  #
  # edgeR-limma pipeline (F1000 paper)
  #
  
  # step 1 :  calculate logCPM units of expression
  dge <- edgeR::DGEList(counts=SummarizedExperiment::assays(object)$counts, 
                        genes=as.data.frame(SummarizedExperiment::rowData(object)))
  SummarizedExperiment::assays(object)$logCPM <- edgeR::cpm(dge, log=TRUE, prior.count=0.5)
  
  # step 2 :  filtering of lowly-expressed genes
  
  if(!is.null(group)){
    if(group %in% varLabelsDS(object)){
      keep <- edgeR::filterByExpr(dge, 
                                  group=object[[group]])
    } else{
      stop('Group [', group, '] is not present on the object. Check "ds.varLabels()"')
    }
    
  } else{
    keep <- edgeR::filterByExpr(dge)
  }
  dge.filt <- dge[keep, ]
  object.filt <- object[keep,]
  
  # step 4: normalization
  dge.filt.noTMM <- dge.filt
  dge.filt <- edgeR::calcNormFactors(dge.filt)
  
  SummarizedExperiment::assays(object.filt)$logCPM <- edgeR::cpm(dge.filt, normalized.lib.sizes=TRUE,
                                                                 log=TRUE, prior.count=0.5)
  
  SummarizedExperiment::assays(object.filt)$noTMM <- dge.filt.noTMM
  
  
  return(object.filt)
}