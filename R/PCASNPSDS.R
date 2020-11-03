#' Title
#'
#' @param gds 
#' @param ld.threshold 
#'
#' @return
#' @export
#'
#' @examples
PCASNPSDS <- function(gds, prune, ld.threshold){
  gdsPCA <- SNPRelate::snpgdsOpen(gds@filename, allow.duplicate = TRUE)
  if(prune){
    snpset <- SNPRelate::snpgdsLDpruning(gdsPCA, ld.threshold=ld.threshold)
    snpset.id <- unlist(unname(snpset))
    pca <- SNPRelate::snpgdsPCA(gdsPCA, snp.id=snpset.id, num.thread=2)
  }
  else{
    pca <- SNPRelate::snpgdsPCA(gdsPCA, num.thread=2)
  }
  
  tab <- data.frame(sample.id = pca$sample.id,
                    EV1 = pca$eigenvect[,1], # the first eigenvector
                    EV2 = pca$eigenvect[,2], # the second eigenvector
                    stringsAsFactors = FALSE)
  return(tab)
}