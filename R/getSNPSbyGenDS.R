#' Title
#'
#' @param gds 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
getSNPSbyGenDS <- function(gds, old_assign, ...){

  if(old_assign){
    gds <- GdsGenotypeReader(gds@filename)
    return(gds)
  }
  
  dots <- unlist(list(...))
  
  df <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = dots, columns = c("SYMBOL", "ENTREZID"), keytype = "SYMBOL")
  gr <- GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, filter=list(gene_id=df$ENTREZID))
  GenomeInfoDb::seqlevelsStyle(gr) <- "NCBI"
  
  snp.ranges <- NULL
  for(i in 1:length(dots)){
    snp.ranges <- c(snp.ranges, seq(gr@ranges[i]@start, length.out = gr@ranges[i]@width))
  }
  
  snp.ids <- GWASTools::getSnpID(gds)[GWASTools::getPosition(gds) %in% snp.ranges]
  
  new_gds <- tempfile(fileext = ".gds")
  GWASTools::close(gds)
  GWASTools::gdsSubset(gds@filename, new_gds, snp.include = snp.ids, allow.fork = TRUE)
  
  gds_new <- GWASTools::GdsGenotypeReader(new_gds)
  
  return(gds_new)
  
}