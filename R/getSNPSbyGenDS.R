#' @title Subset GDS with gene(s)
#'
#' @param gds \code{GdsGenotypeReader} GDS object
#' @param old_assign \code{bool} Internal variable used to re-read an old gds, use at your own risk.
#' @param ... \code{character} Genes to subset the GDS
#'
#' @return GDS object
#' @export

getSNPSbyGenDS <- function(gds, old_assign, ...){
  
  if(old_assign){
    gds <- GdsGenotypeReader(gds@filename)
    return(gds)
  }
  
  dots <- unlist(list(...))
  
  df <- tryCatch({
    AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = dots, columns = c("SYMBOL", "ENTREZID"), keytype = "SYMBOL")
  }, error = function(w){
    stop('Queried gene(s) [', paste(dots, collapse = ","), '] were not found on `org.Hs.eg.db` database.')
  }) 
  gr <- GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, filter=list(gene_id=df$ENTREZID))
  GenomeInfoDb::seqlevelsStyle(gr) <- "NCBI"
  
  snp.ranges <- NULL
  for(i in 1:length(dots)){
    snp.ranges <- c(snp.ranges, seq(gr@ranges[i]@start, length.out = gr@ranges[i]@width))
  }

  snp.ids <- GWASTools::getSnpID(gds)[GWASTools::getPosition(gds) %in% snp.ranges]
  
  if(isEmpty(snp.ids)){
    stop("Queried gene(s) not found on the GDS file")
  }
  
  new_gds <- tempfile(fileext = ".gds")
  
  if(inherits(gds, "GdsGenotypeReader")){
    GWASTools::gdsSubset(gds@filename, new_gds, snp.include = snp.ids, allow.fork = TRUE)
    gds_new <- GWASTools::GdsGenotypeReader(new_gds, allow.fork = TRUE)
    return(gds_new)
  } else if(inherits(gds, "GenotypeData")){
    covars <- gds@scanAnnot@data
    GWASTools::gdsSubset(gds@data@filename, new_gds, snp.include = snp.ids, allow.fork = TRUE)
    gds_new <- GWASTools::GdsGenotypeReader(new_gds, allow.fork = TRUE)
    return(GenotypeDataDS(gds_new, covars, paste0(charToRaw("scanID"), collapse = ""), NULL, NULL, NULL, NULL, NULL, NULL))
  } else{
    stop('Objcect "gds" is not of class "GdsGenotypeReader" or "GenotypeData"')
  }
}