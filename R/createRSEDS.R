#' @title create a RangedSummarizedExperiment from a RNAseq count table and 
#' a phenotypes table
#' 
#' @details The RNAseq count table has to have the following structure: \cr
#' - First column: Entrez identificator (column named "EntrezID") \cr
#' - Second column: Annotated symbol (optional, but if present must be on this position, if not present, 
#' make sure to set annot = FALSE) \cr
#' - All the remaining columns: One col. per individual \cr
#' The phenotypes table has to have the following structure: \cr
#' - First column: ID of the individuals. To match with the column names 
#' of the RNAseq table \cr
#' - All the remaining columns: One col. per phenotype
#'
#' @param rnaseq \code{data.frame} RNAseq count table
#' @param phe \code{data.frame} Phenotypes table
#' @param ... \code{character vector} Columns of the counts table to be used as 
#' annotation when creating the SummarizedExperiment. If \code{NULL}, biomaRt will be used to find 
#' annotation using the first column of the counts table (EntrezID)
#'
#' @return A \code{RangedSummarizedExperiment} object
#' @export

createRSEDS <- function(rnaseq, phe, ...){
  
  annot_cols <- unlist(list(...))
  
  if(is.null(annot_cols)){
    mart <- biomaRt::useEnsembl(dataset = "hsapiens_gene_ensembl", biomart='ensembl',
                                host = "www.ensembl.org")
    annot <- biomaRt::select(mart, 
                             keys = rnaseq$EntrezID, 
                             columns=c('entrezgene_id', 'chromosome_name',
                                       'start_position', 'end_position', 'hgnc_symbol'),
                             keytype='entrezgene_id')
    annot <- subset(annot, chromosome_name%in%c(1:22, "X", "Y"))
    annot <- annot[!duplicated(annot$entrezgene_id),]
  } else {
    if(!all(annot_cols %in% colnames(rnaseq))){
      stop('Not all annotCols columns are in the counts table. Not found columns: [', 
           paste(annot_cols[!(annot_cols %in% colnames(rnaseq))], collapse = ", "), ']')
    }
    annoted_positions <- !is.na(rnaseq[,2])
    annot <- data.frame(entrezgene_id = rnaseq$EntrezID[annoted_positions])
  }

  # select only genes with annotation
  cc <- dplyr::right_join(rnaseq, annot, by = c("EntrezID" = "entrezgene_id"))
  rownames(cc) <- cc$EntrezID
  # create the RSE
  counts <- as.matrix(cc[ , which(colnames(cc) %in% unlist(phe[,1]))]) # Assuming the ID column
  # is the 1st on the phe table
  
  colData <- S4Vectors::DataFrame(phe)
  rownames(colData) <- colData[,1]
  colData[,1] <- NULL
  
  if(is.null(annot_cols)){
    rowRanges <- GenomicRanges::makeGRangesFromDataFrame(cc, 
                                          seqnames.field = 'chromosome_name',
                                          start.field = 'start_position',
                                          end.field = 'end_position')
    GenomicRanges::mcols(rowRanges) <- data.frame(Gene_Symbol = annot$hgnc_symbol)
    
    rse <- SummarizedExperiment::SummarizedExperiment(assays=S4Vectors::SimpleList(counts=counts),
                                                      rowRanges=rowRanges, colData=colData)
  } else{
    rowData <- rnaseq[, annot_cols][annoted_positions,]
    rownames(rowData) <- cc$EntrezID
    rse <- SummarizedExperiment::SummarizedExperiment(assays=S4Vectors::SimpleList(counts=counts),
                                                      rowData=rowData, 
                                                      colData=colData)}
  return(rse)
}

