#' @title reate a RangedSummarizedExperiment from a RNAseq count table and 
#' a phenotypes table
#' 
#' @details The RNAseq count table has to have the following structure: \cr
#' - First column: Entrez identificator \cr
#' - Second column: Annotated symbol \cr
#' - All the remaining columns: One col. per individual \cr
#' The phenotypes table has to have the following structure: \cr
#' - First column: ID of the individuals. To match with the column names 
#' of the RNAseq table \cr
#' - All the remaining columns: One col. per phenotype
#'
#' @param rnaseq \code{data.frame} RNAseq count table
#' @param phe \code{data.frame} Phenotypes table
#'
#' @return A \code{RangedSummarizedExperiment} object
#' @export

createRSEDS <- function(rnaseq, phe){
  mart <- biomaRt::useEnsembl(dataset = "hsapiens_gene_ensembl", biomart='ensembl',
                     host = "www.ensembl.org")
  annot <- biomaRt::select(mart, 
                           keys = rnaseq$EntrezID, 
                           columns=c('entrezgene_id', 'chromosome_name',
                                     'start_position', 'end_position'),
                           keytype='entrezgene_id')
  annot <- subset(annot, chromosome_name%in%c(1:22, "X", "Y"))
  annot <- annot[!duplicated(annot$entrezgene_id),]
  # select only genes with annotation
  cc <- dplyr::right_join(rnaseq, annot, by = c("EntrezID" = "entrezgene_id"))
  rownames(cc) <- cc$EntrezID
  # create the RSE
  counts <- as.matrix(cc[ , which(colnames(cc) %in% unlist(phe[,1]))]) # Assuming the ID column
                                                                       # is the 1st on the phe table
  rowRanges <- makeGRangesFromDataFrame(cc, 
                                        seqnames.field = 'chromosome_name',
                                        start.field = 'start_position',
                                        end.field = 'end_position')
  mcols(rowRanges) <- data.frame(Gene_Symbol = cc[,2])
  colData <- DataFrame(phe)
  rownames(colData) <- colData[,1]
  colData[,1] <- NULL
  rse <- SummarizedExperiment::SummarizedExperiment(assays=SimpleList(counts=counts),
                              rowRanges=rowRanges, colData=colData)
  rse
}

