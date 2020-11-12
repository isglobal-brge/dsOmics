#' @title Variant discovery on a BAM file
#' 
#' @description Discover variants of a BAM file and write them in a VCF file
#'
#' @param bam \code{character} Path to the BAM file
#' @param grange \code{GRange} GRange of the BAM file to search for variants
#' @param destination \code{character} Path to the VCF file
#'
#' @export

BAM2VCF <- function(bam, grange, destination){
  seqlevelsStyle(seqlevels(grange)) <- "UCSC"
  gen <- Biostrings::DNAStringSet(
    BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19[[seqlevels(grange)]]
    )
  names(gen) <- as.character(GenomeInfoDb::seqnames(grange))
  refgenome <- gmapR::GmapGenome(
                        genome = gen, 
                        name = as.character(GenomeInfoDb::seqnames(grange)), 
                        create = T
                      )
  tally.param <- VariantTools::TallyVariantsParam(
                                    refgenome, high_base_quality = 23L,
                                    indels = TRUE)
  BPPARAM <- MulticoreParam(workers = detectCores())
  raw.variants <- tallyVariants(bam, param = tally.param, BPPARAM = BPPARAM)
  qa.variants <- qaVariants(raw.variants)
  called.variants <- callVariants(qa.variants)
  sampleNames(called.variants) <- "sample"
  mcols(called.variants) <- NULL
  vcf <- asVCF(called.variants)
  writeVcf(vcf, destination)
  
}