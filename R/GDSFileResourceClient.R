#' GDS file resource client
#'
#' Build a GDS handler from a resource object describing access to a 
#' GDS file or a VCF file to be converted into a GDS file.
#'
#' @docType class
#' @format A R6 object of class GDSFileResourceClient
#' @import resourcer
#' @export
GDSFileResourceClient <- R6::R6Class(
  "GDSFileResourceClient",
  inherit = FileResourceClient,
  public = list(
    initialize = function(resource) {
      super$initialize(resource)
    },
    getConnection = function(...) {
      conn <- super$getConnection()
      # if (is.null(conn)) {
        path <- super$downloadFile()
        resource <- super$getResource()
        format <- resource$format
        if ("VCF2GDS" == toupper(format)) {
          private$loadSNPRelate() # also loads gdsfmt
          private$.gds.file.tmp <- tempfile(fileext = ".gds")
          url <- super$parseURL()
          method <- "biallelic.only"
          snpfirstdim <- FALSE
          if (!is.null(url$query)) {
            if (!is.null(url$query$method)) {
              method <- url$query$method
            }
            if (!is.null(url$query$snpfirstdim)) {
              snpfirstdim <- as.logical(url$query$snpfirstdim)
            }
          }
          if(!is.null(list(...))){
            # Load necessary packages
            private$loadGenomicRanges()
            private$loadIRanges()
            private$loadVariantAnnotation()
            # Get genomic ranges to subset VCF
            extra_args <- list(...)
            if(all(c("seq", "start.loc", "end.loc") %in% names(extra_args))){
              # Try catch, some files even having gz extension appear to not be compressed
              index_path <- tryCatch({
                if(!file.exists(paste0(path, ".tbi"))){
                  indexTabix(path, format="vcf")
                } else {paste0(path, ".tbi")}
              }, error = function(w){
                path <- bgzip(path, overwrite=TRUE)
                indexTabix(path, format="vcf")
              })
              # Check which chromosomes are present on the VCF, this will be used to subset the extra_args$seq
              # object, if a chromosome is not present on the VCF the readVcf function crashes
              # Get chromosome names from header
              chromosome_names <- as.list(scanVcfHeader(path)@reference)
              valid_chromosomes <- lapply(chromosome_names, function(x){
                tryCatch({
                  scanTabix(path, param = GenomicRanges::GRanges(x, IRanges::IRanges(1, 1)))
                  return(x)
                }, error = function(w){})
              })
              valid_chromosomes <- unlist(valid_chromosomes)
              # Create the g_range object using only the present chromosomes
              valid_indexes <- extra_args$seq %in% valid_chromosomes
              g_range <- GenomicRanges::GRanges(extra_args$seq[valid_indexes], 
                                                IRanges::IRanges(extra_args$start.loc[valid_indexes], 
                                                                 extra_args$end.loc[valid_indexes]))
              params <- VariantAnnotation::ScanVcfParam(which = g_range)
              # Read subsetted VCF and write it to a new file, compress it and 
              # re-write the path variable
              vcf <- readVcf(index_path, "hg19", params)
              writeVcf(vcf, sprintf("%s_2", sub("\\$", "", path)))
              path <- bgzip(sprintf("%s_2", sub("\\$", "", path)), overwrite=TRUE)
            }else if(all("snps" %in% names(extra_args))){
              # Try catch, some files even having gz extension appear to not be compressed
              index_path <- tryCatch({
                if(!file.exists(paste0(path, ".tbi"))){
                  indexTabix(path, format="vcf")
                } else {paste0(path, ".tbi")}
              }, error = function(w){
                path <- bgzip(path, overwrite=TRUE)
                indexTabix(path, format="vcf")
              })
              isSNP <- function(x) {
                grepl(paste(extra_args$snps, collapse = "|"), x)
                }
              prefilters <- FilterRules(list(snp=isSNP))
              filtered_vcf <- filterVcf(path, "hg19", destination = tempfile(),
                        prefilters=prefilters, verbose=TRUE)
              path <- bgzip(filtered_vcf, overwrite=TRUE)
            }else{
              stop()
            }
          }
          snpgdsVCF2GDS(vcf.fn = path, out.fn = private$.gds.file.tmp, method = method, snpfirstdim = snpfirstdim)
          path <- private$.gds.file.tmp
          snpgdsSummary(path)
        } else {
          private$loadGDSFmt()
        }
        private$loadGWASTools()
        conn <- GWASTools::GdsGenotypeReader(path, allow.fork = TRUE)
        super$setConnection(conn)
      # }
      conn
    },
    getValue = function(...) {
      self$getConnection(...)
    },
    close = function() {
      super$close()
      conn <- super$getConnection()
      if (!is.null(conn)) {
        GWASTools::close(conn)
        if (!is.null(private$.gds.file.tmp) && file.exists(private$.gds.file.tmp)) {
          # house keeping
          ignore <- tryCatch(file.remove(private$.gds.file.tmp), error = function(e) {})
        }
      }
    }
  ),
  private = list(
    .gds.file.tmp = NULL,
    loadSNPRelate = function() {
      if (!require("SNPRelate")) {
        install.packages("SNPRelate", repos="http://R-Forge.R-project.org", dependencies = TRUE)
      }
    },
    loadGDSFmt = function() {
      if (!require("gdsfmt")) {
        if (!require("BiocManager")) {
          install.packages("BiocManager", repos = "https://cloud.r-project.org", dependencies = TRUE)
        }
        BiocManager::install("gdsfmt", ask = FALSE)
      }
    },
    loadGWASTools = function() {
      if (!require("GWASTools")) {
        if (!require("BiocManager")) {
          install.packages("BiocManager", repos = "https://cloud.r-project.org", dependencies = TRUE)
        }
        BiocManager::install("GWASTools", ask = FALSE)
      }
    },
    loadGenomicRanges = function() {
      if (!require("GenomicRanges")) {
        if (!require("BiocManager")) {
          install.packages("BiocManager", repos = "https://cloud.r-project.org", dependencies = TRUE)
        }
        BiocManager::install("GenomicRanges", ask = FALSE)
      }
    },
    loadIRanges = function() {
      if (!require("IRanges")) {
        if (!require("BiocManager")) {
          install.packages("BiocManager", repos = "https://cloud.r-project.org", dependencies = TRUE)
        }
        BiocManager::install("IRanges", ask = FALSE)
      }
    },
    loadVariantAnnotation = function() {
      if (!require("VariantAnnotation")) {
        if (!require("BiocManager")) {
          install.packages("BiocManager", repos = "https://cloud.r-project.org", dependencies = TRUE)
        }
        BiocManager::install("VariantAnnotation", ask = FALSE)
      }
    }
  )
)