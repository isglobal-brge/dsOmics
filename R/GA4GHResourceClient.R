#' GDS file resource client
#'
#' Build a GDS handler from a resource object describing access to a 
#' GDS file or a VCF file to be converted into a GDS file.
#'
#' @docType class
#' @format A R6 object of class GDSFileResourceClient
#' @import resourcer
#' @export
GA4GHResourceClient <- R6::R6Class(
  "GA4GHResourceClient",
  inherit = FileResourceClient,
  public = list(
    initialize = function(resource) {
      super$initialize(resource)
    },
    getConnection = function() {
      conn <- super$getConnection()
      if (is.null(conn)) {
        path <- super$downloadFile()
        resource <- super$getResource()
        format <- resource$format
        if ("GA4GH" == toupper(format)) {
          private$loadSNPRelate()
          private$loadRhtsget()
          private$loadGenomicRanges()
          private$loadVariantAnnotation()
          private$.vcf.file.tmp <- tempfile(fileext = ".vcf")
          private$.gds.file.tmp <- tempfile(fileext = ".gds")
          url <- super$parseURL()
          method <- "biallelic.only"
          snpfirstdim <- FALSE
          private$.gr <- GenomicRanges::GRanges(paste0(url$query$referenceName, ":", url$query$start, "-", url$query$end))
          private$.sample_id <- substr(url$path, start = 10, stop = nchar(url$path))

          Rhtsget::htsget_variants(private$.gr, private$.sample_id, url$hostname, destination = private$.vcf.file.tmp)
          # VariantAnnotation::readVcf(private$.vcf.file.tmp)
          SNPRelate::snpgdsVCF2GDS(private$.vcf.file.tmp, private$.gds.file.tmp)
          SNPRelate::snpgdsSummary(private$.gds.file.tmp)
          
        } else {
          NULL
        }
        private$loadGWASTools()
        conn <- GWASTools::GdsGenotypeReader(private$.gds.file.tmp)
        super$setConnection(conn)
      }
      conn
    },
    getValue = function() {
      self$getConnection()
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
    .vcf.file.tmp = NULL,
    .gr = NULL,
    .sample_id = NULL,
    
    loadSNPRelate = function() {
      if (!require("SNPRelate")) {
        install.packages("SNPRelate", repos="http://R-Forge.R-project.org", dependencies = TRUE)
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
    loadVariantAnnotation = function() {
      if (!require("VariantAnnotation")) {
        if (!require("BiocManager")) {
          install.packages("BiocManager", repos = "https://cloud.r-project.org", dependencies = TRUE)
        }
        BiocManager::install("VariantAnnotation", ask = FALSE)
      }
    },
    loadRhtsget = function() {
      if (!require("Rhtsget")) {
        if (!require("BiocManager")) {
          install.packages("BiocManager", repos = "https://cloud.r-project.org", dependencies = TRUE)
        }
        BiocManager::install("Rhtsget", ask = FALSE)
      }
    },
    loadGWASTools = function() {
      if (!require("GWASTools")) {
        if (!require("BiocManager")) {
          install.packages("BiocManager", repos = "https://cloud.r-project.org", dependencies = TRUE)
        }
        BiocManager::install("GWASTools", ask = FALSE)
      }
    }
  )
)