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
    getConnection = function() {
      conn <- super$getConnection()
      if (is.null(conn)) {
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
          snpgdsVCF2GDS(vcf.fn = path, out.fn = private$.gds.file.tmp, method = method, snpfirstdim = snpfirstdim)
          path <- private$.gds.file.tmp
          snpgdsSummary(path)
        } else {
          private$loadGDSFmt()
        }
        private$loadGWASTools()
        conn <- GWASTools::GdsGenotypeReader(path)
        og_snpid <- conn@snpIDvar
        conn@snpIDvar <- "snp.rs.id"
        if(is.null(getSnpID(conn))){conn@snpIDvar <- og_snpid}
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
    }
  )
)