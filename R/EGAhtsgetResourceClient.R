#' EGA htsget resource client
#'
#' Build a GDS handler from a resource object describing access to a 
#' EGA valid credentials, pointing to a VCF or BAM file to be converted into a GDS file.
#'
#' @docType class
#' @format A R6 object of class EGAhtsgetResourceClient
#' @import resourcer
#' @export
EGAhtsgetResourceClient <- R6::R6Class(
  "EGAhtsgetResourceClient",
  inherit = FileResourceClient,
  public = list(
    initialize = function(resource) {
      super$initialize(resource)
    },
    getConnection = function() {
      conn <- super$getConnection()
      if (is.null(conn)) {
        resource <- super$getResource()
        format <- resource$format
        if ("EGAhtsgetVCF" == format) {
          private$loadSNPRelate()
          private$loadRhtsget()
          private$loadGenomicRanges()
          private$loadVariantAnnotation()
          private$.vcf.file.tmp <- tempfile(fileext = ".vcf")
          private$.gds.file.tmp <- tempfile(fileext = ".gds")
          url <- super$parseURL()
          token <- htsget_get_token("https://ega.ebi.ac.uk:8443/ega-openid-connect-server/token",
                                    resource$identity, resource$secret)
          method <- "biallelic.only"
          snpfirstdim <- FALSE
          private$.gr <- GenomicRanges::GRanges(paste0(url$query$referenceName, ":", url$query$start, "-", url$query$end))
          private$.sample_id <- substr(url$path, start = 24, stop = nchar(url$path))
          Rhtsget::htsget_variants(private$.gr, private$.sample_id, "https://ega.ebi.ac.uk:8052/elixir/tickets/tickets", 
                                   token, destination = private$.vcf.file.tmp)
          # VariantAnnotation::readVcf(private$.vcf.file.tmp)
          SNPRelate::snpgdsVCF2GDS(private$.vcf.file.tmp, private$.gds.file.tmp)
          SNPRelate::snpgdsSummary(private$.gds.file.tmp)
        }
        else if("EGAhtsgetBAM" == toupper(format)){
          private$loadgmapR()
          private$loadVariantTools()
          private$loadRhtsget()
          private$.vcf.file.tmp <- tempfile(fileext = ".vcf")
          private$.bam.file.tmp <- tempfile(fileext = ".bam")
          private$.gds.file.tmp <- tempfile(fileext = ".gds")
          token <- htsget_get_token("https://ega.ebi.ac.uk:8443/ega-openid-connect-server/token",
                                    resource$identity, resource$secret)
          private$.gr <- GenomicRanges::GRanges(paste0(url$query$referenceName, ":", url$query$start, "-", url$query$end))
          private$.sample_id <- substr(url$path, start = 24, stop = nchar(url$path))
          bam <- htsget_reads(private$.gr, private$.sample_id, "https://ega.ebi.ac.uk:8052/elixir/tickets/tickets", 
                              ega = TRUE, token, destination = private$.bam.file.tmp)
          sortedFile <- Rsamtools::sortBam(file = bam, destination = paste0(gsub('.{3}$', '', bam), "sorted"), byQname = FALSE)
          indexedFile <- Rsamtools::indexBam(sortedFile)
          bam <- Rsamtools::BamFile(sortedFile)

          BAM2VCF(bam, gr, private$.vcf.file.tmp)
          method <- "biallelic.only"
          snpfirstdim <- FALSE
          SNPRelate::snpgdsVCF2GDS(private$.vcf.file.tmp, private$.gds.file.tmp, method = method, snpfirstdim = snpfirstdim)
          SNPRelate::snpgdsSummary(private$.gds.file.tmp)
        }
        else {
          NULL
        }
        # private$loadGWASTools()
        # conn <- GWASTools::GdsGenotypeReader(private$.gds.file.tmp)
        # super$setConnection(conn)
      }
      conn
    },
    getValue = function() {
      self$getConnection()
    },
    close = function() {
      # super$close()
      # conn <- super$getConnection()
      # if (!is.null(conn)) {
      #   GWASTools::close(conn)
      #   if (!is.null(private$.gds.file.tmp) && file.exists(private$.gds.file.tmp)) {
      #     # house keeping
      #     ignore <- tryCatch(file.remove(private$.gds.file.tmp), error = function(e) {})
      #   }
      # }
    }
  ),
  private = list(
    .gds.file.tmp = NULL,
    .vcf.file.tmp = NULL,
    .bam.file.tmp = NULL,
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
        devtools::install_github("ESCRI11/Rhtsget")
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
    loadgmapR = function() {
      if (!require("gmapR")) {
        if (!require("BiocManager")) {
          install.packages("BiocManager", repos = "https://cloud.r-project.org", dependencies = TRUE)
        }
        BiocManager::install("gmapR", ask = FALSE)
      }
    },
    loadVariantTools = function() {
      if (!require("VariantTools")) {
        if (!require("BiocManager")) {
          install.packages("BiocManager", repos = "https://cloud.r-project.org", dependencies = TRUE)
        }
        BiocManager::install("VariantTools", ask = FALSE)
      }
    }
  )
)