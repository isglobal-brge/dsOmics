#' GA4GH htsget server resource resolver
#'
#' Build a GA4GH resource client from a resource object describing access to a 
#' GA4GH htsget enabled server, pointing to a VCF or BAM file to be converted into a GDS file.
#'
#' @docType class
#' @format A R6 object of class GA4GHResourceResolver
#' @import resourcer
#' @export
GA4GHResourceResolver <- R6::R6Class(
  "GA4GHResourceResolver",
  inherit = ResourceResolver,
  public = list(
    isFor = function(x) {
      if (super$isFor(x)) {
        !is.null(findFileResourceGetter(x)) && toupper(x$format) %in% c("GA4GHVCF", "GA4GHBAM")
      } else {
        FALSE
      }
    },
    newClient = function(x) {
      if (self$isFor(x)) {
        GA4GHResourceClient$new(x)
      } else {
        NULL
      }
    }
  )
)
