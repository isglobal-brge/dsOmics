#' GDS file resource resolver
#'
#' Build a GDS resource client from a resource object describing access to a 
#' GDS file or a VCF file to be converted into a GDS file.
#'
#' @docType class
#' @format A R6 object of class GDSFileResourceResolver
#' @import resourcer
#' @export
GA4GHResourceResolver <- R6::R6Class(
  "GA4GHResourceResolver",
  inherit = ResourceResolver,
  public = list(
    isFor = function(x) {
      if (super$isFor(x)) {
        !is.null(findFileResourceGetter(x)) && toupper(x$format) %in% c("GA4GH-htsget")
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
