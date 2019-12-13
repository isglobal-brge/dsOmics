#' GDS file resource resolver
#'
#' Build a GDS resource client from a resource object describing access to a 
#' GDS file or a VCF file to be converted into a GDS file.
#'
#' @docType class
#' @format A R6 object of class GDSFileResourceResolver
#' @import resourcer
#' @export
GDSFileResourceResolver <- R6::R6Class(
  "GDSFileResourceResolver",
  inherit = ResourceResolver,
  public = list(
    isFor = function(x) {
      if (super$isFor(x)) {
        !is.null(findFileResourceGetter(x)) && toupper(x$format) %in% c("VCF2GDS", "GDS")
      } else {
        FALSE
      }
    },
    newClient = function(x) {
      if (self$isFor(x)) {
        GDSFileResourceClient$new(x)
      } else {
        NULL
      }
    }
  )
)
