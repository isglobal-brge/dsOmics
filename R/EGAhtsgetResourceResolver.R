#' GDS file resource resolver
#'
#' Build a GDS resource client from a resource object describing access to a 
#' GDS file or a VCF file to be converted into a GDS file.
#'
#' @docType class
#' @format A R6 object of class GDSFileResourceResolver
#' @import resourcer
#' @export
EGAhtsgetResourceResolver <- R6::R6Class(
  "EGAhtsgetResourceResolver",
  inherit = ResourceResolver,
  public = list(
    isFor = function(x) {
      if (super$isFor(x)) {
        !is.null(findFileResourceGetter(x)) && x$format %in% c("EGAhtsgetBAM", "EGAhtsgetVCF")
      } else {
        FALSE
      }
    },
    newClient = function(x) {
      if (self$isFor(x)) {
        EGAhtsgetResourceClient$new(x)
      } else {
        NULL
      }
    }
  )
)
