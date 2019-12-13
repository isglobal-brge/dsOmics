.onAttach <- function(libname, pkgname) {
  resourcer::registerResourceResolver(GDSFileResourceResolver$new())
}

.onDetach <- function(libpath) {
  resourcer::unregisterResourceResolver("GDSFileResourceResolver")
}