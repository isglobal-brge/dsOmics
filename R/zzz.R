.onAttach <- function(libname, pkgname) {
  resourcer::registerResourceResolver(GDSFileResourceResolver$new())
  resourcer::registerResourceResolver(GA4GHResourceResolver$new())
}

.onDetach <- function(libpath) {
  resourcer::unregisterResourceResolver("GDSFileResourceResolver")
  resourcer::unregisterResourceResolver("GA4GHResourceResolver")
}