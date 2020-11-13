.onAttach <- function(libname, pkgname) {
  resourcer::registerResourceResolver(GDSFileResourceResolver$new())
  resourcer::registerResourceResolver(GA4GHResourceResolver$new())
  resourcer::registerResourceResolver(EGAhtsgetResourceResolver$new())
}

.onDetach <- function(libpath) {
  resourcer::unregisterResourceResolver("GDSFileResourceResolver")
  resourcer::unregisterResourceResolver("GA4GHResourceResolver")
  resourcer::unregisterResourceResolver("EGAhtsgetResourceResolver")
}