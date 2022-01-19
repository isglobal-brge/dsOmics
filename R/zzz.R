.onAttach <- function(libname, pkgname) {
  resourcer::registerResourceResolver(GDSFileResourceResolver$new())
  resourcer::registerResourceResolver(GA4GHResourceResolver$new())
  resourcer::registerResourceResolver(EGAhtsgetResourceResolver$new())
  options(pillar.sigfig = 3)
  resourcer::registerFileResourceGetter(Opal2FileResourceGetter$new())
  # resourcer::registerResourceResolver(EGAmetadataResourceResolver$new())
}

.onDetach <- function(libpath) {
  resourcer::unregisterResourceResolver("GDSFileResourceResolver")
  resourcer::unregisterResourceResolver("GA4GHResourceResolver")
  resourcer::unregisterResourceResolver("EGAhtsgetResourceResolver")
  # resourcer::unregisterResourceResolver("EGAmetadataResourceResolver")
}