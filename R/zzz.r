.onAttach <- function(libname, pkgname) {
  # packageStartupMessage(
  #   sprintf("\nThis is %s version %s\n", pkgname, packageVersion(pkgname)))
}

.onLoad <- function(libname, pkgname) {
  .check_msdialr_options()
}

.onUnload <- function(libpath) {
  # library.dynam.unload("xcms", libpath)
}
