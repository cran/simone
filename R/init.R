.onLoad   <- function(lib, pkg) {
  # require(mixer, quietly=TRUE)
}
.onAttach <- function(libname, pkgname){
    packageStartupMessage("")
    packageStartupMessage("----------------------------------------------------------------------")
    packageStartupMessage("")
    packageStartupMessage("      'simone' package version 1.0-0")
    packageStartupMessage("      SIMoNe page (http://stat.genopole.cnrs.fr/software/simone)")
    packageStartupMessage("")
    packageStartupMessage("----------------------------------------------------------------------")
    packageStartupMessage("Note that versions >= 1.0-0 are not compatible with versions < 1.0.0. ")
    packageStartupMessage("----------------------------------------------------------------------")

}
