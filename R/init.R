.onLoad   <- function(lib, pkg) {
  require(mixer, quietly=TRUE)
}
.onAttach <- function(libname, pkgname){
  cat(""                                                                      ,
      "----------------------------------------------------------------------",
      ""                                                                      ,
      "      'simone' package version 1.0-0"                                  ,
      "      SIMoNe page (http://stat.genopole.cnrs.fr/software/simone)"      ,
      ""                                                                      ,
      "----------------------------------------------------------------------",
      "Note that versions >= 1.0-0 are not compatible with versions < 1.0.0. ",
      "----------------------------------------------------------------------",
      sep = "\n")
}
