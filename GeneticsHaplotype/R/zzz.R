.onLoad <- function(libname, pkgname){
	#loadRcppModules()
    loadModule("ReconstructorModule", TRUE)
}
