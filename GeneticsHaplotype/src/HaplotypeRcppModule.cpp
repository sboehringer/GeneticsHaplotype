/*
 * HaplotypeRcppModule.cpp
 * Fri Sep  5 15:09:16 CEST 2014
 * 
 */

#include	<Rcpp.h>
#include	"HaplotypeRwrappers.hpp"

RCPP_MODULE(ReconstructorModule) {
	using namespace Rcpp;

	// class_<c++ class name>('R class name')
	class_<Reconstructor>("DiplotypeReconstructor")
		.constructor<SEXP, SEXP>("genotype matrix, list of pedigrees.")

		.method("countMarkers", &Reconstructor::countMarkers, "Number of markers in the reconstruction.")
		.method("drawFromHfs", &Reconstructor::drawFromHfs, "Draw diplotypes for given haplotype frequency distribution for a set of pedigrees.")
		.method("drawFamFromHfs", &Reconstructor::drawFamFromHfs, "Draw diplotypes for given haplotype frequency distribution for a single family from a pedigrees.")
		.method("reconstructionsFam", &Reconstructor::reconstructionsFam, "Return reconstructions for family i. First column is multiplicative constant, second column is inheritance vector, ensuing columns contain pairs of haplotypes.")
		.method("reconstructionsFam", &Reconstructor::reconstructionsFam, "Return reconstructions for family i. First column is multiplicative constant, second column is inheritance vector, ensuing columns contain pairs of haplotypes.")
		.method("reconstructionsFamFull", &Reconstructor::reconstructionsFamFull, "Return reconstructions for family i. First column is multiplicative constant, ensuing columns contain pairs of haplotypes for all family members.")
		.method("reconstructionsAll", &Reconstructor::reconstructionsAll, "Uses reconstructionsFamFull to create a full list of reconstructions for all families.")
	;
}