/*
 * HaplotypeRcppModule.cpp
 * Fri Sep  5 15:09:16 CEST 2014
 * 
 */

#include	<Rcpp.h>
#include	"HaplotypeRwrappers.hpp"

RCPP_MODULE(Reconstructor) {
	using namespace Rcpp;

	class_<Reconstructor>("DiplotypeReconstructor")
	.constructor<SEXP, SEXP>("genotype matrix, list of pedigrees.")

	.method("drawFromHfs", &Reconstructor::drawFromHfs, "Draw diplotypes for given haplotype frequency distribution for a set of pedigrees.")
	;
}