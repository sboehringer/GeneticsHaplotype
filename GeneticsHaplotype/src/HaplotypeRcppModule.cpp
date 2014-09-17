/*
 * HaplotypeRcppModule.cpp
 * Fri Sep  5 15:09:16 CEST 2014
 * 
 */

#include	<Rcpp.h>
#include	"HaplotypeRwrappers.hpp"

RCPP_MODULE(Haplotype){
	using namespace Rcpp ;

	class_<R_DiplotypeReconstructionSNPunordered>("DiplotypeReconstructionSNPunordered")
	// expose the default constructor
	.constructor()
	
	.method("drawFromLogHfs", &R_DiplotypeReconstructionSNPunordered::drawFromLogHfs,
		"Draw diplotypes for given haplotype frequency distribution, with frequencies speicified on the log scale.")
	.method("drawFromHfs", &R_DiplotypeReconstructionSNPunordered::drawFromHfs,
		"Draw diplotypes for given haplotype frequency distribution.")
	;
}


