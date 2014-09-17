/*
 * Copyright (c) 2014, Stefan Boehringer <email>
 * All rights reserved.
 *
 * Thu Jul 10 16:59:54 CEST 2014
 */

#include	"HaplotypeRwrappers.hpp"


IntegerVector	R_DiplotypeReconstructionSNPunordered::drawFromLogHfs(const hfsv_t &lhfs, const random_t lu)
	const {
	
	haplotypes_t	r(0, 0);
	DiplotypeReconstructionSNPunordered::drawFromLogHfs(hfs_t(lhfs), lu, r);
	return IntegerVector(r.begin(), r.end());
}

IntegerVector	R_DiplotypeReconstructionSNPunordered::drawFromHfs(const hfsv_t &hfs, const random_t u)
	const {

	haplotypes_t	r(0, 0);
	DiplotypeReconstructionSNPunordered::drawFromHfs(hfs_t(hfs), u, r);
	return IntegerVector(r.begin(), r.end());
}



R_DiplotypeReconstructionSNPunordered::R_DiplotypeReconstructionSNPunordered()
	:
	DiplotypeReconstructionSNPunordered(*new Pedigree)	// <!> memory leak
	{}

R_DiplotypeReconstructionSNPunordered::R_DiplotypeReconstructionSNPunordered(Pedigree &_pedigree,
	int _bitsFactor, int _bitsHt, iid_t NreconstrBuffer)
	: DiplotypeReconstructionSNPunordered(_pedigree, _bitsFactor, _bitsHt, NreconstrBuffer)
	{}

R_DiplotypeReconstructionSNPunordered::~R_DiplotypeReconstructionSNPunordered() {}
