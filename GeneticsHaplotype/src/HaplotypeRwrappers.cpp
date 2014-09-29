/*
 * Copyright (c) 2014, Stefan Boehringer <email>
 * All rights reserved.
 *
 * Thu Jul 10 16:59:54 CEST 2014
 */

#include	"HaplotypeRwrappers.hpp"

#if 0
R_Pedigree::R_Pedigree(List pedigree)
: Pedigree(Rcpp::as< vector<int> >(rped["founders"]),
		   Rcpp::as< vector< vector<int> > >(rped["itrios"])) {
	
	
	//for (int i = 0; i < _founder.size(); i++) founder[i] = _founder[i];
	//founder = vectorConvert<int, iid_t>(_founder(_founder.begin(), _founder.end()));
	// 	itrio.resize(_itrio.nrow());
	// 	for (int i = 0; i < _itrio.nrow(); i++) {
	// 		itrio[i] = vector<iid_t>(3);
	// 		//itrio[i].swap(vector<iid_t>(3));	// does not work
	// 		for (int j = 0; j < _itrio.ncol(); j++)
	// 			itrio[i][j] = itrio[i][j];
	// 
	// 	}
}
#endif

IntegerVector	R_DiplotypeReconstructionSNPunordered::drawFromLogHfs(const hfsv_t &lhfs, const random_t lu)
	const {
	
	haplotypes_t	r(pedigree.N(), 0);
	DiplotypeReconstructionSNPunordered::drawFromLogHfs(hfs_t(lhfs), lu, r);
	return IntegerVector(r.begin(), r.end());
}

haplotypes_t	R_DiplotypeReconstructionSNPunordered::drawFromHfs(const hfs_t &hfs, const random_t u)
	const {

	haplotypes_t	r(pedigree.N(), 0);
	DiplotypeReconstructionSNPunordered::drawFromHfs(hfs, u, r);
	return r;
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
