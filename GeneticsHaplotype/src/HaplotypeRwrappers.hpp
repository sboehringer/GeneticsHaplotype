/*
 * Copyright (c) 2014, Stefan Boehringer <email>
 * All rights reserved.
 *
 * Thu Jul 10 16:59:54 CEST 2014
 */

#ifndef HAPLOTYPESRWRAPPERS_H
#define HAPLOTYPESRWRAPPERS_H

#include "diplotypereconstruction.h"

class R_DiplotypeReconstructionSNPunordered : public DiplotypeReconstructionSNPunordered
{
public:
	R_DiplotypeReconstructionSNPunordered(Pedigree &_pedigree,
										int _bitsFactor = 6, int _bitsHt = 10, iid_t NreconstrBuffer = 1024);
	R_DiplotypeReconstructionSNPunordered();
	~R_DiplotypeReconstructionSNPunordered();
	
	IntegerVector 	drawFromLogHfs(const hfsv_t &lhfs, const random_t lu) const;
	IntegerVector	drawFromHfs(const hfsv_t &hfs, const random_t u) const;
};

class R_GenotypeFetcher : public GenotypeFetcher
{
	const IntegerMatrix &m;
public:
	R_GenotypeFetcher(const IntegerMatrix &_m) : m(_m) {}
	~R_GenotypeFetcher() {}

	virtual	marker_t	countMarkers(void) const {
		return m.ncol();
	}
	virtual	genotype_t	genotype(iid_t id, marker_t marker) const {
		return (genotype_t)m(id, marker);
	}
};

class PedigreeCollection {
	vector<Pedigree>	peds;
	
public:
	PedigreeCollection(List &pedvector);
	~PedigreeCollection() {}
	
	void	print(void);
	
	//DiplotypeReconstruction		*reconstruct(void);	// ids of all individuals in the pedigrees
};

class GenotypeData {
	const R_GenotypeFetcher		fetcher;
	PedigreeCollection			peds;
public:
	GenotypeData(const R_GenotypeFetcher _fetcher, List &_peds) : fetcher(_fetcher), peds(_peds) {}
};

class Reconstructor : GenotypeData {

	vector<R_DiplotypeReconstructionSNPunordered>	reconstructions;
public:
	Reconstructor(const R_GenotypeFetcher _fetcher, List &_peds) : GenotypeData(_fetcher, _peds), reconstructions() {}
};


#endif //HAPLOTYPESRWRAPPERS_H
