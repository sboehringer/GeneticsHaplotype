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
	List	&peds;
	
public:
	PedigreeCollection(List &pedvector) : peds(pedvector) {}
	~PedigreeCollection() {}
	
	void	print(void);
	
	//DiplotypeReconstruction		*reconstruct(void);	// ids of all individuals in the pedigrees
};

class GenotypeData {
	const R_GenotypeFetcher		&fetcher;
	const List					&peds;
public:
	GenotypeData(const R_GenotypeFetcher &_fetcher, const List &_peds) : fetcher(_fetcher), peds(_peds) {}
};

class Reconstructor {
	const R_GenotypeFetcher		&fetcher;
	const List					&peds;
	iid_t						N;
	vector<R_DiplotypeReconstructionSNPunordered>	reconstructions;

public:
	Reconstructor(const R_GenotypeFetcher &_fetcher, const List &_peds)
	: fetcher(_fetcher), peds(_peds), reconstructions() {}
	Reconstructor(const IntegerMatrix &_m, const List &_peds)
	: GenotypeData(R_GenotypeFetcher(_m), _peds), reconstructions(), N(0) {
		for (iid_t i = 0; i < peds.size(); i++) {
			reconstructions.push(R_DiplotypeReconstructionSNPunordered(peds[i]));
			N += peds[i].size();
		}
	}

	IntegerMatrix	drawFromHfs(const hfs_t &hfs, const vector<random_t> u) const {
		IntegerMatrix	m(N, 2);
		iid_t			Ni = 0;

		for (iid_t i = 0; i <  peds.size(); Ni += peds[i].N(), i++) {
			haplotypes_t	draw;
			reconstructions[i].drawFromHfs(hfs, u[i], draw);

			// <A> not yet brought in correct order
			for (iid_t j = 0; j < draw.size(); j += 2) {
				m(j/2, 0) = draw[j];
				m(j/2, 1) = draw[j + 1];
			}
		}
		return wrap(m);
	}

};

#endif //HAPLOTYPESRWRAPPERS_H
