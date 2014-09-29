/*
 * Copyright (c) 2014, Stefan Boehringer <email>
 * All rights reserved.
 *
 * Thu Jul 10 16:59:54 CEST 2014
 */

#ifndef HAPLOTYPESRWRAPPERS_H
#define HAPLOTYPESRWRAPPERS_H

#include "diplotypereconstruction.h"
#include "pedigree.h"
#include <memory>
using namespace std;

#if 0
class R_Pedigree : Pedigree {
public:
	R_Pedigree(List &pedigree);
	~R_Pedigree();
};
#endif

class R_DiplotypeReconstructionSNPunordered : public DiplotypeReconstructionSNPunordered
{
public:
	R_DiplotypeReconstructionSNPunordered(Pedigree &_pedigree,
										int _bitsFactor = 6, int _bitsHt = 10, iid_t NreconstrBuffer = 1024);
	R_DiplotypeReconstructionSNPunordered();
	~R_DiplotypeReconstructionSNPunordered();
	
	IntegerVector 	drawFromLogHfs(const hfsv_t &lhfs, const random_t lu) const;
	haplotypes_t	drawFromHfs(const hfs_t &hfs, const random_t u) const;
};

class R_GenotypeFetcher : public GenotypeFetcher
{
	IntegerMatrix &m;
public:
	R_GenotypeFetcher(IntegerMatrix &_m) : m(_m) {}
	~R_GenotypeFetcher() {}

	virtual	marker_t	countMarkers(void) const {
		return m.ncol();
	}
	virtual	genotype_t	genotype(iid_t id, marker_t marker) const {
		return (genotype_t)m(id, marker);
	}
};

class PedigreeCollection : public vector< Pedigree > {
	
public:
	PedigreeCollection(const List &pedvector) : vector< Pedigree >() {
		resize(pedvector.size());
		for (int i = 0; i < pedvector.size(); i++) {
			const List	&rped(Rcpp::as<List>(pedvector[i]));
			const Pedigree	ped(
				Rcpp::as< vector<int> >(rped["founders"]),
				Rcpp::as< vector< vector<int> > >(rped["itrios"])
			);
			push_back(ped);
		}
	}

	~PedigreeCollection() {
		
	}
	
	void	print(void);
	
	//DiplotypeReconstruction		*reconstruct(void);	// ids of all individuals in the pedigrees
};

class GenotypeData {
	R_GenotypeFetcher		&fetcher;
	List					&peds;
public:
	GenotypeData(R_GenotypeFetcher &_fetcher, List &_peds) : fetcher(_fetcher), peds(_peds) {}
};

class Reconstructor {
	R_GenotypeFetcher		&fetcher;
	PedigreeCollection		peds;
	iid_t					N;
	vector<R_DiplotypeReconstructionSNPunordered>
								reconstructions;

public:
	Reconstructor(R_GenotypeFetcher &_fetcher, List &_peds)
	: fetcher(_fetcher), peds(_peds), reconstructions() {
		N = 0;
	}

	Reconstructor(IntegerMatrix &_m, List &_peds)
	: fetcher(* new R_GenotypeFetcher(_m)), peds(_peds), reconstructions() {
		N = 0;
		for (iid_t i = 0; i < peds.size(); i++) {
			R_DiplotypeReconstructionSNPunordered reconstruction(peds[i]);
			reconstructions.push_back(reconstruction);
			N += peds[i].N();
		}
	}

	IntegerMatrix	drawFromHfs(const NumericVector &hfsR, const NumericVector &u) const {
		IntegerMatrix	m(N, 2);
		iid_t			Ni = 0;
		hfs_t			hfs(Rcpp::as< hfsv_t >(hfsR));
		
		for (iid_t i = 0; i <  peds.size(); Ni += peds[i].N(), i++) {
			haplotypes_t	draw(reconstructions[i].drawFromHfs(hfs, u[i]));

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
