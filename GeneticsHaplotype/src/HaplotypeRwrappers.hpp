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
#include "genotype.h"
#include <memory>
using namespace std;

#if 0
class R_Pedigree : Pedigree {
public:
	R_Pedigree(List &pedigree);
	~R_Pedigree();
};
#endif


class R_GenotypeFetcher : public GenotypeFetcher
{
	IntegerMatrix m;
public:
	R_GenotypeFetcher(IntegerMatrix _m) : m(_m) {}
	R_GenotypeFetcher(SEXP _m) : m(_m) {
		m = Rcpp::IntegerMatrix(_m);
	}
	~R_GenotypeFetcher() {}

	virtual	marker_t	countMarkers(void) const {
		return m.ncol();
	}
	virtual	genotype_t	genotype(iid_t id, marker_t marker) const {
		return (genotype_t)m(id, marker);
	}
	virtual iid_t		N(void) const {
		return m.nrow();
	}
};

class R_DiplotypeReconstructionSNPunordered : public DiplotypeReconstructionSNPunordered
{
public:
	R_DiplotypeReconstructionSNPunordered(Pedigree &_pedigree,
										int _bitsFactor = 6, int _bitsHt = 10, iid_t NreconstrBuffer = 1024);
	R_DiplotypeReconstructionSNPunordered();
	~R_DiplotypeReconstructionSNPunordered();
	R_DiplotypeReconstructionSNPunordered(const R_DiplotypeReconstructionSNPunordered&) = delete;
	R_DiplotypeReconstructionSNPunordered(R_DiplotypeReconstructionSNPunordered &&other)
	: DiplotypeReconstructionSNPunordered(std::move(other)) {
	}
	
	IntegerVector 	drawFromLogHfs(const hfsv_t &lhfs, const random_t lu) const;
	haplotypes_t	drawFromHfs(const hfs_t &hfs, const random_t u) const;
};

class PedigreeCollection : public vector< Pedigree > {
	
public:
	PedigreeCollection(List &pedvector) : vector< Pedigree >(0) {
		this->addPedigrees(pedvector);
	}
	PedigreeCollection(SEXP pedvector) : vector< Pedigree >(0) {
		Rcpp::List	pedlist(pedvector);
		this->addPedigrees(pedlist);
	}

	~PedigreeCollection() {}
	
	void	print(void);

	void	addPedigrees(List &pedvector) {
		reserve(pedvector.size());
		//cout << "mycapcity: " << capacity() << " pedvector size:" << pedvector.size() << endl;
		for (int i = 0; i < pedvector.size(); i++) {
			const List		&rped(Rcpp::as<List>(pedvector[i]));
			//cout << "pedigree: " << i << endl;
			vector<iid_t>	founders(Rcpp::as< vector<iid_t> >(rped["founders"]));
			IntegerMatrix	itrios(Rcpp::as<IntegerMatrix>(rped["itrios"]));
			//cout << "got list elements" << endl;
			const Pedigree	ped(founders, itrios);
			//ped.print();
			push_back(ped);
		}
	}
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
	vector<R_DiplotypeReconstructionSNPunordered>
								reconstructions;

public:
	Reconstructor(R_GenotypeFetcher &_fetcher, List &_peds)
	: fetcher(_fetcher), peds(_peds), reconstructions() {
	}

	Reconstructor(IntegerMatrix &_m, List &_peds)
	: fetcher(* new R_GenotypeFetcher(_m)), peds(_peds), reconstructions() {
		this->reconstruct();
	}
	Reconstructor(SEXP _m, SEXP _peds)
	: fetcher(* new R_GenotypeFetcher(_m)), peds(_peds), reconstructions() {
		this->reconstruct();
	}
	void	reconstruct(void) {
		//cout << "#peds: " << peds.size() << endl;
		//fetcher.print();
		iid_t	Ncum = 0;
		for (iid_t i = 0; i < peds.size(); i++, Ncum += peds[i].N()) {
			//cout << "Ncum: " << Ncum << endl;
			R_DiplotypeReconstructionSNPunordered reconstruction(peds[i]);
			GenotypeFetcherOffset	pedfetcher((GenotypeFetcher &)fetcher, Ncum);
			//peds[i].print();
			reconstruction.reconstruct((GenotypeFetcher &)pedfetcher);
			//reconstruction.print();
			reconstructions.push_back(std::move(reconstruction));
		}
		//cout << "did reconstruct" << endl;
	}

	IntegerMatrix	drawFromHfs(const NumericVector &hfsR, const NumericVector &u) const {
		IntegerMatrix	m(fetcher.N(), 2);
		iid_t			Ni = 0;
		hfs_t			hfs(Rcpp::as< hfsv_t >(hfsR));
		
		for (iid_t i = 0; i <  peds.size(); Ni += peds[i].N(), i++) {
			haplotypes_t	draw;
			reconstructions[i].DiplotypeReconstructionSNPunordered::drawFromHfs(hfs, u[i], draw);

			//cout << "draw: ";
			//for (iid_t j = 0; j < draw.size(); j++) cout << (j? ", ": "") << draw[j];
			//cout << endl;
			// <A> not yet brought in correct order
			for (iid_t j = 0; j < draw.size(); j ++) m(Ni + j/2, j%2) = draw[j];
		}
		return wrap(m);
	}

};

#endif //HAPLOTYPESRWRAPPERS_H
