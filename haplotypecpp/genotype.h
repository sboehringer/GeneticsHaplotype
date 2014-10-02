/*
 * genotype.h
 * Thu Jul 10 16:20:15 CEST 2014
 */

#ifndef GENOTYPE_H
#define GENOTYPE_H

#include "BitField.hpp"

/*
 * Allele coding: 0/1
 */

typedef int					iid_t;
typedef int					marker_t;
typedef unsigned int		haplotype_t;
typedef vector<haplotype_t>	haplotypes_t;
typedef enum {genotype_0, genotype_1, genotype_2, genotype_NA } genotype_t;
typedef vector<genotype_t>	gt_t;
typedef struct {
	haplotype_t	d1, d2;
} diplotype_t;
// enumerate geneotype combinations
typedef unsigned int	genotypecomb_t;

inline genotypecomb_t	genotypeCombinationFromDiplotype(
	diplotype_t dt, marker_t countMarker, haplotype_t missing = 0) {

	genotypecomb_t	r = 0;
	BitArray<genotypecomb_t, marker_t>	b(&r, 0, 2, countMarker);

	dt.d1 &= ~missing;
	dt.d2 &= ~missing;
	for (marker_t i = 0; i < countMarker; i++) {
		b.set(i, bitAt(dt.d1, i) + bitAt(dt.d2, i));
	}
	return r;
}

class GenotypeFetcher {

	public:
	GenotypeFetcher(void){}
	virtual ~GenotypeFetcher(){}

	virtual	marker_t	countMarkers(void) const {
		throw("Abstract countMarkers fetching method called.");
	}
	virtual	genotype_t	genotype(iid_t id, marker_t marker) const {
		throw("Abstract genotype fetching method called.");
	}
	virtual	int			countMissing(iid_t id) const {
		int	count = 0;
		for (int i = 0; i < countMarkers(); i++) if (genotype(id, i) == genotype_NA) count++;
		return count;
	}
	virtual	int			countHeterozygous(iid_t id) const {
		int	count = 0;
		for (int i = 0; i < countMarkers(); i++) if (genotype(id, i) == genotype_1) count++;
		return count;
	}
	virtual iid_t		N(void) const {
		throw("Abstract method called (# of individuals).");
	}
	// <N> this returns an upper bound, reconstruction algorithm has to check
	// whether d1 <= d2 for reconstructions to be unordered
	int			countReconstructions(iid_t id) const {
		return !countHeterozygous(id)
		// algorithm has to correct special case of only missings: this value is for ordered reconstructions
		// unordered missing diplotypes are not contingent ordinally and reconstruction can more easliy 
		// achieved by oversampling possible diplotypes and filtering for unoreded diplotypes
		? (1 << (2 * countMissing(id)))
		: (1 << (countHeterozygous(id) - 1 + 2 * countMissing(id)));
	}
	haplotype_t		maskHeterozygosity(iid_t id) const {
		haplotype_t	mask = 0;
		for (marker_t i = 0; i < this->countMarkers(); i++)
			if (genotype(id, i) == genotype_1) mask = bitSet<haplotype_t>(mask, i);
		return mask;
	}
	haplotype_t		maskMissing(iid_t id) const {
		haplotype_t	mask = 0;
		for (marker_t i = 0; i < this->countMarkers(); i++)
			if (genotype(id, i) == genotype_NA) mask = bitSet<haplotype_t>(mask, i);
		return mask;
	}
	diplotype_t		diplotypeTemplate(iid_t id) const {
		diplotype_t	dt = { 0, 0 };
		for (marker_t i = 0; i < this->countMarkers(); i++)
			// set for homozygous for allele 1, 0 otherwise
			if (genotype(id, i) == genotype_2) {
				dt.d1 = bitSet<haplotype_t>(dt.d1, i);
				dt.d2 = bitSet<haplotype_t>(dt.d2, i);
			}
		return dt;
	}
	inline genotypecomb_t	genotypeCombination(iid_t id) const {
		genotypecomb_t	r = 0;
		BitArray<genotypecomb_t, marker_t>	b(&r, 0, 2, countMarkers());
		for (marker_t i = 0; i < countMarkers(); i++) {
			if (genotype(id, i) != genotype_NA) b.set(i, genotype(id, i));
		}
		return r;
	}
	void	print(void) {
		cout << "Genotypes:" << endl;
		for (int i = 0; i < N(); i++) {
			cout << "\t";
			for (int j = 0; j < countMarkers(); j++)
				cout << (j? ", ": "") << genotype(i, j);
			cout << endl;	
		}
	}
};

template<typename T>
class GenotypeFetcherMatrix : public GenotypeFetcher {
	vector< vector<T> >	&genotypes;
public:
	GenotypeFetcherMatrix(vector< vector<T> > &gm) : genotypes(gm) {}

	virtual	marker_t	countMarkers(void) const { return (marker_t)genotypes.size(); }
	virtual iid_t		N(void) const { return (iid_t)genotypes[0].size(); }
	virtual	genotype_t	genotype(iid_t id, marker_t marker) const {
		return (genotype_t)genotypes[marker][id];
	}
};


#endif
