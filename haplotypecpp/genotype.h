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

typedef int				iid_t;
typedef int				marker_t;
typedef unsigned int	haplotype_t;
typedef enum {genotype_0, genotype_1, genotype_2, genotype_NA } genotype_t;
typedef struct {
	haplotype_t	d1, d2;
} diplotype_t;
// enumerate geneotype combinations
typedef unsigned int	genotypecomb_t;

inline genotypecomb_t	genotypeCombinationFromDiplotype(diplotype_t dt, marker_t countMarker, haplotype_t missing = 0) {
	genotypecomb_t	r = 0;
	dt.d1 &= ~missing;
	dt.d2 &= ~missing;
	for (marker_t i = 0; i < countMarker; i++) {
		if (bitAt(dt.d1, i)) r = bitSet<haplotype_t>(r, 2*i);
		if (bitAt(dt.d2, i)) r = bitSet<haplotype_t>(r, 2*i + 1);
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
		throw("Abstract genotype fetching method called.");
	}
	virtual	int			countHeterozygous(iid_t id) const {
		throw("Abstract genotype fetching method called.");
	}
	// <N> this returns an upper bound, reconstruction algorithm has to check
	// whether d1 <= d2 for reconstructions to be unordered
	int			countReconstructions(iid_t id) const {
		return !countHeterozygous(id)
		// algorithm has to correct this special case: this value is for ordered reconstructions
		? (1 << (2 * countMissing(id)))
		: (1 << (countHeterozygous(id) - 1 + 2 * countMissing(id)));
	}
	haplotype_t		maskHeterozygosity(iid_t id) const {
		haplotype_t	mask = 0;
		for (marker_t i; i < this->countMarkers(); i++)
			if (genotype(id, i) == genotype_1) mask = bitSet<haplotype_t>(mask, i);
		return mask;
	}
	haplotype_t		maskMissing(iid_t id) const {
		haplotype_t	mask = 0;
		for (marker_t i; i < this->countMarkers(); i++)
			if (genotype(id, i) == genotype_NA) mask = bitSet<haplotype_t>(mask, i);
		return mask;
	}
	diplotype_t		diplotypeTemplate(iid_t id) const {
		diplotype_t	dt = { 0, 0 };
		for (marker_t i; i < this->countMarkers(); i++)
			// set for homozygous for allele 1, 0 otherwise
			if (genotype(id, i) == genotype_2) {
				dt.d1 = bitSet<haplotype_t>(dt.d1, i);
				dt.d2 = bitSet<haplotype_t>(dt.d2, i);
			}
		return dt;
	}
	inline genotypecomb_t	genotypeCombination(iid_t id) const {
		genotypecomb_t	r = 0;
		BitField<genotypecomb_t, marker_t>	bf(&r);
		for (marker_t i = 0; i < countMarkers(); i++) {
			if (genotype(id, i) != genotype_NA) bf.set(r, i, genotype(id, i));
		}
		return r;
	}
};


#endif
