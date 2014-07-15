/*
 * Copyright (c) 2014, Stefan Boehringer <email>
 * All rights reserved.
 *
 * Thu Jul 10 16:59:54 CEST 2014
 */

#ifndef DIPLOTYPERECONSTRUCTION_H
#define DIPLOTYPERECONSTRUCTION_H

#include "genotype.h"
#include "buffer.h"
#include "pedigree.h"
#include "helpers.h"

class DiplotypeReconstruction
{
protected:
	Pedigree	&pedigree;
public:
	/*
	 *	creation / destruction / boilerplate 
	 */
    DiplotypeReconstruction(Pedigree &_pedigree) : pedigree(_pedigree) {}
    ~DiplotypeReconstruction();
    virtual bool operator==(const DiplotypeReconstruction& other);

	/*
	 *	reconstruction methods
	 */
	virtual void	reconstruct(const GenotypeFetcher &fetcher);
};

typedef unsigned long long	buffer_t;

class DiplotypeReconstructionSNPunordered : public DiplotypeReconstruction
{
	Buffer<buffer_t>	reconstruction;
	int					bitsFactor;
	int					bitsHt;
	size_t				reconstructionSize;
public:
	/*
	 *	creation / destruction / boilerplate 
	 */
	DiplotypeReconstructionSNPunordered(Pedigree &_pedigree,
		int _bitsFactor = 6, int _bitsHt = 10, iid_t NreconstrBuffer = 1024);
	~DiplotypeReconstructionSNPunordered();

	virtual void	reconstruct(const GenotypeFetcher &fetcher);
};

/*
 * DiplotypeReconstructionSNPunorderedRaw
 * This is a helper class that buffers certain computations and provides raw reconstructions
 * for founder diplotypes
 */

// comb is current ordinal iteration
// lower bits interpreted to belong to missingness
// higher bits to account for heterozygous positions
inline diplotype_t diplotypeFromTemplate(diplotype_t t, int comb, haplotype_t het, haplotype_t miss) {
	// number of heterozygous loci
	int	countHet = countBitsSet(het);
	// number of missing loci
	int	countMiss = countBitsSet(miss);

	// missing genotypes
	haplotype_t	miss1 = comb & BitMask<haplotype_t, int>(countMiss);
	t.d1 |= bitEmbedValueInOnes(miss, miss1, countMiss);
	haplotype_t	miss2 = (comb >> countMiss) & BitMask<haplotype_t, int>(countMiss);
	t.d2 |= bitEmbedValueInOnes(miss, miss2, countMiss);

	// heterozygous genotypes
	haplotype_t	het1 = (comb >> (2*countMiss)) & BitMask<haplotype_t, int>(countHet);
	t.d1 |= bitEmbedValueInOnes(miss, het1, countMiss);
	t.d2 |= bitEmbedValueInOnes(miss, ~het1, countMiss);

	return t;
}

class DiplotypeReconstructionSNPunorderedRaw
{
	const vector<iid_t>		&founders;
	const GenotypeFetcher	&fetcher;
	vector<haplotype_t>		missing;
	vector<haplotype_t>		heterozygous;
	// diplotype templates, pre-filled for homozygous positions
	vector<diplotype_t>		templates;
	CartesianIterator<int>	*founderIterator;

public:
	/*
	 *	creation / destruction / boilerplate 
	 */
    DiplotypeReconstructionSNPunorderedRaw(Pedigree &_pedigree, const GenotypeFetcher &_fetcher) :
		founders(_pedigree.founders()), fetcher(_fetcher),
		missing(founders.size()), heterozygous(founders.size()), templates(founders.size()), founderIterator(0) {

		vector<int>	countFounders(founders.size());
		for (iid_t i = 0; i < founders.size(); i++) {
			missing[i] = fetcher.maskMissing(founders[i]);
			heterozygous[i] = fetcher.maskHeterozygosity(founders[i]);
			countFounders[i] = fetcher.countReconstructions(founders[i]);
			templates[i] = fetcher.diplotypeTemplate(founders[i]);
		}
		founderIterator = new CartesianIterator<int>(countFounders);
	}
	~DiplotypeReconstructionSNPunorderedRaw() {
		delete founderIterator;
	}

	// returns a possible unordered reconstruction for all founders
	// successively iterates all possible combinations
	// returns false if combinations are exhausted, diplotypes has to have size #founders
	bool	founderReconstruction(vector<diplotype_t> &diplotypes) {
		if (founderIterator->isExhausted()) return false;

		iid_t	i;
		do {
			for (i = 0; i < founders.size(); i++) {
				diplotypes[i] = diplotypeFromTemplate(templates[i],
					(*founderIterator)[i], heterozygous[i], missing[i]);
				// assure unordered diplotypes
				if (diplotypes[i].d1 > diplotypes[i].d2) break;
			}
			++(*founderIterator);
		// tried diplotype was not unordered, still more to try
		} while (!founderIterator->isExhausted() && i < founders.size());
		return !founderIterator->isExhausted();
	}
};

#endif // DIPLOTYPERECONSTRUCTION_H
