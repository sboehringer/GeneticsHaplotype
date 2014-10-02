/*
 * Copyright (c) 2014, Stefan Boehringer <email>
 * All rights reserved.
 *
 * Thu Jul 10 16:59:54 CEST 2014
 */

#ifndef DIPLOTYPERECONSTRUCTION_H
#define DIPLOTYPERECONSTRUCTION_H

#include <valarray>
#include <memory>
#include "valarray_ext.h"
#include "genotype.h"
#include "buffer.h"
#include "pedigree.h"
#include "helpers.h"

typedef double					haplotypefs_t;
typedef Valarray<haplotypefs_t>	hfs_t;
typedef vector<haplotypefs_t>	hfsv_t;
typedef double					random_t;

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
	/*
	 *	reconstruction methods
	 */
	virtual void	reconstruct(GenotypeFetcher &fetcher);
	virtual void	print(void) const;
	inline iid_t	Nfounders(void) const { return pedigree.sizeFounders(); }
	inline iid_t	Nitrios(void) const { return pedigree.sizeItrios(); }
};

typedef unsigned long long	buffer_t;

class DiplotypeReconstructionSNPunordered : public DiplotypeReconstruction
{
	int					bitsFactor_;
	int					bitsHt_;
	size_t				reconstructionSize;
	Buffer<buffer_t>	reconstruction;
public:
	/*
	 *	creation / destruction / boilerplate 
	 */
	DiplotypeReconstructionSNPunordered(Pedigree &_pedigree,
		int _bitsFactor = 6, int _bitsHt = 10, iid_t NreconstrBuffer = 1024);
	DiplotypeReconstructionSNPunordered(const DiplotypeReconstructionSNPunordered&) = delete;
	~DiplotypeReconstructionSNPunordered();
	DiplotypeReconstructionSNPunordered(DiplotypeReconstructionSNPunordered &&other)
	: DiplotypeReconstruction(other.pedigree)
	{
		bitsFactor = other.bitsFactor;
		bitsHt = other.bitsHt;
		reconstructionSize = other.reconstructionSize;
		reconstruction = std::move(other.reconstruction);
	}

	virtual void	reconstruct(GenotypeFetcher &fetcher);
	virtual void	print(void) const;

	void			drawFromLogHfs(const hfs_t &lhfs, const random_t lu, haplotypes_t &draw) const;
	void			drawFromHfs(const hfs_t &hfs, const random_t u, haplotypes_t &draw) const;

	inline int		bitsFactor(void) const { return bitsFactor_; }
	inline int		bitsHt(void) const { return bitsHt_; }
	inline buffer_t	*reconstructionAt(int i) const { return reconstruction.buffer(i); }
};

struct ReconstructionArray {
	ReconstructionArray(DiplotypeReconstructionSNPunordered& d, int i)
	: hts(d.reconstructionAt(i), 0, d.bitsHt(), 2*d.Nfounders()),
	  iv(BitArrayAfter_e, hts, 1, 2*d.Nitrios()),
	  factor(BitArrayAfter_e, iv, d.bitsFactor(), 1)
	{}

	BitArray<buffer_t, iid_t>	hts;
	BitArray<buffer_t, iid_t>	iv;
	BitArray<buffer_t, iid_t>	factor;
};

/*
 * DiplotypeReconstructionSNPunorderedRaw
 * This is a helper class that buffers certain computations and provides raw reconstructions
 * for founder diplotypes
 */

/* Produce valid diplotype from a template + index for a possible iteration
 * This function produces ordered diplotypes
 * 
 * comb is current ordinal iteration
 * lower bits interpreted to belong to missingness
 * higher bits to account for heterozygous positions
 */

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
	t.d1 |= bitEmbedValueInOnes(het, het1, countHet);
	t.d2 |= bitEmbedValueInOnes(het, ~het1, countHet);

	return t;
}

class DiplotypeReconstructionSNPunorderedRaw
{
	vector<iid_t>			&founders;
	GenotypeFetcher			&fetcher;
	haplotypes_t			missing;
	haplotypes_t			heterozygous;
	// diplotype templates, pre-filled for homozygous positions
	vector<diplotype_t>		templates;
	unique_ptr< CartesianIterator<int> >	founderIterator;

public:
	/*
	 *	creation / destruction / boilerplate 
	 */
// 	DiplotypeReconstructionSNPunorderedRaw() :
// 		missing(), heterozygous(),templates()
// 	{}
	DiplotypeReconstructionSNPunorderedRaw(Pedigree &_pedigree, GenotypeFetcher &_fetcher) :
		founders(_pedigree.founders()), fetcher(_fetcher),
		missing(founders.size()),
		heterozygous(founders.size()),
		templates(founders.size()) {

		vector<int>	countFounders(founders.size());
		for (iid_t i = 0; i < founders.size(); i++) {
			missing[i] = fetcher.maskMissing(founders[i]);
			heterozygous[i] = fetcher.maskHeterozygosity(founders[i]);
			countFounders[i] = fetcher.countReconstructions(founders[i]);
			templates[i] = fetcher.diplotypeTemplate(founders[i]);
		}
		founderIterator = unique_ptr< CartesianIterator<int> >(new CartesianIterator<int>(countFounders));
	}
	~DiplotypeReconstructionSNPunorderedRaw() {
	}

	DiplotypeReconstructionSNPunorderedRaw &operator=(DiplotypeReconstructionSNPunorderedRaw &&other) {
		founders = other.founders;
		fetcher = other.fetcher;
		missing = other.missing;
		heterozygous = other.heterozygous;
		templates = other.templates;
		founderIterator.reset(other.founderIterator.get());
		other.founderIterator.reset(nullptr);
		return *this;
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
		// last element was generated and founderIterator is exhausted
		// the latter condition is checked by the first if
		return i >= founders.size();
	}
};

#endif // DIPLOTYPERECONSTRUCTION_H
