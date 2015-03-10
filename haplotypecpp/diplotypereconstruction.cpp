/*
 * Copyright (c) 2014, Stefan Boehringer <email>
 * All rights reserved.
 *
 * Thu Jul 10 17:01:56 CEST 2014
 */

#include <algorithm>
#include <math.h>
#include "diplotypereconstruction.h"

DiplotypeReconstruction::~DiplotypeReconstruction() {}

void	DiplotypeReconstruction::reconstruct(GenotypeFetcher &fetcher) {
	throw("abstract method called");
}

void	DiplotypeReconstruction::print(void) const {
	throw("abstract method called");
}

iid_t	DiplotypeReconstruction::Nreconstruction(void) const {
	throw("abstract method called");
	return 0;
}

#if 0

DiplotypeReconstructionSNPunordered::DiplotypeReconstructionSNPunordered(
	Pedigree &_pedigree, int _bitsFactor, int _bitsHt, iid_t NreconstrBuffer
	) :
	DiplotypeReconstruction(_pedigree), bitsFactor(_bitsFactor), bitsHt(_bitsHt),
	reconstructionSize(roundUp<size_t>(
		bitsHt + (_pedigree.sizeFounders() * 2 * bitsFactor + _pedigree.sizeItrios()), sizeof(size_t) * 8)),
	reconstruction(reconstructionSize * NreconstrBuffer) {
	
}

#endif

/*
 * Haplotype reconstruction algorithm
 * 
 * GenotypeFetcher:
 * provides:
 * - missingness mask for genotype vector
 * - integer coded vector (for SNPs: gt = \sum_i a_i * 3^i
 * 
 * Reconstruction (unordered):
 * - initialize diplotype counts for founders
 *   - counts are 2^(#het - 1) * 4^(#missing) (unorderedness absorbed into het-positions)
 *   - in case of #het == 0: count = 3 * 4^(#missing-1) (unordered case absorbed into missingness) (*)
 * - jointly iterate counts
 * - for each counter combination
 *   - form diplotypes per founder
 *     - seperate counter into ambiguous, missing part (upper bits missing -> modulo =^= bitshift)
 *     - embed missing counter on both diplotypes in missingness part
 *       - because iteration is exhaustive, take lower half of missingness counter and embed into diplotype 1
 *       - embed higher part into diplotype 2
 *       - assume 2 bits per locus and half the number of bits for the split
 *         - case (*) is special: if #het == 0 && locus == 0 && a1 > a2 continue
 *     - embed ambiguous part into het-bitmap, other diplotype gets complement
 *   - iterate inheritance vectors (including direction)
 *     - compare expected transmission genotype against offspring genotype
 * - determine mulitplicity (due to unorderedness during algorithm above)
 *   - #het > 0 contributes factor of two per founder
 *   - ambiguous transmissions do not contribute (already covered by step above)
 *
 *   - of the 4 possiblities a maximum of 2 can be possible be leading to same unordered offsping diplotype (**)
 *      - proof: mother can transmit both haplotypes to offspring => subtraction in offspring determines other
 *        alleles of father => different diplotype constitutions possible for offspring
 */

inline haplotype_t	selectFromDiplotype(diplotype_t dt, bool index) {
	return !index? dt.d1: dt.d2;
}
inline bool	diplotypeIsHet(diplotype_t dt) {
	return dt.d1 != dt.d2;
}
inline bool	diplotypeIsHom(diplotype_t dt) {
	return !diplotypeIsHet(dt);
}

void	DiplotypeReconstructionSNPunordered::reconstruct(GenotypeFetcher &fetcher) {
	DiplotypeReconstructionSNPunorderedRaw	founderReconst(pedigree, fetcher);
	vector<diplotype_t>						dtFounder(Nfounders());
	// all diplotypes
	vector<diplotype_t>						dts(Nfounders() + Nitrios());
	vector<bool>							bothTransm(Nitrios());
	vector<bool>							iv(2*Nitrios());
	marker_t								cm = fetcher.countMarkers();

	while (founderReconst.founderReconstruction(dtFounder)) {
		iid_t	factor = 0;	//factor of multiplicity for reconstruction
		for (int i = 0; i < dtFounder.size(); i++) {
			dts[pedigree.founders()[i]] = dtFounder[i];
			factor += diplotypeIsHet(dtFounder[i]);
		}
		// inheritance vector iterator
		for (iid_t ivI = 0, j = 0; ivI < (1 << (2*pedigree.sizeItrios())); ivI++) {
			// check for consistent transmission
			for (j = 0; j < pedigree.sizeItrios(); j++) {
				iid_t			id = pedigree.trioIid(j), mid = pedigree.trioMid(j), pid = pedigree.trioPid(j);
				haplotype_t		miss = fetcher.maskMissing(id);
				haplotype_t		htm = selectFromDiplotype(dts[mid], bitAt<iid_t>(ivI, 2*j));
				haplotype_t		htp = selectFromDiplotype(dts[pid], bitAt<iid_t>(ivI, 2*j + 1));
				// expected genotype combination
				genotypecomb_t	gtcE = genotypeCombinationFromDiplotype((diplotype_t){htm, htp}, cm, miss);
				genotypecomb_t	gtcO = fetcher.genotypeCombination(id);	//offspring

				// if founder was homozygous, we only consider one transmission (we are unordered)
				if ( (diplotypeIsHom(dts[mid]) && bitAt<iid_t>(ivI, 2*j) > 0)
				  || (diplotypeIsHom(dts[pid]) && bitAt<iid_t>(ivI, 2*j + 1) > 0)) {
					factor += diplotypeIsHom(dts[mid]) + diplotypeIsHom(dts[pid]);
					break;
				}
				if (gtcO != gtcE) break;
				dts[id] = (diplotype_t){htm, htp};
			}
			// founder diplotypes not compatible with offspring genotypes
			if (j < pedigree.sizeItrios()) continue;
			// save reconstruction
			buffer_t	*e = reconstruction.push();
			ReconstructionArrayIvBlock	r(*this, e);
			r.set(dtFounder, ivI);
		}

	}
}

DiplotypeReconstructionSNPunordered::DiplotypeReconstructionSNPunordered(
	Pedigree &_pedigree, int _bitsFactor, int _bitsHt, iid_t NreconstrBuffer)
	:
	DiplotypeReconstruction(_pedigree),
	bitsFactor_(_bitsFactor),
	bitsHt_(_bitsHt),
	reconstructionSize(roundUpBitsToBytes<buffer_t>(
		bitsFactor_ + bitsHt_ * 2 * Nfounders() + 2 * Nitrios())),
	reconstruction(reconstructionSize, NreconstrBuffer) {
#	if 0
	cout << "reconstructionSize: " << reconstructionSize
		<< " pedigreesize: " << Nfounders() << ", " << Nitrios()
		<< " bitsHt: " << bitsHt_ << " bitsFactor: " << bitsFactor_
		<< " Size: " << _bitsFactor + _bitsHt * 2 * Nfounders() + 2 * Nitrios()
		<< endl;

#	endif
}

typedef DiplotypeReconstructionSNPunordered	DRU;

void	DiplotypeReconstructionSNPunordered::print(void) const {
	cout << "\tReconstruction size: " << reconstruction.size() << endl;
	for (int j = 0; j < reconstruction.size(); j++) {
		ReconstructionArray	r(*(DRU *)this, j);

		cout << "\tFounders[";
		for (iid_t i = 0; i < Nfounders(); i++)
			cout << (i? " ": "") << "(" << r.hts[2*i] << ", " << r.hts[2*i + 1] << ")";
		cout << "] IV[";
		for (iid_t i = 0; i < Nitrios(); i++) {
			cout << r.iv[2*i] << r.iv[2*i + 1];
		}
		cout << "] C[" << r.factor[0] << "]" << endl;
	}
}

// inheritance vector
int	DiplotypeReconstructionSNPunordered::ivAt(iid_t i) const {
	ReconstructionArray	r(*(DRU *)this, i);
	int	iv = 0;
	for (iid_t j = 0; j < Nitrios(); j++) iv |= (r.iv[2*i] << (2*i)) | (r.iv[2*i + 1] << (2*i + 1));
	return iv;
}

// multiplicative constant
int	DiplotypeReconstructionSNPunordered::factorAt(iid_t i) const {
	ReconstructionArray	r(*(DRU *)this, i);
	return r.factor[0];
}

// founder diplotype
diplotype_t	DiplotypeReconstructionSNPunordered::diplotypeAt(iid_t i, iid_t j) const {
	ReconstructionArray	r(*(DRU *)this, i);
	return (diplotype_t) { (haplotype_t)r.hts[2*j], (haplotype_t)r.hts[2*j + 1] };
}


//#define __DEBUG_PROB

/*
 * copy diplotypes
 */
void	DiplotypeReconstructionSNPunordered::codeHaplotypesIntoVector(int i, haplotypes_t &draw) const {
	draw.resize(2*Nfounders() + 2*Nitrios());
	ReconstructionArray	r(*(DRU *)this, i);

	for (int j = 0; j < Nfounders(); j++) {
		draw[2*pedigree.founderAt(j)]     = r.hts[2*j];
		draw[2*pedigree.founderAt(j) + 1] = r.hts[2*j + 1];
	}
	for (int j = 0; j < Nitrios(); j++) {
		draw[2*pedigree.trioIid(j)]     = draw[2*pedigree.trioMid(j) + r.iv[2*j]];
		draw[2*pedigree.trioIid(j) + 1] = draw[2*pedigree.trioPid(j) + r.iv[2*j + 1]];
	}
}

void	DiplotypeReconstructionSNPunordered::drawFromLogHfs(const hfs_t &lhfs, const random_t lu,
															haplotypes_t &draw) const {
#	ifdef __DEBUG_PROB
	cout << endl;
	lhfs.print();
	hfs_t	lhfsE(exp(lhfs));
	lhfsE.print();
#	endif
	/*
	 * compute log-prob per reconstruction (up to normalizing factor)
	 */
	Valarray<haplotypefs_t>	famlh(0, (size_t)reconstruction.size());
	for (int i = 0; i < famlh.size(); i++) {
		ReconstructionArray	r(*(DRU *)this, i);

		famlh[i] = r.factor[0] * M_LN2;
		for (int j = 0; j < 2 * Nfounders(); j++) {
			famlh[i] += lhfs[r.hts[j]];
		}

#		ifdef __DEBUG_PROB
		cout << "Reconstruction " << i << ": " << famlh[i] << " [" << exp(famlh[i]) << "]" << endl;
#		endif
	}

	/*
	 * draw reconstruction number
	 */
	Valarray<haplotypefs_t>	cs(famlh.log_cumsum());
	haplotypefs_t			csmx = cs[cs.size() - 1];
	//haplotypefs_t			csmx = cs.back();
	int						i = cs.binary_search(lu + csmx);
	assert(i < cs.size());
#	ifdef __DEBUG_PROB
	cout << "cs: ";
	cs.print();
	cout << "Csmax " << ": " << exp(csmx) << " [" << csmx << "], Search: " << exp(lu + csmx)
		<< " [" << lu + csmx << "]: " << i << endl;
#	endif
	codeHaplotypesIntoVector(i, draw);
}

void	DiplotypeReconstructionSNPunordered::drawFromHfs(const hfs_t &hfs, const random_t u,
														 haplotypes_t &draw) const {
	hfs_t	hfsN(hfs / hfs.sum());
	drawFromLogHfs(log(hfsN), log(u), draw);
}

DiplotypeReconstructionSNPunordered::~DiplotypeReconstructionSNPunordered() {}
