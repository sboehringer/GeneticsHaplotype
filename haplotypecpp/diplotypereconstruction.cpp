/*
 * Copyright (c) 2014, Stefan Boehringer <email>
 * All rights reserved.
 *
 * Thu Jul 10 17:01:56 CEST 2014
 */

#include "diplotypereconstruction.h"

DiplotypeReconstruction::~DiplotypeReconstruction() {}

bool DiplotypeReconstruction::operator==(const DiplotypeReconstruction& other)
{

}

void	DiplotypeReconstruction::reconstruct(const GenotypeFetcher &fetcher) {
	throw("abstract method called");
}

template <typename T>
inline T roundUp(T v, T modulo) {
	return (v + modulo - 1)/modulo;
}

template <typename T>
inline T roundUpBitsToBytes(T bits) {
	return roundUp(bits, typeBits(T)) / bitsPerByte;
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
 *   - iterate inheritance vector
 *     - iterate all four possible transmissions
 *     - form transmitted genotype (with missingness masked to 0)
 *     - compare against offspring genotype
 *     - of the 4 possiblities a maximum of 2 can be possible leading to same unordered offsping diplotype (**)
 *       - proof: mother can transmit both diplotypes to offspring => one diplotype shared between parents for
 *         first transmission, symmetric argument => parents have same diplotype => offspring has same diplotype
 * - determine mulitplicity (due to unorderedness during algorithm above)
 *   - #het > 0 contributes factor of two per founder
 *   - ambiguous transmissions do not contribute (already covered by step above)
 */

void	DiplotypeReconstructionSNPunordered::reconstruct(const GenotypeFetcher &fetcher) {
	DiplotypeReconstructionSNPunorderedRaw	founderReconst(pedigree, fetcher);
	vector<diplotype_t>						dtFounder(pedigree.sizeFounders());
	vector<diplotype_t>						dtOffspring(pedigree.sizeItrios());
	vector<bool>							bothTransm(pedigree.sizeItrios());
	vector<bool>							iv(2*pedigree.sizeItrios());

	while (founderReconst.founderReconstruction(dtFounder)) {
		iid_t j;
		for (j = 0; j < pedigree.sizeItrios(); j++) {
			iid_t			id = pedigree.trioIid(j), mid = pedigree.trioMid(j), pid = pedigree.trioPid(j);
			diplotype_t		dtm = dtFounder[mid], dtp = dtFounder[pid];
			haplotype_t		miss = fetcher.maskMissing(id);
			genotypecomb_t	gtc00 = genotypeCombinationFromDiplotype((diplotype_t){dtm.d1, dtp.d1}, miss);
			genotypecomb_t	gtc01 = genotypeCombinationFromDiplotype((diplotype_t){dtm.d1, dtp.d2}, miss);
			genotypecomb_t	gtc10 = genotypeCombinationFromDiplotype((diplotype_t){dtm.d2, dtp.d1}, miss);
			genotypecomb_t	gtc11 = genotypeCombinationFromDiplotype((diplotype_t){dtm.d2, dtp.d2}, miss);
			genotypecomb_t	gtcO = fetcher.genotypeCombination(id);	//offspring

			if (gtcO != gtc00 && gtcO != gtc01 && gtcO != gtc10 && gtcO != gtc11) break;
			// compute possible inheritance vector
			iv[2*j]		= gtcO == gtc10 || gtcO == gtc11;
			iv[2*j + 1]	= gtcO == gtc01 || gtcO == gtc11;
			// are both directions of transmission possible?
			bothTransm[j] = (gtcO == gtc00 && gtcO == gtc11) || (gtcO == gtc01 && gtcO == gtc10);
		}
		// founder diplotypes not compatible with offspring genotypes
		if (j < pedigree.sizeItrios()) continue;

		// encode reconstruction
		buffer_t	*e = reconstruction.push();
		int			factor = 0;
		BitArray<buffer_t, iid_t>	Ahts(e, bitsFactor, bitsHt);
		BitArray<buffer_t, iid_t>	Aiv(e, bitsFactor + 2*bitsHt, 1);
		BitArray<buffer_t, int>		Afactor(e, 0, bitsFactor);

		for (iid_t i = 0; i < dtFounder.size(); i++) {
			Ahts.set(2*i,     dtFounder[i].d1);
			Ahts.set(2*i + 1, dtFounder[i].d2);
			if (dtFounder[i].d1 != dtFounder[i].d2) factor++;
		}
		for (iid_t i = 0; i < dtOffspring.size(); i++) {
			Aiv.set(i, iv[2*i] + 2*iv[2*i + 1]);
			if (bothTransm[i]) factor++;
		}
		Afactor.set(0, factor);
	}
// 	//CartesianIterator<int>	&i = *(this->founderIterator(fetcher));
// 
// 	
}

DiplotypeReconstructionSNPunordered::DiplotypeReconstructionSNPunordered(
	Pedigree &_pedigree, int _bitsFactor, int _bitsHt, iid_t NreconstrBuffer)
	:
	bitsFactor(_bitsFactor),
	bitsHt(_bitsHt),
	reconstructionSize(bitsFactor + bitsHt * 2 * pedigree.sizeFounders() + 2 * pedigree.sizeItrios()),
	reconstruction(reconstructionSize, NreconstrBuffer),
	DiplotypeReconstruction(_pedigree) {
}