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

void	DiplotypeReconstruction::print(void){
	throw("abstract method called");
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

inline haplotype_t	selectFromDiplotype(diplotype_t dt, bool index) {
	return !index? dt.d1: dt.d2;
}

void	DiplotypeReconstructionSNPunordered::reconstruct(const GenotypeFetcher &fetcher) {
	DiplotypeReconstructionSNPunorderedRaw	founderReconst(pedigree, fetcher);
	vector<diplotype_t>						dtFounder(pedigree.sizeFounders());
	// all diplotypes
	vector<diplotype_t>						dts(pedigree.sizeFounders() + pedigree.sizeItrios());
	vector<bool>							bothTransm(pedigree.sizeItrios());
	vector<bool>							iv(2*pedigree.sizeItrios());
	marker_t								cm = fetcher.countMarkers();

	while (founderReconst.founderReconstruction(dtFounder)) {
		for (int i = 0; i < dtFounder.size(); i++)
			dts[pedigree.founders()[i]] = dtFounder[i];

		iid_t j;
		for (j = 0; j < pedigree.sizeItrios(); j++) {
			iid_t			id = pedigree.trioIid(j), mid = pedigree.trioMid(j), pid = pedigree.trioPid(j);
			diplotype_t		dtm = dts[mid], dtp = dts[pid];
			haplotype_t		miss = fetcher.maskMissing(id);
			// read index gtcIJ right-to-left
			genotypecomb_t	gtc00 = genotypeCombinationFromDiplotype((diplotype_t){dtm.d1, dtp.d1}, cm, miss);
			genotypecomb_t	gtc10 = genotypeCombinationFromDiplotype((diplotype_t){dtm.d1, dtp.d2}, cm, miss);
			genotypecomb_t	gtc01 = genotypeCombinationFromDiplotype((diplotype_t){dtm.d2, dtp.d1}, cm, miss);
			genotypecomb_t	gtc11 = genotypeCombinationFromDiplotype((diplotype_t){dtm.d2, dtp.d2}, cm, miss);
			genotypecomb_t	gtcO = fetcher.genotypeCombination(id);	//offspring

			if (gtcO != gtc00 && gtcO != gtc01 && gtcO != gtc10 && gtcO != gtc11) break;
			// compute possible inheritance vector: first haplotype from grandmother =^= 0
			// if needed to make both entries of iv consistent (<N> special case: bothTransm)
			if ((gtcO == gtc00 || gtcO == gtc10)) {
				iv[2*j] = 0;
				iv[2*j + 1] = gtcO == gtc10;
			} else {
				iv[2*j] = 1;
				iv[2*j + 1]	= gtcO == gtc11;
			}
			// save offspring diplotype for reference for other I-trios
			dts[id] = (diplotype_t){
				selectFromDiplotype(dtm, iv[2*j]),
				selectFromDiplotype(dtp, iv[2*j + 1]) };
			// are both directions of transmission possible?
			bothTransm[j] = (gtcO == gtc00 && gtcO == gtc11) || (gtcO == gtc01 && gtcO == gtc10);
		}
		// founder diplotypes not compatible with offspring genotypes
		if (j < pedigree.sizeItrios()) continue;

		// encode reconstruction
		buffer_t	*e = reconstruction.push();
		int			factor = 0;
		BitArray<buffer_t, iid_t>	Ahts(e, 0, bitsHt, 2*dtFounder.size());
		BitArray<buffer_t, iid_t>	Aiv(BitArrayAfter_e, Ahts, 1, 2*pedigree.sizeItrios());
		BitArray<buffer_t, iid_t>	Afactor(BitArrayAfter_e, Aiv, bitsFactor, 1);

		for (iid_t i = 0; i < dtFounder.size(); i++) {
			Ahts.set(2*i,     dtFounder[i].d1);
			Ahts.set(2*i + 1, dtFounder[i].d2);
			if (dtFounder[i].d1 != dtFounder[i].d2) factor++;
		}
		for (iid_t i = 0; i < pedigree.sizeItrios(); i++) {
			Aiv.set(2*i,	iv[2*i]);
			Aiv.set(2*i + 1,iv[2*i + 1]);
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
	DiplotypeReconstruction(_pedigree),
	bitsFactor(_bitsFactor),
	bitsHt(_bitsHt),
	reconstructionSize(roundUpBitsToBytes<buffer_t>(
		bitsFactor + bitsHt * 2 * pedigree.sizeFounders() + 2 * pedigree.sizeItrios())),
	reconstruction(reconstructionSize, NreconstrBuffer) {
	cout << "reconstructionSize: " << reconstructionSize
		<< " pedigreesize: " << pedigree.sizeFounders() << ", " << pedigree.sizeItrios()
		<< " bitsHt: " << bitsHt << " bitsFactor: " << bitsFactor
		<< " Size: " << _bitsFactor + _bitsHt * 2 * _pedigree.sizeFounders() + 2 * _pedigree.sizeItrios()
		<< endl;
}

void	DiplotypeReconstructionSNPunordered::print(void) {
	cout << "\tReconstruction size: " << reconstruction.size() << endl;
	for (int j = 0; j < reconstruction.size(); j++) {
		buffer_t					*e = reconstruction.buffer(j);
		BitArray<buffer_t, iid_t>	Ahts(e, 0, bitsHt, 2*pedigree.sizeFounders());
		BitArray<buffer_t, iid_t>	Aiv(BitArrayAfter_e, Ahts, 1, 2*pedigree.sizeItrios());
		BitArray<buffer_t, iid_t>	Afactor(BitArrayAfter_e, Aiv, bitsFactor, 1);

		cout << "\tFounders[";
		for (iid_t i = 0; i < pedigree.sizeFounders(); i++)
			cout << (i? " ": "") << "(" << Ahts[2*i] << ", " << Ahts[2*i + 1] << ")";
		cout << "] IV[";
		for (iid_t i = 0; i < pedigree.sizeItrios(); i++) {
			cout << Aiv[2*i] << Aiv[2*i + 1];
		}
		cout << "] C[" << Afactor[0] << "]" << endl;
	}
}

DiplotypeReconstructionSNPunordered::~DiplotypeReconstructionSNPunordered() {}
