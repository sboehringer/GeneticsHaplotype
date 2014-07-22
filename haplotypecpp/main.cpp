#include <iostream>
#include <vector>
using namespace std;
#include "pedigree.h"
#include "diplotypereconstruction.h"

typedef vector<int>	g;
typedef struct {
	vector<iid_t>			founder;
	vector< vector<iid_t> >	trios;
} pedigree_t;

pedigree_t	ped1 {
	{ 0, 1, 4 },
	{ 	vector<iid_t> { 2, 0, 1 },
		vector<iid_t> { 3, 0, 1 },
		vector<iid_t> { 5, 3, 4 } }
};
vector<g>	ped1Genotypes1 {	// genotypes by row, individuals by column 0..5
	g { 1, 1, 1, 1, 1, 1 },
	g { 1, 1, 1, 1, 1, 1 },
	g { 1, 1, 1, 1, 1, 1 }
};
vector<g>	ped1Genotypes2 {	// genotypes by row, individuals by column 0..5
	g { 2, 1, 1, 1, 1, 1 },
	g { 1, 1, 1, 1, 1, 1 },
	g { 1, 1, 1, 1, 1, 1 }
};
vector<g>	ped1Genotypes3 {	// genotypes by row, individuals by column 0..5
	g { 3, 1, 1, 1, 1, 1 },
	g { 1, 1, 1, 1, 1, 1 },
	g { 1, 1, 1, 1, 1, 1 }
};
vector<g>	ped1Genotypes4 {	// genotypes by row, individuals by column 0..5
	g { 2, 1, 2, 1, 1, 1 },
	g { 1, 1, 1, 1, 1, 1 },
	g { 1, 1, 1, 1, 1, 1 }
};

void	reconstruct(pedigree_t &ped, vector<g> &gts) {
	Pedigree					pedigree(ped.founder, ped.trios);
	GenotypeFetcherMatrix<int>	gtF(gts);

	DiplotypeReconstructionSNPunordered	dts(pedigree);
	dts.reconstruct(gtF);
	dts.print();
}

int main(int argc, char **argv) {
	reconstruct(ped1, ped1Genotypes1);
	reconstruct(ped1, ped1Genotypes2);
	reconstruct(ped1, ped1Genotypes3);
	reconstruct(ped1, ped1Genotypes4);
	return 0;
}
