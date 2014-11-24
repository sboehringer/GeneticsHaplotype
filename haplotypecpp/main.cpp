#include <iostream>
#include <vector>
using namespace std;
#include "pedigree.h"
#include "diplotypereconstruction.h"

double	randL(void) {
	return ((double) rand())/((double)RAND_MAX + 1);
}

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
vector<g>	ped1Genotypes5 {	// genotypes by row, individuals by column 0..5
	g { 2, 2, 2, 2, 2, 2 },
	g { 0, 0, 0, 0, 0, 0 },
	g { 1, 1, 1, 1, 1, 1 }
};


pedigree_t	ped2 {
	{ 0, 1, 3 },
	{ 	vector<iid_t> { 2, 0, 1 },
		vector<iid_t> { 4, 2, 3 } }
};

vector<g>	ped2Genotypes1 {	// genotypes by row, individuals by column 0..5
	g { 1, 2, 2, 2, 2 },
	g { 0, 1, 1, 0, 1 },
	g { 1, 1, 1, 1, 1 }
};

pedigree_t	ped3 {
	{ 0, 1 },
	{ 	vector<iid_t> { 2, 0, 1 } }
};

vector<g>	ped3Genotypes1 {	// genotypes by row, individuals by column 0..2
	g { 1, 2, 2 },
	g { 0, 1, 1 }
};
vector<g>	ped3Genotypes2 {	// genotypes by row, individuals by column 0..2
	g { 1, 2, 2 },
	g { 1, 1, 1 }
};
vector<g>	ped3Genotypes3 {	// genotypes by row, individuals by column 0..2
	g { 1, 1, 1 }
};
vector<g>	ped3Genotypes4 {	// genotypes by row, individuals by column 0..2
	g { 2, 2, 2 }
};
vector<g>	ped3Genotypes5 {	// genotypes by row, individuals by column 0..2
	g { 2, 1, 1 }
};

void	reconstruct(pedigree_t &ped, vector<g> &gts) {
	Pedigree					pedigree(ped.founder, ped.trios);
	GenotypeFetcherMatrix<int>	gtF(gts);

	DiplotypeReconstructionSNPunordered	dts(pedigree);
	dts.reconstruct(gtF);
	dts.print();
}

void	draw(pedigree_t &ped, vector<g> &gts, const random_t u, const hfs_t &hfs) {
	Pedigree					pedigree(ped.founder, ped.trios);
	GenotypeFetcherMatrix<int>	gtF(gts);
	GenotypeFetcherOffset		gfo(gtF, 0);

	DiplotypeReconstructionSNPunordered	dts(pedigree);
	dts.reconstruct(gfo);
// 	dts.print();

	haplotypes_t hts(0);
	dts.drawFromHfs(hfs, u, hts);

	for (int i = 0; i < hts.size(); i++) cout << (i > 0? (i%2? "|": ", "): "") << hts[i];
	cout << endl;
}

void	drawN(int Ndraws, pedigree_t &ped, vector<g> &gts, const hfs_t &hfs) {
	srand(time(NULL));
	reconstruct(ped, gts);
	for (int i = 0; i < Ndraws; i++) draw(ped, gts, randL(), hfs);
}

int main(int argc, char **argv) {
#	if 0
	reconstruct(ped1, ped1Genotypes1);
	reconstruct(ped1, ped1Genotypes2);
	reconstruct(ped1, ped1Genotypes3);
	reconstruct(ped1, ped1Genotypes4);
#	endif

#	if 0
#	define	Ndraws	4
	//const hfs_t		hfs1(vector<haplotypefs_t> { 1, 1, 1, 1, 5, 6, 7, 8 });
	//const hfs_t		hfs1(vector<haplotypefs_t> { 8, 7, 6, 5, 1, 1, 1, 1 });
	const hfs_t		hfs1(vector<haplotypefs_t> { 1, 1, 1, 1, 1, 1, 1, 1 });
	//drawN(Ndraws, ped1, ped1Genotypes4, hfs1);
	drawN(Ndraws, ped2, ped2Genotypes1, hfs1);
#	endif
#	if 0
#	define	Ndraws	4
	const hfs_t		hfs1(vector<haplotypefs_t> { 8, 4, 1, 1 });
	//drawN(Ndraws, ped1, ped1Genotypes4, hfs1);
	drawN(Ndraws, ped3, ped3Genotypes2, hfs1);
#	endif
#	if 1
#	define	Ndraws	4
	const hfs_t		hfs1(vector<haplotypefs_t> { 8, 4 });
	//drawN(Ndraws, ped1, ped1Genotypes4, hfs1);
	drawN(Ndraws, ped3, ped3Genotypes5, hfs1);
#	endif
	return 0;
}
