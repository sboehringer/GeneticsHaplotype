#
#	GeneticsHaplotype.R
#Wed Sep  3 18:18:06 CEST 2014

library('devtools');
library('Rcpp');
source('haplotype.R');

if (F) {
	system('rm GeneticsHaplotype/src/*.o GeneticsHaplotype/src/*.so');
	Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
	install('GeneticsHaplotype', threads = 6);
	require('GeneticsHaplotype');
}

if (0) {
	M = Module('Reconstructor', PACKAGE = 'GeneticsHaplotype');
}

if (0) {
	ped1 = list(
		founders = c( 0, 1, 4 ),
		itrios = matrix(c(c(2, 0, 1 ), c(3, 0, 1), c(5, 3, 4)), ncol = 3, byrow = T)
	);
	peds = list(ped1);
	gts = t(matrix(c(c( 1, 1, 1, 1, 1, 1), c(1, 1, 1, 1, 1, 1), c(1, 1, 1, 1, 1, 1)), ncol = 6, byrow = T));
	R = new(M$DiplotypeReconstructor, gts, peds);
}

if (0) {
print(hts);
}

freqs = function(v)table(v)/sum(table(v))

if (1) {
	N = 1e4;
	hts = sapply(1:N, function(i, ped, hfs) {
		hts = R$drawFromHfs(hfs, runif(1));
		as.vector(hts[ped$founders + 1, ])
	}, ped = ped1, hfs = 1:8);
	print(freqs(as.vector(hts)));
}