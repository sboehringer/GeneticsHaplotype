#
#	GeneticsHaplotype.R
#Wed Sep  3 18:18:06 CEST 2014

library('devtools');
library('Rcpp');
source('GeneticsHaplotype/R/pedigree.R');
source('GeneticsHaplotype/R/simulation.R');
source('RgenericAll.R');

if (F) {
	system('rm GeneticsHaplotype/src/*.o GeneticsHaplotype/src/*.so');
	Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
	install('GeneticsHaplotype', threads = 6);
}
if (T) {
	require('GeneticsHaplotype');
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

if (0) {
	N = 1e4;
	hts = sapply(1:N, function(i, ped, hfs) {
		hts = R$drawFromHfs(hfs, runif(1));
		as.vector(hts[ped$founders + 1, ])
	}, ped = ped1, hfs = 1:8);
	print(freqs(as.vector(hts)));
}

if (0) {
	ped = Df(names = c('iid', 'fid', 'mid', 'pid'), matrix(
		c(	10, 1, NA, NA,
			20, 1, NA, NA,
			30, 1, 10, 20,
			40, 1, NA, NA,
			55, 1, 30, 40,

			10, 2, NA, NA,
			20, 2, NA, NA,
			30, 2, 10, 20,
			40, 2, NA, NA,
			55, 2, 30, 40
		)
	, byrow = T, ncol = 4));
	ped2 = pedSplit2ivTrios(ped);
	s = simulateDiplotypes(ped2, hfs = 1:8);
}

if (0) {
	pedP = ped[sample(nrow(ped), nrow(ped)), ];
}
if (0) {
	print(data.frame(ped, sex = pedInferSex(ped)));
	r = data.frame(pedP, sex = pedInferSex(pedP));
	print(r[order(r$fid, r$iid), ]);
}

if (0) {
	plotPedigree(ped);
}

if (1) {
	N = 25;
	pedTemplate = Df(names = c('iid', 'mid', 'pid'), matrix(
		c(	1, NA, NA,
			2, NA, NA,
			3, 1, 2,
			4, NA, NA,
			5, 3, 4
	), byrow = T, ncol = 3));
	ped = familiesFromTemplate(pedTemplate, N = N);
	#print(data.frame(ped, sex = pedInferSex(ped)));

	dts = simulateDiplotypesPed(ped, 1:8);
	#plotPedigrees(ped, tag = apply(dts, 1, function(r)paste(r, collapse = '|')));
}

if (0) {
	gts = diplotypes2gts(dts, summarize = sum);
	dtfs = diplotypeFs(ped, dts);
	print(dtfs);
	print(diplotypes2gts(dts));
	print(diplotypes2gts(dts, summarize = sum));

	dtTags = apply(s, 1, function(r)paste(r, collapse = '/'));
	gtTags = apply(gts, 1, function(r)paste(r, collapse = ':'));
	plotPedigrees(ped, tag = paste(dtTags, gtTags, sep = "\n"));

}

simulateFromTemplate = function(pedTemplate, N = 25, hfs = 1:8) {
	ped = familiesFromTemplate(pedTemplate, N = N);
	dts = simulateDiplotypesPed(ped, hfs);
	gts = diplotypes2gts(dts, summarize = sum);
	dtfs = diplotypeFs(ped, dts);
	o = pedForwardOrder(ped);
	pedO = ped[o, ];

	r = list(ped = pedO, dts = dts, gts = gts, dtfs = dtfs, peds = pedSplit2ivTrios(ped));
	r
}

if (1) {
	Nsims = c(1e2, 1e3, 1e4);
	Ns = c(25, 50, 100, 400);
	Nsim = 1e3;
	mses = sapply(Ns, function(N) {
		sim = simulateFromTemplate(pedTemplate, N, hfs = 1:8);
		R = new(M$DiplotypeReconstructor, sim$gts, pedsItrios2rcpp(sim$peds));
		dts = sapply(1:Nsim, function(i) {
			dtsS = R$drawFromHfs(sim$dtfs, runif(length(sim$peds)));
			diplotypeFs(sim$ped, dtsS)
		});
		dtsMean = apply(dts, 1, mean);
		dtsSd = apply(dts, 1, sd);
		mse = mean((dtsMean - sim$dtfs)^2);
		cat(Sprintf('N: %{N}d MSE:%{mse}.2e\n'));
		list(Nsim = Nsim, mse = mse, dtsMean = dtsMean, dtsSd = dtsSd)
	});
	
}
