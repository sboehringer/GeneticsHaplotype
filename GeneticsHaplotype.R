#
#	GeneticsHaplotype.R
#Wed Sep  3 18:18:06 CEST 2014

library('devtools');
source('~/.Rprofile');

if (F) {
	system('rm GeneticsHaplotype/src/*.o GeneticsHaplotype/src/*.so');
	Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
	install('GeneticsHaplotype', threads = 6);
}

library('Rcpp');
require('GeneticsHaplotype');
source('GeneticsHaplotype/R/Rdata.R');
source('GeneticsHaplotype/R/mcmc.R');
source('GeneticsHaplotype/R/mcmcRegression.R');
source('GeneticsHaplotype/R/mcmcLinear.R');
source('GeneticsHaplotype/R/mcmcBinomial.R');
source('GeneticsHaplotype/R/pedigree.R');
source('GeneticsHaplotype/R/simulation.R');

if (F) {
	library('GeneticsHaplotype');
	#M = Module('Reconstructor', PACKAGE = 'GeneticsHaplotype');
}
if ( T) {
	pedTemplate = Df(names = c('iid', 'mid', 'pid'), matrix(
		c(	1, NA, NA,
			2, NA, NA,
			3, 1, 2,
			4, NA, NA,
			5, 3, 4
	), byrow = T, ncol = 3));
	pedTemplate1 = Df(names = c('iid', 'mid', 'pid'), matrix(
		c(	1, NA, NA,
			2, NA, NA,
			3, 1, 2
	), byrow = T, ncol = 3));
}
if (F) {
	hfs = rev(vector.std(1:8));
	N = 5e2;
	d = simulateFromTemplate(pedTemplate, N = N, hfs = hfs);
}


if (0) {
	ped1 = list(
		founders = c( 0, 1, 4 ),
		itrios = matrix(c(c(2, 0, 1 ), c(3, 0, 1), c(5, 3, 4)), ncol = 3, byrow = T)
	);
	peds = list(ped1);
	gts = t(matrix(c(c( 1, 1, 1, 1, 1, 1), c(1, 1, 1, 1, 1, 1), c(1, 1, 1, 1, 1, 1)), ncol = 6, byrow = T));
	#R = new(M$DiplotypeReconstructor, gts, peds);
	R = new(DiplotypeReconstructor, gts, peds);
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

if (0) {
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

if (0) {
	Nsims = c(1e2, 1e3, 1e4);
	Ns = c(25, 50, 100, 400);
	Nsim = 1e3;
	hfs = 1:8;
	mses = lapply(Ns, function(N) {
		r = simulationFromPed(familiesFromTemplate(pedTemplate, N), hfs, Nsim);
		mse = mean(r$hfsB^2 + r$hfsSd^2);
		cat(Sprintf('N: %{N}d MSE:%{mse}.2e\n'));
		list(Nsim = Nsim, mse = mse, hfsE = r$hfsE, hfsSd = r$hfsSd)
	});
}

if (0) {
	r = simulationsFromPed(familiesFromTemplate(pedTemplate, 25), 1:8, 1e2, 1e1);
}

# simulations requiring parallelization resources <p>
if (0) {
	#source('RgenericAll.R');
	source('RcomputeResources.R');
	parallelize_initialize(Parallelize_config__, backend = 'snow', force_rerun = T, parallel_count = 8,
		libraries = c('Rcpp', 'GeneticsHaplotype'));
	#parallelize_setEnable(F);
	modelList = list(
		NsimPerPed = c(1e3),
		Npeds = c(25, 50, 100, 400),
		#Npeds = c(25, 50),
		Nsim = 5e2,
		hfs = list(1:8, rep(1, 8), 1/(1:8))
	);
	parF = function(modelList, template) {
		iterateModels(modelList, function(Npeds, Nsim, NsimPerPed, hfs, pedTemplate) {
			M = Module('Reconstructor', PACKAGE = 'GeneticsHaplotype');
			r0 = simulationsFromPed(familiesFromTemplate(pedTemplate, Npeds), hfs, NsimPerPed, Nsim,
				module = M);
			r1 = c(list(Npeds = Npeds), r0);
			unlist(r1)
		}, pedTemplate = template)
	};
	r = parallelize_call(parF(modelList, template = pedTemplate));
	r1 = t(simplify2array(r$results));
	r1 = r1[order(r1[,'hfs1']),];
	print(r1, digits = 3);
	write.csv(r1, file = 'simulations/2014-10-haplotypedrawing.csv');
}

# testing
if (0) {
	Npeds = 1;
	hfs = rep(1, 8);
	ped = familiesFromTemplate(pedTemplate, Npeds);
	test_drawCompatibility(ped, hfs, Nsim = 1e4, M);
}
if (0) {
	source('GeneticsHaplotype/R/pedigree.R');
	ped = Df(names = c('iid', 'fid', 'mid', 'pid'), matrix(
		c(	10, 1, NA, NA,
			20, 1, NA, NA,
			30, 1, 10, 20,
			40, 1, NA, NA,
			55, 1, 30, 40
		)
	, byrow = T, ncol = 4));
	# genotypes by row, individuals by column 0..5
	gtsRaw = c(
		2, 2, 2, 2, 2,
		0, 0, 0, 0, 0,
		1, 1, 1, 1, 1);
	gtsRaw = c(
		1, 2, 2, 2, 2,
		0, 1, 1, 0, 1,
		1, 1, 1, 1, 1);

	gts = t(matrix(gtsRaw, byrow = T, ncol = 5));
	#sim = simulateFromPed(ped, hfs = 1:8);
	peds = pedsItrios2rcpp(pedSplit2ivTrios(ped));
	M = Module('Reconstructor', PACKAGE = 'GeneticsHaplotype');
	R = new(M$DiplotypeReconstructor, gts, peds);
	sapply(seq(.01, .99, length.out = 5e0), function(u) {
		dtsS = R$drawFromHfs(rep(1, 8), u);
		print(dtsS);
		gtsR = diplotypes2gtsDose(dtsS);
		print(all(gts == gtsR));
		if (any(gts != gtsR)) {
			print(gtsR);
		}
	});
}


#
#	<p> test MCMC genotype imputation
#

if (0) {
	#source('GeneticsHaplotype/R/mcmc.R');
	source('GeneticsHaplotype/R/pedigree.R');
}

freqHat = function(dts, j) {
	htfs = table.n.freq(dts[-j, ], min = 0, n = 7);
	htfsF
}

redrawFamily = function(dts, j, ped) {
	dtsJ = R$drawFamFromHfs(j, freqHat(dts, j), runif(1));
	r = dtsJ[ped$founders, ];
	r
}

if (0) {
	i = pedFounderIdcsForward(d$ped);
	Ns = c(1, cumsum(pedFounderSizes(ped)) + 1);
	dts = R$drawFromHfs(1:8, runif(length(d$peds)))[i, ];
	dtTab = table.n(as.vector(dts), min = 0, n = 7);
	j = 1;
	print(freqHat(dts, j));
	dtsI = R$drawFamFromHfs(j, freqHat(dts, j), runif(1));
	print(freqHat(dts, 2));
	dtsJ = redrawFamily(dts, j, d$peds[[j]]);
	dts[Ns[j]:(Ns[j + 1] -1), ] = dtsJ;
}

# debug reconstructions
if (0) {
	r = sapply(1:length(d$peds), function(i)R$reconstructionsFam(i - 1));
	reconstGts = cbind(t(r), matrix(as.vector(d$gts), byrow = T, ncol = 3));
}

if (0) {
	R = new(DiplotypeReconstructor, d$gts, pedsItrios2rcpp(d$peds));
	mcmc = MCMCimputationClass$new(
		reconstruction = R,
		peds = pedSplit2ivTrios(d$ped),
		Nburnin = 0L,
		Nchain = 1e5L,
		NsampleSpacing = length(d$peds),
		prior = list(haplotypes = 1:4)
	);
	mcmc$run();
	chain = sapply(mcmc$chain, identity);
}
if (0) {
	require(ggplot2);
	require(grid);
	require(gridExtra);

	ps = lapply(1:nrow(chain), function(i, hfs = NULL) {
		df = data.frame(par = chain[i, ], it = 1:ncol(chain));
		p = ggplot() + geom_line(data = df, aes(it, par)) + 
			scale_y_reverse() + 
			theme_bw();
		if (!is.null(hfs)) p = p + geom_hline(aes(yintercept = hfs), alpha = .5);
		ggplot_gtable(ggplot_build(p))
	}, hfs = hfs);
	do.call(grid.arrange, c(ps, list(ncol = 1, nrow = length(ps))));

}

if (0) {
	#source('GeneticsHaplotype/R/mcmc.R');
	source('GeneticsHaplotype/R/simulation.R');
	gts = c(0:2, 2:0);
	pa = simulatePhenotypesBinRaw(gts, c(-3, 4), scoreGt = scoresL$additive);
	print(pa);
	pg = simulatePhenotypesBinRaw(gts, c(-3, 4, 2), scoreGt = scoresL$genotype);
	print(pg);
}

if (0) {
	source('GeneticsHaplotype/R/Rdata.R');
	source('GeneticsHaplotype/R/mcmc.R');
	b = MCMCBlockClass$new();
	MCMCBlockClass$methods();
	#b$hello();	# <!> mb called before b[['hello']]() works
	activateMethods(b, 'hello');
	b[['hello']]();
	
}

if (0) {
	source('GeneticsHaplotype/R/Rdata.R');
	source('GeneticsHaplotype/R/mcmc.R');
	Nhts = 4;
	hts = listKeyValue(rep('hts', Nhts), splitN(15, Nhts));
	linB = list(beta = 1);
	linS = list(sigma = 1);
	lol = list(hts, linB, linS);
	mesh = matrix(c(1:length(hts), rep(1, Nhts), rep(1, Nhts)), ncol = length(lol));
	r = meshLists(lol, mesh);
	print(r);
}

if (0) {
	source('GeneticsHaplotype/R/Rdata.R');
	source('GeneticsHaplotype/R/mcmc.R');

	Nhts = 4;
	hts = listKeyValue(rep('hts', Nhts), splitN(15, Nhts));
	linB = list(beta = 1);
	lol = list(hts, linB);
	mesh = matrix(c(1:length(hts), rep(1, Nhts)), ncol = length(lol));
	r = meshLists(lol, mesh);

	b = MCMCBlockClass$new(blocking = r, Nburnin = 1e1L, Nchain = 4e1L);
	b$run();
}

if (0) {
	require('GeneticsHaplotype');
	R = new(DiplotypeReconstructor, d$gts, pedsItrios2rcpp(d$peds));
	reconstructions = R$reconstructionsAll();
	gts = lapply(reconstructions, function(m) {
		# alleles
		as = apply(m[, -1, drop = F], c(1, 2), function(e)e%%2);
		t(apply(as, 1, function(r)apply(matrix(r, byrow = T, ncol = 2), 1, sum)))
	});
}

if (0) {
	source('GeneticsHaplotype/R/Rdata.R');
	source('GeneticsHaplotype/R/simulation.R');
	require('GeneticsHaplotype');
	# get reconstructions for debuggin
	R = new(DiplotypeReconstructor, d$gts, pedsItrios2rcpp(d$peds));
	reconstructions = R$reconstructionsAll();
	test_reconstruction(d, reconstructions);

	# simulate
	y = simulatePhenotypesLinear(d$gts[, 1], c(0, 7), sd = 2)[, 1];
	#X = model.matrix(~ gts, data.frame(gts = d$gts[, 1]));
	X = model.matrix(~ 1, data.frame(dummy = rep(1, length(y))));
	R = new(DiplotypeReconstructor, d$gts, pedsItrios2rcpp(d$peds));
	reconstructions = R$reconstructionsAll();

	# chain
	mcmcLin = new('MCMCLinear', y = y, X = X, peds = d$peds, reconstructor = R,
		Nburnin = 1e3L, Nchain = 1e5L);
	mcmcLin$run();
}

# chain is matrix with parameter per row
plotChain = function(chain, parSim = NULL) {
	require(ggplot2);
	require(grid);
	require(gridExtra);

	ps = lapply(1:nrow(chain), function(i) {
		df = data.frame(par = chain[i, ], it = 1:ncol(chain), yinter = parSim[i]);
		p = ggplot() + geom_line(data = df, aes(it, par)) + 
			#scale_y_reverse() + 
			 scale_y_continuous(name = row.names(chain)[i]) + theme_bw();
		if (!is.null(parSim)) {
			p = p + geom_line(data = df, aes(it, yinter), alpha = .5);
		}
		ggplot_gtable(ggplot_build(p))
	});
	p = do.call(grid.arrange, c(ps, list(ncol = 2, nrow = ceiling(length(ps)/2))));
	p
}

if (0) {
	pars = sapply(mcmcLin$chain, unlist)
	pars1 = pars[, -(1:5)];
	p = plotChain(pars, parSim = c(vector.std(d$dtfs), c(0, 7), 4));
	p1 = plotChain(pars1, parSim = c(vector.std(d$dtfs), c(0, 7), 4));
}


# 9.3.2015: debugging
if (F) save(d, file = 'debugging/reconstructionError.Rdata');
if (0) {
	#load('debugging/reconstructionError.Rdata')
	#reconstructions = R$reconstructionsAll();
	gts64 = d$gts[pedsIdcs(d$peds)[[64]], ];
	plotPedigree(d$ped[pedsIdcs(d$peds)[[64]],], tag = apply(gts64, 1, function(r)paste(r, collapse = ',')))
}



if (0) {
	print(matrixM12(diag(rep(3, 3))))
}

#
#	<p> RE model
#

if (0) {
	# get reconstructions for debuggin
	R = new(DiplotypeReconstructor, d$gts, pedsItrios2rcpp(d$peds));
	reconstructions = R$reconstructionsAll();
	test_reconstruction(d, reconstructions);

	# simulate
	y = simulatePhenotypesLinearReFam(pedsIdcs(d$peds), d$gts[, 1], c(0, 7), sd = 2, sdRe = 1)[, 1];
	#X = model.matrix(~ gts, data.frame(gts = d$gts[, 1]));
	X = model.matrix(~ 1, data.frame(dummy = rep(1, length(y))));

	# chain
	mcmcLinRe = new('MCMCLinearReFam',
		y = y, X = X, peds = d$peds, reconstructor = R,
		Nburnin = 1e3L, Nchain = 1e5L
	);
	mcmcLinRe$run();
}
if (0) {
	pars = sapply(mcmcLinRe$chain, unlist)
	pars1 = pars[, -(1:5)];
	p1 = plotChain(pars1, parSim = c(vector.std(d$dtfs), c(0, 7), 4, 1));
}

if (0) {
	ped = peds[[1]];
	print(pedAncestry(ped$itrios));
	print(pedCoeffOfRel(ped));
}
# test cases for pedCoeffOfRel
if (0) {
	tests = list(
		# three generations, stacked trios
		c(	10, 1, NA, NA,
			20, 1, NA, NA,
			30, 1, 10, 20,
			40, 1, NA, NA,
			55, 1, 30, 40),
		# sibs
		c(	10, 1, NA, NA,
			20, 1, NA, NA,
			30, 1, 10, 20,
			40, 1, 10, 20,
			55, 1, 10, 20),
		# cousins
		c(	10, 1, NA, NA,
			11, 1, NA, NA,
			20, 1, 10, 11,
			21, 1, NA, NA,
			22, 1, 10, 11,
			23, 1, NA, NA,
			30, 1, 20, 21,
			31, 1, 22, 23)
	);
	cors = lapply(tests[2:3], function(test) {
		ped = Df(names = c('iid', 'fid', 'mid', 'pid'), matrix(test, byrow = T, ncol = 4));
		ped2 = pedSplit2ivTrios(ped);
		print(ped2[[1]]);
		pedCoeffOfRel(ped2[[1]])
	});
	print(cors);
}

if (0) {
	# Multiple inheritance conclusion: initialize must not expect well defined instance variables
	#	-> limit initialize to basic initialization
	#	-> other initialization need to be put into other methods
AClass = setRefClass('A',
	fields = list( a = 'integer' ),
	methods = list(
	initialize = function(...) {
		print('called A');
		callSuper(...);
		.self
	}
	)
);
BClass = setRefClass('B',
	fields = list( b = 'integer' ),
	methods = list(
	initialize = function(...) {
		print('called B');
		callSuper(...);
		.self
	}
	)
);
CClass = setRefClass('C', contains = c('A', 'B'),
	methods = list(
	initialize = function(...) {
		callSuper(...);
		.self
	}
	)
);
DClass = setRefClass('D', contains = 'C',
	methods = list(
	initialize = function(...) {
		callSuper(...);
		.self
	}
	)
);
	Dinst = new('D', a = 10L, b = 20L);
}

#
#	<p> RE model coefficients of relationship
#

if (0) {
	# get reconstructions for debuggin
	if (F) {
	R = new(DiplotypeReconstructor, d$gts, pedsItrios2rcpp(d$peds));
	reconstructions = R$reconstructionsAll();
	test_reconstruction(d, reconstructions);
	cors = pedsCoeffOfRel(d$peds);
	}

	# simulate
	y = simulatePhenotypesLinearReRel(pedsIdcs(d$peds), d$gts[, 1], c(0, 7), sd = 2,
		sdRe = 1, cors = cors)[, 1];
	#X = model.matrix(~ gts, data.frame(gts = d$gts[, 1]));
	X = model.matrix(~ 1, data.frame(dummy = rep(1, length(y))));
	# chain
	mcmcLinReRel = new('MCMCLinearReRel',
		y = y, X = X, peds = d$peds, reconstructor = R,
		Nburnin = 1e3L, Nchain = 5e5L, cors = cors
	);
	mcmcLinReRel$run();
}

if (0) {
	pars = sapply(mcmcLinReRel$chain, unlist)
	pars1 = pars[, -(1:5)];
	p1 = plotChain(pars1, parSim = c(vector.std(d$dtfs), c(0, 7), 4, 1));
}

#
#	<p> logistic regression
#

if (0) {
	# get reconstructions for debuggin
	if (F) {
	R = new(DiplotypeReconstructor, d$gts, pedsItrios2rcpp(d$peds));
	reconstructions = R$reconstructionsAll();
	}

	# simulate
	y = simulatePhenotypesBin(d$gts[, 1], c(-.5, .5));
	#X = model.matrix(~ gts, data.frame(gts = d$gts[, 1]));
	X = model.matrix(~ 1, data.frame(dummy = rep(1, length(y))));
	# chain
	mcmcBin = new('MCMCBinomial',
		y = y, X = X, peds = d$peds, reconstructor = R,
		Nburnin = 1e3L, Nchain = 5e5L
	);
	mcmcBin$run();
}

qqDist = function(Nqts = 1e2, qdist, ...) {
	qtls = (1:Nqtls)/(Nqtls + 1);
	qtlsExp = qdist(qtls, ...);
	qtlsObs = quantile(rtn, qtls);
	qq = qplot(qtlsExp, qtlsObs) + theme_bw();
	qq
}

qqSim = function(Nsim, dist = 'truncnorm', Nqts = Nsim/10, ...) {
	rdist = get(Sprintf('r%{dist}s'));
	r = rdist(Nsim, ...);
	qdist = get(Sprintf('q%{dist}s'));
	qq = qqDist(Nqts, qdist, ...);
	qq
}

if (0) {
	#rtn = rtruncnorm(1e4, 0, 1, 0, 1);
	#qqRtn = qqDist(1e3, qtruncnorm, mean = 0, sd = 1, lower = 0, upper = 1);
	qqRtn = qqSim(1e4, 'truncnorm', mean = 0, sd = 1, lower = 0, upper = 1);
}

#
#	<p> probit regression
#

if (0) {
	# get reconstructions for debuggin
	if (F) {
	R = new(DiplotypeReconstructor, d$gts, pedsItrios2rcpp(d$peds));
	reconstructions = R$reconstructionsAll();
	}

	# simulate
	y = simulatePhenotypesBin(d$gts[, 1], c(-1, 1), link = pnorm);
	#X = model.matrix(~ gts, data.frame(gts = d$gts[, 1]));
	X = model.matrix(~ 1, data.frame(dummy = rep(1, length(y))));
	# chain
	mcmcBinProbit = new('MCMCBinProbit',
		y = y, X = X, peds = d$peds, reconstructor = R,
		Nburnin = 1e3L, Nchain = 5e5L
	);
	mcmcBinProbit$run();
}

if (1) {
	# get reconstructions for debuggin
	if (F) {
	R = new(DiplotypeReconstructor, d$gts, pedsItrios2rcpp(d$peds));
	reconstructions = R$reconstructionsAll();
	}

	# simulate
	y = simulatePhenotypesBinReFam(pedsIdcs(d$peds), d$gts[, 1], c(-1, 1), sdRe = 2, link = pnorm);
	#X = model.matrix(~ gts, data.frame(gts = d$gts[, 1]));
	X = model.matrix(~ 1, data.frame(dummy = rep(1, length(y))));
	# chain
	mcmcBinProbitReFam = new('MCMCBinProbitReFam',
		y = y, X = X, peds = d$peds, reconstructor = R,
		Nburnin = 1e3L, Nchain = 5e5L
	);
	mcmcBinProbitReFam$run();
}

