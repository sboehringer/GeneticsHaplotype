#
#	GeneticsHaplotype.R
#Wed Sep  3 18:18:06 CEST 2014

library('devtools');
library('Rcpp');

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

if (T) {
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
	source('GeneticsHaplotype/R/simulation.R');
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
#	<p> imputation
#

MCMCClass = setRefClass('MCMC',
	fields = list(
		# list of numerics specifying prior distributions
		prior = 'list',
		# start chain with these parameters
		start = 'numeric',
		# compute that many cycles before starting to sample the chain
		Nburnin = 'integer',
		# run the chain for this many cycles
		Nchain = 'integer',
		# store samples of the chain at these intervals
		NsampleSpacing = 'integer',
		# samples from the chain
		chain = 'list'
	),
	methods = list(
	#
	#	<p> methods
	#
	initialize = function(...) {
		chain <<- list();
		.self$initFields(...);
		.self
	},
	update = function() {
		stop('Abstract method called');
	},
	getParameter = function() {
		stop('Abstract method called');
	},
	sample = function() {
		chain <<- c(chain, list(getParameter()));
	},
	run = function() {
		N = Nburnin + Nchain;
		for (i in 1:N) {
			.self$update(i);
			if (i > Nburnin && ((i - Nburnin - 1) %% NsampleSpacing) == 0) .self$sample();
		}
	}
	#
	#	</p> methods
	#
	)
);
MCMCClass$accessors(names(MCMCClass$fields()));

MCMCimputationClass = setRefClass('MCMCimputation', contains = 'MCMC',
	fields = list(
		# Diplotype recontstruction to use
		reconstruction = 'Rcpp_DiplotypeReconstructor',
		peds = 'list',
		state = 'matrix',
		# <p> pre-computed values
		N = 'integer',
		Nhts = 'integer',
		Ncum = 'integer',
		Ifounders = 'integer',
		IfoundersPerFamily = 'list'
	),
	methods = list(
	#
	#	<p> methods
	#
	initialize = function(...) {
		.self$initFields(...);
		# <p> pre-compute
		Nhts <<- as.integer(2^reconstruction$countMarkers());
		# add default prior
		if (is.null(prior$haplotypes)) prior$haplotypes <<- rep(1, Nhts);
		# draw first state
		htfs = rep(1, Nhts);	# <i> draw from Dirichlet
		state <<- reconstruction$drawFromHfs(htfs, runif(length(peds)));

		# <p> pre-compute
		N <<- length(peds);
		Ncum <<- as.integer(c(0L, cumsum(pedsFamilySizes(peds))) + 1L);
		IfoundersPerFamily <<- pedsFounderIdcs(peds);
		Ifounders <<- unlist(IfoundersPerFamily);
		#assign('state0', state[Ifounders, ], envir = .GlobalEnv);
		.self
	},
	#
	# <p> helpers
	#
	freqHat = function(j) {
		htfs = table.n(state[setdiff(Ifounders, IfoundersPerFamily[[j]]), ], min = 0, n = Nhts - 1);
		htfs
	},
	redrawFamily = function(j) {
		# posterior distribution of haplotypes
		htfsPost = freqHat(j) + prior$haplotypes;
		#print(round(vector.std(htfsPost)*36, 1));
		# <A> module indexes from 0
		dtsJ = R$drawFamFromHfs(j - 1, htfsPost, runif(1));
		#if (any(state[Ncum[j]:(Ncum[j + 1] - 1), ] - dtsJ != 0)) browser();

		state[Ncum[j]:(Ncum[j + 1] - 1), ] <<- dtsJ;
		NULL
	},
	getParameter = function() {
		#if (any(state[Ifounders, ] - state0 != 0)) print(which(state[Ifounders, ] - state0 != 0));
		table.n.freq(state[Ifounders, ], min = 0, n = Nhts - 1)
	},
	update = function(i) {
		# iteratively update families
		famI = ((i - 1) %% N) + 1;
		redrawFamily(famI);
	}
	#
	#	</p> methods
	#
	)
);
MCMCimputationClass$accessors(names(MCMCimputationClass$fields()));

#
#	<p> test MCMC genotype imputation
#

if (0) {
	#source('GeneticsHaplotype/R/mcmc.R');
	source('GeneticsHaplotype/R/pedigree.R');
}

freqHat = function(dts, j) {
	htfs = table.n.freq(dts[-j, ], min = 0, n = 7);
	htfs
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

if (1) {
	hfs = vector.std(1:8);
	N = 250;
	d = simulateFromTemplate(pedTemplate1, N = N, hfs = hfs);
	R = new(M$DiplotypeReconstructor, d$gts, pedsItrios2rcpp(d$peds));
}
# debug reconstructions
if (0) {
	r = sapply(1:length(d$peds), function(i)R$reconstructionsFam(i - 1));
	reconstGts = cbind(t(r), matrix(as.vector(d$gts), byrow = T, ncol = 3));
}

if (1) {
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
if (1) {
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
