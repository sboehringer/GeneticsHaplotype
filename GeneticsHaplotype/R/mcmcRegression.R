#
#	mcmcRegression.R
#Wed Mar 25 15:57:28 CET 2015

# base class for linear and logistic regression
# implements haplotype reconstruction and precompuation of data based on reconstructions

MCMCRegressionClass = setRefClass('MCMCRegression', contains = c('MCMCBlock', 'HaplotypeHelper'),
	fields = list(
		# Rcpp module to compute reconstructions
		reconstructor = 'envRefClass',
		# Diplotype recontstruction to use
		reconstructions = 'list',
		# risk genotypes for each reconstruction
		reconstructionsGts = 'list',
		peds = 'list',
		state = 'list',
		Nloci = 'integer',
		gtScores = 'numeric',
		NpedSplit = 'integer'
	),
	methods = list(
	#
	#	<p> methods
	#
	# precompute genotypes for all reconstructions -> haplotype drawing
	genotypesPrecompute = function() {
		reconstructionsGts <<- lapply(reconstructions, function(m) {
			# alleles
			as = apply(m[, -1, drop = F], c(1, 2), function(e)e%%2);
			t(apply(as, 1, function(r)apply(matrix(r, byrow = T, ncol = 2), 1, sum)))
		});
		#apply(state$hts, 1, function(hts)(hts[1] %% 2 + hts[2] %% 2))
		NULL
	},
	initialize = function(..., NpedSplit = 4L) {
		callSuper(..., NpedSplit = NpedSplit);
		.self
	},
	runInitialize = function() {
		callSuper();
		# determine number of loci, reconstructions
		reconstructions <<- R$reconstructionsAll();
		Nloci <<- as.integer(log2(max(unlist(reconstructions)) + 1));
		Npeds = length(peds);

		# Haplotype Helper and other initialization
		initialize_cache();
		genotypesPrecompute();
		
		# priors
		if (is.null(prior$hts)) prior$hts <<- rep(1, 2^Nloci);
		
		# initial state (chain)
		htfs = rep(1, Nhts);	# <i> draw from Dirichlet
		gtScores <<- scoresL$additive;
		NULL
	},
	drawFromPrior = function() {
		htfs = rep(1, Nhts);	# <i> draw from Dirichlet
		state$hts <<- R$drawFromHfs(htfs, runif(length(peds)));
		NULL
	},
	getCountMarkers = function()Nloci,
	# by convention we regress on the locus 0, corresponding to index 1 in R
	# loci have to be rearranged in advance to follow this convention if marker order differs
	# locus number corresponds to binary position in hts
	genotypes = function(i = NULL) {
		myHts = if (is.null(i)) state$hts else state$hts[Ifams[[i]], , drop = F];
		apply(myHts, 1, function(hts)(hts[1] %% 2 + hts[2] %% 2))
	},
	#
	# update family i %% N
	#	compute conditional outcome likelihood
	#	compute likelihood of complete families
	#	draw from joint family distribution
	# 
	update_hts = function(i) {
		# preparation
		N <<- length(peds);
		Ncum <<- as.integer(c(0L, cumsum(pedsFamilySizes(peds))) + 1L);
		iF = (i - 1) %% N + 1;
		famSel = Ncum[iF]:(Ncum[iF + 1] - 1);

		# haplotype frequencies
		NhtsI = freqHat(state$hts, iF);
		logPhts = log(vector.std(NhtsI + prior$hts));
		Ifdrs = peds[[iF]]$founders;	#relative indeces

		# log-probs reconstructions
		Preconstructions = sapply(1:nrow(reconstructions[[iF]]), function(k) {
			# linear predictor
			E = cbind(X[famSel, , drop = F], gtScores[reconstructionsGts[[iF]][k, ] + 1]) %*%
				state$beta;
			# likelihood outcome
			llPts = sum(llOutcome(iF, y[famSel], E[, 1], famSel));
			# likelihood haplotyeps
			logFactor = reconstructions[[iF]][k, 1] * log(2);
			htsFounders = matrix(reconstructions[[iF]][k, -1], byrow = T, nrow = 2)[, Ifdrs, drop = F];
			llHts = logFactor + sum(logPhts[as.vector(htsFounders + 1)]);
			ll = llPts + llHts;
			ll
		});
		# draw family haplotypes as a block
		draw = 1:nrow(reconstructions[[iF]]) %*% rmultinomLog(1, 1, Preconstructions);
		htsI = matrix(reconstructions[[iF]][draw, -1], byrow = T, ncol = 2);
		state$hts[Ifams[[iF]], ] <<- htsI;
		NULL
	},
	llOutcome = function(i, famSel, yFam, lpredFam) {
		stop('abstract method llOutcome@MCMCregression called');
	},
	# assume state$beta exists <N>
	getParameter = function() {
		#if (any(state[Ifounders, ] - state0 != 0)) print(which(state[Ifounders, ] - state0 != 0));
		list(htfs = table.n.freq(state$hts[Ifounders, ], min = 0, n = Nhts - 1),
			beta = state$beta
		)
	}

	#
	#	</p> methods
	#
	)
);
MCMCRegressionClass$accessors(names(MCMCRegressionClass$fields()));