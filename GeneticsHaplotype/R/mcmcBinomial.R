#
#	mcmcBionmial.R
#Mon Mar 23 15:15:40 CET 2015

library('Matrix');

#
#	<p> MCMC logistic regression
#

MCMCBinomialClass = setRefClass('MCMCBinomial', contains = c('MCMCBlock', 'HaplotypeHelper'),
	fields = list(
		# Rcpp module to compute reconstructions
		reconstructor = 'envRefClass',
		# Diplotype recontstruction to use
		reconstructions = 'list',
		# risk genotypes for each reconstruction
		reconstructionsGts = 'list',
		peds = 'list',
		state = 'list',
		X = 'matrix',	# design matrix
		y = 'integer',	# response vector
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
		state$hts <<- R$drawFromHfs(htfs, runif(length(peds)));
		state$beta <<- rnorm(ncol(X) + 1, 0, 1);	# use prior parameters <!>
		state$sigma <<- 1;	# use prior parameters <!>, sigma is variance <!>
		prior$betaMu <<- rep(0, ncol(X) + 1);			# use prior parameters <!>
		prior$betaVar <<- diag(rep(5, ncol(X) + 1));	# use prior parameters <!>
		prior$giscale <<- 2;
		prior$gishape <<- 5;

		# <p> precompute prior-derivates
		prior$betaVarM12 <<- matrixM12(prior$betaVar);	# use prior parameters <!>
		prior$betaVarInv <<- solve(prior$betaVar);	# use prior parameters <!>
		prior$betaMuScaled <<- prior$betaVarInv %*% prior$betaMu;	# use prior parameters <!>
		NULL
	},
	drawFromPrior = function() {
		htfs = rep(1, Nhts);	# <i> draw from Dirichlet
		state$hts <<- R$drawFromHfs(htfs, runif(length(peds)));
		state$beta <<- rnorm(ncol(X) + 1, 0, 1);	# use prior parameters <!>
		state$scaling <<- rep(1, length(y));	# use prior parameters <!>
		state$liability <<- rep(0, length(y));	# use prior parameters <!>
		NULL
	},
	blocking = function() {
		# determine blocking
		blockHts = listKeyValue(rep('hts', NpedSplit), splitN(N, NpedSplit));
		blockBinB = list(beta = 1);
		blockBinLiab = list(liability = 1);
		blockBinScal = list(scaling = 1);
		lol = list(blockHts, blockBinB, blockBinLiab, blockBinScal);
		mesh = cbind(1:NpedSplit, matrix(1, ncol = length(lol) - 1, nrow = NpedSplit));
		blocking = meshLists(lol, mesh);
	},
	getCountMarkers = function()Nloci,
	# by convention we regress on the locus 0, corresponding to index 1 in R
	# loci have to be rearranged in advance to follow this convention if marker order differs
	# locus number corresponds to binary position in hts
	genotypes = function(i = NULL) {
		myHts = if (is.null(i)) state$hts else state$hts[Ifams[[i]], , drop = F];
		apply(myHts, 1, function(hts)(hts[1] %% 2 + hts[2] %% 2))
	},
	update_beta = function(i) {
		# build design matrix
		Xg <- cbind(X, gtScores[genotypes() + 1]);
		W = Diagonal(x = state$scaling);
		# scaled X
		Xs = t(Xg) %*% W;
		# compute posterior distribution (mean, cov-mat of MVN)
		S1inv = (prior$betaVarInv + Xs %*% Xg);
		S1 = solve(S1inv);
		mu = S1 %*% (prior$betaMuScaled + Xs %*% state$liability);
		# draw new state
		state$beta <<- mvrnorm(1, mu, S1)[, 1];
		print(sprintf('Beta: %.2f', state$beta))
		NULL
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
			# likelihood phenotypes
			p = expit(E);
			llPts = sum(log(ifelse(y[famSel], p, 1 - p)));
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
	# latent liability
	update_liability = function(i) {
		lpred = cbind(X, gtScores[genotypes() + 1]) %*% state$beta;
		liab = rnorm(length(lpred), lpred, 1/state$scaling);
		# truncate normal according to outcome
		#liab = ifelse(y == 1, abs(liab), -abs(liab));
		liab = ifelse(y == 1, ifelse(liab > 0, liab, 0), ifelse(liab < 0, liab, 0));
stem(liab);
		state$liability <<- liab;
		NULL
	},
	update_scaling = function(i) {
		nu = 7.3;	# t-distribution based approximation of the logistic (Kinney & Dunson 2006)
		lpred = cbind(X, gtScores[genotypes() + 1]) %*% state$beta;
		shape = (nu + 1)/2;
		scale = 2/(nu + (state$liability - lpred)^2);
		state$scaling <<- rgamma(length(lpred), shape, scale = scale);
		#stem(state$scaling);
		NULL
	},
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
MCMCBinomialClass$accessors(names(MCMCBinomialClass$fields()));

