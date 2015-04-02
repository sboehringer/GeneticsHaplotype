#
#	mcmcLinear.R
#Mon Mar 23 15:06:19 CET 2015


#
#	<p> MCMC linear
#

MCMCLinearClass = setRefClass('MCMCLinear', contains = c('MCMCRegression'),
	fields = list(
		X = 'matrix',	# design matrix
		y = 'numeric'	# response vector
	),
	methods = list(
	#
	#	<p> methods
	#
	initialize = function(..., NpedSplit = 4L) {
		callSuper(..., NpedSplit = NpedSplit);
		.self
	},
	runInitialize = function() {
		callSuper();
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
		callSuper();
		state$beta <<- rnorm(ncol(X) + 1, 0, 1);	# use prior parameters <!>
		state$sigma <<- 1;	# use prior parameters <!>, sigma is variance <!>
		NULL
	},
	blocking = function() {
		# determine blocking
		blockHts = listKeyValue(rep('hts', NpedSplit), splitN(N, NpedSplit));
		blockLinB = list(beta = 1);
		blockLinS = list(sigma = 1);
		lol = list(blockHts, blockLinB, blockLinS);
		mesh = cbind(1:NpedSplit, matrix(1, ncol = length(lol) - 1, nrow = NpedSplit));
		blocking = meshLists(lol, mesh);
	},
	update_beta = function(i) {
		# build design matrix
		Xg <- cbind(X, gtScores[genotypes() + 1]);
		# compute posterior distribution (mean, cov-mat of MVN)
		S1inv = prior$betaVarInv + t(Xg) %*% Xg / state$sigma;
		S1 = solve(S1inv);
		mu = S1 %*% ((prior$betaMuScaled + t(Xg) %*% y) / state$sigma);
		# draw new state
		state$beta <<- mvrnorm(1, mu, S1);
		#print(sprintf('Beta: %.2f', state$beta))
		NULL
	},
	update_sigma = function(i) {
		Xg <- cbind(X, gtScores[genotypes() + 1]);
		res = as.vector(y) - Xg %*% state$beta;
		gishape = prior$gishape + nrow(Xg)/2;
		giscale = prior$giscale + (t(res) %*% res)/2;
		state$sigma <<- 1/rgamma(1, shape = gishape, scale = 1/as.numeric(giscale));
		#print(sprintf('Sigma: %.2f', state$sigma));
		NULL
	},
	llOutcome = function(i, yFam, lpredFam, famSel) {
		ll = dnorm(yFam, lpredFam, sd = state$sigma, log = TRUE);
		ll
	},
	getParameter = function() {
		#if (any(state[Ifounders, ] - state0 != 0)) print(which(state[Ifounders, ] - state0 != 0));
		list(htfs = table.n.freq(state$hts[Ifounders, ], min = 0, n = Nhts - 1),
			beta = state$beta, sigma = state$sigma
		)
	}

	#
	#	</p> methods
	#
	)
);
MCMCLinearClass$accessors(names(MCMCLinearClass$fields()));

MCMCLinearReFamClass = setRefClass('MCMCLinearReFam', contains = c('MCMCRegression'),
	fields = list(
		X = 'matrix',	# design matrix
		y = 'numeric'	# response vector
	),
	methods = list(
	#
	#	<p> methods
	#
	initialize = function(..., NpedSplit = 4L) {
		callSuper(..., NpedSplit = NpedSplit);
		.self
	},
	runInitialize = function() {
		callSuper();
		prior$betaMu <<- rep(0, ncol(X) + 1);			# use prior parameters <!>
		prior$betaVar <<- diag(rep(5, ncol(X) + 1));	# use prior parameters <!>
		prior$giscale <<- 2;
		prior$gishape <<- 5;
		prior$giscaleRe <<- 2;
		prior$gishapeRe <<- 5;

		# <p> precompute prior-derivates
		prior$betaVarM12 <<- matrixM12(prior$betaVar);	# use prior parameters <!>
		prior$betaVarInv <<- solve(prior$betaVar);	# use prior parameters <!>
		prior$betaMuScaled <<- prior$betaVarInv %*% prior$betaMu;	# use prior parameters <!>
		NULL
	},
	drawFromPrior = function() {
		callSuper();
		state$beta <<- rnorm(ncol(X) + 1, 0, 1);	# use prior parameters <!>
		state$sigma <<- 1;	# use prior parameters <!>, sigma is variance <!>
		state$sigmaRe <<- 1;	# use prior parameters <!>, sigmaRe is variance <!>
		update_re(0);	# draw realization of re
		NULL
	},
	blocking = function() {
		Npeds = length(peds);

		# determine blocking
		blockHts = listKeyValue(rep('hts', NpedSplit), splitN(Npeds, NpedSplit));
		blockLinB = list(beta = 1);
		blockLinS = list(sigma = 1);
		blockLinRe = list(re = 1);
		blockLinReS = list(sigmaRe = 1);
		lol = list(blockHts, blockLinB, blockLinS, blockLinRe, blockLinReS);
		mesh = cbind(1:NpedSplit, matrix(1, ncol = length(lol) - 1, nrow = NpedSplit));
		blocking = meshLists(lol, mesh);
	},
	update_beta = function(i) {
		# build design matrix
		Xg <- cbind(X, gtScores[genotypes() + 1]);
		# compute posterior distribution (mean, cov-mat of MVN)
		S1inv = prior$betaVarInv + t(Xg) %*% Xg / state$sigma;
		S1 = solve(S1inv);
		mu = S1 %*% ((prior$betaMuScaled + t(Xg) %*% (y - state$re)) / state$sigma);
		# draw new state
		state$beta <<- mvrnorm(1, mu, S1);
		#print(sprintf('Beta: %.2f', state$beta))
		NULL
	},
	update_sigma = function(i) {
		Xg <- cbind(X, gtScores[genotypes() + 1]);
		res = as.vector(y) - (Xg %*% state$beta + state$re);
		gishape = prior$gishape + nrow(Xg)/2;
		giscale = prior$giscale + (t(res) %*% res)/2;
		state$sigma <<- 1/rgamma(1, shape = gishape, scale = 1/as.numeric(giscale));
		#print(sprintf('Sigma: %.2f', state$sigma));
		NULL
	},
	update_re = function(i) {
		Xg <- cbind(X, gtScores[genotypes() + 1]);
		res = as.vector(y) - (Xg %*% state$beta);
		# <p> compute posterior distribution
		sigmas = sapply(1:N, function(i)1/(1/state$sigmaRe + Nfams[i]/state$sigma));
		means = sapply(1:N, function(i)(sigmas[i] * sum(res[Ifams[[i]]])/state$sigma));
		# <p> draw realizations of random effect per family
		reFam = rnorm(N, means, sqrt(sigmas));
		scoreReNO = unlist(lapply(1:N, function(i)rep(reFam[i], Nfams[i])));
		state$re <<- scoreReNO[IfamsOinv];
		NULL
	},
	update_sigmaRe = function(i) {
		gishapeRe = prior$gishapeRe + length(y)/2;
		giscaleRe = prior$giscaleRe + (state$re %*% state$re)[1, 1]/2;
		state$sigmaRe <<- 1/rgamma(1, shape = gishapeRe, scale = 1/as.numeric(giscaleRe));
		#print(sprintf('SigmaRe: %.2f', state$sigmaRe));
		NULL
	},
	llOutcome = function(i, yFam, lpredFam, famSel) {
		ll = dnorm(yFam, lpredFam + state$re[famSel], sd = state$sigma, log = TRUE);
		ll
	},
	getParameter = function() {
		#if (any(state[Ifounders, ] - state0 != 0)) print(which(state[Ifounders, ] - state0 != 0));
		list(htfs = table.n.freq(state$hts[Ifounders, ], min = 0, n = Nhts - 1),
			beta = state$beta, sigma = state$sigma, sigmaRe = state$sigmaRe
		)
	}

	#
	#	</p> methods
	#
	)
);
MCMCLinearReFamClass$accessors(names(MCMCLinearReFamClass$fields()));

#
#	<p> MCMC based on random effect with correlation structure proportional to coefficients of relationship
#

MCMCLinearReRelClass = setRefClass('MCMCLinearReRel', contains = 'MCMCLinearReFam',
	fields = list(
		# coefficients of relationships
		cors = 'list',
		# coefficients of relationships, precomputed inverses
		corsInv = 'list',
		# coefficients of relationships, precomputed determinants of inverse square
		corsDetM12 = 'numeric'
	),
	methods = list(
	#
	#	<p> methods
	#
	initialize = function(..., cors = NULL) {
		callSuper(..., cors = cors);
		.self
	},
	runInitialize = function() {
		callSuper();
		corsInv <<- lapply(cors, solve);
		corsDetM12 <<- sapply(cors, function(cor)det(matrixM12(cor)));
		NULL
	},
	update_re = function(i) {
		Xg <- cbind(X, gtScores[genotypes() + 1]);
		res = as.vector(y) - (Xg %*% state$beta);
		# <p> compute posterior distribution
		scoreReNO = unlist(lapply(1:N, function(i) {
			Sigma = solve(corsInv[[i]]/state$sigmaRe + diag(rep(1, Nfams[i]))/state$sigma);
			Mu = as.vector(res[Ifams[[i]]] %*% Sigma)/state$sigma;
			mvrnorm(mu = Mu, Sigma = Sigma)
		}));
		state$re <<- scoreReNO[IfamsOinv];
		NULL
	}
	#	</p> methods
	#
	)
);
MCMCLinearReRelClass$accessors(names(MCMCLinearReRelClass$fields()));
