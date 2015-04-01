#
#	mcmcBionmial.R
#Mon Mar 23 15:15:40 CET 2015

library('Matrix');

#
#	<p> helper functions
#

qtruncnorm = function(p, mean = 0, sd = 1, lower = -Inf, upper = Inf) {
	plow = pnorm(lower, mean, sd);
	pupp = pnorm(upper, mean, sd);
	qnorm(p * (pupp - plow) + plow, mean, sd = sd)
}
rtruncnorm = function(N, mean = 0, sd = 1, lower = -Inf, upper = Inf) {
	qtruncnorm(runif(N), mean, sd = sd, lower, upper)
}


#
#	<p> MCMC logistic regression
#

MCMCBinomialClass = setRefClass('MCMCBinomial', contains = c('MCMCRegression'),
	fields = list(
		X = 'matrix',	# design matrix
		y = 'integer'	# response vector
	),
	methods = list(
	#
	#	<p> methods
	#
	initialize = function(...) {
		callSuper(...);
		.self
	},
	runInitialize = function() {
		callSuper();
		# initial state (chain)
		prior$betaMu <<- rep(0, ncol(X) + 1);			# use prior parameters <!>
		prior$betaVar <<- diag(rep(5, ncol(X) + 1));	# use prior parameters <!>
		prior$giscale <<- 2;
		prior$gishape <<- 5;

		# <p> precompute prior-derivates
		prior$betaVarM12 <<- matrixM12(prior$betaVar);	# use prior parameters <!>
		prior$betaVarInv <<- solve(prior$betaVar);	# use prior parameters <!>
		prior$betaMuScaled <<- prior$betaVarInv %*% prior$betaMu;	# use prior parameters <!>
		prior$nuRange <<- 1:10;
		prior$nu <<- rep(1, length(prior$nuRange));
		prior$nuLog <<- log(prior$nu);
		NULL
	},
	drawFromPrior = function() {
		callSuper();
		state$beta <<- rnorm(ncol(X) + 1, 0, 1);	# use prior parameters <!>
		state$liabilityVar <<- rep(1, length(y));	# use prior parameters <!>
		state$liability <<- rep(0, length(y));	# use prior parameters <!>
		state$nu <<- 7;	# t-distribution based approximation of the logistic (Kinney & Dunson 2006)
		NULL
	},
	blocking = function() {
		# determine blocking
		blockHts = listKeyValue(rep('hts', NpedSplit), splitN(N, NpedSplit));
		blockBinLiab = list(liability = 1);
		blockBinLiabVar = list(liabilityVar = 1);
		blockBinB = list(beta = 1);
		blockBinNu = list(nu = 1);
		lol = list(blockHts, blockBinLiab, blockBinLiabVar, blockBinB, blockBinNu);
		mesh = cbind(1:NpedSplit, matrix(1, ncol = length(lol) - 1, nrow = NpedSplit));
		blocking = meshLists(lol, mesh);
		blocking
	},
	update_beta = function(i) {
		# build design matrix
		Xg <- cbind(X, gtScores[genotypes() + 1]);
		W = Diagonal(x = 1/state$liabilityVar);
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
	llOutcome = function(i, yFam, lpredFam, famSel) {
		p = expit(lpredFam);
		ll = log(ifelse(yFam, p, 1 - p));
		ll
	},
	# latent liability
	update_liability = function(i) {
		lpred = cbind(X, gtScores[genotypes() + 1]) %*% state$beta;
		liab = rtruncnorm(length(lpred), lpred, sqrt(state$liabilityVar),
			lower = ifelse(y == 1, 0, -Inf), upper = ifelse(y == 1, Inf, 0)
		);
		# truncate normal according to outcome
		#liab = ifelse(y == 1, ifelse(liab > 0, liab, 0), ifelse(liab < 0, liab, 0));
#print('LiabilityVar');
#stem(state$liabilityVar);
#stem(liab);
#by(liab, y, stem)


		state$liability <<- liab;
		NULL
	},
	update_liabilityVar = function(i) {
		lpred = (cbind(X, gtScores[genotypes() + 1]) %*% state$beta)[, 1];
		shape = (state$nu + 1)/2;
		scaleInv = (state$nu + (state$liability - lpred)^2)/2;
#by(scaleInv, y, stem)
		state$liabilityVar <<- 1/rgamma(length(lpred), shape, scale = 1/scaleInv);
		NULL
	},
	update_nu = function(i) {
		psLog = prior$nuLog + sapply(prior$nuRange, function(nu)
			sum(-lgamma(nu/2) -log(nu/2)*(nu/2) -log(state$liabilityVar)*(nu/2 - 1) - nu * state$liabilityVar/2));
		state$nu <<- (prior$nuRange %*% rmultinom(1, 1, exp(psLog - max(psLog))))[1, 1];
		print(sprintf('Nu: %.5f', state$nu))
		#print( exp(psLog - max(psLog)));
		NULL
	}
	#
	#	</p> methods
	#
	)
);
MCMCBinomialClass$accessors(names(MCMCBinomialClass$fields()));


#
#	<p> MCMC probit regression
#

MCMCBinProbitClass = setRefClass('MCMCBinProbit', contains = c('MCMCRegression'),
	fields = list(
		X = 'matrix',	# design matrix
		y = 'integer'	# response vector
	),
	methods = list(
	#
	#	<p> methods
	#
	initialize = function(...) {
		callSuper(...);
		.self
	},
	runInitialize = function() {
		callSuper();
		# initial state (chain)
		prior$betaMu <<- rep(0, ncol(X) + 1);			# use prior parameters <!>
		prior$betaVar <<- diag(rep(5, ncol(X) + 1));	# use prior parameters <!>

		# <p> precompute prior-derivates
		prior$betaVarM12 <<- matrixM12(prior$betaVar);	# use prior parameters <!>
		prior$betaVarInv <<- solve(prior$betaVar);	# use prior parameters <!>
		prior$betaMuScaled <<- prior$betaVarInv %*% prior$betaMu;	# use prior parameters <!>
		NULL
	},
	drawFromPrior = function() {
		callSuper();
		state$beta <<- rnorm(ncol(X) + 1, 0, 1);	# use prior parameters <!>
		state$liability <<- rep(0, length(y));	# use prior parameters <!>
		NULL
	},
	blocking = function() {
		# determine blocking
		blockHts = listKeyValue(rep('hts', NpedSplit), splitN(N, NpedSplit));
		blockBinLiab = list(liability = 1);
		blockBinB = list(beta = 1);
		lol = list(blockHts, blockBinLiab, blockBinB);
		mesh = cbind(1:NpedSplit, matrix(1, ncol = length(lol) - 1, nrow = NpedSplit));
		blocking = meshLists(lol, mesh);
		blocking
	},
	update_beta = function(i) {
		# build design matrix
		Xg <- cbind(X, gtScores[genotypes() + 1]);
		# compute posterior distribution (mean, cov-mat of MVN)
		S1inv = (prior$betaVarInv + t(Xg) %*% Xg);
		S1 = solve(S1inv);
		mu = S1 %*% (prior$betaMuScaled + t(Xg) %*% state$liability);
		# draw new state
		state$beta <<- mvrnorm(1, mu, S1);
		print(sprintf('Beta: %.2f', state$beta))
		NULL
	},
	llOutcome = function(i, yFam, lpredFam, famSel) {
		p = pnorm(lpredFam);
		ll = log(ifelse(yFam, p, 1 - p));
		ll
	},
	# latent liability
	update_liability = function(i) {
		lpred = cbind(X, gtScores[genotypes() + 1]) %*% state$beta;
		liab = rtruncnorm(length(lpred), lpred, 1,
			lower = ifelse(y == 1, 0, -Inf), upper = ifelse(y == 1, Inf, 0)
		);
		state$liability <<- liab;
		NULL
	}
	#
	#	</p> methods
	#
	)
);
MCMCBinProbitClass$accessors(names(MCMCBinProbitClass$fields()));

MCMCBinProbitReFamClass = setRefClass('MCMCBinProbitReFam', contains = c('MCMCRegression'),
	fields = list(
		X = 'matrix',	# design matrix
		y = 'integer'	# response vector
	),
	methods = list(
	#
	#	<p> methods
	#
	initialize = function(...) {
		callSuper(...);
		.self
	},
	runInitialize = function() {
		callSuper();
		# initial state (chain)
		prior$betaMu <<- rep(0, ncol(X) + 1);			# use prior parameters <!>
		prior$betaVar <<- diag(rep(5, ncol(X) + 1));	# use prior parameters <!>

		# <p> precompute prior-derivates
		prior$betaVarM12 <<- matrixM12(prior$betaVar);	# use prior parameters <!>
		prior$betaVarInv <<- solve(prior$betaVar);	# use prior parameters <!>
		prior$betaMuScaled <<- prior$betaVarInv %*% prior$betaMu;	# use prior parameters <!>
		prior$giscaleRe <<- 2;
		prior$gishapeRe <<- 5;
		NULL
	},
	drawFromPrior = function() {
		callSuper();
		state$beta <<- rnorm(ncol(X) + 1, 0, 1);	# use prior parameters <!>
		state$liability <<- rep(0, length(y));	# use prior parameters <!>
		state$sigmaRe <<- 1;	# use prior parameters <!>, sigmaRe is variance <!>
		update_re(0);	# draw realization of re
		NULL
	},
	blocking = function() {
		# determine blocking
		blockHts = listKeyValue(rep('hts', NpedSplit), splitN(N, NpedSplit));
		blockBinLiab = list(liability = 1);
		blockBinB = list(beta = 1);
		blockBinRe = list(re = 1);
		blockBinReS = list(sigmaRe = 1);
		lol = list(blockHts, blockBinLiab, blockBinB, blockBinRe, blockBinReS);
		mesh = cbind(1:NpedSplit, matrix(1, ncol = length(lol) - 1, nrow = NpedSplit));
		blocking = meshLists(lol, mesh);
		blocking
	},
	update_beta = function(i) {
		# build design matrix
		Xg <- cbind(X, gtScores[genotypes() + 1]);
		# compute posterior distribution (mean, cov-mat of MVN)
		S1inv = (prior$betaVarInv + t(Xg) %*% Xg);
		S1 = solve(S1inv);
		mu = S1 %*% (prior$betaMuScaled + t(Xg) %*% (state$liability + state$re));
		# draw new state
		state$beta <<- mvrnorm(1, mu, S1);
		print(sprintf('Beta: %.2f', state$beta))
		NULL
	},
	llOutcome = function(i, yFam, lpredFam, famSel) {
		p = pnorm(lpredFam);
		ll = as.numeric(na.omit(log(ifelse(yFam, p, 1 - p))));
		ll
	},
	# latent liability
	update_liability = function(i) {
		lpred = cbind(X, gtScores[genotypes() + 1]) %*% state$beta;
		liab = rtruncnorm(length(lpred), lpred + state$re, 1,
			lower = ifelse(y == 1, 0, -Inf), upper = ifelse(y == 1, Inf, 0)
		);
		state$liability <<- liab;
		NULL
	},
	update_re = function(i) {
		lpred = cbind(X, gtScores[genotypes() + 1]) %*% state$beta;
		res = state$liability - lpred;
		# <p> compute posterior distribution
		sigmas = sapply(1:N, function(i)1/(1/state$sigmaRe + Nfams[i]));
		means = sapply(1:N, function(i)(sigmas[i] * sum(res[Ifams[[i]]])));
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
		print(sprintf('SigmaRe: %.2f', state$sigmaRe));
		NULL
	}
	#
	#	</p> methods
	#
	)
);
MCMCBinProbitReFamClass$accessors(names(MCMCBinProbitReFamClass$fields()));

#
#	<p> MCMC based on random effect with correlation structure proportional to coefficients of relationship
#		copied from MCMCLinearReRelClass
#

MCMCBinProbitReRelClass = setRefClass('MCMCBinProbitReRel', contains = 'MCMCBinProbitReFam',
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
		lpred = cbind(X, gtScores[genotypes() + 1]) %*% state$beta;
		res = state$liability - lpred;
		# <p> compute posterior distribution
		scoreReNO = unlist(lapply(1:N, function(i) {
			Sigma = solve(corsInv[[i]]/state$sigmaRe + diag(rep(1, Nfams[i])));
			Mu = as.vector(res[Ifams[[i]]] %*% Sigma);
			mvrnorm(mu = Mu, Sigma = Sigma)
		}));
		state$re <<- scoreReNO[IfamsOinv];
		NULL
	}
	#	</p> methods
	#
	)
);
MCMCBinProbitReRelClass$accessors(names(MCMCBinProbitReRelClass$fields()));
