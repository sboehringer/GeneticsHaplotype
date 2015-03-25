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
		NULL
	},
	drawFromPrior = function() {
		callSuper();
		state$beta <<- rnorm(ncol(X) + 1, 0, 1);	# use prior parameters <!>
		state$liabilityVar <<- rep(1, length(y));	# use prior parameters <!>
		state$liability <<- rep(0, length(y));	# use prior parameters <!>
		NULL
	},
	blocking = function() {
		# determine blocking
		blockHts = listKeyValue(rep('hts', NpedSplit), splitN(N, NpedSplit));
		blockBinLiab = list(liability = 1);
		blockBinLiabVar = list(liabilityVar = 1);
		blockBinB = list(beta = 1);
		lol = list(blockHts, blockBinLiab, blockBinLiabVar, blockBinB);
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
		nu = 7.3;	# t-distribution based approximation of the logistic (Kinney & Dunson 2006)
		lpred = (cbind(X, gtScores[genotypes() + 1]) %*% state$beta)[, 1];
		shape = (nu + 1)/2;
		scaleInv = (nu + (state$liability - lpred)^2)/2;
#by(scaleInv, y, stem)
		state$liabilityVar <<- 1/rgamma(length(lpred), shape, scale = 1/scaleInv);
		NULL
	}
	#
	#	</p> methods
	#
	)
);
MCMCBinomialClass$accessors(names(MCMCBinomialClass$fields()));

