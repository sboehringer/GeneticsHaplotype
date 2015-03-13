#
#	mcmc.R
#Tue Nov 25 16:05:28 CET 2014

#
#	<p> MCMC base class
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

#
#	<p> helper function
#

# hack to vivify method lookup via list-syntax
activateMethods = function(class, methods) {
	sapply(methods, function(m)eval(parse(text = sprintf('class$%s', m))))
	NULL
}

last = function(v)v[length(v)]
#
#	<p> block updating class
#
#	Inherits from MCMC and adds the ability to block updating across the parameter components.
#	A list of integers is provided that determines the order and count of updating steps.
#	Each parameter component is garuanteed to be called with a sequential number (w/o gaps)
#	for updating steps. The following pseudo-code applies:
#		N = sum(unlist(blocking))
#		Ncum = cumsum(unlist(blocking))
#		for (i in 1:Ncycles) {
#			Iblock = (i - 1) %/% N + 1;
#			# update within block
#			j = i %% N;
#			# component
#			Ic = which(j <= Ncum)[1];
#			# index within component
#			Iwi = blocking[[Ic]] * (j - 1) + (j - Ncum[Ic]) + 1;
#			.self[[updateComponentName]](Iwi);
#		}


MCMCBlockClass = setRefClass('MCMCBlock', contains = 'MCMC',
	fields = list(
		.blocking = 'list',
		# <p> pre-computed vars
		updateMethods = 'character'
	),
	methods = list(
	#
	#	<p> methods
	#
	initialize = function(..., blocking = list()) {
		callSuper(...);
		.blocking <<- blocking;
		updateMethods <<- as.character(sapply(names(.blocking), function(e)sprintf('update_%s', e)));
		activateMethods(.self, updateMethods);	# <A> hack

		.self
	},
	run = function() {
		# number of MCMC cycles
		Ncycles = Nburnin + Nchain;
		# cycles to spend in parameter components "blocks"
		N = sum(unlist(.blocking));
		Ncum = cumsum(unlist(.blocking));
		# padded version to subtract #updates taken so far
		NcumPad = c(0, Ncum);
		ns = names(.blocking);
		nsU = unique(names(.blocking));
		NperComp = nlapply(nsU, function(n)sum(as.integer(.blocking[which(ns == n)])));
		NcumPerComp = nlapply(nsU, function(n)c(0, cumsum(as.integer(.blocking[which(ns == n)]))));
		# cumulative component index: the how-manieth update is performed for the current block
		IcumComp = lapply(1:length(ns), function(i)sum(ns[1:i] == ns[i]));

		for (i in 1:Ncycles) {
			# how many full rounds of updating (indexed from 0)?
			Ifull = (i - 1) %/% N;
			# update within block (one iteration of all updates in .blocking), indexed from 0
			j = (i - 1) %% N;
			# component to update (indexed from 1)
			Ic = which(j + 1 <= Ncum)[1];
			#
			# index within component
			#
			# name of current component
			n = ns[Ic];
			# start with offset from previous rounds
			Iwi = last(NperComp[[n]]) * Ifull +
				# add #updates taking from begin of current component
				j - NcumPad[Ic] +
				# add cumulative #updates taken for the current component in the current round
				NcumPerComp[[n]][IcumComp[[Ic]]] + 1;
			method_ = .self[[updateMethods[Ic]]];
			method_(Iwi);
		}
	},
	update_beta = function(i) {
		print(sprintf('beta: #%d', i));
	},
	update_hts = function(i) {
		print(sprintf('hts: #%d', i));
	}
	#
	#	</p> methods
	#
	)
);
MCMCBlockClass$accessors(names(MCMCBlockClass$fields()));


#
#	<p> imputation
#

HaplotypeHelperClass = setRefClass('HaplotypeHelper',
	fields = list(
		# <p> pre-computed values
		N = 'integer',
		Nhts = 'integer',
		Ncum = 'integer',
		Ifams = 'list',
		Ifounders = 'integer',
		IfoundersPerFamily = 'list'
	),
	methods = list(
	#
	#	<p> methods
	#
	initialize_cache = function() {
		peds_ = getPeds();
		countMarkers = getCountMarkers();
		# <p> pre-compute
		Nhts <<- as.integer(2^countMarkers);
		N <<- length(peds_);
		Ncum <<- as.integer(c(0L, cumsum(pedsFamilySizes(peds_))) + 1L);
		IfoundersPerFamily <<- pedsFounderIdcs(peds_);
		Ifounders <<- as.integer(unlist(IfoundersPerFamily));
		Ifams <<- pedsIdcs(peds);
		#assign('state0', state[Ifounders, ], envir = .GlobalEnv);
		.self
	},
	getPeds = function()stop('abstract methdod'),
	getCountMarkers = function()stop('abstract methods'),
	#
	# <p> helpers
	#
	freqHat = function(state, j) {
		htfs = table.n(state[setdiff(Ifounders, IfoundersPerFamily[[j]]), ], min = 0, n = Nhts - 1);
		htfs
	}
	#
	#	</p> methods
	#
	)
);
HaplotypeHelperClass$accessors(names(HaplotypeHelperClass$fields()));

# forward declaration
#setRefClass('DiplotypeReconstructor');

MCMCimputationClass = setRefClass('MCMCimputation', contains = c('MCMC', 'HaplotypeHelper'),
	fields = list(
		peds = 'list',
		# Diplotype recontstruction to use
		#reconstruction = 'DiplotypeReconstructor',
		reconstruction = 'envRefClass',
		state = 'matrix'
	),
	methods = list(
	#
	#	<p> methods
	#
	initialize = function(...) {
		.self$initFields(...);
		# Haplotype Helper
		initialize_cache();
		# add default prior
		if (is.null(prior$haplotypes)) prior$haplotypes <<- rep(1, Nhts);
		# draw first state
		htfs = rep(1, Nhts);	# <i> draw from Dirichlet
		state <<- reconstruction$drawFromHfs(htfs, runif(length(peds)));
		.self
	},
	getCountMarkers = function()reconstruction$countMarkers(),
	#
	# <p> helpers
	#
	redrawFamily = function(j) {
		# posterior distribution of haplotypes
		htfsPost = freqHat(state, j) + prior$haplotypes;
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
#	<p> MCMC linear
#

pop = function(v)(v[-length(v)])
# differences between successive elements, first diff is first element with start
vectorLag = function(v, start = 0)pop(c(v, start) - c(start, v))
splitN = function(N, by = 4) vectorLag(round(cumsum(rep(N/by, by))));

# @par: lolRaw: list of list containing lists to be meshed
# @par: lolMesh: containing the mesh: matrix of elements to be picked row-wise
#	column indicates list index from lolRaw value is the list position, NA's cause component to be skipped
meshLists = function(lolRaw, lolMesh) {
	idcs = merge.multi.list(list(list(1:nrow(lolMesh)), list(1:ncol(lolMesh))));
	r = lapply(1:nrow(idcs), function(i) {
		Iel = lolMesh[idcs[i, 1], idcs[i, 2]];
		print(lolRaw[[idcs[i, 2]]][[Iel]]);
		if (is.na(Iel)) NULL else lolRaw[[idcs[i, 2]]][Iel]
	});
	List_(unlist.n(r, 1), rm.null = T);
}

rmultinomLog = function(n, size = 1, logProb)rmultinom(n, size, vector.std(exp(logProb)))

require('MASS');

sqrtMatrix = function(V) {
	Vsvd = svd(V);
	Vsqrt = Vsvd$u %*% diag(sqrt(Vsvd$d)) %*% t(Vsvd$v);
}
matrixM12 = function(V) solve(sqrtMatrix(V))

MCMCLinearClass = setRefClass('MCMCLinear', contains = c('MCMCBlock', 'HaplotypeHelper'),
	fields = list(
		# Diplotype recontstruction to use
		reconstructions = 'list',
		# risk genotypes for each reconstruction
		reconstructionsGts = 'list',
		peds = 'list',
		state = 'list',
		X = 'matrix',	# design matrix
		y = 'numeric',	# response vector
		Nloci = 'integer'
	),
	methods = list(
	#
	#	<p> methods
	#
	# <!> obsolete
	genotypesPrecompute = function() {
		reconstructionsGts <<- lapply(reconstructions, function(m) {
			# alleles
			as = apply(m[, -1, drop = F], c(1, 2), function(e)e%%2);
			t(apply(as, 1, function(r)apply(matrix(r, byrow = T, ncol = 2), 1, sum)))
		});
		#apply(state$hts, 1, function(hts)(hts[1] %% 2 + hts[2] %% 2))
		NULL
	},
	initialize = function(peds = NULL, reconstructor = NULL, ..., NpedSplit = 4) {
		# determine number of loci, reconstructions
		reconstructions <<- R$reconstructionsAll();
		Nloci <<- as.integer(log2(max(unlist(reconstructions)) + 1));
		Npeds = length(peds);

		# determine blocking
		blockHts = listKeyValue(rep('hts', NpedSplit), splitN(Npeds, NpedSplit));
		blockLinB = list(beta = 1);
		blockLinS = list(sigma = 1);
		lol = list(blockHts, blockLinB, blockLinS);
		mesh = matrix(c(1:NpedSplit, rep(1, NpedSplit), rep(1, NpedSplit)), ncol = 3);
		blocking = meshLists(lol, mesh);
		callSuper(blocking = blocking, peds = peds, reconstructions = reconstructions, ...);

		# Haplotype Helper and other initialization
		initialize_cache();
		genotypesPrecompute();

		# priors
		if (is.null(prior$hts)) prior$hts <<- rep(1, 2^Nloci);
		
		# initial state (chain)
		htfs = rep(1, Nhts);	# <i> draw from Dirichlet
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
		.self
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
		Xg <- cbind(X, genotypes());
		# compute posterior distribution (mean, cov-mat of MVN)
		S1inv = prior$betaVarInv + t(Xg) %*% Xg / state$sigma;
		S1 = solve(S1inv);
		mu = S1 %*% (prior$betaMuScaled + t(Xg) %*% y) / state$sigma;
		# draw new state
print(Xg);
print(mu);
print(S1);
		state$beta <<- mvrnorm(1, mu, S1);
		print(sprintf('Beta: %.2f', state$beta))
		NULL
	},
	update_sigma = function(i) {
		Xg <- cbind(X, genotypes());
		res = as.vector(y) - Xg %*% state$beta;
		gishape = prior$gishape + nrow(Xg)/2;
		giscale = prior$giscale + (t(res) %*% res)/2;
		state$sigma <<- 1/rgamma(1, shape = gishape, scale = 1/as.numeric(giscale));
		print(sprintf('Sigma: %.2f', state$sigma));
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
			E = cbind(X[famSel, , drop = F], reconstructionsGts[[iF]][k, ]) %*%
				state$beta;
			# likelihood phenotypes
			llPts = sum(dnorm(y[famSel], E, sd = state$sigma, log = TRUE));
			# likelihood haplotyeps
			logFactor = reconstructions[[iF]][k, 1] * log(2);
			htsFounders = matrix(reconstructions[[iF]][k, -1], byrow = T, nrow = 2)[, Ifdrs, drop = F];
			llHts = logFactor + sum(logPhts[as.vector(htsFounders + 1)]);
			ll = llPts + llHts;
			ll
		});
		# draw family haplotypes as a block
		draw = 1:nrow(reconstructions[[iF]]) %*% rmultinomLog(1, 1, Preconstructions);
		htsI = t(matrix(reconstructions[[iF]][draw, -1], byrow = T, nrow = 2));
		state$hts[Ifams[[iF]], ] <<- htsI;
		NULL
	}
	#
	#	</p> methods
	#
	)
);
MCMCLinearClass$accessors(names(MCMCLinearClass$fields()));
