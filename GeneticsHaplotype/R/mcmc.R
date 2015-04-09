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
	runInitialize = function() {
		NULL
	},
	drawFromPrior = function() {
		NULL
	},
	run = function() {
		runInitialize();
		drawFromPrior();
		N = Nburnin + Nchain;
		for (i in 1:N) {
			.self$update(i);
			if (i > Nburnin && ((i - Nburnin - 1) %% NsampleSpacing) == 0) .self$sample();
		}
	},
	summary = function(names = NULL, credibleAlpha = .05) {
		pars = do.call(rbind, lapply(chain, unlist));
		summaries = apply(pars, 2, function(par) {
			list(
				mean = mean(par), median = median(par), sd = sd(par),
				credible = quantile(par, probs = c(credibleAlpha/2, 1 - credibleAlpha/2))
			)
		});
		if (!is.null(names)) names(summaries) = paste(names, 1:(length(summaries)/length(names)), sep = '');
		summaries
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
	initialize = function(...) {
		callSuper(...);
		.self
	},
	blocking = function() {
		stop('abstract methdod: blocking scheme needs to be returned by subclass')
	},
	run = function() {
		runInitialize();
		drawFromPrior();
		# cache blocking scheme
		.blocking <<- blocking();
		# make sure methods exist
		updateMethods <<- as.character(sapply(names(.blocking), function(e)sprintf('update_%s', e)));
		activateMethods(.self, updateMethods);	# <A> hack
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
			if ((length(.self$NsampleSpacing) &&
				(i > Nburnin && ((i - Nburnin - 1) %% NsampleSpacing) == 0)) || j == 0) {
				.self$sample();
			}
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
		Nfams = 'integer',
		Ncum = 'integer',
		Ifams = 'list',
		Ifounders = 'integer',
		IfoundersPerFamily = 'list',
		IfamsOinv = 'integer'
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
		Nfams <<- pedsFamilySizes(peds_);
		Ncum <<- as.integer(c(0L, cumsum(Nfams)) + 1L);
		IfoundersPerFamily <<- pedsFounderIdcs(peds_);
		Ifounders <<- as.integer(unlist(IfoundersPerFamily));
		Ifams <<- pedsIdcs(peds);
		IfamsOinv <<- inverseOrder(unlist(Ifams));
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
	},
	R2 = function(gtsReal) {
		dosage = imputation %*% 0:2;
		R2 = (cor(gtsReal, dosage, use = 'pairwise.complete.obs')^2)[1, 1];
		R2
	}
	#
	#	</p> methods
	#
	)
);
HaplotypeHelperClass$accessors(names(HaplotypeHelperClass$fields()));

#
#	<p> helper functions
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
#		print(lolRaw[[idcs[i, 2]]][[Iel]]);
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


