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
		.blocking = 'list'
	),
	methods = list(
	#
	#	<p> methods
	#
	initialize = function(..., blocking = list()) {
		callSuper(...);
		.blocking <<- blocking
		.self
	},
	run = function() {
		# number of MCMC cycles
		Ncycles = Nburnin + Nchain;
		# cycles to spend in parameter components "blocks"
		N = sum(unlist(blocking));
		Ncum = cumsum(unlist(blocking));

		for (i in 1:Ncycles) {
			Iblock = (i - 1) %/% N + 1;
			# update within block
			j = i %% N;
			# component
			Ic = which(j <= Ncum)[1];
			# index within component
			Iwi = blocking[[Ic]] * (j - 1) + (j - Ncum[Ic]) + 1;
			.self[[updateComponentName]](Iwi);
		}
	},
	hello = function() {
		print("hello world");
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
#	<p> MCMC linear
#

pop = function(v)(v[-length(v)])
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

MCMCLinearClass = setRefClass('MCMCLinear', contains = 'MCMCBlock',
	fields = list(
		# Diplotype recontstruction to use
		#reconstruction = 'Rcpp_DiplotypeReconstructor',
		peds = 'list',
		state = 'matrix'
	),
	methods = list(
	#
	#	<p> methods
	#
	initialize = function(peds = NULL, ..., NpedSplit = 4) {
		Npeds = length(peds);
		blockHts = listKeyValue(rep('hts', NpedSplit), splitN(Npeds, NpedSplit));
		blockLinB = list(beta = 1);
		blockLinS = list(sigma = 1);
		lol = list(blockHts, blockLinB, blockLinS);
		mesh = matrix(c(1:NpedSplit, rep(1, NpedSplit), rep(1, NpedSplit)), ncol = 3);
		blocking = meshLists(lol, mesh);
		callSuper(peds, ..., blocking = blocking);
		.self
	},
	updateBeta = function(i) {
	}
	#
	#	</p> methods
	#
	)
);
MCMCLinearClass$accessors(names(MCMCLinearClass$fields()));
