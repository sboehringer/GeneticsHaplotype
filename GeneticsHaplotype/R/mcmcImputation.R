#
#	mcmcImputation.R
#Wed Apr  1 17:28:07 CEST 2015

# forward declaration
#setRefClass('DiplotypeReconstructor');

MCMCimputationClass = setRefClass('MCMCimputation', contains = c('MCMC', 'HaplotypeHelper'),
	fields = list(
		peds = 'list',
		# Diplotype recontstruction to use
		#reconstructor = 'DiplotypeReconstructor',
		reconstructor = 'envRefClass',
		state = 'matrix',
		imputation = 'matrix'
	),
	methods = list(
	#
	#	<p> methods
	#
	initialize = function(peds, ...) {
		callSuper(NsampleSpacing = length(peds), peds = peds, ...);
	},
	getCountMarkers = function()reconstructor$countMarkers(),
	genotypes = function() {
		apply(state, 1, function(hts)(hts[1] %% 2 + hts[2] %% 2))
	},
	#
	# <p> helpers
	#
	redrawFamily = function(j) {
		# posterior distribution of haplotypes
		htfsPost = freqHat(state, j) + prior$haplotypes;
		#print(round(vector.std(htfsPost)*36, 1));
		# <A> module indexes from 0
		dtsJ = reconstructor$drawFamFromHfs(j - 1, htfsPost, runif(1));
		#if (any(state[Ncum[j]:(Ncum[j + 1] - 1), ] - dtsJ != 0)) browser();

		state[Ncum[j]:(Ncum[j + 1] - 1), ] <<- dtsJ;
		NULL
	},
	getParameter = function() {
		#if (any(state[Ifounders, ] - state0 != 0)) print(which(state[Ifounders, ] - state0 != 0));
		table.n.freq(state[Ifounders, ], min = 0, n = Nhts - 1)
	},
	#
	# <p> chain implementation
	#
	runInitialize = function() {
		# Haplotype Helper
		initialize_cache();
		# add default prior
		if (is.null(prior$haplotypes)) prior$haplotypes <<- rep(1, Nhts);
		NULL
	},
	drawFromPrior = function() {
		# draw first state
		htfs = rep(1, Nhts);	# <i> draw from Dirichlet
		state <<- reconstructor$drawFromHfs(htfs, runif(length(peds)));
		NULL
	},
	update = function(i) {
		# iteratively update families
		famI = ((i - 1) %% N) + 1;
		redrawFamily(famI);
	},
	sample = function() {
		callSuper();
		updateImputation();
	},
	updateImputation = function() {
		dfg = data.frame(g = as.factor(genotypes()));
		gmat = model.matrix(model.frame(~ 0 + g, dfg), dfg);
		Isample = length(chain);
		imputation <<- if (Isample == 1) gmat else ((gmat + (Isample - 1)*imputation)/Isample);
		NULL
	}


	#
	#	</p> methods
	#
	)
);
MCMCimputationClass$accessors(names(MCMCimputationClass$fields()));

