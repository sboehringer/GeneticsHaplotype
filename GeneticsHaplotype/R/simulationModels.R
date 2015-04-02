#
#	simulationModels.R
#Wed Apr  1 18:16:36 CEST 2015

simulationModelComparisonSingle = function(
	# genetic parameters
	htfs = rev(vector.std(1:8)), N = 2e2, pedTemplate,
	# phenotypes
	simulatePhenotypes, useCors = F, usePedIdcs = F,
	# models to compare
	models,	#list(mcmcClass, passPhenotype, coerceFounders)
	# other parameters
	Nburnin = 1e3L, Nchain = 1e4L, missingness = 0, ...
	) {
	# diplotypes
	d = simulateFromTemplate(pedTemplate, N = N, hfs = htfs);
	# phenotypes
	phenotypePars = list(...);
	cors = NULL;
	if (usePedIdcs) phenotypePars = c(phenotypePars, list(pedIdcs = pedsIdcs(d$peds)));
	if (useCors) {
		cors = pedsCoeffOfRel(d$peds);
		phenotypePars = c(phenotypePars, list(cors = cors));
	}
	y = do.call(simulatePhenotypes, c(list(gts = d$gts[, 1]), phenotypePars));
	if (is.matrix(y)) y = y[, 1];
	# Covariates
	X = model.matrix(~ 1, data.frame(dummy = rep(1, length(y))));

	# induce missingness
	dMiss = if (missingness > 0) createMissing(d, missing = missingness) else d;

	r = lapply(models, function(model) with(model, {
		Logs('Analyzing with class: %{mcmcClass}s', 4);
		# <p> data set
		d0 = dMiss;
		cors0 = cors;
		if (coerceFounders) {
			d0 = createIndependent(dMiss);
			cors0 = pedsCoeffOfRel(d0$peds);
		}
		# <p> reconstructor, instantiate object
		R = new(DiplotypeReconstructor, d0$gts, pedsItrios2rcpp(d0$peds));
		# <p> arguments for mcmc instantiation
		mcmcArgs = list(peds = d0$peds, reconstructor = R, Nburnin = Nburnin, Nchain = Nchain);
		if (passPhenotype)			mcmcArgs = c(mcmcArgs, list(y = y, X = X));
		if (nif(model$passCors))	mcmcArgs = c(mcmcArgs, list(cors = cors0));
		mcmc = do.call(new, c(list(Class = mcmcClass), mcmcArgs));
		# <p> run chain
		mcmc$run();
		# <p> summary
		r = c(mcmc$summary(), R2 = mcmc$R2(d$gts[, 1]));
		r
	}));
	r
}

simulationModelComparison = function(..., Nrepetition = 1L, lapply__ = lapply) {
	print(date());
	r = Lapply__(1:Nrepetition, function(i, ...) simulationModelComparisonSingle(...), ...);
	print(date());
	r
}

simulationModelComparisonSpec = function(spec, htfs, N, pedTemplate, beta,
	Nburnin = 1e2L, Nchain = 1e3L, Nrepetition = 1L, missingness = 0, lapply__ = lapply) {
	args = c(list(htfs = htfs, N = N, pedTemplate = pedTemplate, missingness = missingness,
		simulatePhenotype = get(spec$phenotypeFunction), beta = beta,
		models = spec$models
	), spec$phenotype);
	sim = do.call(simulationModelComparison, args);
}