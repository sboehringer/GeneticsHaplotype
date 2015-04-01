#
#	simulationModels.R
#Wed Apr  1 18:16:36 CEST 2015

simulationModelComparison = function(
	# genetic parameters
	htfs = rev(vector.std(1:8)), N = 2e2, pedTemplate,
	# phenotypes
	simulatePhenotypes, useCors = F, usePedIdcs = F,
	# models to compare
	models,	#list(mcmcClass, usePhenotype, coerceFounders)
	# other parameters
	Nburnin = 1e3L, Nchain = 1e4L, missingness = 0, ...
	) {
	# diplotypes
	d = simulateFromTemplate(pedTemplate, N = N, hfs = htfs);
	# phenotypes
	phenotypePars = list(...);
	if (usePedIdcs) phenotypePars = c(phenotypePars, list(pedIdcs = pedsIdcs(d$peds)));
	if (useCors) phenotypePars = c(phenotypePars, list(useCors = pedsCoeffOfRel(d$peds)));
	y = do.call(simulatePhenotypes, c(list(gts = d$gts[, 1]), phenotypePars))[, 1];
	# Covariates
	X = model.matrix(~ 1, data.frame(dummy = rep(1, length(y))));

	# induce missingness
	dMiss = if (missingness > 0) createMissing(d, missing = missingness) else d;

	r = lapply(models, function(model) with(model, {
		Logs('Analyzing with class: %{mcmcClass}s', 4);
		d0 = if (coerceFounders) createIndependent(dMiss) else dMiss;
		R = new(DiplotypeReconstructor, d0$gts, pedsItrios2rcpp(d0$peds));
		# instantiate object
		mcmcArgs = list(peds = d0$peds, reconstructor = R, Nburnin = Nburnin, Nchain = Nchain);
		if (usePhenotype) mcmcArgs = c(mcmcArgs, list(y = y, X = X));
		mcmc = do.call(new, c(list(Class = mcmcClass), mcmcArgs));
		# run chain
		mcmc$run();
		# summary
		r = c(mcmc$summary(), R2 = mcmc$R2(d$gts[, 1]));
		r
	}));
	r
}
