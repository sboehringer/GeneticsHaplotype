#
#	simulationModelsExample.R
#Thu Apr  2 21:31:43 CEST 2015

modelListMCMCsimulationsPedTemplate = Df(names = c('iid', 'mid', 'pid'), matrix(
	c(	1, NA, NA,
		2, NA, NA,
		3, 1, 2,
		4, NA, NA,
		5, 3, 4
	), byrow = T, ncol = 3)
);

modelListMCMCsimulations = list(
	pars = list(list(
		Nburnin = 1e4L, Nchain = 1e5L, Nrepetition = 5e2, N = 2e2,
		pedTemplate = modelListMCMCsimulationsPedTemplate
	)),
	htfs = list(
		8:1,
		c(1, .1, 1, 1, .1, .1, .1, .1),
		c(1, .05, .05, .05, .05, .05, .05, .05)
	),
	beta = list(
		c(0, 0),
		c(0, .5),
		c(0, 1)
	),
	missingness = list(0, .2, .5),
	mcmcSpec = inlist(list(list(
		phenotypeFunction = 'simulatePhenotypesBinReRel',
		phenotype = list(usePedIdcs = T, useCors = T),
		models = list(
			list(mcmcClass = 'MCMCimputation', passPhenotype = F, coerceFounders = F),
			list(mcmcClass = 'MCMCimputation', passPhenotype = F, coerceFounders = T),
			list(mcmcClass = 'MCMCBinProbitReRel', passPhenotype = T, coerceFounders = F, passCors = T),
			list(mcmcClass = 'MCMCBinProbitReRel', passPhenotype = T, coerceFounders = T, passCors = T)
		)
	), list(
		phenotypeFunction = 'simulatePhenotypesBinReFam',
		phenotype = list(usePedIdcs = T),
		models = list(
			list(mcmcClass = 'MCMCimputation', passPhenotype = F, coerceFounders = F),
			list(mcmcClass = 'MCMCimputation', passPhenotype = F, coerceFounders = T),
			list(mcmcClass = 'MCMCBinProbitReFam', passPhenotype = T, coerceFounders = F),
			list(mcmcClass = 'MCMCBinProbitReFam', passPhenotype = T, coerceFounders = T)
		)
	), list(
		phenotypeFunction = 'simulatePhenotypesLinearReRel',
		phenotype = list(usePedIdcs = T, useCors = T, sd = 1),
		models = list(
			list(mcmcClass = 'MCMCimputation', passPhenotype = F, coerceFounders = F),
			list(mcmcClass = 'MCMCimputation', passPhenotype = F, coerceFounders = T),
			list(mcmcClass = 'MCMCLinearReRel', passPhenotype = T, coerceFounders = F, passCors = T),
			list(mcmcClass = 'MCMCLinearReRel', passPhenotype = T, coerceFounders = T, passCors = T)
		)
	), list(
		phenotypeFunction = 'simulatePhenotypesLinearReFam',
		phenotype = list(usePedIdcs = T, sd = 1),
		models = list(
			list(mcmcClass = 'MCMCimputation', passPhenotype = F, coerceFounders = F),
			list(mcmcClass = 'MCMCimputation', passPhenotype = F, coerceFounders = T),
			list(mcmcClass = 'MCMCLinearReFam', passPhenotype = T, coerceFounders = F),
			list(mcmcClass = 'MCMCLinearReFam', passPhenotype = T, coerceFounders = T)
		)
	), list(
		phenotypeFunction = 'simulatePhenotypesLinear',
		phenotype = list(sd = 1),
		models = list(
			list(mcmcClass = 'MCMCimputation', passPhenotype = F, coerceFounders = F),
			list(mcmcClass = 'MCMCimputation', passPhenotype = F, coerceFounders = T),
			list(mcmcClass = 'MCMCLinear', passPhenotype = T, coerceFounders = F),
			list(mcmcClass = 'MCMCLinear', passPhenotype = T, coerceFounders = T)
		)
	)))
);

modelListMCMCsimulationsDebug = merge.lists(modelListMCMCsimulations, list(
	pars = list(list(
		Nburnin = 1e1L, Nchain = 1e2L, Nrepetition = 2, N = 2e2,
		pedTemplate = modelListMCMCsimulationsPedTemplate
	)), htfs = list(8:1), beta = list(c(0, 0)), missingness = list(0)
));

runSimulation = function(Nburnin, Nchain, Nrepetition, N, pedTemplate, htfs, beta, missingness, mcmcSpec) {
	simulationModelComparisonSpec(mcmcSpec,
		htfs = vector.std(htfs), N = N, pedTemplate = pedTemplate, beta = beta,
		Nburnin = Nburinin, Nchain = Nchain, Nrepetition = Nrepetition, missingness = missingness)
}
