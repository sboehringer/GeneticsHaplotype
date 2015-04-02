#
#	simulationModelsExample.R
#Thu Apr  2 21:31:43 CEST 2015

modelListMCMCsimulations = list(
	htfs = list(
		8:1,
		c(1, .1, 1, 1, .1, .1, .1, .1),
		c(1, .05, .05, .05, .05, .05, .05, .05)
	),
	missing = list(0, .2, .5),
	mcmcSpec = list(
		phenotypeFunction = 'simulatePhenotypesBinReRel',
		phenotype = list(usePedIdcs = T, useCors = T),
		models = list(
		list(mcmcClass = 'MCMCimputation', passPhenotype = F, coerceFounders = F),
		list(mcmcClass = 'MCMCimputation', passPhenotype = F, coerceFounders = T),
		list(mcmcClass = 'MCMCBinProbitReRel', passPhenotype = T, coerceFounders = F, passCors = T),
		list(mcmcClass = 'MCMCBinProbitReRel', passPhenotype = T, coerceFounders = T, passCors = T)
	));

);
