#
#	simulation.R
#Mon Oct  6 13:54:25 2014


simulateDiplotypes = function(peds, hfs) {
	Nfounder = sapply(peds, function(e)length(e$founder));
	Nfcum = c(0, cumsum(Nfounder));
	Nfounders = sum(Nfounder);
	Nitrio = sapply(peds, function(e)nrow(e$itrios));
	Nicum = c(0, cumsum(Nitrio));
	Nitrios = sum(Nitrio);

	hts = as.vector(t(rmultinom(2 * Nfounders, 1, hfs)) %*% t(t(0:(length(hfs) - 1))));
	dts = matrix(hts, byrow = T, ncol = 2);
	# draw inheritance vectors
	iv = rbinom(2 * Nitrios, 1, .5) + 1;
	dts = lapply(seq_along(peds), function(i) {
		dtsF = matrix(0, ncol = 2, nrow = Nfounder[i] + Nitrio[i]);
		dtsF[peds[[i]]$founders, ] = dts[(Nfcum[i] + 1):(Nfcum[i] + Nfounder[i]), ];
		for (j in 1:nrow(peds[[i]]$itrios)) {
			trio = peds[[i]]$itrios[j, ];
			dtsF[trio$iid, ] = c(
				dtsF[trio$mid, ][iv[Nicum[i] + 2*j - 1]],
				dtsF[trio$pid, ][iv[Nicum[i] + 2*j]]);
		}
		dtsF
		
	});
	r = do.call(rbind, dts);
	r
}

familiesFromTemplate = function(template, N) {
	peds = lapply(1:N, function(i)data.frame(fid = i, template));
	do.call(rbind, peds)
}

simulateDiplotypesPed = function(ped, hfs) {
	pedSplit = pedSplit2ivTrios(ped);
	dts = simulateDiplotypes(pedSplit, hfs = hfs);
	dts
}

diplotypeFs = function(ped, dts, countLoci = ceiling(log2(max(dts) + 1))) {
	dts = dts[pedFounders(ped), ];
	t0 = table.n(as.vector(dts), min = 0, n = 2^countLoci - 1);
	r = t0 / sum(t0);
	r
}

# convert diplotypes to genotypes
# diplotypes as from simulateDiplotypes: ncol = 2, nrow = N
# genotypes are vector of alleles, 2 alleles per genotype, countLoci times per individual
# summarize is applied to pairs of alleles, sum applies for SNPs

diplotypes2gts = function(dts, countLoci = ceiling(log2(max(dts) + 1)), summarize = NULL, ...) {
	gts = apply(dts, 1, function(d){
		g1 = ord2bin(d[1], digits = countLoci);
		g2 = ord2bin(d[2], digits = countLoci);
		as.vector(rbind(g1, g2));	# cog alleles
	});
	if (!is.null(summarize)) gts = matrix(as.vector(t(
		apply(matrix(as.vector(gts), ncol = 2, byrow = T), 1, summarize, ...)
	)), nrow = countLoci);
	t(gts)
}
diplotypes2gtsDose = function(dts, countLoci = ceiling(log2(max(dts) + 1)), summarize = sum, ...) {
	diplotypes2gts(dts, countLoci, summarize, ...)
}

#
#	<p> simulation functions
#

simulateFromPed = function(ped, hfs = 1:8) {
	dts = simulateDiplotypesPed(ped, hfs);
	gts = diplotypes2gts(dts, summarize = sum);
	dtfs = diplotypeFs(ped, dts);
	o = pedForwardOrder(ped);
	pedO = ped[o, ];

	r = list(ped = pedO, dts = dts, gts = gts, dtfs = dtfs, peds = pedSplit2ivTrios(ped));
	r
}
simulateFromTemplate = function(pedTemplate, N = 25, hfs = 1:8) {
	simulateFromPed(familiesFromTemplate(pedTemplate, N = N), hfs)
}

vector.std = function(v, C = 1)(C * v / sum(v))
simulationFromPed = function(ped, hfs, Nsim, module) {
	hfs = vector.std(hfs);
	sim = simulateFromPed(ped, hfs = hfs);
	R = new(module$DiplotypeReconstructor, sim$gts, pedsItrios2rcpp(sim$peds));
	dts = sapply(1:Nsim, function(i) {
		dtsS = R$drawFromHfs(sim$dtfs, runif(length(sim$peds)));
		diplotypeFs(sim$ped, dtsS)
	});
	dtsMean = apply(dts, 1, mean);
	dtsSd = apply(dts, 1, sd);
	r = list(Nsim = Nsim, hfs = hfs, hfsE = dtsMean, hfsSd = dtsSd, hfsB = abs(dtsMean - hfs));
	r
}
simplify2a = simplify2array;
simulationsFromPed = function(ped, hfs, NsimPped, Nsim, module) {
	r = lapply(1:Nsim, function(i)simulationFromPed(ped, hfs, NsimPped, module));
	hfsB = apply(simplify2a(list.kp(r, 'hfsB')), 1, mean);
	hfsSd = apply(simplify2a(list.kp(r, 'hfsSd')), 1, mean);
	r = list(Nsim = Nsim, NsimPped = NsimPped, hfs = vector.std(hfs), hfsB = hfsB, hfsSd = hfsSd,
		mse = hfsB^2 + hfsSd^2
	);
	r
}

#
#	<p> missingness
#

createMissing = function(d, missing = .05) {
	d$gts[runif(length(d$gts)) < missing] = NA;
	d
}

#
#	<p> phenotype helpers
#

af2hwe = function(q) { c((1-q)^2, 2*q*(1-q), q^2 ) }

logitI = expit = function(x, min = 0, max = 1) { (max - min)/(1 + exp(-x)) + min }
logit = function(x, min = 0, max = 1) { log((x - min)/(max - x)) }

scoresL = scoresS = scoresStd = list(
	dom = c(0, 1, 1), add =  c(0, .5, 1), rec = c(0, 0, 1), gen = factor(0:2)
);
names(scoresL) = c("dominant", "additive", "recessive", 'genotype');
names(scoresS) = c("D", "A", "R", 'G');

# @param gts: vector of genotypes
# @param beta: vector of effects (including baseline)
# @param score: vector of risk scores per individual
# @param scoresGt: scoring vector for gts

simulatePhenotypesBinRaw = function(gts, beta, score = 0, scoreGt = scoresL$additive, link = expit) {
	d = data.frame(scores = scoreGt[gts + 1]);
	mm = model.matrix(model.frame(~ scores, d), d);
	risk = link(mm %*% beta + score);
	y = as.integer(risk > runif(length(risk)));
	#qplot(risk, y, position = position_jitter(w = 0.01, h = 0.1), alpha = .1) + theme_bw()
	list(mm = mm, risk = risk, y = y)
}

simulatePhenotypesBin = function(gts, beta, score = 0, scoreGt = scoresL$additive, link = plogis) {
	simulatePhenotypesBinRaw(gts, beta, score, scoreGt, link)$y
}

simulatePhenotypesBinReFam = function(pedIdcs, gts, beta, score = 0,
	scoreGt = scoresL$additive, sdRe = 1, link) {

	reFam = rnorm(length(pedIdcs), 0, sdRe);
	scoreReNO = unlist(lapply(1:length(pedIdcs), function(i)rep(reFam[i], length(pedIdcs[[i]]))));
	scoreRe = scoreReNO[inverseOrder(unlist(pedIdcs))];
	simulatePhenotypesBinRaw(gts, beta, score + scoreRe, scoreGt = scoreGt, link = link)$y
}
simulatePhenotypesBinReRel = function(pedIdcs, gts, beta, score = 0,
	scoreGt = scoresL$additive, sdRe = 1, link = plogis, cors) {

	scoreReNO = unlist(lapply(1:length(pedIdcs),
		function(i)mvrnorm(mu = rep(0, nrow(cors[[i]])), Sigma = sdRe^2 * cors[[i]])
	));
	scoreRe = scoreReNO[inverseOrder(unlist(pedIdcs))];
	simulatePhenotypesBinRaw(gts, beta, score + scoreRe, scoreGt = scoreGt, link = link)$y
}

simulatePhenotypesLinearRaw = function(gts, beta, sd = 1, score = 0, scoreGt = scoresL$additive) {
	d = data.frame(scores = scoreGt[gts + 1]);
	mm = model.matrix(model.frame(~ scores, d), d);
	linpred = mm %*% beta + score;
	y = linpred + rnorm(nrow(mm), sd = sd);
	list(mm = mm, linpred = linpred, y = y)
}

simulatePhenotypesLinear = function(gts, beta, sd = 1, score = 0, scoreGt = scoresL$additive) {
	simulatePhenotypesLinearRaw(gts, beta, sd, score, scoreGt)$y
}

simulatePhenotypesLinearReFam = function(pedIdcs, gts, beta, sd = 1, score = 0,
	scoreGt = scoresL$additive, sdRe = 1) {

	reFam = rnorm(length(pedIdcs), 0, sdRe);
	scoreReNO = unlist(lapply(1:length(pedIdcs), function(i)rep(reFam[i], length(pedIdcs[[i]]))));
	scoreRe = scoreReNO[inverseOrder(unlist(pedIdcs))];
	simulatePhenotypesLinearRaw(gts, beta, sd, score + scoreRe, scoreGt)$y
}

simulatePhenotypesLinearReRel = function(pedIdcs, gts, beta, sd = 1, score = 0,
	scoreGt = scoresL$additive, sdRe = 1, cors) {

	reFam = rnorm(length(pedIdcs), 0, sdRe);
	scoreReNO = unlist(lapply(1:length(pedIdcs),
		function(i)mvrnorm(mu = rep(0, nrow(cors[[i]])), Sigma = sdRe^2 * cors[[i]])
	));
	scoreRe = scoreReNO[inverseOrder(unlist(pedIdcs))];
	simulatePhenotypesLinearRaw(gts, beta, sd, score + scoreRe, scoreGt)$y
}

#
#	unit tests
#

# test wether drawn diplotypes are compatible with genotype input
test_drawCompatibility = function(ped, hfs, Nsim = 1e2, module) {
	hfs = vector.std(hfs);
	o = pedForwardOrder(ped);
	sim = simulateFromPed(ped, hfs = hfs);
	R = new(module$DiplotypeReconstructor, sim$gts, pedsItrios2rcpp(sim$peds));
	dts = sapply(1:Nsim, function(i) {
		dtsS = R$drawFromHfs(sim$dtfs, runif(length(sim$peds)));
		gts = diplotypes2gtsDose(dtsS);
		if (any(sim$gts[o, ] != gts)) {
			print(o);
			print(cbind(sim$gts[o, ], gts));
			stop('drawn diplotypes incompatible with observed genotypes.');
		}
	});
	NULL
}

test_reconstruction = function(d, reconstructions) {
	pedI = pedsIdcs(d$peds);
	countLoci = ceiling(log2(max(d$dts) + 1));

	sapply(1:length(reconstructions), function(i) {
		sapply(1:nrow(reconstructions[[i]]), function(j) {
			dts = matrix(reconstructions[[i]][j, -1], byrow = T, ncol = 2);
			if (any(d$gts[pedI[[i]], ] != diplotypes2gtsDose(dts, countLoci))) {
				stop('diplotypes incompatible with genotypes.');
			}
		})
	})
}
