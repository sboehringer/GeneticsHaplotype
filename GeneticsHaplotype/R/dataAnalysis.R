#
#	dataAnalysis.R
#Wed Apr  1 14:35:38 CEST 2015

getSnps = function(plink, snps) {
	gtRaw = as.data.frame(plink$genotypes[, snps]);
	gtRaw = Df_(gtRaw, as_integer = names(gtRaw));
	gtRaw[gtRaw == 0] = NA;
	gts = gtRaw - 1;
	gts
}
getSnpsWindow = function(plink, snp, windowSize = 5) {
	snpI = which(row.names(plink$map) == snp);
	window = max(snpI - windowSize, 1):min(snpI + windowSize, nrow(plink$map));
	gts = getSnps(plink, row.names(plink$map)[window]);
	gts
}

gtsPruneR2 = function(gts, taboo = NULL, threshold = .95, NsnpsMin = 3) {
	corW = cor(gts, use = 'pairwise.complete.obs')^2;
	while (max(corW[lower.tri(corW)]) > threshold && nrow(corW) > NsnpsMin) {
		isGreater = sapply(1:nrow(corW), function(i)(max(corW[i, -i]) > threshold)) & 
			!(row.names(corW) %in% taboo);
		Iremove = min(which(isGreater));
		corW = corW[-Iremove, -Iremove];
	}
	corW
}
snpsSelectR2 = function(corW, snp, N = 2) {
	#cors = sapply(1:nrow(corW), function(i)max(corW[i, -i]));
	snpI = which(row.names(corW) == snp);
	cors = corW[snpI, -snpI];
	#names(cors) = row.names(corW);
	r = names(sort(cors, decreasing = T))[1:N];
	r
}

prepareMCMC = function(plinkGts, snp, Ntag = 2, Nwindow = 5, NfidThreshold = 4) {
	# <p> get genotypes
	gts = getSnpsWindow(plinkGts, snp, windowSize = Nwindow);
	corW = gtsPruneR2(gts, taboo = snp);
	snpsW = snpsSelectR2(corW, snp, N = Ntag);
	gtsA = getSnps(plinkGts, c(snp, snpsW));

	# <p> prepare ped
	ped = Df_(plinkGts$fam,
		headerMap = list(pedigree = 'fid', member = 'iid', father = 'pid', mother = 'mid'),
		min_ = 'affected');
	peduO = pedu = ped2uniqueId(ped);
	# <p> exclude large families
	fids = unique(pedu$fid);
	Nfid = sapply(fids, function(fid)sum(pedu$fid == fid));
	fidExcl = fids[Nfid > NfidThreshold];
	# incestual families
	fidExcl = c(fidExcl, 3);
	idcsExcl = which(pedu$fid %in% fidExcl);
	pedu = peduO[-idcsExcl, ];
	# <p> add missing founders
	pedS = pedSanitize(pedu);
	# <p> embedding map
	mapE = which.indeces(pedu$iid, pedS$iid);

	# <p> prepare final data set
	gtsMCMC0 = matrix(as.integer(NA), ncol = ncol(gtsA), nrow = nrow(pedS));
	gtsMCMC = matrix.assign(gtsMCMC0, mapE, as.matrix(gtsA)[-idcsExcl, ]);
	y = as.integer(vector.assign(rep(NA, nrow(pedS)), mapE, plinkGts$fam$affected[-idcsExcl] - 1));
	peds = pedSplit2ivTrios(pedS);
	R = new(DiplotypeReconstructor, gtsMCMC, pedsItrios2rcpp(peds));
	reconstructions = R$reconstructionsAll();

	# <p> prune empty reconstructions
	Nreconstructions = sapply(reconstructions, nrow)
	famIempty = which(Nreconstructions == 0 | Nreconstructions > 1e3);
	rowRemove = unlist(pedsIdcs(peds)[famIempty]);
	pedSP = pedS[-rowRemove, ];
	pedsP = pedSplit2ivTrios(pedSP);
	gtsMCMCP = gtsMCMC[-rowRemove, ];
	yP = y[-rowRemove];

	# <p> re-reconstruct
	RP = new(DiplotypeReconstructor, gtsMCMCP, pedsItrios2rcpp(pedsP));
	print('check pruning');
	print(which(sapply(RP$reconstructionsAll(), nrow) == 0));
	#pedsCoeffOfRel = function(peds)lapply(1:length(peds), function(i){print(i);print(peds[[i]]); pedCoeffOfRel(peds[[i]]);})
	cors = pedsCoeffOfRel(pedsP);
	X = model.matrix(~ 1, data.frame(dummy = rep(1, length(yP))));
	r = list(y = yP, X = X, R = RP, peds = pedsP, cors = cors);
	r
}
