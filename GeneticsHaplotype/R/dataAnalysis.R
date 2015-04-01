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
