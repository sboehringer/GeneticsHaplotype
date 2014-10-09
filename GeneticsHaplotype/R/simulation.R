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
diplotypes2gtsDose = function(dt, countLoci = ceiling(log2(max(dts) + 1)), summarize = sum, ...) {
	diplotypes2gts(dt, countLoci, summarize, ...)
}
