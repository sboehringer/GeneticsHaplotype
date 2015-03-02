#
#	pedigree.R
#Mon Oct  6 13:53:45 2014

#	By definition the unique ids used 
#
#

uniqueId = function(fid, iid, mx = max(iid)) {
	# form iid by fid * 10^exponent + iid; exponent chosen to encompass all iids
	shift = 10^ceiling(log10(mx));
	#new_iid = function(iid)(shift * fiid + iid);
	iidU = as.character(shift * fid + iid);	# unique iid
	iidU
}
pedUniqueId = function(ped, mx = max(ped$iid)) with(ped, uniqueId(fid, iid, mx = mx))
ped2uniqueId = function(ped, mx = max(ped$iid)) with(ped, {
	iidU = pedUniqueId(ped, mx = mx);
	idMap = listKeyValue(iidU, 1:nrow(ped));
	r = data.frame(fid = fid,
		iid = avu(idMap[iidU]),
		mid = avu(idMap[uniqueId(fid, mid, mx)]),
		pid = avu(idMap[uniqueId(fid, pid, mx)]));
	r
})

pedDistFromRoot = function(pedu, col = 'mid', maxDepth = 5) {
	# use maternal line
	cols = c('iid', col);
	d0 = pedu[cols];
	pedp = pedu[cols];
	i = 0;
	while (any(!is.na(d0[[col]])) && i < maxDepth) {
		i = i + 1;
		nnew = Sprintf('%{col}s%{i}02d')
		d0 = Df_(d0, names = listKeyValue(col, nnew));
		d0 = merge(d0, pedp, by.x = nnew, by.y = 'iid', all.x = T);
	}
	if (i == maxDepth) stop(Sprintf('pedigree depth too large [depth >= %{maxDepth}d]'));
	depth = apply(Df_(d0, min_ = cols), 1, function(r)sum(is.na(r)));
	depth = max(depth) - depth;	# distance from root maternally
	depthO = depth[order_align(pedu$iid, d0$iid)];
	depthO
}

ped2ivTrios = function(ped) {
	pedu = ped2uniqueId(ped);
	depthM = pedDistFromRoot(pedu, 'mid');
	depthP = pedDistFromRoot(pedu, 'pid');
	depth = apply(cbind(depthM, depthP), 1, max);

	# original ped sorted according to depth (max - number of ancestors in tree)
	peds = data.frame(pedu, depth = depth);
	trios = peds[peds$depth > 0, ];
	r = list(
		itrios = trios[order(trios$depth), c('iid', 'mid', 'pid')],
		founders = peds$iid[peds$depth == 0]);
	r
}

pedSplit2ivTrios = function(ped) {
	r = lapply(unique(ped$fid), function(fid)ped2ivTrios(ped[ped$fid == fid, ]));
}

pedsItrios2rcpp = function(peds) {
	lapply(peds, function(ped) {
		list(itrios = as.matrix(ped$itrios) - 1, founders = ped$founders -1 )
	})
}

# sex: 0: female, 1: male, assume normalized id 1:N
ivTriosInferSex = function(ped) {
	N = length(ped$founders) + nrow(ped$itrios);
	mapping = apply(ped$itrios, 1, function(r) with(as.list(r), c(mid, 0, pid, 1)));
	mapping = matrix(as.vector(mapping), byrow = T, ncol = 2);
	# add as yet unmapped ids
	mapping = unique(rbind(mapping, as.matrix(data.frame(iid = setdiff(1:N, mapping[, 1]), sex = NA))));
	if (nrow(mapping) != N) {
		print(mapping);
		stop(Sprintf('Sex mapping inconsistent in pedigree.'));
	}
	sex = mapping[order(mapping[, 'iid']), 'sex'];
	sex
}

# compute order that maps back a vector with values for iids 1:N_1, 1:N_2, ... back to the original ped-order
pedInverseOrder = function(ped) {
	idOrig = pedUniqueId(ped);
	idUniq = unlist(lapply(unique(ped$fid), function(fid)pedUniqueId(ped[ped$fid == fid, ], max(ped$iid))));
	o = order_align(idOrig, idUniq);
	o
}
pedForwardOrder = function(ped) {
	idOrig = pedUniqueId(ped);
	idUniq = unlist(lapply(unique(ped$fid), function(fid)pedUniqueId(ped[ped$fid == fid, ], max(ped$iid))));
	o = order_align(idUniq, idOrig);
	o
}

pedInferSex = function(ped) {
	peds = pedSplit2ivTrios(ped);
	sexes = unlist(lapply(peds, ivTriosInferSex));
	sexes[pedInverseOrder(ped)]
}

pedFounders = function(ped) {
	apply(ped[, c('mid', 'pid')], 1, function(i)all(is.na(i)))
}

# founder indeces on uniquified version of ped
pedFounderIdcsForward = function(ped) {
	o = pedForwardOrder(ped);
	r = which(pedFounders(ped)[o]);
	r
}

pedsFamilySizes = function(peds)
	sapply(peds, function(ped)as.integer(length(ped$founders) + nrow(ped$itrios)));
pedFamilySizes = function(ped) pedsFamilySizes(pedSplit2ivTrios(ped));
pedsFounderSizes = function(peds) sapply(peds, function(ped)length(ped$founders));
pedFounderSizes = function(ped) pedsFounderSizes(pedSplit2ivTrios(ped));
pedsFounderIdcs = function(peds) {
	if (!length(peds)) return(list());
	Ns = as.integer(pop(c(0, cumsum(pedsFamilySizes(peds)))));
	r = lapply(1:length(peds), function(i)peds[[i]]$founders + Ns[i]);
	r
}
pedsIdcs = function(peds) {
	if (!length(peds)) return(list());
	Nfams = pedsFamilySizes(peds);
	Ns = as.integer(pop(c(0, cumsum(Nfams))));
	r = lapply(1:length(peds), function(i)(1:Nfams[i]) + Ns[i]);
	r
}

plotPedigree = plotPedigree_kinship2 = function(ped, tag = '') {
	pedu = ped2uniqueId(ped);
	require('kinship2');
	# <p> infer and recode sex
	sexu = pedInferSex(pedu);
	sex = 2 - sexu;	# recode
	sex[is.na(sex)] = 3;			# 3 value for missing
	# <p> plot using function pedigree
	pedPlot = with(pedu, pedigree(iid, pid, mid, sex, affected = rep(0, nrow(pedu))));
	par(xpd = T);
	plot(pedPlot, id = paste(paste(ped$iid, ped$fid, sep = '*'), tag, sep = "\n"));
}

plotPedigree_gap = function(ped) {
	pedu = ped2uniqueId(ped);
	require('gap');
	# ... pedtodot
}

# NfamPerRow: #families per row / #rows
plotPedigrees = function(ped, tag = '', NfamPerRow = 5) {
	o = order(ped$fid);
	# re-order ped by families
	ped = ped[o, ];
	tag = tag[o];
	fids = unique(ped$fid);

	# split families into rows
	Nrows = ceiling(length(fids)/NfamPerRow);
	layout(matrix(1:Nrows, ncol = 1));
	idcs = splitListIndcs(length(fids), Nrows);

	# plot
	apply(idcs, 1, function(i) {
		fidI = fids[i[1]:i[2]];
		pedI = which(ped$fid %in% fidI);
		plotPedigree(ped[pedI, ], tag[pedI]);
		NULL
	});
}

# c++ code orders correctly as of 6.11.2014
# # drawn haplotypes are ordered by: Founders, IV-trios
# reorderDraw = function(peds) {
# 	Ns = c(0, cumsum(sapply(peds, function(p)(length(p$founders) + nrow(p$itrios))))) + 1;
# 	o = unlist(lapply(1:length(peds), function(i) {
# 		c(peds[[i]]$founders, peds[[i]]$itrios[, 'iid']) + Ns[i]
# 	}));
# 	o
# }
