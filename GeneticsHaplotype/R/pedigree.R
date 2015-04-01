#
#	pedigree.R
#Mon Oct  6 13:53:45 2014

library('sets');

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
	if (nrow(ped$itrios) == 0) {
		return(cbind(ped$founders, NA));
	}
	N = length(ped$founders) + nrow(ped$itrios);
	mapping = apply(ped$itrios, 1, function(r) with(as.list(r), c(mid, 0, pid, 1)));
	mapping = matrix(as.vector(mapping), byrow = T, ncol = 2);
	# add as yet unmapped ids
	mapping = unique(rbind(mapping, as.matrix(data.frame(iid = setdiff(1:N, mapping[, 1]), sex = NA))));
	if (nrow(mapping) != N) {
		print(mapping);
		stop(sprintf('Sex mapping inconsistent in pedigree.'));
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
	for (p in peds) {
		if (any(is.na(p$itrios[, c('mid', 'pid')]))) browser();
	}
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

plotPedigree = plotPedigree_kinship2 = function(ped, tag = '', sex = NULL) {
	pedu = ped2uniqueId(ped);
	require('kinship2');
	# <p> infer and recode sex
	sexu = if (is.null(sex)) pedInferSex(pedu) else sex;
	sex = 2 - sexu;	# recode
	sex[is.na(sex)] = 3;	# 3 value for missing
	# <p> plot using function pedigree
	pedPlot = with(pedu, {par(cex = .3); pedigree(iid, pid, mid, sex, affected = rep(0, nrow(pedu)))});
	par(xpd = T);
	plot(pedPlot, id = paste(paste(ped$iid, ped$fid, sep = '*'), tag, sep = "\n"));
}

plotPedigree_gap = function(ped) {
	pedu = ped2uniqueId(ped);
	require('gap');
	# ... pedtodot
}

# NfamPerRow: #families per row / #rows
plotPedigrees = function(ped, tag = '', NfamPerRow = 5, sex = NULL) {
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
		plotPedigree(ped[pedI, ], tag[pedI], sex[pedI]);
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

#
#	<p> derive coefficients of relationships matrix from ped matrix
#

pedIdParents = function(id, itrios) with(itrios, {
	pedParentsFollow = function(id, dist) {
		if (!(id %in% iid)) return(NULL);
		idI = which(iid == id);
		rbind(
			matrix(c(mid[idI], dist + 1, pid[idI], dist + 1), byrow = T, ncol = 2),
			pedParentsFollow(mid[idI], dist + 1),
			pedParentsFollow(pid[idI], dist + 1)
		)
	};
	pedParentsFollow(id, 0)
})
# ped is in founders/itrios format
pedAncestry = function(itrios) {
	ancestry = nlapply(itrios$iid, pedIdParents, itrios = itrios);
	ancestry
}

#
#	coefficients of relationship between non-founders
#
pedCoeffOfRelNonFounders = function(ped) {
	ancest = pedAncestry(ped$itrios);
	iids = ped$itrios$iid;
	if (length(iids) == 1) return(Df_(as.data.frame(t(rep(0, 3))), names = c('id1', 'id2', 'cor'))[-1, ]);
	cor = lapply(as.list(set_combn(iids, 2)), function(pair) {
		pair = as.integer(as.list(pair));
		# add pair itself to ancestors as non-founders could be parents of each other
		anc1 = rbind(ancest[[as.character(pair[1])]], c(pair[1], NA));
		anc2 = rbind(ancest[[as.character(pair[2])]], c(pair[2], NA));
		shared = intersect(anc1[, 1], anc2[, 1]);
		meioses = cbind(
			anc1[which.indeces(shared, anc1[, 1]), 2],
			anc2[which.indeces(shared, anc2[, 1]), 2]);
		meiosesDists = apply(meioses, 1, sum, na.rm = T);
		meiosesDist = min(meiosesDists);
		# sum(meiosesDists == meiosesDist) count how many minimum paths exist: e.g. sibs, cousins
		# we assume no loops <!><N>
		cor = 2^(-min(meiosesDist)) * sum(meiosesDists == meiosesDist);
		c(pair, cor)
	});
	Df_(do.call(rbind, cor), names = c('id1', 'id2', 'cor'));
}
#
#	coefficients of relationship between founders and non-founders
#
pedCoeffOfRelFounderNonFounders = function(ped) {
	ancest = pedAncestry(ped$itrios);
	pairs = Df_(merge.multi(ped$founders, ped$itrios$iid), names = c('f', 'nf'));
	iids = ped$itrios$iid;
	cor = apply(pairs, 1, function(pair) {
		ancNf = ancest[[as.character(pair['nf'])]];
		meiosesDist = ancNf[which(ancNf[, 1] == pair['f']), 2];
		cor = if (length(meiosesDist) > 0) 2^(-meiosesDist) else 0;
		c(pair, cor);
	});
	Df_(t(cor), names = c('id1', 'id2', 'cor'));
}

pedCoeffOfRel = function(ped) {
	if (length(ped$founders) == 1) return(matrix(1, nrow = 1, ncol = 1));
	corNf = pedCoeffOfRelNonFounders(ped);
	corFNf = pedCoeffOfRelFounderNonFounders(ped);
	cor = rbind(corNf, corFNf);
	corMatRaw = matrixFromIndexedDf(cor, 'id1', 'id2', 'cor', 1:(nrow(ped$itrios) + length(ped$founders)));
	corMat = symmetrizeMatrix(corMatRaw);
	corMat[is.na(corMat)] = 0;
	diag(corMat) = 1;
	corMat
}
pedsCoeffOfRel = function(peds)lapply(peds, pedCoeffOfRel)

#
#	pedigree cleaning
#

C = as.character;
U = unlist;
# split a pedigree data frame (collection of families) into a list of families
pedigreeSeparate = function(ped) {
	# <p> build clusters
	clusters = list();		# holds members of clusters
	n2c = list();			# maps members to clusters; name2cluster
	# <p> founders
	iidFdr = ped$iid[is.na(ped$mid) & is.na(ped$pid)];
	for (i in iidFdr) {
		it = ped[which(ped$iid == i), ];			# inheritance trio
		# open new cluster for founders
		clusters[[C(i)]] = it;
		# i is garuanteed to be a unique new cluster name
		n2c[[C(it$iid)]] = C(i);
	}
	# <p> offspring
	# <!> only process inheritance trio if both parents are already known
	# -> forward references
	stack = setdiff(ped$iid, iidFdr);
	while (length(stack) > 0) {
		it = ped[which(ped$iid == stack[1]), ];			# inheritance trio
		# is any parents not yet processed?
		if (is.null(n2c[[C(it$mid)]]) || is.null(n2c[[C(it$pid)]])) {
			# rotate stack
			stack = c(stack[-1], stack[1]);
			next;
		}
		# join clusters if needed
		if (n2c[[C(it$mid)]] != n2c[[C(it$pid)]]) {
			clusters[[n2c[[C(it$mid)]]]] = rbind(clusters[[n2c[[C(it$mid)]]]], clusters[[n2c[[C(it$pid)]]]]);
			clusters[[n2c[[C(it$pid)]]]] = NULL;
			n2c[which(n2c == n2c[[C(it$pid)]])] = n2c[[C(it$mid)]];
		}
		n2c[[C(it$iid)]] = n2c[[C(it$mid)]];
		clusters[[n2c[[C(it$iid)]]]] = rbind(clusters[[n2c[[C(it$iid)]]]], it);
		stack = stack[-1];
	}
	clusters
}

pedSanitizeSingle = function(ped, pedFactor = 1e3L) {
	pedO = ped;
	# <p> missing parents
	newParents = list();
	myFid = ped$fid[1];
	iidCounter = pedFactor * myFid;
	for(i in 1:nrow(ped)) {
		if (sum(is.na(ped[i, c('mid', 'pid')])) == 1) {
			if (is.na(ped$mid[i])) {
				if (!is.null(newParents[[C(ped$pid[i])]])) {
					ped$mid[i] = newParents[[C(ped$pid[i])]];
				} else {
					newParents[[C(ped$pid[i])]] = ped$mid[i] = iidCounter;
					iidCounter = iidCounter + 1;
				}
			} else {
				if (!is.null(newParents[[C(ped$mid[i])]])) {
					ped$pid[i] = newParents[[C(ped$mid[i])]];
				} else {
					newParents[[C(ped$mid[i])]] = ped$pid[i] = iidCounter;
					iidCounter = iidCounter + 1;
				}
			}
		}
	}
	if (length(newParents) > 0) {
		pedNew = Df(fid = myFid, iid = unique(unlist(newParents)), mid = NA, pid = NA);
		ped = rbind(ped, pedNew);
	}
	# <p> separate into components
	cls = pedigreeSeparate(ped);
	if (length(cls) > 1) {
		clsN = lapply(1:length(cls), function(i) {
			ped = cls[[i]];
			ped$fid = pedFactor * myFid + i;
			ped
		});
		ped = do.call(rbind, clsN);
	}

# 	# <p> split unreachable members
# 	iidParents = unique(na.omit(c(ped$mid, ped$pid)));
# 	fidCounter = pedFactor * ped$fid[1];
# 	for(i in 1:nrow(ped)) {
# 		# founder && not reachable
# 		if (sum(is.na(ped[i, c('mid', 'pid')])) == 2 && !(ped$iid[i] %in% iidParents)) {
# 			ped$fid[i] = fidCounter;
# 			fidCounter = fidCounter + 1;
# 		}
# 	}

	ped
}
# assume ped to have been run through ped2uniqueId
pedSanitize = function(pedu) {
	peds = lapply(unique(pedu$fid), function(fid)pedSanitizeSingle(pedu[pedu$fid == fid, ]));
	pedS = do.call(rbind, peds);
	pedS
}

pedRemoveSingletons = function(ped) {
	fids = unique(ped$fid);
	pedSingletons = fids[which(sapply(fids, function(fid)length(ped$iid[ped$fid == fid])) == 1)];
	pedFam = ped[!(ped$fid %in% pedSingletons), ];
	pedFam
}
