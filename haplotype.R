#
#	haplotype.R
#Tue Jul  8 17:23:59 CEST 2014


ped2uniqueId = function(ped) with(ped, {
	# form id by fid * 10^exponent + iid; exponent chosen to encompass all iids
	shift = 10^ceiling(log10(max(ped$iid)));
	#new_id = function(id)(shift * fid + id);
	idMap = listKeyValue(shift * fid + iid, 1:nrow(ped));
	r = data.frame(fid = fid, id = avu(idMap[iid]), m =  avu(idMap[mid]), p =  avu(idMap[pid]));
	r
})

pedDistFromRoot = function(pedu, col = 'm') {
	# use maternal line
	cols = c('id', col);
	d0 = pedu[cols];
	pedp = pedu[cols];
	i = 0;
	while (any(!is.na(d0[[col]]))) {
		i = i + 1;
		nnew = Sprintf('%{col}s%{i}02d')
		d0 = Df_(d0, names = listKeyValue(col, nnew));
		d0 = merge(d0, pedp, by.x = nnew, by.y = 'id', all.x = T);
	}
	depth = apply(Df_(d0, min_ = cols), 1, function(r)sum(is.na(r)));
	depth = max(depth) - depth;	# distance from root maternally
}

ped2ivTrios = function(ped) {
	pedu = ped2uniqueId(ped);
	depthM = pedDistFromRoot(pedu, 'm');
	depthP = pedDistFromRoot(pedu, 'p');
	depth = apply(cbind(depthM, depthP), 1, max);

	# original ped sorted according to depth (max - number of ancestors in tree)
	depthO = depth[order_align(pedu$id, d0$id)];
	peds = data.frame(pedu, depth = depthO)[order(depthO), ];
	r = list(trios = peds[peds$depth > 0, c('id', 'm', 'p')], founders = peds$id[peds$depth == 0]);
	r
}

pedSplit2ivTrios = function(ped) {
	r = lapply(unique(ped$fid), function(fid)ped2ivTrios(ped[ped$fid == fid, ]));
}
