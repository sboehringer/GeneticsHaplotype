#
#	haplotypeExp.R
#Tue Jul  8 17:21:11 CEST 2014

# project to implement c++-based code for haplotype reconstruction based on earlier code

source('haplotype.R');
system("cd ~/src/Rprivate ; ./exportR.sh"); source("RgenericAll.R"); source("Rgenetics.R"); loadLibraries();

if (T) {
	ped = data.frame(fid = 1, iid = 1:6, mid = c(NA, NA, 1, 1, NA, 4), pid = c(NA, NA, 2, 2, NA, 5));
	ped = ped[c(1,6,3,4,5,2), ];
}

if (0) {
	print(ped);
	pedu = ped2uniqueId(ped);
	print(pedu);
	#print(pedigreeInheritanceTrios(pedu));
}

if (0) {
	# use maternal line
	cols = c('id', 'm');
	d0 = pedu[cols];
	pedp = pedu[cols];
	i = 0;
	while (any(!is.na(d0$m))) {
		i = i + 1;
		mnew = Sprintf('m%{i}02d')
		d0 = Df_(d0, names = list(m = mnew));
		d0 = merge(d0, pedp, by.x = mnew, by.y = 'id', all.x = T);
print(d0);
	}
	depth = apply(Df_(d0, min_ = cols), 1, function(r)sum(is.na(r)));
	peds = pedu[order_align(depth, pedu$id), ]
	peds = pedu[order(depth[order_align(pedu$id, d0$id)]), ]
}

if (0) {
	peds = ped2ivTrios(ped)
}

if (1) {
	ivTrios = pedSplit2ivTrios(ped);
	print(ivTrios);
}
