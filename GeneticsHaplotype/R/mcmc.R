#
#	mcmc.R
#Mon Jul 17 16:48:30 EDT 2006

dec.to.bin = function(number = NA, digits = 4) {
	ret = rep(0, digits);
	for (i in 1:digits) {
		ret[i] = number %% 2;
		number = number %/% 2;
	}
	ret
}
bin.to.dec = function(bin = NA) {
	sum(sapply(1:length(bin), function(i)bin[i]*2^(i-1)))
}
instantiate.list = function(l, n = 1) {
	for (nm in names(l)) {
		if (is.integer(l[[nm]])) {
			eval.parent(parse(file = "", text = sprintf("%s = %d", nm, l[[nm]])), n = n);
		} else if (is.numeric(l[[nm]])) {
			eval.parent(parse(file = "", text = sprintf("%s = %f", nm, l[[nm]])), n = n);
		} else {
			eval.parent(parse(file = "", text = sprintf("%s = \"%s\"", nm, l[[nm]])), n = n);
		};
	}
}

gt.std = function(g) {
	gt = matrix(g, nrow = 2);
	apply(gt, 2, sort)
}

countHets = function(g) {
	pos = which(sapply(1:dim(g)[2], function(i){ g[1, i] != g[2, i]}));
	list(positions = pos, count = length(pos))
}

# pairs of alleles in each line == individual
geno2haplo = function(g) {
	hc = countHets(g);
	hetBits = if (hc$count <= 1) 0 else (hc$count - 1);
	countDts = 2^hetBits;

	dts = NULL;
	for (i in 1:countDts) {
		dt = g;	# pre-initialize homocygous positions
				# (and last heterozygous position)
		dt.bin = dec.to.bin(i - 1);
		for (j in 1:hetBits) {
			dt[1, hc$positions[j]] = dt.bin[j];
			dt[2, hc$positions[j]] = 1 - dt.bin[j];
		}
		dts = rbind(dts, dt);
	}
	apply(dts, 1, bin.to.dec);
}

count.haplotypes = function(dts, countLoci) {
	counts = rep(0, 2^countLoci);
	for (i in 1:length(dts)) counts[dts[i] + 1] = counts[dts[i] + 1] + 1;
	counts;
}

# <!> implementation only true for observed genotypes
draw.from = function(freqs, comp.dts) {
	sub.freqs = freqs[comp.dts + 1] / sum(freqs[comp.dts + 1]);
	#print(list(freqs = freqs, sub.freqs = sub.freqs, dts = comp.dts));
	ht = which(rmultinom(1, 1, sub.freqs) == 1);
	# one ht determines the other
	if (ht %% 2 != 0) c(comp.dts[ht], comp.dts[ht + 1])
	else c(comp.dts[ht - 1], comp.dts[ht])
}

# <N> remedies draw.from for missing data
# comp.dts is matrix with dt in row, 2 columns
# comp.dts is unique in an unordered sense
draw.from.constraint = function(freqs, comp.dts) {
	# draw complete dts
	dt.freqs = apply(comp.dts, 1, function(dt){
		freqs[dt[1] + 1] * freqs[dt[2] + 1] * ifelse(dt[1] == dt[2], 1, 2)
	});
	dt = comp.dts[which(rmultinom(1, 1, dt.freqs) == 1), ];
	dt
}

ll.from.constraint = function(dt, freqs, comp.dts) {
	dt.freqs = apply(comp.dts, 1, function(dt){
		freqs[dt[1] + 1] * freqs[dt[2] + 1] * ifelse(dt[1] == dt[2], 1, 2)
	});
	dt.freq = dt.freqs[which.row(comp.dts, dt)] / sum(dt.freqs);
	log(dt.freq)
}

draw.from.constraint.ll = function(freqs, comp.dts) {
	realization = draw.from.constraint(freqs.comp.dts);
	ll = ll.from.constraint(realization, freqs, comp.dts);
	list(realization = realization, ll = ll)
}

draw.from.counts = function(counts, comp.dts, prior) {
	freqs = (prior + counts) / sum(prior + counts);
	draw.from(freqs, comp.dts)
}

# assume matrix with n rows and 2 * m cols == # loci
mcmc.iid.initialize = function(gts, prior) {
	dts = list();
	n = dim(gts)[1];
	countLoci = dim(gts)[2] / 2;
	for (i in 1:n) {
		dts[[i]] = geno2haplo(gt.std(gts[i, ]));
	}
	counts = count.haplotypes(unlist(dts), countLoci);
	multi = rdirichlet(1, prior / sum(prior));
	imputations = NULL;
	for (i in 1:n) {
		imputations = c(imputations, draw.from.constraint(multi, dts[[i]]));
	}

	list(diplotypes = dts, counts = counts, imputations = imputations,
		countLoci = countLoci, prior = prior
	)
}

mcmc.iid.update = function(state) {
	counts.imp = state$counts;
	imputations = NULL;
	for (i in 1:length(state$diplotypes)) {
		counts.imp = counts;
		# remove current observation
		c.i = state$imputations[c(2*i - 1, 2*i)] + 1;
		counts.imp[c.i[1]] = counts.imp[c.i[1]] - 1;
		counts.imp[c.i[2]] = counts.imp[c.i[2]] - 1;

		draw.new = draw.from.counts.missing(counts.imp, state$diplotypes[[i]], state$prior);
		imputations = c(imputations, draw.new);

		counts.imp[draw.new[1] + 1] = counts.imp[draw.new[1] + 1] + 1;
		counts.imp[draw.new[2] + 1] = counts.imp[draw.new[2] + 1] + 1;
	}
	# summarize imputations
	freqs.new = counts.imp / sum(counts.imp);
	merge.lists(state, list(imputations = imputations, freqs = freqs.new, counts = counts.imp));
}

expand.block = function(count, blocksize, indeces) {
	as.vector(apply(to.col(1:count), 1,
		function(i){ (i - 1) * blocksize + t(to.col(indeces)) }
	));
}

runif.int = function(n = 1, min = 1, max = 2) {
	r = runif(n, min, max + 1);
	sapply(r, as.integer)
}

search.block = function(l, s) {
	b.sz = length(s);
	which(sapply(
		1:(length(l)/b.sz), function(i){all(l[((i - 1) * b.sz + 1):(i * b.sz)] == s)}
	));
}
which.row = function(m, row) {
	rows = 1:(dim(m)[1]);
	rows.found = rows[sapply(rows, function(i){ all(m[i, ] == row) })];
	rows.found
}

OFF.IDC = function(i) { (3 + 2 * i):(4 + 2 * i) }
# <i> subsort diplotypes
ll.for.fam.from = function(fam, dist, configs, famSize) {
	# P(H_m | G_i), P(H_p | H_m, G_i), P(H_j | H_m, H_p, G_i)
	cfg.size = 4 + (famSize - 2) * 8;
	c.cfgs = length(configs) / cfg.size;
	c.offspring = famSize - 2;

	mat.dts = configs[expand.block(c.cfgs, cfg.size, 1:2)];
	mat.dts.u = unique(matrix(mat.dts, ncol = 2, byrow = T));
	ll.m = ll.from.constraint(fam[1:2], dist, mat.dts.u);	# ll for maternal dt
	ll = ll.m;

	pat.dts = configs[expand.block(c.cfgs, cfg.size, 3:4)];
	# pat.dts are conditional on maternal dts
	mat.i = search.block(mat.dts, fam[1:2]);
	cnd.pat.dts = matrix(pat.dts, ncol = 2, byrow = T)[mat.i, ];
	cnd.pat.dts.u = unique(cnd.pat.dts);
	ll.p = ll.from.constraint(fam[3:4], dist, cnd.pat.dts.u);	# cond ll for paternal dt
	ll = ll + ll.p;

	# select unique parental diplotype combination
	par.dts = configs[expand.block(c.cfgs, cfg.size, 1:4)];
	p.cfg = search.block(par.dts, fam[1:4]);
	if (0) {
		print(col.frame(list(fam = fam, mat.dts = mat.dts, cnd.pat = as.vector(cnd.pat.dts),
		ll.m = ll.m, ll.p = ll.p, ll = ll), do.paste = " ", digits = c(0, 0, 0, 3, 3, 3)));
	}
	if (0) {
		print(col.frame(list(p.cfg = p.cfg, par.dts = par.dts, fam = fam[1:4]),
		do.paste = " ", digits = c(0, 0, 0)));
	}
	fam.cfg = configs[((p.cfg - 1) * cfg.size + 1):(p.cfg * cfg.size)];
	off.cfg = fam.cfg[-(1:4)];
	new.off.dts = c();
	for (i in 1:c.offspring) {
		ll = ll + log(mendel(fam[1:2], fam[3:4], fam[OFF.IDC(i)]));
	}
	ll
}

draw.fam.from = function(dist, configs, famSize, debug = F) {
	# P(H_m | G_i), P(H_p | H_m, G_i), P(H_j | H_m, H_p, G_i)
	cfg.size = 4 + (famSize - 2) * 8;
	c.cfgs = length(configs) / cfg.size;
	c.offspring = famSize - 2;
#print(matrix(configs, ncol = cfg.size, byrow = T));
#print(expand.block(c.cfgs, cfg.size, 1:2));

	mat.dts = configs[expand.block(c.cfgs, cfg.size, 1:2)];
#print(col.frame(list(mat.dts = mat.dts), digits = 0, do.paste = ","));
	mat.dts.u = unique(matrix(mat.dts, ncol = 2, byrow = T));
	new.mat.dt = draw.from.constraint(dist, mat.dts.u);
	pat.dts = configs[expand.block(c.cfgs, cfg.size, 3:4)];
	# pat.dts are conditional on maternal dts
	mat.i = search.block(mat.dts, new.mat.dt);
	cnd.pat.dts = matrix(pat.dts, ncol = 2, byrow = T)[mat.i, ];
	if (debug) print(col.frame(list(
		mat.dts = as.vector(t(mat.dts.u)), new.mat.dt= new.mat.dt,
		pat.dts = pat.dts, cnd.pat.dts = as.vector(t(cnd.pat.dts))), digits = 0, do.paste = ","
	));
	if (length(cnd.pat.dts) == 0) print(list(cfgs = matrix(configs, ncol = cfg.size, byrow = T), 
		new.mat.dt = new.mat.dt));
#print(list(
#	mat.dts.u = unique(matrix(mat.dts, ncol = 2, byrow = T)), new.dt = new.mat.dt,
#	mat.i = mat.i,
#	pat.dts.u = unique(matrix(cnd.pat.dts, ncol = 2, byrow = T)),
#	pat.dts = pat.dts,
#	cfgs = configs,
#));
	cnd.pat.dts.u = unique(cnd.pat.dts);
	new.pat.dt = draw.from.constraint(dist, cnd.pat.dts.u);
	# select unique parental diplotype combination
	par.dts = configs[expand.block(c.cfgs, cfg.size, 1:4)];
	# trick to compute config in which the parental dts are present <!>; seems bogus
	#	p.cfg = max(which(par.dts == rep(c(new.mat.dt, new.pat.dt), c.cfgs)) %/% 4);
	p.cfg = search.block(par.dts, c(new.mat.dt, new.pat.dt));
	fam.cfg = configs[((p.cfg - 1) * cfg.size + 1):(p.cfg * cfg.size)];
	off.cfg = fam.cfg[-(1:4)];
	new.off.dts = c();
	for (i in 1:c.offspring) {
		this.off.cfg = off.cfg[((i - 1)*8 + 1):(i * 8)];
		# select valid dts
		this.off.cfg = this.off.cfg[this.off.cfg >= 0];
		#print(unlist(list(mat = new.mat.dt, pat = new.pat.dt, off = this.off.cfg)));
		# uniformly select from possible dts
		dt.no = runif.int(1, 1, length(this.off.cfg) / 2);
		new.off.dt = this.off.cfg[(dt.no * 2 - 1):(dt.no * 2)];
		new.off.dts = c(new.off.dts, new.off.dt)
	}
#	# draw offspring dts from parents Mendelian-ish
#	# draw inheritance vector
#	inh.v = matrix(as.integer(runif(2 * c.offspring) > .5) + 1, ncol = 2);
#	# draw maternal alleles
#	m.dts = new.mat.dt[inh.v[,1]];
#	# draw paternal alleles
#	p.dts = new.pat.dt[inh.v[,2]];
#	new.off.dts = as.vector(t(cbind(m.dts, p.dts)));
	c(new.mat.dt, new.pat.dt, new.off.dts)
}


# just count parents
count.haplotypes.fam = function(dts, famCounts, countLoci) {
	counts = rep(0, 2^countLoci);
	mat.dts.i = cumsum(2 * famCounts) - 2 * famCounts + 1;
	par.dts = dts[sort(c(mat.dts.i, mat.dts.i + 1, mat.dts.i + 2, mat.dts.i + 3))];
	for (i in 1:length(par.dts)) counts[par.dts[i] + 1] = counts[par.dts[i] + 1] + 1;
	counts
}

mcmc.nucfam.initialize = function(gts, famCounts, countLoci, prior) {
	dts = geno2haploFam(gts, famCounts, countLoci, 0);
#print(list(gts = gts, dts = dts));
	configCounts = dts$configCounts;
	multi = rdirichlet(1, prior / sum(prior));
	imputations = NULL;
	n = length(dts$configCounts);

	mat.dts.i = cumsum(2 * famCounts) - 2 * famCounts + 1;
	cfgs.s = configCounts * sapply(famCounts, function(i)(4 + (i - 2) * 8));
	cfgs.i = cumsum(cfgs.s) - cfgs.s + 1;
	cfgs.i = c(cfgs.i, last(cfgs.i) + last(cfgs.s));
#print(list(cfg.c = configCounts, cfgs.i = cfgs.i));
#print(list(reconstruction = matrix(
#	dts$diplotypes,
#	byrow = T, ncol = 4 + (famCounts[1] - 2) * 8
#)[, 1:4]));

	for (i in 1:n) {
		new.fam = draw.fam.from(multi,
			dts$diplotypes[cfgs.i[i]:(cfgs.i[i + 1] - 1)], famCounts[i]);
		imputations = c(imputations, new.fam);
	}
 	counts = count.haplotypes.fam(imputations, famCounts, countLoci);
	list(
		diplotypes = dts, famCounts = famCounts, countLoci = countLoci, prior = prior,
		imputations = imputations, counts = counts, cfgs.i = cfgs.i, mat.dts.i = mat.dts.i,
	)
}

last = function(l)l[length(l)];

mcmc.nucfam.update = function(state) {
	attach(state, warn.conflicts = F);
	configCounts = diplotypes$configCounts;
	dts = diplotypes$diplotypes;
	imp.new = NULL;

	n = length(state$diplotypes$configCounts);
	for (i in 1:n) {
		# remove current observation (parents)
		par.dts.c = table.n(imputations[mat.dts.i[i]:(mat.dts.i[i] + 3)] + 1, 2^countLoci);
#		print(list(counts = counts, par.dts.c = par.dts.c));
		counts = counts - par.dts.c;
#		print(list(counts.post = counts));

		cfg = dts[cfgs.i[i]:(cfgs.i[i + 1] - 1)];
		multi = (prior + counts) / sum(prior + counts);
		draw.new = draw.fam.from(multi, cfg, famCounts[i]);
		imp.new = c(imp.new, draw.new);

		# add imputation
		draw.new.c = table.n(draw.new[1:4] + 1, 2^countLoci);
		counts = counts + draw.new.c;
#		print(list(draw.new = draw.new, counts.post.post = counts));
	}
	# summarize imputations
	counts.new = counts;
	freqs.new = counts / sum(counts);
	detach();
	merge.lists(state, list(imputations = imp.new, freqs = freqs.new, counts = counts.new));
}

# all proposal distributions are symmetric here.
# pass NULL for lh.prop if lh.prop is uniform and cancels out
metropolis = function(new, ll.new, old, ll.old) {
	# <!> hack to avoid starting problems
	if (is.na(ll.old)) return(list(value = new, rejection = 0));
	if (is.na(ll.new)) print("ll.new is na!");
	alpha = log(runif(1));
	alpha.p = if (ll.new == -Inf) -Inf else { if (ll.old > -Inf) ll.new - ll.old else 0};
#print(list(alpha = alpha, alpha.p = alpha.p, ll.old = ll.old, ll.new = ll.new));
	if (alpha > alpha.p) list(value = old, rejection = 1, lh = ll.old)
	else list(value = new, rejection = 0, lh = ll.new)
}
metropolis.hastings = function(new, ll.new, old, ll.old, ll.prop.new, ll.prop.old) {
	# <!> hack to avoid starting problems
	if (is.na(ll.old)) return(list(value = new, rejection = 0));
	if (is.na(ll.new)) print("ll.new is na!");
	alpha = log(runif(1));
	# assume ll.prop.new, ll.prop.old to be real
	alpha.p = if (ll.new == -Inf) -Inf else {
		if (ll.old > -Inf) ll.new - ll.prop.new - (ll.old - ll.prop.old) else 0
	};
#print(list(alpha = alpha, alpha.p = alpha.p, ll.old = ll.old, ll.new = ll.new));
	if (alpha > alpha.p) list(value = old, rejection = 1, ll = ll.old, ll.prop = ll.prop.old)
	else list(value = new, rejection = 0, ll = ll.new, ll.prop = ll.prop.new)
}

get.lh.hstar = function(countLoci, baseline, penetrance = penetrance.add) {
	function(pars, hstar, phenotypes) ll.hts.f.R(pars, hstar, phenotypes, countLoci, baseline, penetrance);
}

get.lh.beta = function(countLoci, baseline, penetrance = penetrance.add, phenotypes) {
	function(pars, hts) ll.hts.R(pars, hts, phenotypes, countLoci, baseline, penetrance);
}

# update H*
mcmc.full.update.hstar = function(s, cycle.no) {
	#instantiate.list(state);
	# draw haplotypes from a gibbs sampler (cmp mcmc.nucfam.update)
	configCounts = s$diplotypes$configCounts;
	dts = s$diplotypes$diplotypes;
	imp.new = NULL;
	counts.old = counts = s$counts;
	freqs.old = counts / sum(counts);

	rejections = 0;
	n = length(s$diplotypes$configCounts);
	ll = 0;
	for (i in 1:n) {
		# predictive updating step
		# remove current observation (parents)
		m.i = s$mat.dts.i[i];
		fam.old.par = s$imputations[m.i:(m.i + 3)];	#parents
		par.dts.c = table.n(fam.old.par + 1, 2^s$countLoci);
#		print(list(counts = counts, par.dts.c = par.dts.c));
		counts = counts - par.dts.c;
#		print(list(counts.post = counts));

		cfg = dts[s$cfgs.i[i]:(s$cfgs.i[i + 1] - 1)];
		#multi.real.prop = counts + 1;	# avoid zeros
		multi.real.prop = sapply(counts, function(e)max(e,1))
		multi.real = multi.real.prop / sum(multi.real.prop);
		# choose a proposal distribution
		multi.prop =  multi.real;
		#multi.prop = counts + 2 * rep(sum(counts) / 2^s$countLoci, 2^s$countLoci);
		#multi.prop = rep(sum(counts) / 2^s$countLoci, 2^s$countLoci);
		multi = multi.prop / sum(multi.prop);

#		multi = (s$prior$haplotypes + counts) / sum(s$prior$haplotypes + counts);
# hack up multi
#		norm.f = sum(multi[-(1:2)]) / (1 - sum(s$free.pars[1:2]));
#		multi = c(s$free.pars[1:2], multi[-(1:2)]/norm.f);
#print(c(multi, sum(multi)));
		draw.new = draw.fam.from(multi, cfg, s$famCounts[i]);

		phenotypes = s$phenotypes[s$pts.i[i]: (s$pts.i[i + 1] - 1)];
		# compute old likelihood
		fam.old.dts = s$imputations[m.i:(m.i + 4 + 2 * (s$famCounts[i] - 2) - 1)];
#		lh.old = ll.hts.f.c.R(c(multi.real, s$beta), fam.old.dts, phenotypes, s$countLoci,
#			s$baseline, s$penetrance);
#		lh.new = ll.hts.f.c.R(c(multi.real, s$beta), draw.new, phenotypes, s$countLoci,
#			s$baseline, s$penetrance);
#if (cycle.no == 6) print(col.frame(list(lh.old = lh.old, lh.new = lh.new, multi = multi, draw.new = draw.new, old = s$imputations[s$mat.dts.i[i]:(s$mat.dts.i[i + 1] - 1)], p = phenotypes), do.paste = ", "));
		lh.prop.old = ll.for.fam.from(fam.old.dts, multi.real, cfg, length(fam.old.dts) / 2);
		lh.prop.new = ll.for.fam.from(draw.new, multi.real, cfg, length(fam.old.dts) / 2);
		lh.old = lh.prop.old + ll.fam.pts.R(s$beta, fam.old.dts, phenotypes, s$countLoci,
			s$baseline, s$penetrance);
		lh.new = lh.prop.new + ll.fam.pts.R(s$beta, draw.new, phenotypes, s$countLoci,
			s$baseline, s$penetrance);

		# metropolis-hastings sampling
		new.m = metropolis.hastings(draw.new, lh.new, fam.old.dts, lh.old, lh.prop.new, lh.prop.old);
		rejections = rejections + new.m$rejection;	# count rejections
		ll = ll + new.m$ll;

		# add imputation
		draw.new.c = table.n(new.m$value[1:4] + 1, 2^s$countLoci);
		counts = counts + draw.new.c;
#		print(list(draw.new = draw.new, counts.post.post = counts));
		imp.new = c(imp.new, new.m$value);
	}
#print(col.frame(list(cycle = cycle.no, freqs.old = freqs.old, freqs.new = counts / sum(counts), rejections = rejections, ll = ll), do.paste = " "));
	# summarize imputations
	counts.new = counts;
	freqs.new = counts / sum(counts);

	merge.lists(s, list(imputations = imp.new, hfreqs = freqs.new, counts = counts.new,
		h.rej = s$h.rej + rejections)
	)
}

# we expect the ll for haplotypes to be passed as argument, the missing phenotype part is caculated
mcmc.full.update.beta = function(s, cycle.no = NA, count.updates = 1) {
	# we sample to betas per cycle, the haplotype likelihood keeps constant
	# therefore it is precalculated for the beta-updates
	# ll.hts cancels out in the metropolis step, therefore doesn't have to be computed
#	print(unlist(list(hfs = s$hfreqs)));
#	ll.hts = ll.hts.c.constraint.R(s$hfreqs, s$imputations, s$famCounts, s$countLoci, s$baseline,
#			s$penetrance, s$diplotypes$diplotypes, s$cfgs.i);
	ll.hts = 0;

	new.rej = 0;
	old = s$beta;
	lh.old = ll.hts + ll.pts.c.constraint.R(old, s$imputations, s$phenotypes,
			s$famCounts, s$countLoci, s$baseline, s$penetrance);
	new.m = list(value = s$beta, rejection = 0, lh = lh.old);

	for (i in 1:s$beta.prop$updates) {
		# draw form a proposal with mean prior value -> is symmetric <N>
		new = rnorm(1, new.m$value, s$beta.prop$var);
		lh.new = ll.hts + ll.pts.c.constraint.R(new, s$imputations, s$phenotypes,
			s$famCounts, s$countLoci, s$baseline, s$penetrance);
#print(t(col.frame(list(beta.old = new.m$value, beta.new = new, lh.old = new.m$lh, lh.new = lh.new, rej = new.m$rejection))));
		# metropolis sampling; new.m$lh is likelihood from last round or initialization
		new.m = metropolis(new, lh.new, new.m$value, new.m$lh);
		new.rej = new.rej + new.m$rejection;
	}
	merge.lists(s, list(beta = new.m$value, beta.rej = s$beta.rej + new.rej))
}

rnorm.trunc = function(n, m, sd, a, b) {
	prop = a;
	while (prop <= a || prop >= b) { prop = rnorm(1, m, sd); }
	prop
}

dnorm.trunc = function(x, m = 0, sd = 1, a = -Inf, b = Inf) {
	ret = dnorm(x, m, sd) / (1 - (pnorm(a, m, sd) + (1 - pnorm(b, m, sd))));
	ret
}

norm.multi = function(p) {
	p = c(p, max(0, 1 - sum(p)));
	ret = p / sum(p);
	ret
}

mcmc.full.update = function(state, cycle.no = NA) {
	state.new = mcmc.full.update.hstar(state, cycle.no);
	state.new = mcmc.full.update.beta(state.new, cycle.no, state.new$beta.prop$updates);
	if (cycle.no == state.new$beta.prop$switch.cycle) state.new$beta.prop$var = state.new$beta.prop$var.late;
	state.new
}

sort.pairs = function(l) {
	n = length(l)/2;
	l.ret = as.vector(sapply(1:n, function(i)sort(l[(2*i-1):(2*i)])));
	l.ret
}

# count loci are all loci included in the model
mcmc.full.initialize = function(fixed.pars, free.pars,
	gts, pts, famCounts, countLoci, prior, dont.update.beta = F) {
	instantiate.list(fixed.pars);
	# compute haplotype distribution from ld pars
	hfreqs = ld.f2hFreq(free.pars);

	# create fist state for haploytpes from prior
	dts = geno2haploFam(gts, famCounts, countLoci - 1, 1);
	# we need ordered pairs for quick indexing
	dts = merge.lists(dts, list(diplotypes = sort.pairs(dts$diplotypes)));
	configCounts = dts$configCounts;
	multi = rdirichlet(1, prior$haplotypes / sum(prior$haplotypes));
	# start with true distribution
	#	multi = hfreqs;
	imputations = NULL;
	n = length(dts$configCounts);

	pts.i = c(cumsum(famCounts) - famCounts + 1, sum(famCounts) + 1);
	mat.dts.i = cumsum(2 * famCounts) - 2 * famCounts + 1;
	cfgs.s = configCounts * sapply(famCounts, function(i)(4 + (i - 2) * 8));
	cfgs.i = cumsum(cfgs.s) - cfgs.s + 1;
	cfgs.i = c(cfgs.i, last(cfgs.i) + last(cfgs.s));

	for (i in 1:n) {
		new.fam = draw.fam.from(multi,
			dts$diplotypes[cfgs.i[i]:(cfgs.i[i + 1] - 1)], famCounts[i]);
		imputations = c(imputations, new.fam);
	}
 	counts = count.haplotypes.fam(imputations, famCounts, countLoci);

	# create first instance of beta
	if (dont.update.beta) beta = free.pars$beta
	else beta = rnorm(1, prior$beta$mean, prior$beta$variance);

	# prepare likelihood
	penetrance = get(sprintf("penetrance.%s", model.est));

	state.0 = list(
		diplotypes = dts, phenotypes = pts,
		famCounts = famCounts, countLoci = countLoci, prior = prior,
		imputations = imputations, counts = counts, cfgs.i = cfgs.i, mat.dts.i = mat.dts.i,
		pts.i = pts.i,
		penetrance = penetrance,
		baseline = fixed.pars$baseline,
		countLoci = countLoci,
		beta = beta,
		beta.prop = merge.lists(list(switch.cycle = 20, var = 2, var.late = .1,
			updates = ifelse(dont.update.beta, 0, 1)
		), fixed.pars$beta.prop),

		dont.update.beta = dont.update.beta,
		lh.hstar = get.lh.hstar(countLoci, baseline, penetrance),
		free.pars = hfreqs,
		cycle = 0,

		# rejection counts
		beta.rej = 0,
		h.rej = 0
	);
	state.0
}

run.chain.full = function(d, free.pars, fixed.pars,
	prior = list(haplotypes = rep(1e-3, 2^d$pars$fixed.pars$countLoci),
				 beta = list(mean = 1.5, variance = 3)),
	dont.update.beta = F) {
	f = fixed.pars;
	state.0 = mcmc.full.initialize(fixed.pars, free.pars,
		d$sim$genotypes, d$sim$phenotypes, d$sim$familySizes,
		f$countLoci, prior = prior, dont.update.beta
	);
	state = state.0;

	results = list();
	m = (f$burn.in + f$mcmc.cycles);
	for (i in 1:m) {
		#print(list(cycle = i));
		state = mcmc.full.update(state, i);
		results[[i]] = list(beta = state$beta, hfreqs = state$hfreqs);
	}
	betas = sapply(1:m, function(i){ results[[i]]$beta });
	hfreqs = t(apply(t(t(1:m)), 1, function(i){results[[i]]$hfreqs }));
	ret = list(
		betas = betas, rej.betas = state$beta.rej,
		hfreqs = hfreqs, rej.hfreqs = state$h.rej,
		beta.prop = state.0$beta.prop
	);
	ret
}


estimator.mcmc = function(ctxt, d, free.pars, fixed.pars, no) {
	estimates = list();
	for (i in 1:d$cycles) {
		data = (d = extractDataSetFromDataSets(d, i))$data;
		estimates[[i]] = run.chain.full(data, free.pars, fixed.pars);
	}
	estimates
}

plot.estimates.mcmc = function(pars, names = paste("pars", 1:dim(pars)[2], sep = ""),
	path = "/tmp/mcmc-test-1.pdf", interval = 10) {
	pdf(path, height = 12);
	n.p = dim(pars)[2];
	layout(matrix(1:n.p));
	picks = seq(from = interval, to = dim(pars)[1], by = interval);
	for (i in 1:n.p) {
		plot(picks, pars[picks, i], type = "l", xlab = "Iteration", ylab = names[i]);
	}
	dev.off();

}

print.estimates.mcmc = function(pars, names, path, pars.true, pars.emp, fixed.pars) {
	pars.real = pars[(fixed.pars$burn.in + 1):(fixed.pars$burn.in + fixed.pars$mcmc.cycles - 1), ];
#	print(pars.real);
	pars.mcmc = apply(pars.real, 2,
		function(r){ c(mean(r, na.rm = T), as.vector(quantile(r, probs = c(0.025, 0.975), na.rm = T))) }
	);
#	print(pars.mcmc[, 1]);
	if (!is.null(path)) {
		print(col.frame(list(
			True = as.vector(pars.true), Sample = pars.emp, "Estimates" = pars.mcmc[1, ],
			"KI-lower" = pars.mcmc[2, ], "KI-upper" = pars.mcmc[3, ],
			Difference = pars.mcmc[1, ] - pars.true,
		), col.name = names));
	}
	ret = as.vector(rbind(as.vector(pars.true), pars.emp, pars.mcmc[1, ], pars.mcmc[2, ], pars.mcmc[3, ]));
}

print.and.plot.mcmc = function(pars, names = paste("pars", 1:dim(pars)[2], sep = ""),
	path.prefix = "/tmp/mcmc-test-1", interval = 1, pars.true, pars.emp, fixed.pars) {

	# <!> reactivate
	plot.estimates.mcmc(pars, names, path = sprintf("%s-plot.pdf", path.prefix), interval);
	print.estimates.mcmc(pars, names, path = NULL, pars.true, pars.emp, fixed.pars);
}

summary.mcmc = function(ctxt, d, free.pars, fixed.pars, no) {
	d = (d = extractDataSetFromDataSets(d, 1))$data;
	syms = try(load(ctxt$output.name));	# expect estimates
	if (class(syms) == "try-error") return(NULL);
	r = simulation[[1]];
	#rep = Rreporter(path = sprintf("%s/%s.rr", output.dir, output.prefix), to.screen = F);

	print(col.frame(fixed.pars, minus = c("models")));
	print(unlist(free.pars));
	mcmc.count = (fixed.pars$burn.in + fixed.pars$mcmc.cycles);
	rej.fracs = list(
		rejections.dts = r$rej.hfreqs / (mcmc.count * fixed.pars$count),
		rejections.beta = r$rej.betas / (mcmc.count * r$beta.prop$updates)
	);
	print(col.frame(rej.fracs));

	# print likelihood at true parameters
	true.pars = c(ld.f2hFreq(free.pars), free.pars$beta);
	pars.emp = c(table(d$sim$diplotypes) / sum(table(d$sim$diplotypes)), free.pars$beta);

	# interval: count of beta updates per hfreq-update
	parsHap = print.and.plot.mcmc(cbind(r$hfreqs, r$betas),
		names = c(sapply(1:dim(r$hfreqs)[2],
			function(i)sprintf("eta^s_%d", i)), "beta"),
		path.prefix = splitPath(ctxt$output.name)$fullbase,
		interval = r$beta.prop$updates,
		c(ld2hFreq(complete.ld(free.pars)), free.pars$beta),
		pars.emp, fixed.pars
	);

	ld.freqs = t(apply(r$hfreqs, 1, hFreq2r2.free));
	ld.pars.emp = hFreq2r2.free(pars.emp[-length(pars.emp)]);
	ld.names = c(sapply(1:(dim(r$hfreqs)[2]/2 - 1),
		function(i)sprintf("eta_%d", i)),
		sapply(1:(dim(r$hfreqs)[2]/2 - 1),
			function(i)sprintf("r^2_%d", i)), "p");
	parsLd = print.and.plot.mcmc(cbind(ld.freqs, r$betas),
		names = c(ld.names, "beta"),
		path.prefix = sprintf("%s%s", splitPath(ctxt$output.name)$fullbase, "-ld"),
		interval = r$beta.prop$updates,
		unlist(free.pars), c(ld.pars.emp, free.pars$beta), fixed.pars
	);
	#r$finalize();
	ret = list(parsHap = parsHap, parsLd = parsLd);
	ret
}

summaryOld.mcmc.finalize = function(ctxt, free.pars, fixed.pars, final.results) {
	h = ctxt$hints;
	#pars = parCombinations(fixed.pars, free.pars, h$fixed.keys, h$names.fixed.keys);
	countLoci = fixed.pars[[1]]$countLoci;

	r = Rreporter(path = sprintf("%s/%s.rr", ctxt$output.dir, ctxt$prefix), to.screen = F);
	for (i in 1:length(final.results)) {

		df = final.results[[i]]$parsHap;
		dfNames = c(paste("\\eta_", 1:(2^countLoci), sep = ""), "p^\\star", "\\beta");
		dfNames = sapply(dfNames, function(n)sprintf("$%s$", n));
		names(df) = dfNames;
		r$report.data.frame(df,
			row.formatters = c(function(v, dg=1){if (is.na(to.numeric(v))) v else sprintf("%.*f", dg, v)}),
			digits = 2,
			caption = sprintf("Simulation no %d", i)
		);

		df = final.results[[i]]$parsLd;
		dfNames = c(paste("\\eta_", 1:(2^(countLoci - 1) - 1), sep = ""),
			paste("R^2_", 1:(2^(countLoci - 1) - 1), sep = ""), "p^\\star", "\\beta");
		dfNames = sapply(dfNames, function(n)sprintf("$%s$", n));
		names(df) = dfNames;
		r$report.data.frame(df,
			row.formatters = c(function(v, dg=1){if (is.na(to.numeric(v))) v else sprintf("%.*f", dg, v)}),
			digits = 2,
			caption = sprintf("Simulation no %d", i)
		);
	}
	r$finalize(final.path = sprintf("%s/%s.pdf", ctxt$output.dir, ctxt$prefix));
	
}

summary.mcmc.finalize = function(ctxt, free.pars, fixed.pars, final.results) {
	parsHap = lapply(final.results, function(l)l$parsHap);
	print(parsHap);
	parameterNames = c(paste("\\eta^\\star_", 1:length(unlist(free.pars[[1]])), sep = ""), "\\beta");

	summarizeEstimates(ctxt, free.pars, fixed.pars, parsHap,
		# 4 cols for each par: mean, median, bias, se
		rowFormatter = function(r){
			n.r = sapply(5*(0:(length(r)/5 - 1)), function(j) {
				sprintf("%.2f (%.2f)", r[j + 3], r[j + 2]) });
			n.r1 = sapply(5*(0:(length(r)/5 - 1)), function(j) { 
				sprintf("(%.2f, %.2f)", r[j + 4], r[j + 5]) });
			n.r = c(n.r, n.r1);
			n.r
	}, splitRow = 2, parameterNames, reParam = function(p){c(ld.f2hFreq(p), p$beta)});
}
