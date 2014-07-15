/*
 * pedigree.h
 * Thu Jul 10 16:07:16 CEST 2014
 * 
 * classes for collections of pedigrees (families)
 * individual pedigrees
 * 
 * interface to reconstruction
 */

#ifndef PEDIGREE_H
#define PEDIGREE_H

#include <vector>
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

#include "genotype.h"

template<class TF, class TT> TT	identity(TF from) { return (TT)from; }
// convert vector type from class TF to class TT
template<class TF, class TT> vector<TT>	*vectorConvert(const vector<TF> &vf) {
	vector<TT>	&vt = * new vector<TT>();
	vt.resize( vf.size() );
	transform( vf.begin(), vf.end(), vt.begin(), identity<TF, TT> );
	return &vt;
}

typedef vector<iid_t>	iidVector_t;

template<class T> class vectorR : public vector<T> {
public:
	vectorR<T>() : vector<T>() {}
	vectorR<T>(IntegerVector &v) {
		this->reserve(v.length());
		for (int i = 0; i < v.length(); i++) (*this)[i] = v[i];
	}
};

template<class T> class matrixR : public vector< vectorR<T> > {
public:
	matrixR<T>() : vector< vectorR<T> >() {}
	matrixR<T>(IntegerMatrix &m) : vector< vectorR<T> >(m.nrow()) {
		this->reserve(m.nrow());
		for (int i = 0; i < m.nrow(); i++) {
			IntegerVector	iv = m.row(i);
			(*this)[i] = iv;
		}
	}
};

class Pedigree {
	const vectorR<iid_t>	founder;	// iids of founders
	const matrixR<iid_t>	itrio;		// inheritance trios: vector[0][0 .. 2] iid, mid, pid of trio 0, ...

	public:
	/*
	 * initialization, boilerplate	
	 */
	Pedigree() : founder(), itrio() {}
	Pedigree(IntegerVector &_founder, IntegerMatrix &_itrio);
	~Pedigree();

	/*
	 * pedigree methods
	 */
	iid_t	sizeFounders(void) { return founder.size(); }
	iid_t	sizeItrios(void) { return itrio.size(); }
	const vector<iid_t>	&founders(void) { return founder; }
	inline iid_t	trioIid(iid_t i) { return itrio[i][0]; }
	inline iid_t	trioMid(iid_t i) { return itrio[i][1]; }
	inline iid_t	trioPid(iid_t i) { return itrio[i][2]; }
};

class PedigreeCollection {
	vector<Pedigree>	peds;
	
	public:
	PedigreeCollection(List &pedvector);
	~PedigreeCollection() {}

	void	print(void);

	//DiplotypeReconstruction		*reconstruct(void);	// ids of all individuals in the pedigrees
};

#endif
