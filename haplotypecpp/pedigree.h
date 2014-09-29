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
	vectorR<T>(const vector<T> &other) : vector<T>(other) {}
	vectorR<T>(const vectorR<T> &other) : vector<T>((vector<T>)other) {}
	vectorR<T>(const IntegerVector &v) : vector<T>(Rcpp::as<vector< T> >) {}
	//inline T	operator[](int i) { return (*this)[i]; }
};
#if 0
template<class T> class matrixR : public vector< vectorR<T> > {
public:
	matrixR<T>() : vector< vectorR<T> >() {}
	matrixR<T>(vector< vector<iid_t> > &other) : vector< vectorR<T> >(other.size()) {
		this->reserve(other.size());
		for (int i = 0; i < other.size(); i++) {
			(*this)[i] = other[i];
		}
	}
	matrixR<T>(IntegerMatrix &m) : vector< vectorR<T> >(m.nrow()) {
		this->reserve(m.nrow());
		for (int i = 0; i < m.nrow(); i++) {
			IntegerVector	iv = m.row(i);
			(*this)[i] = iv;
		}
	}
};
#endif
template<class T> class matrixR : public vector< T > {
	int	nrow_;
	int	ncol_;
public:
	matrixR<T>() : nrow_(0), ncol_(0), vector<T>() {}
	matrixR<T>(int _nrow, vector< iid_t > &other)
	: vector< T >(other), nrow_(_nrow), ncol_(other.size() / nrow_) { }
	matrixR<T>(const vector< vector<T> > &other)
	: vector< T >(other.size()? other.size() * other[0].size(): 0),
	  nrow_(other.size()),
	  ncol_(nrow_? other[0].size(): 0) 
	{
		this->reserve(other.size());
		for (int r = 0; r < nrow_; r++)
			for (int c = 0; c < ncol_; c++)
				(*this)[r + c*nrow_] = other[r][c];
	}
	matrixR<T>(const IntegerMatrix &m)
	: vector<T>(Rcpp::as<vector< T > >(m)), nrow_(m.nrow()), ncol_(m.ncol()) { }
	matrixR<T>(const matrixR<T> &m)
	: vector<T>( (vector< T >)m ), nrow_(m.nrow()), ncol_(m.ncol()) { }

	inline T	at(int r, int c) const { return (*this)[r + nrow_*c]; }
	inline int	nrow(void) const { return this->nrow_; }
	inline int	ncol(void) const { return this->ncol_; }

	inline matrixR<T>& operator=(const matrixR<T> &other) {
		((vector<T>)(*this)) = (vector<T>)other;
		nrow_ = other.nrow_;
		ncol_ = other.ncol_;
		return *this;
	}
};

class Pedigree {
	vectorR<iid_t>	founder;	// iids of founders
	matrixR<iid_t>	itrio;		// inheritance trios: vector[0][0 .. 2] iid, mid, pid of trio 0, ...

	public:
	/*
	 * initialization, boilerplate	
	 */
	Pedigree() : founder(), itrio() {}
	Pedigree(vectorR<iid_t> &_founder, matrixR<iid_t> &_itrio);
	Pedigree(vector<iid_t> &_founder, vector< vector<iid_t> > &_itrio);
	Pedigree(vector<iid_t> _founder, vector< vector<iid_t> > _itrio);
	Pedigree(const Pedigree &other) : founder(other.founder), itrio(other.itrio) {}
	~Pedigree() {}

	/*
	 * pedigree methods
	 */
	Pedigree &operator=(const Pedigree &other) {
		founder = other.founder;
		itrio = other.itrio;
		return *this;
	}
	iid_t	sizeFounders(void) const { return founder.size(); }
	iid_t	sizeItrios(void) const { return itrio.nrow(); }
	inline iid_t	N(void) const { return sizeFounders() + sizeItrios(); }
	vector<iid_t>	&founders(void) const { return *(vector<iid_t> *)&founder; }
	inline iid_t	trioIid(iid_t i) const { return itrio.at(i, 0); }
	inline iid_t	trioMid(iid_t i) const { return itrio.at(i, 1); }
	inline iid_t	trioPid(iid_t i) const { return itrio.at(i, 2); }
};

#endif
