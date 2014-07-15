/*
 * pedigree.cpp
 * Thu Jul 10 16:19:32 CEST 2014
 */

#include "pedigree.h"

Pedigree::Pedigree(IntegerVector &_founder, IntegerMatrix &_itrio) : founder(_founder), itrio(_itrio) {
	//for (int i = 0; i < _founder.size(); i++) founder[i] = _founder[i];
	//founder = vectorConvert<int, iid_t>(_founder(_founder.begin(), _founder.end()));
// 	itrio.resize(_itrio.nrow());
// 	for (int i = 0; i < _itrio.nrow(); i++) {
// 		itrio[i] = vector<iid_t>(3);
// 		//itrio[i].swap(vector<iid_t>(3));	// does not work
// 		for (int j = 0; j < _itrio.ncol(); j++)
// 			itrio[i][j] = itrio[i][j];
// 
// 	}
}

