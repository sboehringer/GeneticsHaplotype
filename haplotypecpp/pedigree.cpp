/*
 * pedigree.cpp
 * Thu Jul 10 16:19:32 CEST 2014
 */

#include "pedigree.h"

Pedigree::Pedigree(vectorR<iid_t> &_founder, matrixR<iid_t> &_itrio)
: founder(_founder), itrio(_itrio) {
	
}
Pedigree::Pedigree(vector<iid_t> &_founder, IntegerMatrix &_itrio)
: founder(_founder), itrio(_itrio) {
	
}

Pedigree::Pedigree(vector<iid_t> &_founder, vector< vector<iid_t> > &_itrio)
	: founder(_founder), itrio(_itrio) {
	
}
