/*
 * Valarray.h
 * Wed Aug 13 13:59:34 CEST 2014
 */

#ifndef VALARRAY_EXT_H
#define VALARRAY_EXT_H

#include <valarray>
#include <math.h>
#include <vector>
#include <iostream>
#include <iomanip>
using namespace std;


template<class T> class Valarray : public valarray<T>
{
public:
	Valarray(const valarray<T> other) : valarray<T>(other) {}
// 	Valarray(const vector<T> other) : valarray<T>((T)0, other.size()) {
// 		for (int i = 0; i < other.size(); i++) (*this)[i] = other[i];
// 	}
	Valarray(const vector<T> &other) : valarray<T>((T)0, other.size()) {
		for (int i = 0; i < other.size(); i++) (*this)[i] = other[i];
	}
	Valarray(T val, int size) : valarray<T>(val, size) {}
	// R code
	//logSumExpRaw = function(v, pivot = median(v))(log(sum(exp(v - pivot))) + pivot)
	// <o>
	// remove max, then r = pivot + log1p((((*this) - pivot).exp()).sum())
	inline T	logSumExpRaw(const T pivot) const {
		// <o> optimize: component-wise computation
		T	r = pivot + log((((*this) - pivot).exp()).sum());
		return r;
	}
	// assume array is ordered in ascending order
	inline T	logSumExpSorted(void) const {
		// assert(this->sorted());	// <i>
		return logSumExpRaw(this->back());	// <!> last element is biggest
	}
	inline T	logSumExpSortedInv(void) const {
		// assert(this->rev().sorted());	// <i>
		return logSumExpRaw(this->front());	// <!> last element is biggest
	}
	// garuantees sorted output if all entries positive
	inline Valarray<T>	cumsum(void) const {
		Valarray<T>	r(0, this->size());
		T			acc = 0;

		for (int i = 0 ; i < this->size(); i++) r[i] = (acc += (*this)[i]);
		return r;
	}
	// garuantees sorted output if all entries positive
	inline Valarray<T>	log_cumsum(void) const {
		T			pivot = this->max();
		Valarray<T>	r(exp(*this - pivot));
		T			acc = 0;

		for (int i = 0 ; i < this->size(); i++) r[i] = log(acc += r[i]);
		r += pivot;
		return r;
	}
	/*
	 * return index for which element in valarray is infimum of elements bigger than value
	 * might point at end if value > all elements in valarray
	 * alternative interpretation: where to insert value to keep valarray sorted
	 * insertion into sequences of identical values is unstable: the position is undefined
	 */
	inline int	binary_search(const T value) const {
		int	min = 0, max = this->size() - 1;

		while (max >= min) {
			int	mid = min + (max - min)/2;
			if ((*this)[mid] == value) return mid;
			if ((*this)[mid] < value) min = mid + 1;
			else max = mid - 1;
		}
		return min;
	}
	inline bool	exists(const T value) const {
		int	pos = binary_search(value);
		return pos < this->size() && this[pos] == value;
	}
	inline T	back(void) const { return this->back(); }
	inline void	print(void) const {
		cout << setprecision(4);
		for (int i = 0 ; i < this->size(); i++) cout << (i? ", ": "") << (*this)[i];
		cout << endl;
	}
};

template<class T> inline Valarray<T> log(const Valarray<T>& a) { return Valarray<T>(log((valarray<T>)a)); }
template<class T> inline void print(const Valarray<T>& a) { a.print(); }

#endif // VALARRAY_EXT_H
