/*
 * helpers.h
 * Fri Jul 11 17:35:19 CEST 2014
 */

#ifndef CARTESIAN_H
#define CARTESIAN_H

#include <vector>
using namespace std;

template<typename T> class Cartesian;

template<typename T>
class CartesianIterator  : public vector<T> {

protected:
	Cartesian<T>		counts;

public:
	CartesianIterator() : Cartesian<T>(), vector<T>() {}
	CartesianIterator(Cartesian<T> &_counts, vector<T> &_indeces) : counts(_counts), vector<T>(_indeces) {}
	CartesianIterator(Cartesian<T> &_counts) : counts(_counts), vector<T>(0, _counts.size()) {}
	CartesianIterator(vector<T> &_counts) : counts(_counts), vector<T>(0, _counts.size()) {}
    CartesianIterator(const CartesianIterator& other) : counts(other.counts), vector<T>(other) {}
    ~CartesianIterator() {}
    CartesianIterator& operator=(const CartesianIterator& other) {
		*(vector<T> *)this = (vector<T>)other;
		return *this;
	}

	bool operator==(const CartesianIterator &other) {
		*(vector<T> *)this == (vector<T>)other;
	}
	bool operator++(void) {
		int i = 0;
		while (i < this->size() && (*this)[i] == counts[i] - 1) i++;
		if (i == this->size()) return false;	// exhaustion of iterator
		(*this)[i]++;
		while (--i >= 0) (*this)[i] = 0;
		return true;
	}
	bool	isExhausted(void) {
		for (int i = 0; i < this->size(); i++) if ((*this)[i] < counts[i] - 1) return false;
		return true;
	}
};

template<typename T>
class Cartesian : public vector<T>
{
public:
    Cartesian(vector<T> &_counts) : vector<T>(_counts) {}
    Cartesian(const Cartesian& other) : vector<T>(other) {}
    ~Cartesian() {}

	CartesianIterator<T>	*begin(void) {
		vector<T>	start(0, this->size());
		return new CartesianIterator<T>(*this, start);
	}
	CartesianIterator<T>	*end(void) {
		vector<T>	stop(*this);
		return new CartesianIterator<T>(*this, stop);
	}
};

#endif // CARTESIAN_H
