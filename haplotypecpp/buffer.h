/*
 * Copyright (c) 2014, Stefan Boehringer <email>
 * All rights reserved.
 *
 * buffer.h
 * Fri Jul 11 15:48:18 CEST 2014
 *
 */

#ifndef BUFFER_H
#define BUFFER_H

#include <cstdlib>
#include <cerrno>
#include <iostream>
using namespace std;

extern "C" {
#include <memory.h>
}

template <typename T>
class Buffer
{
	size_t	N;
	size_t	capacity;
	size_t	elementSize;
	T		*_buffer;
public:
    Buffer() : elementSize(1), capacity(0), _buffer(0), N(0) {}
	Buffer(size_t _elementSize, size_t _capacity)
	: elementSize(_elementSize), capacity(_capacity), N(0), _buffer((T *)calloc(elementSize, capacity)) {
		//cout << "ElementSize: " << elementSize << " Capacity: " << capacity << endl;
		if (!_buffer) throw("Out of memory.");
	}
    Buffer(const Buffer& other) {
		this->capacity = other.capacity;
		this->elementSize = other.elementSize;
		this->N = other.N;

		if (_buffer) free(buffer);
		_buffer = (unsigned char *)calloc(elementSize, capacity);
		if (!_buffer) throw("Out of memory.");
		memcpy(_buffer, other.buffer(), elementSize * capacity);
	}
    ~Buffer() {
		if (_buffer) free(_buffer);
	}

	size_t	size(void) const { return(N); }
	T 		*buffer(size_t i = 0) const { return((T *) ((unsigned char *)_buffer + i * elementSize)); }
	T		*operator[](int i) {
		return (T *)((char *)_buffer + elementSize * i);
	}
	T		*resize(size_t _N) {
		if (_N <= capacity) {
			N = _N;
			return (*this)[N - 1];
		}
		capacity <<= 1;
		_buffer = (T *)realloc(_buffer, capacity * elementSize);
		if (!_buffer) throw("Out of memory.");
		N = _N;
		return (*this)[N - 1];
	}	
	void	shrink(void) {
		capacity = N;
		_buffer = (T *)realloc(_buffer, capacity * elementSize);
		if (ENOMEM == errno) throw("Reallocation exception.");
	}
	T	*push(void) {
		return this->resize(N + 1);
	}
	size_t	size(void) { return N; }

};

#endif // BUFFER_H
