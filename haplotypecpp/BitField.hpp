/*
	BitField.hpp
	Wed Apr 18 14:20:50 CEST 2012
 */

#if !defined(__BITFIELD_HPP__)
#define __BITFIELD_HPP__

#include	<vector>
#include	<iostream>
#include	<stdlib.h>
#include	<assert.h>
using namespace std;
//#include	"bitArray.hpp"
/*
 * bit-mask functions
 */

#define	bitsPerByte	8
#define	typeBits(T)	(sizeof(T) * bitsPerByte)
// bit mask setting to 1 bit 0 through bit n - 1 and to 0 all other bits
//	takes into account the special case MaxBits that is not naturally handled on the x86 architecture
#define safeMask(n, T)		( (n) <= 0? 0: ( (n)>=typeBits(T)? (T)-1: (((T)1<<(n)) - (T)1) ) )
// shortcut for safeMask
#define	Mask(n,T)			safeMask(n, T)
// natural version working on 680x0 architectures
//	define	Mask(n)		( (1<<(n))-1 )
//	same as mask
#define	RightMask(n, T)		Mask(n, T)
//	the bits 0...n-1 are zero and the rest is 1
#define	LeftMask(n, T)			(~Mask(n, T))

template<typename T, typename I> inline T	BitMask(I n) {
	return (n) <= 0? 0: ( (n)>=typeBits(T)? (T)-1: (((T)1 << n ) - (T)1) );
}

// round value up to the next greatest power of 2
//	__builtin_clz: count leading zeros (left 0 bits)
template <class T, class I> inline T	ceilP2(T value, I power2) {
	// count of leading zeros
	I	clz;
	if (sizeof(T) == 4)			clz = __builtin_clz(value);
	else if (sizeof(T) == 8)	clz = __builtin_clzl(value);
	else if (sizeof(T) == 16)	clz = __builtin_clzll(value);
	if (clz == typeBits(T)) return 0;
	// floor log2 of value
	I	pv = bitsPerByte * sizeof(T) - clz;
	
	return (value + RightMask(pv - 1, T)) & LeftMask(pv, T);
}

template <typename T>
inline T roundUp(T v, T modulo) {
//	std::cout << "roundUp: " << v << " modulo: " << modulo << " == " << v + modulo - v % modulo << std::endl;
	return (v + modulo - v % modulo);
}

template <typename T>
inline T roundUpBitsToBytes(T bits) {
	return roundUp<T>(bits, typeBits(T)) / bitsPerByte;
}


/*
 * The bitfield is interpreted as having units, i.e. a certain basic type (int8 == char, int16 == short, ..., int128),
 * within which bits are numbered 0, ..., typeBits(T) -1 and read from right to left (least to most significant).
 * The LeftMask therefore has the most significant bits set, whereas the RightMask least significant bits
 * 
 * class I: type of index variables: has to be signed
 * class B: type of the buffer: values 
 * 
 * ToDo: signed values using <limits>, numeric_limits<B>::is_signed;
 */
// B: buffer type, I index type
template <class B, class I>
class BitField {
protected:
	B	*buffer;
	// buffer size in bits
	I	bufferSize;

	public:
	BitField(B *_buffer, I _bufferSize) : bufferSize(_bufferSize), buffer(_buffer) { }
	BitField(B *_buffer) : bufferSize(1), buffer(_buffer) { }
	inline void	set(I from, I size, B value) {
		assert(size <= typeBits(B) && (from + size) <= bufferSize);
		I	u = from / typeBits(B), b = from % typeBits(B);	// number of unit, number of bit
		// number of least significant bits in the first part
		I	sl = size >= (typeBits(B) - b)? (typeBits(B) - b): size;
		// number of most significant bits in the second part
		I	sm = size - sl;
		// set least significant part
		buffer[u] = (buffer[u] & ~(RightMask(sl, B) << b)) | ((value & RightMask(sl, B)) << b);
		// set most significant part
		if (sm > 0) buffer[u + 1] = (buffer[u + 1] & LeftMask(sm, B)) | ((value >> sl) & RightMask(sm, B));
	}
	inline B	get(I from, I size) {
		assert(size <= typeBits(B) && (from + size) <= bufferSize);
		I	u = from / typeBits(B), b = from % typeBits(B);	// number of unit, number of bit
		// number of least significant bits in the first part
		I	sl = size >= (typeBits(B) - b)? (typeBits(B) - b): size;
		// number of most significant bits in the second part
		I	sm = size - sl;
		// get least significant part
		B	value = ((buffer[u] >> b) & RightMask(sl, B));
		// set most significant part
		if (sm > 0) value |= (buffer[u + 1] & RightMask(sm, B)) << sl;
		return value;
	}
};

template <class B, class I,  class Allocator = allocator<B> >
class BitFieldBuffer : public BitField<B, I> {
	public:
	//BitFieldBuffer(I _bufferSize) : bufferSize(_bufferSize), buffer(Allocator(ceilP2<B, I>(bufferSize) / typeBits(B))) { }
	BitFieldBuffer(I _bufferSize) {
		BitField<B, I>(Allocator(ceilP2<B, I>(_bufferSize) / typeBits(B)), _bufferSize);
	}
};

/*
 * Bitarray allows to access bitfields of equal size arranged in tandem in the buffer
 * offset allows to shift that array along the buffer
 */
typedef enum { BitArrayAfter_e } BitArrayAfter_t;
template <class B, class I>
class BitArray : public BitField<B, I> {
	I	offset;
	I	elementSize;
	I	count;

	public:
	/*
	 * construction/destruction
	 */
	BitArray(B *_buffer, I _bufferSize, I _elementSize) :
		BitField<B, I>(_buffer, _bufferSize),
		offset(0),
		elementSize(_elementSize),
		count(_bufferSize/_elementSize) { }
	BitArray(B *_buffer, I _offset, I _elementSize, I _count) :
		BitField<B, I>(_buffer, _offset + _elementSize * _count),
		offset(_offset),
		elementSize(_elementSize),
		count(_count) { }
	BitArray(BitArrayAfter_t _t, BitArray<B, I> &o, I _elementSize, I _count) :
		BitField<B, I>(o.buffer, o.offset + o.elementSize * o.count + _elementSize * _count),
		offset(o.offset + o.elementSize * o.count),
		elementSize(_elementSize),
		count(_count) { }
	~BitArray() {}

	/*
	 * methods
	 */
	inline B	operator[](const I i) { return this->get(offset + i * elementSize, elementSize); }
	inline void	set(const I i, B v) { BitField<B, I>::set(offset + i * elementSize, elementSize, v); }
};

/*
 * BitField helpers
 */

typedef	int	bitidx_t;

// http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
inline int countBitsSet(unsigned int i) {
    i = i - ((i >> 1) & 0x55555555U);
    i = (i & 0x33333333U) + ((i >> 2) & 0x33333333U);
    return (((i + (i >> 4)) & 0x0F0F0F0FU) * 0x01010101U) >> 24;
}
inline int	bitCountSet(unsigned int i, bitidx_t maxBits) {
	return countBitsSet(i & ( (1 << maxBits) - 1));
}
inline int bitCountUnset(unsigned int i, bitidx_t maxBits) {
	bitidx_t		bitsInT = sizeof(int) * 8;
	unsigned int	mask = ~(unsigned int)0 - ( (1 << maxBits) - 1);
	return bitsInT - countBitsSet(i | mask);
}

template<class T> inline T	bitSet(T v, bitidx_t i) {
	return v | (((T)1) << i);
}
template<class T> inline T	bitClear(T v, bitidx_t i) {
	return v & (~ (((T)1) << i));
}
template<class T> inline T	bitAt(T v, bitidx_t i) {
	return !!(v & (((T)1) << i));
}
template<class T> inline T	bitEmbedValueInZeros(T d, T v, bitidx_t maxBits) {
	//bitidx_t	bitsToSet = bitCountUnset(d, maxBits);

	for (bitidx_t i = 0, j = 0; j < maxBits; i++) {
		if (!bitAt<T>(d, i)) {
			if (bitAt<T>(v, j)) d = bitSet<T>(d, i);
			j++;
		}
	}
	return d;
}
template<class T> inline T	bitEmbedValueInOnes(T d, T v, bitidx_t maxBits) {
	//bitidx_t	bitsToSet = bitCountSet(d, maxBits);

	for (bitidx_t i = 0, j = 0; j < maxBits; i++) {
		if (bitAt<T>(d, i)) {
			if (!bitAt<T>(v, j)) d = bitClear<T>(d, i);
			j++;
		}
	}
	return d;
}

#endif // __BITFIELD_HPP__
