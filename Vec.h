/*******************************************************************************************[Vec.h]
Copyright (c) 2003-2007, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/
/* Copyright (c) 2020, Muhammad Osama */
/**************************************/

#ifndef __VEC_
#define __VEC_

class MEMOUTEXCEPTION {};
inline void* mem_out(void* ptr, size_t size)
{
	void* mem = realloc(ptr, size);
	if (mem == NULL && errno == ENOMEM)
		throw MEMOUTEXCEPTION();
	else
		return mem;
}

template<class T, class SZ_T = int>
class Vec {
private:
	T* data;
	SZ_T sz;
	SZ_T cap;

	Vec<T>& operator=(Vec<T>& rhs);
	Vec(Vec<T>& rhs);

	inline SZ_T max(SZ_T x, SZ_T y) {
		return (x > y) ? x : y;
	}

public:
	Vec() {
		data = NULL;
		sz = 0;
		cap = 0;
	}
	explicit Vec(SZ_T size) {
		data = NULL;
		sz = 0;
		cap = 0;
		incMem(size);
	}
	Vec(SZ_T size, const T& val) {
		data = NULL;
		sz = 0;
		cap = 0;
		incMem(size, val);
	}
	~Vec() { clear(true); }

	// pointer to raw data
	operator T* (void) {
		return data;
	}
	T* d_ptr(void) {
		return data;
	}
	// indexing operators:
	const T& operator [] (SZ_T index) const {
		assert(index < sz);
		return data[index];
	}
	T& operator [] (SZ_T index) {
		assert(index < sz);
		return data[index];
	}

	// memory 
	void capacity(SZ_T min_cap);
	void incMem(SZ_T size);
	void incMem(SZ_T size, const T& val);
	void clear(bool dealloc = false);
	bool empty() { return sz == 0; }
	SZ_T size(void) const { return sz; }
	int capacity(void) const { return cap; }
	void resize(SZ_T n) {
		if (n == sz) return;
		if (n < sz) shrink(sz - n);
		else incMem(n);
	}
	void shrink(SZ_T n) {
		assert(n <= sz);
		for (SZ_T i = 0; i < n; i++) {
			sz--;
			data[sz].~T();
		}
	}

	// new elements
	void push(void) {
		if (sz == cap) capacity(sz + 1);
		new (&data[sz]) T();
		sz++;
	}

	void push(const T& val) {
		if (sz == cap) capacity(sz + 1);
		new (&data[sz++]) T(val);
	}

	void pop(void) {
		assert(sz > 0);
		sz--;
		data[sz].~T();
	}

	void copyFrom(Vec<T>& copy) const {
		assert(copy.size() <= sz);
		for (SZ_T i = 0; i < sz; i++) data[i] = copy[i];
	}

	template<class CSZ_T>
	void copyFrom(T* copy, const CSZ_T& copy_sz) const {
		assert(copy_sz <= sz);
		for (SZ_T i = 0; i < sz; i++) data[i] = copy[i];
	}
};

template<class T, class SZ_T>
void Vec<T, SZ_T>::capacity(SZ_T min_cap) {
	if (cap >= min_cap) return;
	SZ_T add = max((min_cap - cap + 1) & ~1, ((cap >> 1) + 2) & ~1);
	const SZ_T size_max = std::numeric_limits<SZ_T>::max();
	if (((size_max <= std::numeric_limits<int>::max()) && (add > size_max - cap))
		|| (((data = (T*)::realloc(data, (cap += add) * sizeof(T))) == NULL) && errno == ENOMEM))
		throw MEMOUTEXCEPTION();
}

template<class T, class SZ_T>
void Vec<T, SZ_T>::incMem(SZ_T size, const T& val) {
	if (sz >= size) return;
	capacity(size);
	for (SZ_T i = sz; i < size; i++) data[i] = val;
	sz = size;
}

template<class T, class SZ_T>
void Vec<T, SZ_T>::incMem(SZ_T size) {
	if (sz >= size) return;
	capacity(size);
	for (SZ_T i = sz; i < size; i++) new (&data[i]) T();
	sz = size;
}

template<class T, class SZ_T>
void Vec<T, SZ_T>::clear(bool dealloc) {
	if (data != NULL) {
		for (SZ_T i = 0; i < sz; i++) 
			data[i].~T();
		sz = 0;
		if (dealloc) {
			free(data);
			data = NULL;
			cap = 0;
		}
	}
}

#endif
