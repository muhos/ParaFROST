/******************************************************************************************[Vec.h]
 ParaFROST -- Copyright (c) 2020, Muhammad Osama - Anton Wijs

ParaFROST "Vec.h" header file is a modified version, based on MiniSat (see below MiniSat copyrights).
Therefore, permissions and copyrights of this file ONLY (Vec.h) are exactly the same as Minisat.
ParaFROST other sources are based on another copyright. This program is free software: 
you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  
If not, see <https://www.gnu.org/licenses/>. 

- The above and below copyrights notices and this permission notice shall be included in all
copies or substantial portions of the Software.

- Original Minisat Copyrights:

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
	SZ_T sz, cap;

	Vec<T>& operator=(Vec<T>& rhs);
	Vec(Vec<T>& rhs);
	inline SZ_T max(SZ_T x, SZ_T y) { return (x > y) ? x : y; }

public:
	Vec() { data = NULL, sz = 0, cap = 0; }
	explicit Vec(SZ_T size) { data = NULL, sz = 0, cap = 0; resize(size); }
	Vec(SZ_T size, const T& val) { data = NULL, sz = 0, cap = 0; resize(size, val); }
	~Vec() { clear(true); }

	// pointer to raw data
	operator T* (void) { return data; }
	T* d_ptr(void) { return data; }
	// indexing operators:
	const T& operator [] (SZ_T index) const { assert(index < sz); return data[index]; }
	T& operator [] (SZ_T index) { assert(index < sz); return data[index]; }

	// memory 
	void reallocMem(SZ_T min_cap);
	void expand(SZ_T size);
	void expand(SZ_T size, const T& val);
	void clear(bool dealloc = false);
	bool empty() { return sz == 0; }
	SZ_T size(void) const { return sz; }
	void resize(SZ_T n) {
		if (n == sz) return;
		if (n < sz) shrink(sz - n);
		else expand(n);
	}
	void resize(SZ_T n, const T& val) {
		if (n == sz) {
			for (SZ_T i = 0; i < sz; i++) data[i] = val;
			return;
		}
		if (n < sz) shrink(sz - n, val);
		else expand(n, val);
	}
	void shrink(SZ_T n) {
		assert(n <= sz);
		for (SZ_T i = 0; i < n; i++) data[--sz].~T();
	}
	void shrink(SZ_T n, const T& val) {
		assert(n <= sz);
		for (SZ_T i = 0; i < n; i++) data[--sz].~T();
		for (SZ_T i = 0; i < sz; i++) data[i] = val;
	}

	// new elements
	void push(void) { if (sz == cap) reallocMem(sz + 1); new (&data[sz++]) T(); }
	void push(const T& val) { if (sz == cap) reallocMem(sz + 1); new (&data[sz++]) T(val); }
	void pop(void) { assert(sz > 0); data[--sz].~T(); }
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
void Vec<T, SZ_T>::reallocMem(SZ_T min_cap) {
	if (cap >= min_cap) return;
	SZ_T add = max((min_cap - cap + 1) & ~1, ((cap >> 1) + 2) & ~1);
	const SZ_T size_max = std::numeric_limits<SZ_T>::max();
	if (((size_max <= std::numeric_limits<int>::max()) && (add > size_max - cap))
		|| (((data = (T*)::realloc(data, (cap += add) * sizeof(T))) == NULL) && errno == ENOMEM))
		throw MEMOUTEXCEPTION();
}

template<class T, class SZ_T>
void Vec<T, SZ_T>::expand(SZ_T size, const T& val) {
	if (sz >= size) return;
	reallocMem(size);
	for (SZ_T i = sz; i < size; i++) data[i] = val;
	sz = size;
}

template<class T, class SZ_T>
void Vec<T, SZ_T>::expand(SZ_T size) {
	if (sz >= size) return;
	reallocMem(size);
	for (SZ_T i = sz; i < size; i++) new (&data[i]) T();
	sz = size;
}

template<class T, class SZ_T>
void Vec<T, SZ_T>::clear(bool dealloc) {
	if (data != NULL) {
		for (SZ_T i = 0; i < sz; i++) data[i].~T();
		sz = 0;
		if (dealloc) {
			free(data);
			data = NULL;
			cap = 0;
		}
	}
}

#endif
