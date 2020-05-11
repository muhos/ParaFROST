/***********************************************************************[pfvec.h]
Copyright(c) 2020, Muhammad Osama - Anton Wijs,
Technische Universiteit Eindhoven (TU/e).

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
**********************************************************************************/

#ifndef __VEC_
#define __VEC_

#include <cstdlib>
#include <limits>
#include <cassert>
#include <typeinfo>
#include <iostream>

class MEMOUTEXCEPTION {};

template<class T, class S = int>
class Vec {
	T* _mem;
	S sz, cap;
	bool check(const S& idx) const {
		if (idx >= sz) {
			std::cout << "ERROR - index is out of vector boundary:" << std::endl;
			std::cout << " - type  : " << typeid(T).name() << std::endl <<
				         " - index : " << (long long)idx << std::endl <<
				         " - size  :" << (long long)sz << std::endl;
			return false;
		}
		return true;
	}
public:
	inline Vec() { _mem = NULL, sz = 0, cap = 0; }
	inline explicit Vec(const S& size) { _mem = NULL, sz = 0, cap = 0; resize(size); }
	inline Vec(const S&  size, const T& val) { _mem = NULL, sz = 0, cap = 0; resize(size, val); }
	inline ~Vec() { clear(true); }
	inline Vec<T>& operator=(Vec<T>& rhs) { assert(0); }
	inline const T& operator [] (const S& index) const { assert(check(index)); return _mem[index]; }
	inline T& operator [] (const S& index) { assert(check(index)); return _mem[index]; }
	inline const T& back(void) const { assert(sz); return _mem[sz - 1]; }
	inline operator T* (void) { return _mem; }
	inline T* data(void) { return _mem; }
	inline bool empty(void) const { return sz == 0; }
	inline S size(void) const { return sz; }
	inline S capacity(void) const { return cap; }
	inline void push(const T& val) { if (sz == cap) reserve(sz + 1); new (_mem + sz) T(val); sz++; }
	inline void pop(void) { assert(sz > 0); _mem[--sz].~T(); }
	inline void reserve(const S& min_cap, const S& size) { reserve(min_cap); sz = size; }
	inline void resize(const S& n) {
		if (n == sz) return;
		if (n < sz) shrink(sz - n);
		else expand(n);
	}
	inline void resize(const S& n, const T& val) {
		if (n == sz) {
			for (S i = 0; i < sz; i++) _mem[i] = val;
			return;
		}
		if (n < sz) shrink(sz - n, val);
		else expand(n, val);
	}
	inline void shrink(const S& n) {
		assert(n <= sz);
		for (S i = 0; i < n; i++) _mem[--sz].~T();
	}
	inline void shrink(const S& n, const T& val) {
		assert(n <= sz);
		for (S i = 0; i < n; i++) _mem[--sz].~T();
		for (S i = 0; i < sz; i++) _mem[i] = val;
	}
	inline void appendFrom(T* copy, const S& copy_sz, const bool& _new = 0) {
		if (_new) for (S i = 0; i < copy_sz; i++) push(copy[i]);
		else {
			assert(copy_sz <= cap);
			assert(sz + copy_sz <= cap);
			for (S i = 0; i < copy_sz; i++) _mem[sz++] = copy[i];
		}
	}
	inline void copyFrom(Vec<T>& copy) const {
		assert(copy.size() <= sz);
		for (S i = 0; i < sz; i++) _mem[i] = copy[i];
	}
	template<class SS> inline void copyFrom(T* copy, const SS& copy_sz) const {
		assert(copy_sz <= sz);
		for (S i = 0; i < sz; i++) _mem[i] = copy[i];
	}
	inline void expand(const S& size) {
		if (sz >= size) return;
		reserve(size);
		for (S i = sz; i < size; i++) new (&_mem[i]) T();
		sz = size;
	}
	inline void expand(const S& size, const T& val) {
		if (sz >= size) return;
		reserve(size);
		for (S i = sz; i < size; i++) _mem[i] = val;
		sz = size;
	}
	inline void reserve(const S& min_cap) {
		if (cap >= min_cap) return;
		if (cap > (std::numeric_limits<S>::max() - cap)) cap = min_cap;
		else { cap <<= 1; if (cap < min_cap) cap = min_cap; }
		T* __mem = (T*)std::realloc(_mem, cap * sizeof(T));
		if (__mem == NULL) throw MEMOUTEXCEPTION();
		_mem = __mem;
	}
	inline void shrinkCap(void) {
		if (sz == 0) { clear(true); return; }
		T* __mem = NULL;
		__mem = (T*)std::realloc(__mem, cap * sizeof(T));
		if (__mem == NULL) throw MEMOUTEXCEPTION();
		for (S i = 0; i < sz; i++) __mem[i] = _mem[i];
		std::free(_mem);
		_mem = __mem;
		cap = sz;
	}
	inline void clear(const bool& _free = false) {
		if (_mem != NULL) {
			for (S i = 0; i < sz; i++) _mem[i].~T();
			sz = 0;
			if (_free) { std::free(_mem); _mem = NULL; cap = 0; }
		}
	}
};
// vector types
typedef Vec<int> Vec1D;
typedef Vec<unsigned int> uVec1D;

#endif
