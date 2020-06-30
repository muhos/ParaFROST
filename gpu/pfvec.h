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

#include "pfralloc.h"
#include <cstdlib>
#include <limits>
#include <cassert>
#include <typeinfo>
#include <iostream>

namespace pFROST {

#ifdef __GNUC__
#define __forceinline __attribute__((always_inline))
#endif

	template<class T, class S = uint32>
	class Vec {
		T* _mem;
		S sz, cap, maxCap;
		bool check(const S& idx) const {
			if (idx >= sz) {
				std::cout << "Error - index is out of vector boundary (type: " << typeid(T).name() <<
					", idx: " << (long long)idx << ", sz:" << (long long)sz << std::endl;
				return false;
			}
			return true;
		}
	public:
		__forceinline			~Vec		() { clear(true); }
		__forceinline			Vec			() { maxCap = std::numeric_limits<S>::max(), _mem = NULL, sz = 0, cap = 0; }
		__forceinline explicit	Vec			(const S& size) {
			maxCap = std::numeric_limits<S>::max();
			_mem = NULL, sz = 0, cap = 0; resize(size);
		}
		__forceinline			Vec			(const S& size, const T& val) {
			maxCap = std::numeric_limits<S>::max();
			_mem = NULL, sz = 0, cap = 0; resize(size, val);
		}
		__forceinline Vec<T>&	operator=	(Vec<T>& rhs) { return *this; }
		__forceinline const T&	operator[]	(const S& index) const { assert(check(index)); return _mem[index]; }
		__forceinline T&		operator[]	(const S& index) { assert(check(index)); return _mem[index]; }
		__forceinline const T&	back		() const { assert(sz); return _mem[sz - 1]; }
		__forceinline T&		back		() { assert(sz); return _mem[sz - 1]; }
		__forceinline			operator T* () { return _mem; }
		__forceinline T*		data		() { return _mem; }
		__forceinline T*		end			() { return _mem + sz; }
		__forceinline bool		empty		() const { return sz == 0; }
		__forceinline S			size		() const { return sz; }
		__forceinline S			capacity	() const { return cap; }
		__forceinline void		pop			() { assert(sz > 0); _mem[--sz].~T(); }
		__forceinline void		insert		(const T& val) { assert(cap > sz);  _mem[sz++] = val; }
		__forceinline void		push		(const T& val) { if (sz == cap) reserve(sz + 1); new (_mem + sz) T(val); sz++; }
		__forceinline void		reserve		(const S& min_cap, const S& size) { reserve(min_cap); sz = size; }
		__forceinline void		init		(const T& val) {
			for (S i = 0; i < sz; i++) _mem[i] = val;
		}
		__forceinline void		resize		(const S& n) {
			if (n == sz) return;
			if (n < sz) shrink(sz - n);
			else expand(n);
		}
		__forceinline void		resize		(const S& n, const T& val) {
			if (n == sz) { init(val); return; }
			if (n < sz) shrink(sz - n, val);
			else expand(n, val);
		}
		__forceinline void		shrink		(const S& n) {
			assert(n <= sz);
			for (S i = 0; i < n; i++) _mem[--sz].~T();
		}
		__forceinline void		shrink		(const S& n, const T& val) {
			assert(n <= sz);
			for (S i = 0; i < n; i++) _mem[--sz].~T();
			init(val);
		}
		__forceinline void		copyFrom	(Vec<T, S>& copy) {
			resize(copy.size());
			std::memcpy(_mem, copy, sz * sizeof(T));
		}
		template<class SS = uint32>
		__forceinline void		copyFrom	(T* copy, const SS& copy_sz) {
			assert(copy_sz <= sz);
			std::memcpy(_mem, copy, sz * sizeof(T));
		}
		__forceinline void		expand		(const S& size) {
			if (sz >= size) return;
			reserve(size);
			for (S i = sz; i < size; i++) new (&_mem[i]) T();
			sz = size;
		}
		__forceinline void		expand		(const S& size, const T& val) {
			if (sz >= size) return;
			reserve(size);
			for (S i = sz; i < size; i++) _mem[i] = val;
			sz = size;
		}
		__forceinline void		reserve		(const S& min_cap) {
			if (cap >= min_cap) return;
			if (cap > (maxCap - cap)) cap = min_cap;
			else { cap <<= 1; if (cap < min_cap) cap = min_cap; }
			pfalloc(_mem, sizeof(T) * cap);
		}
		__forceinline void		shrinkCap	() {
			if (sz == 0) { clear(true); return; }
			pfshrinkAlloc(_mem, sizeof(T) * sz);
			cap = sz;
		}
		__forceinline void		clear		(const bool& _free = false) {
			if (_mem != NULL) {
				for (S i = 0; i < sz; i++) _mem[i].~T();
				sz = 0;
				if (_free) { std::free(_mem); _mem = NULL; cap = 0; }
			}
		}
	};
	// vector types
	typedef Vec<int> Vec1D;
	typedef Vec<uint32> uVec1D;
	typedef Vec<uint32, int> Lits_t;
}

#endif
