/***********************************************************************[vector.cuh]
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

#ifndef __CU_VECTOR_
#define __CU_VECTOR_

#include "definitions.cuh"
#include "const.h"

namespace pFROST {

	template<typename T>
	class cuVec {
		T* _mem;
		uint32 sz, cap;

		_PFROST_D_ bool checkAtomicBound(const uint32 idx, const uint32 cap)
		{
			if (idx >= cap) {
				printf("c ERROR: vector atomic returned index (%d) exceeding allocated capacity (%d)\n", idx, cap);
				return false;
			}
			return true;
		}

	public:
		_PFROST_H_D_			cuVec		() : _mem(NULL), sz(0), cap(0) {}
		_PFROST_H_D_			~cuVec		() { clear(true); }
		_PFROST_H_D_ void		alloc		(T* head) { _mem = head; }
		_PFROST_H_D_ void		alloc		(T* head, const uint32& cap) { _mem = head, this->cap = cap; }
		_PFROST_H_D_ void		alloc		(const uint32& cap) { _mem = (T*)(this + 1), this->cap = cap; }
		_PFROST_IN_D_ T*        jump		(const uint32&);
		_PFROST_IN_D_ void		insert		(const T&);
		_PFROST_IN_D_ void		push		(const T&);
		_PFROST_H_D_ void		_pop		() { sz--; }
		_PFROST_H_D_ void		_shrink		(const uint32& n) { sz -= n; }
		_PFROST_H_D_ void		_push		(const T& val) { _mem[sz++] = val; }
		_PFROST_H_D_ cuVec<T>&  operator=	(cuVec<T>& rhs) { return *this; }
		_PFROST_H_D_ const T&   operator [] (const uint32& idx) const { assert(idx < sz); return _mem[idx]; }
		_PFROST_H_D_ T&         operator [] (const uint32& idx)		  { assert(idx < sz); return _mem[idx]; }
		_PFROST_H_D_ T&         at			(const uint32& idx)		  { return _mem[idx]; }
		_PFROST_H_D_ const T&   at			(const uint32& idx) const { return _mem[idx]; }
		_PFROST_H_D_			operator T* () { return _mem; }
		_PFROST_H_D_ T*         data		() { return _mem; }
		_PFROST_H_D_ T*         end			() { return _mem + sz; }
		_PFROST_H_D_ T&         back		() { assert(sz); return _mem[sz - 1]; }
		_PFROST_H_D_ bool		empty		() const { return !sz; }
		_PFROST_H_D_ uint32		size		() const { return sz; }
		_PFROST_H_D_ uint32		capacity	() const { return cap; }
		_PFROST_H_D_ void		resize		(const uint32& n) { assert(n <= cap); sz = n; }
		_PFROST_H_D_ void		shareTo		(T* dest) {
			assert(sz && _mem);
			T* s = _mem, * e = s + sz;
			while (s != e) *dest++ = *s++;
		}
		_PFROST_H_D_ void		copyFrom	(T* src) {
			assert(sz);
			T* d = _mem, * e = d + sz;
			while (d != e) *d++ = *src++;
		}
		_PFROST_H_D_ void		copyFrom	(T* src, const uint32& n) {
			assert(n <= sz);
			T* d = _mem, * e = d + n;
			while (d != e) *d++ = *src++;
		}
		_PFROST_H_D_ void		clear		(const bool& _free = false) {
			if (_free) _mem = NULL, cap = 0;
			sz = 0;
		}
		_PFROST_H_D_ void		print		(const bool& litType = false, const bool& clauseType = false) {
			printf("->(size = %d)[", sz);
			if (litType) {
				for (uint32 i = 0; i < sz; i++) {
					printf("%2d  ", SIGN(_mem[i]) ? -int(ABS(_mem[i])) : int(ABS(_mem[i])));
					if (clauseType && _mem[i] < 2) printf("\nc \t\t");
					else if (!clauseType && i && i < sz - 1 && i % 10 == 0) printf("\nc \t\t");
				}
			}
			else {
				for (uint32 i = 0; i < sz; i++) {
					printf("%2lld  ", uint64(_mem[i]));
					if (i && i < sz - 1 && i % 10 == 0) printf("\nc \t\t");
				}
			}
			printf("]\n");
		}
	};
	
}

#endif