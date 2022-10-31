/***********************************************************************[heap.h]
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

#ifndef __HEAP_
#define __HEAP_

#include "vector.h"
#include "definitions.h"

namespace ParaFROST {

#define ILLEGAL_POS UINT32_MAX

	template <class CMP>
	class HEAP {
		uVec1D	heap;
		uVec1D	pos;
		CMP		cmp;
		__forceinline uint32	L			(const uint32& x) { return (pos[x] << 1) + 1; }
		__forceinline uint32	R			(const uint32& x) { return (pos[x] + 1) << 1; }
		__forceinline uint32	left		(const uint32& x) { return heap[L(x)]; }
		__forceinline uint32	right		(const uint32& x) { return heap[R(x)]; }
		__forceinline uint32	parent		(const uint32& x) { return heap[(pos[x] - 1) >> 1]; }
		__forceinline bool		hasLeft		(const uint32& x) { return L(x) < heap.size(); }
		__forceinline bool		hasRight	(const uint32& x) { return R(x) < heap.size(); }
		__forceinline bool		hasParent	(const uint32& x) { return pos[x] > 0; }
		__forceinline void		exch		(const uint32& x, const uint32& y) {
			uint32& i = pos[x], & j = pos[y];
			std::swap(heap[i], heap[j]);
			std::swap(i, j);
		}
		__forceinline void		bubbleUp	(uint32 x) {
			assert(x < ILLEGAL_POS);
			uint32 p;
			while (hasParent(x) && cmp((p = parent(x)), x)) exch(p, x);
		}
		__forceinline void		bubbleDown	(uint32 x) {
			assert(x < ILLEGAL_POS);
			while (hasLeft(x)) {
				uint32 child = left(x);
				if (hasRight(x)) {
					uint32 r_x = right(x);
					if (cmp(child, r_x)) child = r_x;
				}
				if (!cmp(x, child)) break;
				exch(x, child);
			}
		}
	public:
								HEAP		(const CMP& _cmp) : cmp(_cmp) {}
								HEAP		() {}
								~HEAP		() { destroy(); }
		__forceinline uint32*	data		() { return heap; }
		__forceinline uint32	top			() { assert(!empty()); return *heap; }
		__forceinline uint32	size		() const { return heap.size(); }
		__forceinline bool		empty		() const { return !heap.size(); }
		__forceinline uint32	pop			() {
			uint32 top_x = heap[0], last = heap.back();
			if (heap.size() > 1) exch(top_x, last);
			pos[top_x] = ILLEGAL_POS;
			heap.pop();
			if (heap.size() > 1) bubbleDown(last);
			return top_x;
		}
		__forceinline void		clear		() {
			for (uint32 i = 0; i < heap.size(); i++) pos[heap[i]] = ILLEGAL_POS;
			heap.clear();
		}
		__forceinline void		destroy		() {
			heap.clear(true);
			pos.clear(true);
		}
		__forceinline uint32&	operator[]	(const uint32& i) { assert(i < heap.size()); return heap[i]; }
		__forceinline bool		has			(const uint32& x) const { return pos.size() > x && pos[x] != ILLEGAL_POS; }
		__forceinline void		bump		(const uint32& x) { if (has(x)) update(x); }
		__forceinline void		update		(uint32 x) {
			assert(has(x));
			bubbleUp(x);
			bubbleDown(x);
		}
		__forceinline void		insert		(const uint32& x) {
			pos.expand(x + 1, ILLEGAL_POS);
			assert(!has(x));
			pos[x] = heap.size();
			heap.push(x);
			bubbleUp(x);
			bubbleDown(x);
		}
		__forceinline void		rebuild		(uVec1D& vars) {
			destroy();
			for (uint32 i = 0; i < vars.size(); i++)
				insert(vars[i]);
		}
		__forceinline void		init		(const CMP& _cmp) { cmp = _cmp; }
		
	};

}

#endif