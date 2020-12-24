/***********************************************************************[pfqueue.h]
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

#ifndef __QUEUE_
#define __QUEUE_

#include "pfvec.h"
#include "pfdefs.h"

namespace pFROST {

	struct LINK {
		uint32 prev, next;
		LINK() : prev(0), next(0) {}
	};
	typedef Vec<LINK> Links;
	class QUEUE {
		Links links;
		int64 _bumped;
		uint32 _first, _last, _free;

		__forceinline void		inQue		(const uint32& v) {
			LINK& link = links[v];
			link.prev = _last;
			if (_last) links[_last].next = v;
			else _first = v;
			_last = v;
			link.next = 0;
		}
		__forceinline void		outQue		(const uint32& v) {
			LINK& link = links[v];
			if (link.prev) links[link.prev].next = link.next;
			else _first = link.next;
			if (link.next) links[link.next].prev = link.prev;
			else _last = link.prev;
		}

	public:
								QUEUE		() : _first(0), _last(0), _free(0), _bumped(0) {}
		__forceinline void		init		(const uint32& v) {
			links.expand(v + 1);
			LINK& link = links[v];
			link.next = 0;
			if (_last) {
				assert(!links[_last].next);
				links[_last].next = v;
			}
			else {
				assert(!_first);
				_first = v;
			}
			link.prev = _last;
			_last = v;
		}
		__forceinline void		map			(uint32* mapped, const uint32& firstDL0) {
			uint32 q, _next, _prev = 0, _mPrev = 0;
			for (q = _first; q; q = _next) {
				_next = next(q);
				if (q == firstDL0) continue;
				uint32 mVar = mapped[q];
				if (mVar) {
					if (_prev) links[_prev].next = mVar;
					else _first = mVar;
					links[q].prev = _mPrev;
					_mPrev = mVar;
					_prev = q;
				}
			}
			if (_prev) links[_prev].next = 0;
			else _first = 0;
			_free = _last = _mPrev;

		}
		__forceinline void		update		(const uint32& v, const int64& bump) { _free = v, _bumped = bump; PFLOG2(4, "  - Queue free updated to (v: %d, bump: %lld)", _free, _bumped); }
		__forceinline void		toFront		(const uint32& v) { outQue(v), inQue(v); }
		__forceinline uint32	previous	(const uint32& v) { return links[v].prev; }
		__forceinline uint32	next		(const uint32& v) { return links[v].next; }
		__forceinline Links&	data		() { return links; }
		__forceinline uint32	free		() { return _free; }
		__forceinline uint32	first		() { return _first; }
		__forceinline uint32	last		() { return _last; }
		__forceinline int64		bumped		() { return _bumped; }
					  void		print		() {
			PFLOG1(" Queue (first: %d, last: %d, free: %d, bumped: %lld):", _first, _last, _free, _bumped);
			for (uint32 i = 0; i < links.size(); i++) 
				PFLOG1(" Q(%d)->(p: %d, n: %d)", i, links[i].prev, links[i].next);
		}
	};

	template<class T>
	class BQUEUE {
		Vec<T>  _mem;
		uint64 _sum;
		int h, t, _bound, _sz;

	public:
		BQUEUE() { h = 0, t = 0, _bound = 0, _sz = 0, _sum = 0ULL; }
		~BQUEUE() { _mem.clear(true); }
		inline void		clear(bool _free = false) { _mem.clear(_free), reset(); }
		inline void		alloc(int cap) {
			_bound = cap;
			_mem.resize(_bound, 0);
			reset();
		}
		inline void		push(T x) {
			if (_sz == _bound) {
				assert(t == h);
				_sum -= _mem[t++];
				if (t == _bound) t = 0;
			}
			else _sz++;
			_sum += x;
			_mem[h++] = x;
			if (h == _bound) { h = 0, t = 0; }
		}
		inline void		pop() {
			_sum -= _mem[t++];
			if (t == _bound) t = 0;
			_sz--;
		}
		inline int		bound() const { return _bound; }
		inline bool		full() const { return _sz == _bound; }
		inline T		head() const { assert(_sz > 0); return *_mem; }
		inline T		tail() const { assert(t < _sz); return _mem[t]; }
		inline int		size() const { return _sz; }
		inline void		reset() { h = 0, t = 0, _sz = 0, _sum = 0ULL; }
		inline uint64	sum() { return _sum; }
		inline uint32	average() { return (uint32)(_sum / (uint64)_sz); }
	};

}

#endif 