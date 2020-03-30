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

#include "pfdefs.h"

template<class T>
class BQUEUE {
	Vec<T>  data;
	uint64 _sum;
	int h, t, _bound, _sz;

public:
	BQUEUE(){ h = 0, t = 0, _bound = 0, _sz = 0, _sum = 0ULL; }
	~BQUEUE() { h = 0, t = 0, _bound = 0, _sz = 0, _sum = 0ULL; }

	inline void alloc(int cap) {
		_bound = cap;
		data.resize(_bound, 0);
		reset();
	}

	inline void push(T x) {
		if (_sz == _bound) {
			assert(t == h); 
			_sum -= data[t++];
			if (t == _bound) t = 0;
		}
		else _sz++;
		_sum += x;
		data[h++] = x;
		if (h == _bound) { h = 0, t = 0; }
	}

	inline void pop() {
		_sum -= data[t++]; 
		if (t == _bound) t = 0;
		_sz--; 
	}

	inline int bound() { return _bound; }
	inline bool full() { return _sz == _bound; }
	inline T tail() { assert(_sz > 0); return data[t]; }
	inline int size() { return _sz; }
	inline void reset() { h = 0, t = 0, _sz = 0, _sum = 0ULL; }
	inline void clear(bool dealloc = false) { data.clear(dealloc), reset(); }
	// stats calculation
	inline uint64 sum() { return _sum; }
	inline uint32 average() { return (uint32)(_sum / (uint64)_sz); }
};

#endif // __BQ_