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
	BQUEUE(void){
		h = 0, t = 0, _bound = 0, _sz = 0, _sum = 0ULL;
	}
	~BQUEUE() {}

	void alloc(int cap) {
		_bound = cap;
		data.resize(_bound, 0);
		reset();
	}

	void push(T x) {
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

	void pop() { 
		_sum -= data[t++]; 
		if (t == _bound) t = 0;
		_sz--; 
	}

	int bound() { return _bound; }

	bool isFull() { return _sz == _bound; }

	T tail() { assert(_sz > 0); return data[t]; }

	int size() { return _sz; }

	void reset() { h = 0, t = 0, _sz = 0, _sum = 0ULL; }

	void clear(bool dealloc = false) { 
		data.clear(dealloc); 
		reset();
	}
	
	// stats calculation
	uint64 sum() { return _sum; }

	uint32 average() { return (uint32)(_sum / (uint64)_sz); }
};

#endif // __BQ_