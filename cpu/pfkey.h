/***********************************************************************[pfkey.h]
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

#ifndef __SORT_KEY_
#define __SORT_KEY_

#include "pfvec.h"
#include "pfdefs.h"

namespace pFROST {

	//============================//
	//  Default Comparators       //
	//============================//
	template <class T>
	struct DEFAULT_RANK {
		T operator () (T val) { return val; }
	};

	struct PTR_RANK {
		size_t operator () (void* ptr) { return (size_t)ptr; }
	};

	template<class T>
	struct LESS {
		bool operator () (const T& x, const T& y) const {
			return x < y;
		}
	};

	template<class T>
	struct GREATER {
		bool operator () (const T& x, const T& y) const {
			return x > y;
		}
	};
	//============================//
	//  Custom Comparators        //
	//============================//
	struct LCV_CMP {
		uint32* scores;
		LCV_CMP(uint32* _scores) {
			assert(_scores != NULL);
			scores = _scores;
		}
		bool operator () (const uint32& a, const uint32& b) const {
			uint32 x = scores[a], y = scores[b];
			if (x < y) return true;
			if (x > y) return false;
			return a < b;
		}
	};
	struct MCV_CMP {
		uint32* scores;
		MCV_CMP(uint32* _scores) {
			assert(_scores != NULL);
			scores = _scores;
		}
		bool operator () (const uint32& a, const uint32& b) const {
			uint32 x = scores[a], y = scores[b];
			if (x > y) return true;
			if (x < y) return false;
			return a > b;
		}
	};
	struct HEAP_CMP {
		const Vec<double>& act;
		HEAP_CMP(const Vec<double>& _act) : act(_act) {}
		bool operator () (const uint32& a, const uint32& b) const {
			const double xact = act[a], yact = act[b];
			if (xact < yact) return true;
			if (xact > yact) return false;
			return a > b;
		}
	};
	struct ANALYZE_CMP {
		Vec<int64>& bumped;
		ANALYZE_CMP(Vec<int64>& _bumped) : bumped(_bumped) {}
		bool operator () (const uint32& a, const uint32& b) const {
			return bumped[a] < bumped[b];
		}
	};
	struct SUBSUME_OCCURS_CMP {
		uVec1D& hist;
		SUBSUME_OCCURS_CMP(uVec1D& _hist) : hist(_hist) {}
		bool operator () (const uint32& a, const uint32& b) const {
			uint32 xh = hist[a], yh = hist[b];
			if (xh < yh) return true;
			if (xh > yh) return false;
			return ABS(a) < ABS(b);
		}
	};
	struct KEY_CMP_ACTIVITY {
		double* acts;
		uint32* scores;
		KEY_CMP_ACTIVITY(double* _acts, uint32* _scores) {
			assert(_acts != NULL);
			assert(_scores != NULL);
			acts = _acts;
			scores = _scores;
		}
		bool operator()(const uint32& a, const uint32& b) const {
			double dx = acts[a], dy = acts[b];
			if (dx > dy) return true;
			if (dx < dy) return false;
			uint32 x = scores[a], y = scores[b];
			if (x > y) return true;
			if (x < y) return false;
			return a > b;
		}
	};
	struct KEY_CMP_BUMP {
		int64* bumps;
		KEY_CMP_BUMP(int64* _bumps) {
			assert(_bumps != NULL);
			bumps = _bumps;
		}
		bool operator()(const uint32& x, const uint32& y) const {
			return bumps[x] > bumps[y];
		}
	};
	//============================//
	//  Custom Rankers	          //
	//============================//
	struct ANALYZE_RANK {
		Vec<int64>& bumped;
		ANALYZE_RANK(Vec<int64>& _bumped) : bumped(_bumped) {}
		uint64 operator () (const uint32& a) const {
			return bumped[a];
		}
	};
	struct KEY_RANK_BUMP {
		int64* bumps;
		KEY_RANK_BUMP(int64* _bumps) {
			assert(_bumps != NULL);
			bumps = _bumps;
		}
		uint64 operator()(const uint32& x) const {
			uint64 b = bumps[x];
			return ~b;
		}
	};
	struct LCV_RANK {
		uint32* scores;
		LCV_RANK(uint32* _scores) {
			assert(_scores != NULL);
			scores = _scores;
		}
		uint32 operator () (const uint32& a) const {
			return scores[a];
		}
	};
	struct MCV_RANK {
		uint32* scores;
		MCV_RANK(uint32* _scores) {
			assert(_scores != NULL);
			scores = _scores;
		}
		uint32 operator () (const uint32& a) const {
			return ~scores[a];
		}
	};

}

#endif 