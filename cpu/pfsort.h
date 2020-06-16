/***********************************************************************[pfsort.h]
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

#ifndef __SORT_
#define __SORT_

#include "pfvec.h"
#include "pfdefs.h"

namespace pFROST {

#define INSORT_THR 20

	//============================//
	//  Template Comparators      //
	//============================//
	template<class T>
	struct LESS {
		bool operator () (T x, T y) {
			return x < y;
		}
	};

	template<class T>
	struct GREATER {
		bool operator () (T x, T y) {
			return x > y;
		}
	};

	template<class A, class S>
	struct KEY_CMP_ACTIVITY {
		A* acts;
		S* scores;
		KEY_CMP_ACTIVITY(A* _acts, S* _scores) {
			assert(_acts != NULL);
			assert(_scores != NULL);
			acts = _acts;
			scores = _scores;
		}
		bool operator()(uint32& x, uint32& y) {
			if (acts[x] != acts[y]) return acts[x] > acts[y];
			else return scores[x] != scores[y] ? scores[x] > scores[y] : x > y;
		}
	};

	template<class B>
	struct KEY_CMP_BUMP {
		B* bumps;
		KEY_CMP_BUMP(B* _bumps) {
			assert(_bumps != NULL);
			bumps = _bumps;
		}
		bool operator()(uint32& x, uint32& y) {
			return bumps[x] > bumps[y];
		}
	};

	struct LCV_CMP {
		uint32* scores;
		LCV_CMP(uint32* _scores) {
			assert(_scores != NULL);
			scores = _scores;
		}
		bool operator () (uint32& x, uint32& y) {
			return (scores[x] != scores[y]) ? scores[x] < scores[y] : x < y;
		}
	};

	struct MCV_CMP {
		uint32* scores;
		MCV_CMP(uint32* _scores) {
			assert(_scores != NULL);
			scores = _scores;
		}
		bool operator () (uint32& x, uint32& y) {
			return (scores[x] != scores[y]) ? scores[x] > scores[y] : x > y;
		}
	};

	//============================//
	//  Sorting Functions         //
	//============================//
	template<class T, class SZ, class CMP>
	inline bool isSorted(T* d, const SZ& sz, CMP cmp) {
		for (SZ i = 1; i < sz; i++)
			if (cmp(d[i], d[i - 1])) return false;
		return true;
	}

	template<class T, class SZ, class CMP>
	void insertion_sort(T* d, const SZ& sz, CMP cmp)
	{
		if (sz == 2 && cmp(d[1], d[0]))
			swap(d[1], d[0]);
		else if (sz > 2) {
			SZ i, j;
			for (i = 1; i < sz; i++) {
				register T tmp = d[i];
				for (j = i; j > 0 && cmp(tmp, d[j - 1]); j--)
					d[j] = d[j - 1];
				d[j] = tmp;
			}
		}
	}

	template<class T, class SZ>
	void Sort(T* d, const SZ& sz) {
		assert(d != NULL);
		assert(sz > 0);
		if (sz <= INSORT_THR)
			insertion_sort(d, sz, LESS<T>());
		else
			std::sort(d, d + sz, LESS<T>());
		assert(isSorted(d, sz, LESS<T>()));
	}

	template<class T, class SZ, class CMP>
	void Sort(T* d, const SZ& sz, CMP cmp) {
		assert(d != NULL);
		assert(sz > 0);
		if (sz <= INSORT_THR)
			insertion_sort(d, sz, cmp);
		else
			std::sort(d, d + sz, cmp);
		assert(isSorted(d, sz, cmp));
	}

	template<class T, class CMP>
	void Sort(Vec<T>& d, CMP cmp) {
		assert(d.data() != NULL);
		assert(d.size() > 0);
		if (d.size() <= INSORT_THR)
			insertion_sort(d.data(), d.size(), cmp);
		else
			std::sort(d.data(), d.data() + d.size(), cmp);
		assert(isSorted(d.data(), d.size(), cmp));
	}

	template<class T, class SZ, class CMP>
	void Sort(Vec<T, SZ>& d, CMP cmp) {
		assert(d.data() != NULL);
		assert(d.size() > 0);
		if (d.size() <= INSORT_THR)
			insertion_sort(d.data(), d.size(), cmp);
		else
			std::sort(d.data(), d.data() + d.size(), cmp);
		assert(isSorted(d.data(), d.size(), cmp));
	}

}

#endif // __SORT_