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

template<class T>
struct KEY_CMP_LESS1 {
	T* key;
	KEY_CMP_LESS1(T* key) {
		assert(key != NULL);
		this->key = key;
	}
	bool operator()(int x, int y) {
		return key[x] < key[y];
	}
};

template<class T>
struct KEY_CMP_GREATER1 {
	T* key;
	KEY_CMP_GREATER1(T* key) {
		assert(key != NULL);
		this->key = key;
	}
	bool operator()(int x, int y) {
		return key[x] > key[y];
	}
};

template<class A>
struct KEY_CMP_ACTIVITY{
	A* acts;
	KEY_CMP_ACTIVITY(A* acts) {
		assert(acts != NULL);
		this->acts = acts;
	}
	bool operator()(SCORE& x, SCORE& y) {
		if (acts[x.v] != acts[y.v]) return (acts[x.v] > acts[y.v]);
		else return (x.sc != y.sc) ? (x.sc > y.sc) : (x.v > y.v);
	}
};

struct VAR_CMP {
	bool operator () (SCORE& a, SCORE& b) {
		return a.sc != b.sc ? a.sc < b.sc : a.v < b.v;
	}
};

//============================//
//  Sorting Functions         //
//============================//
template<class T>
inline void Swap(T* a, T* b)
{
	T t = *a;
	*a = *b;
	*b = t;
}

template<class T, class SZ, class CMP>
inline bool isSorted(T* d, const SZ& sz, CMP cmp) {
	for (int i = 1; i < sz; i++)
		if (cmp(d[i], d[i - 1])) return false;
	return true;
}

template<class T, class SZ, class CMP>
void insertion_sort(T* d, const SZ& sz, CMP cmp)
{
	if (sz == 2 && cmp(*(d + 1), *d))
		Swap(d + 1, d);
	else if (sz > 2) {
		int i, j;
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

#endif // __SORT_