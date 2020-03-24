/***********************************************************************
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
************************************************************************/
#ifndef __SORT_
#define __SORT_

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

template<class A, class H>
struct KEY_CMP_ACTIVITY{
	A* acts;
	H* scores;
	KEY_CMP_ACTIVITY(A* acts, H* scores) {
		assert(acts != NULL);
		assert(scores != NULL);
		this->acts = acts;
		this->scores = scores;
	}
	bool operator()(int x, int y) {
		if (acts[x] != acts[y]) return (acts[x] > acts[y]);
		else return (scores[x] != scores[y]) ? (scores[x] > scores[y]) : (x > y);
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
inline bool isSorted(T* data, const SZ& sz, CMP cmp) {
	for (int i = 1; i < sz; i++)
		if (cmp(data[i], data[i - 1])) return false;
	return true;
}

template<class T, class SZ, class CMP>
void insertion_sort(T* data, const SZ& sz, CMP cmp)
{
	if (sz == 2 && cmp(*(data + 1), *data))
		Swap(data + 1, data);
	else if (sz > 2) {
		int i, j;
		for (i = 1; i < sz; i++) {
			register T tmp = data[i];
			for (j = i; j > 0 && cmp(tmp, data[j - 1]); j--)
				data[j] = data[j - 1];
			data[j] = tmp;
		}
	}
}

template<class T, class SZ>
void Sort(T* data, const SZ& sz) {
	assert(data != NULL);
	assert(sz > 0);
	if (sz <= INSORT_THR)
		insertion_sort(data, sz, LESS<T>());
	else
		std::sort(data, data + sz, LESS<T>());
	assert(isSorted(data, sz, LESS<T>()));
}

template<class T, class SZ, class CMP>
void Sort(T* data, const SZ& sz, CMP cmp) {
	assert(data != NULL);
	assert(sz > 0);
	if (sz <= INSORT_THR)
		insertion_sort(data, sz, cmp);
	else 
		std::sort(data, data + sz, cmp);
	assert(isSorted(data, sz, cmp));
}

template<class T, class CMP>
void Sort(Vec<T>& data, CMP cmp) {
	assert(data.d_ptr() != NULL);
	assert(data.size() > 0);
	if (data.size() <= INSORT_THR)
		insertion_sort(data.d_ptr(), data.size(), cmp);
	else
		std::sort(data.d_ptr(), data.d_ptr() + data.size(), cmp);
	assert(isSorted(data.d_ptr(), data.size(), cmp));
}

#endif // __SORT_