/***********************************************************************[sort.h]
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

#include "key.h"
#include "pdqsort.h"
#include "wolfsort.h"
#include "radixsort.h"

namespace ParaFROST {

#define INSORT_THR 20
#define RSORT_THR 800

	//============================//
	//  Sorting Functions         //
	//============================//
	template<class T, class SZ, class CMP>
	inline bool isSorted(T* d, const SZ& sz, CMP cmp) 
	{
		for (SZ i = 1; i < sz; i++)
			if (cmp(d[i], d[i - 1])) return false;
		return true;
	}

	template<class T, class SZ, class CMP>
	inline void insertionSort(T* d, const SZ& sz, CMP cmp)
	{
		if (sz == 2 && cmp(d[1], d[0]))
			swap(d[1], d[0]);
		else if (sz > 2) {
			SZ i, j;
			for (i = 1; i < sz; i++) {
				T tmp = d[i];
				for (j = i; j > 0 && cmp(tmp, d[j - 1]); j--)
					d[j] = d[j - 1];
				d[j] = tmp;
			}
		}
	}

	template<class T, class SZ>
	void Sort(T* d, const SZ& sz) 
	{
		assert(d != NULL);
		if (!sz) return;
		if (sz <= INSORT_THR)
			insertionSort(d, sz, LESS<T>());
		else
			std::sort(d, d + sz, LESS<T>());
		assert(isSorted(d, sz, LESS<T>()));
	}

	template<class T, class SZ, class CMP>
	void Sort(T* d, const SZ& sz, CMP cmp)
	{
		assert(d != NULL);
		if (!sz) return;
		if (sz <= INSORT_THR)
			insertionSort(d, sz, cmp);
		else
			std::sort(d, d + sz, cmp);
		assert(isSorted(d, sz, cmp));
	}

	template<class T, class CMP>
	void Sort(Vec<T>& d, CMP cmp) 
	{
		const auto size = d.size();
		T* data = d.data();
		assert(data != NULL);
		if (!size) return;
		if (size <= INSORT_THR)
			insertionSort(data, size, cmp);
		else
			std::sort(data, data + size, cmp);
		assert(isSorted(d.data(), size, cmp));
	}

	template<class T, class SZ, class CMP>
	void Sort(Vec<T, SZ>& d, CMP cmp)
	{
		const auto size = d.size();
		T* data = d.data();
		assert(data != NULL);
		if (!size) return;
		if (size <= INSORT_THR)
			insertionSort(data, size, cmp);
		else
			std::sort(data, data + size, cmp);
		assert(isSorted(d.data(), size, cmp));
	}

	template<class T, class SZ>
	void rSort(T* d, const SZ& sz) {
		assert(d != NULL);
		if (!sz) return;
		if (sz <= RSORT_THR) {
			if (sz <= INSORT_THR) insertionSort(d, sz, LESS<T>());
			else std::sort(d, d + sz, LESS<T>());
			assert(isSorted(d, sz, LESS<T>()));
		}
		else radixSort(d, d + sz, DEFAULT_RANK<T>());
	}

	template<class T, class SZ, class CMP, class RANK>
	void rSort(T* d, const SZ& sz, CMP cmp, RANK rank) {
		assert(d != NULL);
		if (!sz) return;
		if (sz <= RSORT_THR) {
			if (sz <= INSORT_THR) insertionSort(d, sz, cmp);
			else std::sort(d, d + sz, cmp);
			assert(isSorted(d, sz, cmp));
		}
		else radixSort(d, d + sz, rank);
	}

	template<class T, class CMP, class RANK>
	void rSort(Vec<T>& d, CMP cmp, RANK rank) {
		assert(d.data() != NULL);
		const auto size = d.size();
		if (!size) return;
		if (size <= RSORT_THR) {
			if (size <= INSORT_THR) insertionSort(d.data(), size, cmp);
			else std::sort(d.data(), d.end(), cmp);
			assert(isSorted(d.data(), size, cmp));
		}
		else radixSort(d.data(), d.end(), rank);
	}

	template<class T, class CMP, class RANK>
	void qrSort(Vec<T>& d, CMP cmp, RANK rank) {
		assert(d.data() != NULL);
		const auto size = d.size();
		if (!size) return;
		if (size <= RSORT_THR) {
			if (size <= INSORT_THR) insertionSort(d.data(), size, cmp);
			else pdqsort(d.data(), d.end(), cmp);
			assert(isSorted(d.data(), size, cmp));
		}
		else radixSort(d.data(), d.end(), rank);
	}

}

#endif // __SORT_