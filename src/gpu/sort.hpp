/***********************************************************************[sort.hpp]
Copyright(c) 2020, Muhammad Osama - Anton Wijs,
Copyright(c) 2022-present, Muhammad Osama.

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

#include "key.hpp"

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
	inline void insertion_sort(T* d, const SZ& sz, CMP cmp)
	{
		if (sz == 2 && cmp(d[1], d[0]))
			std::swap(d[1], d[0]);
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
			insertion_sort(d, sz, LESS<T>());
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
			insertion_sort(d, sz, cmp);
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
			insertion_sort(data, size, cmp);
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
			insertion_sort(data, size, cmp);
		else
			std::sort(data, data + size, cmp);
		assert(isSorted(d.data(), size, cmp));
	}

	//==================//
	//	radix sort		//
	//==================//
	#define RADIXBITS        8
	#define RADIXWIDTH       (1 << RADIXBITS)
    #define RADIXWIDTHBYTES  (RADIXWIDTH * 8)
	#define RADIXMASK        (RADIXWIDTH - 1)
	extern size_t  RADIXBUFFER[RADIXWIDTH]; // defined in 'pfinput.cpp'

	template<class T, class RANK>
	inline bool isSortedRadix(T* d, T* e, RANK rank)
	{
		for (T* p = d; p + 1 != e; p++)
			if (rank(p[0]) > rank(p[1])) return false;
		return true;
	}

	template<class T, class RANK>
	void radixSort(T* data, T* end, RANK rank)
	{
		assert(data <= end);
		size_t n = end - data;
		if (n < 2) return;
		const size_t dsize = (sizeof(rank(*data)) << 3);
		Vec<T, size_t> b(n);
		T* a = data, * c = a;
		for (size_t i = 0; i < dsize; i += RADIXBITS) {
			memset(RADIXBUFFER, 0, RADIXWIDTHBYTES);
			T* end = c + n;
			size_t upper = 0, lower = SIZE_MAX;
			for (T* p = c; p != end; p++) {
				auto s = rank(*p) >> i;
				auto m = s & RADIXMASK;
				lower &= s;
				upper |= s;
				RADIXBUFFER[m]++;
			}
			if (lower == upper) break;
			size_t pos = 0;
			for (size_t j = 0; j < RADIXWIDTH; j++) {
				size_t delta = RADIXBUFFER[j];
				RADIXBUFFER[j] = pos;
				pos += delta;
			}
			T* d = (c == a) ? b : a;
			for (T* p = c; p != end; p++) {
				auto s = rank(*p) >> i;
				auto m = s & RADIXMASK;
				d[RADIXBUFFER[m]++] = *p;
			}
			c = d;
		}
		if (c == b) {
			for (size_t i = 0; i < n; i++)
				a[i] = b[i];
		}
		assert(isSortedRadix(data, end, rank));
	}

	template<class T, class SZ>
	void rSort(T* d, const SZ& sz) {
		assert(d != NULL);
		if (!sz) return;
		if (sz <= RSORT_THR) {
			if (sz <= INSORT_THR) insertion_sort(d, sz, LESS<T>());
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
			if (sz <= INSORT_THR) insertion_sort(d, sz, cmp);
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
			if (size <= INSORT_THR) insertion_sort(d.data(), size, cmp);
			else std::sort(d.data(), d.end(), cmp);
			assert(isSorted(d.data(), size, cmp));
		}
		else radixSort(d.data(), d.end(), rank);
	}

	template<class T, class SZ, class CMP, class RANK>
	void rSort(Vec<T, SZ>& d, CMP cmp, RANK rank) {
		assert(d.data() != NULL);
		const auto size = d.size();
		if (!size) return;
		if (size <= RSORT_THR) {
			if (size <= INSORT_THR) insertion_sort(d.data(), size, cmp);
			else std::sort(d.data(), d.end(), cmp);
			assert(isSorted(d.data(), size, cmp));
		}
		else radixSort(d.data(), d.end(), rank);
	}

}

#endif