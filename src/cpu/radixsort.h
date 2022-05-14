/***********************************************************************[radixsort.hpp]
Copyright(c) 2022, Muhammad Osama - Anton Wijs,
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

#ifndef __RADIXSORT_
#define __RADIXSORT_

#include "malloc.h"

namespace pFROST {

	constexpr size_t RADIXBITS			= 8;
	constexpr size_t RADIXWIDTH			= (1 << RADIXBITS);
	constexpr size_t RADIXBYTES			= (RADIXWIDTH * sizeof(size_t));
	constexpr size_t RADIXMASK			= (RADIXWIDTH - 1);

	static size_t RADIXBUFFER[RADIXWIDTH];

	template<class T, class RANK>
	inline bool isSortedRadix(T* d, T* e, RANK rank)
	{
		for (T* p = d; p + 1 != e; ++p)
			if (rank(p[0]) > rank(p[1])) 
				return false;
		return true;
	}

	template<class T, class RANK>
	void radixSort(T* data, T* end, RANK rank)
	{
		assert(data <= end);
		size_t n = end - data;

		if (n < 2) return;

		const size_t dsize = (sizeof(rank(*data)) * RADIXBITS);

		T* b = pfmalloc<T>(n);

		T* a = data, * c = a;
		for (size_t i = 0; i < dsize; i += RADIXBITS) {
			std::memset(RADIXBUFFER, 0, RADIXBYTES);
			T* end = c + n;
			size_t upper = 0, lower = SIZE_MAX;
			for (T* p = c; p != end; ++p) {
				auto s = rank(*p) >> i;
				auto m = s & RADIXMASK;
				lower &= s;
				upper |= s;
				RADIXBUFFER[m]++;
			}
			if (lower == upper) break;
			size_t pos = 0;
			for (size_t j = 0; j < RADIXWIDTH; ++j) {
				size_t delta = RADIXBUFFER[j];
				RADIXBUFFER[j] = pos;
				pos += delta;
			}
			T* d = (c == a) ? b : a;
			for (T* p = c; p != end; ++p) {
				auto s = rank(*p) >> i;
				auto m = s & RADIXMASK;
				d[RADIXBUFFER[m]++] = *p;
			}
			c = d;
		}

		if (c == b)
			std::memcpy(a, b, n * sizeof(T));

		std::free(b);

		assert(isSortedRadix(data, end, rank));
	}

}

#endif