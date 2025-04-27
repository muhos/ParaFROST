/***********************************************************************[primitives.cuh]
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

#ifndef __GPU_PRIMITIVE_
#define __GPU_PRIMITIVE_

#include "key.cuh"
#include "cnf.cuh"
#include "grid.cuh"
#include "vector.cuh"
#include "sclause.cuh"
#include "options.cuh"
#include "constants.cuh"
#include "histogram.cuh"
#include "definitions.hpp"

namespace ParaFROST {

	#define ORIGINIZELIT(ORGLIT, LIT)							\
		const uint32* VORGPTR = DC_PTRS->d_vorg;				\
		assert(VORGPTR);										\
		assert(LIT > 1 && LIT < NOVAR);							\
		uint32 ORGLIT = V2DEC(VORGPTR[ABS(LIT)], SIGN(LIT));	\
		assert(ORGLIT > 1 && ORGLIT < NOVAR);					\

	template<class T, class CMP>
	_PFROST_H_D_ bool devIsSorted(T* d, const int& sz, CMP cmp) {
		for (int i = 1; i < sz; i++)
			if (cmp(d[i], d[i - 1])) return false;
		return true;
	}

	template<typename T, typename CMP>
	_PFROST_IN_D_ void cswap(T& x, T& y, CMP cmp)
	{
		const bool which = cmp(x, y);
		const T ta = which ? x : y;
		const T tb = which ? y : x;
		x = ta, y = tb;
	}

	template<typename T, typename CMP>
	_PFROST_IN_D_ void sort3(T& x, T& y, T& z, CMP cmp)
	{
		cswap(y, z, cmp);
		cswap(x, z, cmp);
		cswap(x, y, cmp);
	}

	template<typename T, typename CMP>
	_PFROST_IN_D_ void sort4(T* d, CMP cmp)
	{
		cswap(d[0], d[1], cmp);
		cswap(d[2], d[3], cmp);
		cswap(d[0], d[2], cmp);
		cswap(d[1], d[3], cmp);
		cswap(d[1], d[2], cmp);
	}

	template<typename T, typename CMP>
	_PFROST_IN_D_ void sort5(T* d, CMP cmp)
	{
		cswap(d[0], d[1], cmp);
		cswap(d[3], d[4], cmp);
		cswap(d[2], d[4], cmp);
		cswap(d[2], d[3], cmp);
		cswap(d[0], d[3], cmp);
		cswap(d[0], d[2], cmp);
		cswap(d[1], d[4], cmp);
		cswap(d[1], d[3], cmp);
		cswap(d[1], d[2], cmp);
	}

	template<typename T, typename CMP>
	_PFROST_IN_D_ void sort6(T* d, CMP cmp)
	{
		cswap(d[1], d[2], cmp);
		cswap(d[0], d[2], cmp);
		cswap(d[0], d[1], cmp);
		cswap(d[4], d[5], cmp);
		cswap(d[3], d[5], cmp);
		cswap(d[3], d[4], cmp);
		cswap(d[0], d[3], cmp);
		cswap(d[1], d[4], cmp);
		cswap(d[2], d[5], cmp);
		cswap(d[2], d[4], cmp);
		cswap(d[1], d[3], cmp);
		cswap(d[2], d[3], cmp);
	}

	template<typename T, typename S, typename CMP>
	_PFROST_D_ void devInsertionSort(T* data, const S& size, CMP cmp)
	{
		int i, j;
#pragma unroll
		for (i = 1; i < size; i++) {
			T tmp = data[i];
			for (j = i; j > 0 && cmp(tmp, data[j - 1]); j--) data[j] = data[j - 1];
			data[j] = tmp;
		}
	}

	template<typename T, typename S, typename CMP>
	_PFROST_D_ void devSort(T* data, const S& size, CMP cmp)
	{
		if (size <= 1) return;
		assert(data != NULL);
		if (size == 2) cswap(data[0], data[1], cmp);
		else if (size == 3) sort3(data[0], data[1], data[2], cmp);
		else if (size == 4) sort4(data, cmp);
		else if (size == 5) sort5(data, cmp);
		else if (size == 6) sort6(data, cmp);
		else devInsertionSort(data, size, cmp);
		assert(devIsSorted(data, size, cmp));
	}

	template<typename T>
	_PFROST_D_ void devSort(T* data, const int& size)
	{
		if (size <= 1) return;
		assert(data != NULL);
		devSort(data, size, GPU_DEFAULT_CMP<T>());
		assert(devIsSorted(data, size, GPU_DEFAULT_CMP<T>()));
	}

	template<typename T, typename CMP>
	_PFROST_D_ void devMergeSort(T* data, T* aux, const int& size, const int& from, const int& mid, const int& to, CMP cmp)
	{
		int k = from, i = from, j = mid + 1;
		while (i <= mid && j <= to) {
			if (cmp(data[i], data[j])) aux[k++] = data[i++];
			else					   aux[k++] = data[j++];
		}
		while (i < size && i <= mid) aux[k++] = data[i++];
		// copy sorted segment to input data
#pragma unroll
		for (k = from; k <= to; k++) data[k] = aux[k];
	}

	template <typename T, class CMP>
	_PFROST_D_ void devMergeSort(T* data, T* aux, const int& size, CMP cmp)
	{
		int high = size - 1;
		// divide & conquer to power-of-2 segments
#pragma unroll
		for (int s = 1; s < size; s <<= 1) {
			int ds = s << 1;
#pragma unroll
			for (int i = 0; i < high; i += ds) {
				int mid = i + s - 1;
				int to = MIN(i + ds - 1, high);
				devMergeSort(data, aux, size, i, mid, to, cmp);
			}
		}
		assert(devIsSorted(data, size, cmp));
	}

	template<typename T>
	_PFROST_IN_D_ void devSwap(T& a, T& b) { T c = a; a = b, b = c; }


	_PFROST_IN_D_ void calcSig(SCLAUSE& c)
	{
		if (c.size() <= 1) return;
		uint32 sig = 0;
		forall_clause(c, lit) {
			sig |= MAPHASH(*lit);
		}
		c.set_sig(sig);
	}

	_PFROST_IN_D_ void calcSig(uint32* data, const int& size, uint32& _sig)
	{
		if (size <= 1) return;
		assert(_sig == 0);
		uint32* end = data + size;
		while (data != end)
			_sig |= MAPHASH(*data++);
	}



} 

#endif