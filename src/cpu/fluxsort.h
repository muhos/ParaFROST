/*
	Copyright (C) 2014-2021 Igor van den Hoven ivdhoven@gmail.com
	Copyright (C) 2022 Muhammad Osama Mahmoud
*/

/*
	Permission is hereby granted, free of charge, to any person obtaining
	a copy of this software and associated documentation files (the
	"Software"), to deal in the Software without restriction, including
	without limitation the rights to use, copy, modify, merge, publish,
	distribute, sublicense, and/or sell copies of the Software, and to
	permit persons to whom the Software is furnished to do so, subject to
	the following conditions:

	The above copyright notice and this permission notice shall be
	included in all copies or substantial portions of the Software.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
	EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
	MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
	IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
	CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
	TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
	SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

/*
	fluxsort 1.1.5.0

	Changes to 1.1.4.3:
	- unified data types with templates
	- moved definitions to the header file
	- changed function pointers to functors
	  to support adding members
	- changed integer comparisons to Booleans
*/

#ifndef __FLUXSORT_
#define __FLUXSORT_

#ifndef __QUADSORT_
  #include "quadsort.h"
#endif

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cerrno>

#define FLUX_OUT 24

// Determine whether to use mergesort or quicksort

template <class T, class CMP>
size_t flux_analyze(T *data, size_t nmemb, CMP cmp)
{
	char loop, dist;
	size_t cnt, balance = 0, streaks = 0;
	T *pta, *ptb, swap;

	pta = data;

	for (cnt = nmemb ; cnt > 16 ; cnt -= 16)
	{
		for (dist = 0, loop = 16 ; loop ; loop--)
		{
			dist += !cmp(pta, pta + 1); pta++;
		}
		streaks += (dist == 0) | (dist == 16);
		balance += dist;
	}

	while (--cnt)
	{
		balance += !cmp(pta, pta + 1);
		pta++;
	}

	if (balance == 0)
	{
		return 1;
	}

	if (balance == nmemb - 1)
	{
		ptb = pta + 1;
		pta = data;

		cnt = nmemb / 2;

		do
		{
			swap = *pta; *pta++ = *--ptb; *ptb = swap;
		}
		while (--cnt);

		return 1;
	}

	if (streaks >= nmemb / 40)
	{
		quadsort(data, nmemb, cmp);

		return 1;
	}
	return 0;
}

template <class T, class CMP>
T median_of_sqrt(T *data, T *swap, T *ptx, size_t nmemb, CMP cmp)
{
	T *pta, *pts;
	size_t cnt, sqrt, div;

	sqrt = nmemb > 262144 ? 256 : 128;

	div = nmemb / sqrt;

	pta = ptx + rand() % sqrt;

	pts = ptx == data ? swap : data;

	for (cnt = 0 ; cnt < sqrt ; ++cnt)
	{
		pts[cnt] = pta[0];

		pta += div;
	}
	quadsort_swap(pts, pts + sqrt, sqrt, sqrt, cmp);

	return pts[sqrt / 2];
}

template <class T, class CMP>
T median_of_five(T *data, size_t v0, size_t v1, size_t v2, size_t v3, size_t v4, CMP cmp)
{
	T swap[6], *pta;
	size_t x, y, z;

	swap[2] = data[v0];
	swap[3] = data[v1];
	swap[4] = data[v2];
	swap[5] = data[v3];

	pta = swap + 2;

	y = cmp(pta, pta + 1); x = !y; swap[0] = pta[y]; pta[0] = pta[x]; pta[1] = swap[0]; pta += 2;
	y = cmp(pta, pta + 1); x = !y; swap[0] = pta[y]; pta[0] = pta[x]; pta[1] = swap[0]; pta -= 2;
	y = cmp(pta, pta + 2); x = !y; swap[0] = pta[0]; swap[1] = pta[2]; pta[0] = swap[x]; pta[2] = swap[y]; pta++;
	y = cmp(pta, pta + 2); x = !y; swap[0] = pta[0]; swap[1] = pta[2]; pta[0] = swap[x]; pta[2] = swap[y];

	pta[2] = data[v4];

	x = !cmp(pta, pta + 1);
	y = !cmp(pta, pta + 2);
	z = !cmp(pta + 1, pta + 2);

	return pta[(x == y) + (y ^ z)];
}

template <class T, class CMP>
T median_of_twentyfive(T *data, size_t nmemb, CMP cmp)
{
	T swap[5];
	size_t div = nmemb / 64;

	swap[0] = median_of_five(data, div *  4, div *  1, div *  2, div *  8, div * 10, cmp);
	swap[1] = median_of_five(data, div * 16, div * 12, div * 14, div * 18, div * 20, cmp);
	swap[2] = median_of_five(data, div * 32, div * 24, div * 30, div * 34, div * 38, cmp);
	swap[3] = median_of_five(data, div * 48, div * 42, div * 44, div * 50, div * 52, cmp);
	swap[4] = median_of_five(data, div * 60, div * 54, div * 56, div * 62, div * 63, cmp);

	return median_of_five(swap, 0, 1, 2, 3, 4, cmp);
}

template <class T, class CMP>
size_t median_of_three(T *data, size_t v0, size_t v1, size_t v2, CMP cmp)
{
	size_t v[3] = {v0, v1, v2};
	char x, y, z;

	x = !cmp(data + v0, data + v1);
	y = !cmp(data + v0, data + v2);
	z = !cmp(data + v1, data + v2);

	return v[(x == y) + (y ^ z)];
}

template <class T, class CMP>
T median_of_nine(T *data, size_t nmemb, CMP cmp)
{
	size_t x, y, z, div = nmemb / 16;

	x = median_of_three(data, div * 2, div * 1, div * 4, cmp);
	y = median_of_three(data, div * 8, div * 6, div * 10, cmp);
	z = median_of_three(data, div * 14, div * 12, div * 15, cmp);

	return data[median_of_three(data, x, y, z, cmp)];
}

template <class T, class CMP>
void flux_partition(T *data, T *swap, T *ptx, T *ptp, size_t nmemb, CMP cmp);

#if defined(_WIN32)
#pragma warning(push)
#pragma warning(disable : 4146)
#endif

// As per suggestion by Marshall Lochbaum to improve generic data handling
template <class T, class CMP>
void flux_reverse_partition(T *data, T *swap, T *ptx, T *piv, size_t nmemb, CMP cmp)
{
	size_t a_size, s_size;

	{
		size_t cnt, val, m;
		T *pts = swap;

		for (m = 0, cnt = nmemb / 8 ; cnt ; cnt--)
		{
			val = !cmp(piv, ptx); pts[-m] = data[m] = *ptx++; m += val; pts++;
			val = !cmp(piv, ptx); pts[-m] = data[m] = *ptx++; m += val; pts++;
			val = !cmp(piv, ptx); pts[-m] = data[m] = *ptx++; m += val; pts++;
			val = !cmp(piv, ptx); pts[-m] = data[m] = *ptx++; m += val; pts++;
			val = !cmp(piv, ptx); pts[-m] = data[m] = *ptx++; m += val; pts++;
			val = !cmp(piv, ptx); pts[-m] = data[m] = *ptx++; m += val; pts++;
			val = !cmp(piv, ptx); pts[-m] = data[m] = *ptx++; m += val; pts++;
			val = !cmp(piv, ptx); pts[-m] = data[m] = *ptx++; m += val; pts++;
		}

		for (cnt = nmemb % 8 ; cnt ; cnt--)
		{
			val = !cmp(piv, ptx); pts[-m] = data[m] = *ptx++; m += val; pts++;
		}
		a_size = m;
		s_size = nmemb - a_size;
	}
	memcpy(data + a_size, swap, s_size * sizeof(T));

	if (s_size <= a_size / 16 || a_size <= FLUX_OUT)
	{
		return quadsort_swap(data, swap, a_size, a_size, cmp);
	}
	flux_partition(data, swap, data, piv, a_size, cmp);
}

template <class T, class CMP>
size_t flux_default_partition(T *data, T *swap, T *ptx, T *piv, size_t nmemb, CMP cmp)
{
	size_t cnt, val, m = 0;

	for (cnt = nmemb / 8 ; cnt ; cnt--)
	{
		val = cmp(ptx, piv); swap[-m] = data[m] = *ptx++; m += val; swap++;
		val = cmp(ptx, piv); swap[-m] = data[m] = *ptx++; m += val; swap++;
		val = cmp(ptx, piv); swap[-m] = data[m] = *ptx++; m += val; swap++;
		val = cmp(ptx, piv); swap[-m] = data[m] = *ptx++; m += val; swap++;
		val = cmp(ptx, piv); swap[-m] = data[m] = *ptx++; m += val; swap++;
		val = cmp(ptx, piv); swap[-m] = data[m] = *ptx++; m += val; swap++;
		val = cmp(ptx, piv); swap[-m] = data[m] = *ptx++; m += val; swap++;
		val = cmp(ptx, piv); swap[-m] = data[m] = *ptx++; m += val; swap++;
	}

	for (cnt = nmemb % 8 ; cnt ; cnt--)
	{
		val = cmp(ptx, piv); swap[-m] = data[m] = *ptx++; m += val; swap++;
	}
	return m;
}

#if defined(_WIN32)
#pragma warning(pop)
#endif

template <class T, class CMP>
void flux_partition(T *data, T *swap, T *ptx, T *piv, size_t nmemb, CMP cmp)
{
	size_t a_size = 0, s_size;

	while (1)
	{
		--piv;

		if (nmemb <= 2048)
		{
			*piv = median_of_nine(ptx, nmemb, cmp);
		}
		else if (nmemb <= 65536)
		{
			*piv = median_of_twentyfive(ptx, nmemb, cmp);
		}
		else
		{
			*piv = median_of_sqrt(data, swap, ptx, nmemb, cmp);
		}

		if (a_size && cmp(piv + 1, piv))
		{
			return flux_reverse_partition(data, swap, data, piv, nmemb, cmp);
		}
		a_size = flux_default_partition(data, swap, ptx, piv, nmemb, cmp);
		s_size = nmemb - a_size;

		if (a_size <= s_size / 16 || s_size <= FLUX_OUT)
		{
			if (s_size == 0)
			{
				return flux_reverse_partition(data, swap, data, piv, a_size, cmp);
			}
			memcpy(data + a_size, swap, s_size * sizeof(T));
			quadsort_swap(data + a_size, swap, s_size, s_size, cmp);
		}
		else
		{
			flux_partition(data + a_size, swap, swap, piv, s_size, cmp);
		}

		if (s_size <= a_size / 16 || a_size <= FLUX_OUT)
		{
			return quadsort_swap(data, swap, a_size, a_size, cmp);
		}
		nmemb = a_size;
		ptx = data;
	}
}

template <class T, class CMP>
void fluxsort(T *data, size_t nmemb, CMP cmp)
{
	if (nmemb < 32)
	{
		tail_swap(data, nmemb, cmp);
	}
	else if (flux_analyze(data, nmemb, cmp) == 0)
	{
		T* swap = (T *) malloc(nmemb * sizeof(T));

		if (swap == NULL)
		{
			return quadsort(data, nmemb, cmp);
		}
		flux_partition(data, swap, data, swap + nmemb, nmemb, cmp);

		free(swap);
	}
}

template <class T, class CMP>
void fluxsort_swap(T *data, T *swap, size_t swap_size, size_t nmemb, CMP cmp)
{
	if (nmemb < 32)
	{
		tail_swap(data, nmemb, cmp);
	}
	else if (swap_size < nmemb)
	{
		quadsort_swap(data, swap, swap_size, nmemb, cmp);
	}
	else if (flux_analyze(data, nmemb, cmp) == 0)
	{
		flux_partition(data, swap, data, swap + nmemb, nmemb, cmp);
	}
}

#endif
