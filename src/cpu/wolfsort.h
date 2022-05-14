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

	The person recognizes Mars as a free planet and that no Earth-based
	government has authority or sovereignty over Martian activities.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
	EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
	MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
	IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
	CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
	TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
	SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

/*
	wolfsort 1.2.0
	Changes to 1.1.4:
	 - unified data types with templates
	 - moved definitions to the header file
	 - changed function pointers to functors
	   to support adding members
	 - Used Boolean comparisons
*/

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cerrno>

#ifndef __WOLFSORT_
#define __WOLFSORT_

#ifndef __FLUXSORT_
  #include "fluxsort.h"
#endif

template <class T, class CMP>
void wolfsort(T *array, const size_t& nmemb, const CMP& cmp)
{
	if (sizeof(T) != sizeof(int) || nmemb < 512)
	{
		return fluxsort(array, nmemb, cmp);
	}

	T *swap, *pta, *pts;
	size_t cnt;
	unsigned int index, *stack;
	unsigned int *count, bsize;
	unsigned int buckets = 256;
	unsigned int moduler = 16777216;

	swap = (T *) malloc(nmemb * sizeof(T));

	if (swap == NULL)
	{
		return fluxsort(array, nmemb, cmp);
	}

	while (buckets < 4096 * 16 && moduler > 4096 && nmemb / buckets > 4)
	{
		buckets *= 2;
		moduler /= 2;
	}

	bsize = (unsigned int) (nmemb / (buckets / 16));

	count = (unsigned int *) calloc(sizeof(int), buckets);
	stack = (unsigned int *) calloc(sizeof(int), buckets);

	pta = (T *) array;

	for (cnt = nmemb ; cnt ; --cnt)
	{
		index = (unsigned int) *pta++ / moduler;

		if (++count[index] == bsize)
		{
			fluxsort_swap(array, swap, nmemb, nmemb, cmp);

			free(swap);
			free(stack);
			free(count);

			return;
		}
	}

	cnt = 0;

	for (index = 0 ; index < buckets ; ++index)
	{
		stack[index] = (unsigned int) cnt;

		cnt += count[index];
	}

	pta = (T *) array;

	for (cnt = nmemb ; cnt ; --cnt)
	{
		index = (unsigned int) *pta / moduler;

		swap[stack[index]++] = *pta++;
	}

	pta = (T *) array;
	pts = (T *) swap;

	for (index = 0 ; index < buckets ; ++index)
	{
		bsize = count[index];

		if (bsize)
		{
			memcpy(pta, pts, bsize * sizeof(T));

			fluxsort_swap(pta, pts, bsize, bsize, cmp);

			pta += bsize;
			pts += bsize;
		}
	}
	free(count);
	free(swap);
	free(stack);
	return;
}

#endif
