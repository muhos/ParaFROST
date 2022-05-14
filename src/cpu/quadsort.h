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
quadsort 1.1.6.0

Changes to 1.1.5.2:
- unified data types with templates
- moved definitions to the header file
- changed function pointers to functors
to support adding members
- changed integer comparisons to Booleans
*/

#ifndef __QUADSORT_
#define __QUADSORT_

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cerrno>
#include <cstring>

#define parity_merge_two(data, swap, x, y, ptl, ptr, pts, cmp)  \
{  \
	ptl = data + 0; ptr = data + 2; pts = swap + 0;  \
	x = cmp(ptl, ptr); y = !x; pts[x] = *ptr; ptr += y; pts[y] = *ptl; ptl += x; pts++;  \
	*pts = cmp(ptl, ptr) ? *ptl : *ptr;  \
  \
	ptl = data + 1; ptr = data + 3; pts = swap + 3;  \
	x = cmp(ptl, ptr); y = !x; pts--; pts[x] = *ptr; ptr -= x; pts[y] = *ptl; ptl -= y;  \
	*pts = !cmp(ptl, ptr) ? *ptl : *ptr;  \
}

#define parity_merge_four(data, swap, x, y, ptl, ptr, pts, cmp)  \
{  \
	ptl = data + 0; ptr = data + 4; pts = swap;  \
	x = cmp(ptl, ptr); y = !x; pts[x] = *ptr; ptr += y; pts[y] = *ptl; ptl += x; pts++;  \
	x = cmp(ptl, ptr); y = !x; pts[x] = *ptr; ptr += y; pts[y] = *ptl; ptl += x; pts++;  \
	x = cmp(ptl, ptr); y = !x; pts[x] = *ptr; ptr += y; pts[y] = *ptl; ptl += x; pts++;  \
	*pts = cmp(ptl, ptr) ? *ptl : *ptr;  \
  \
	ptl = data + 3; ptr = data + 7; pts = swap + 7;  \
	x = cmp(ptl, ptr); y = !x; pts--; pts[x] = *ptr; ptr -= x; pts[y] = *ptl; ptl -= y;  \
	x = cmp(ptl, ptr); y = !x; pts--; pts[x] = *ptr; ptr -= x; pts[y] = *ptl; ptl -= y;  \
	x = cmp(ptl, ptr); y = !x; pts--; pts[x] = *ptr; ptr -= x; pts[y] = *ptl; ptl -= y;  \
	*pts = !cmp(ptl, ptr) ? *ptl : *ptr;  \
}

template <class T, class CMP>
void unguarded_insert(T *data, size_t offset, size_t nmemb, CMP cmp)
{
	T key, *pta, *end;
	size_t i, top;

	for (i = offset ; i < nmemb ; ++i)
	{
		pta = end = data + i;

		if (cmp(--pta, end))
		{
			continue;
		}

		key = *end;

		if (!cmp(data, &key))
		{
			top = i;

			do
			{
				*end-- = *pta--;
			}
			while (--top);

			*end = key;
		}
		else
		{
			do
			{
				*end-- = *pta--;
			}
			while (!cmp(pta, &key));

			*end = key;
		}
	}
}

template <class T, class CMP>
void bubble_sort(T *data, size_t nmemb, CMP cmp)
{
	T swap, *pta;
	size_t x, y;

	if (nmemb > 1)
	{
		pta = data;

		if (nmemb > 2)
		{
			y = cmp(pta, pta + 1); x = !y; swap = pta[y]; pta[0] = pta[x]; pta[1] = swap; pta++;
			y = cmp(pta, pta + 1); x = !y; swap = pta[y]; pta[0] = pta[x]; pta[1] = swap; pta--;
		}
		y = cmp(pta, pta + 1); x = !y; swap = pta[y]; pta[0] = pta[x]; pta[1] = swap;
	}
}

template <class T, class CMP>
void quad_swap_four(T *data, CMP cmp)
{
	T *pta, swap;
	size_t x, y;

	pta = data;
	y = cmp(pta, pta + 1); x = !y; swap = pta[y]; pta[0] = pta[x]; pta[1] = swap; pta += 2;
	y = cmp(pta, pta + 1); x = !y; swap = pta[y]; pta[0] = pta[x]; pta[1] = swap; pta--;

	if (!cmp(pta, pta + 1))
	{
		swap = pta[0]; pta[0] = pta[1]; pta[1] = swap; pta--;

		y = cmp(pta, pta + 1); x = !y; swap = pta[y]; pta[0] = pta[x]; pta[1] = swap; pta += 2;
		y = cmp(pta, pta + 1); x = !y; swap = pta[y]; pta[0] = pta[x]; pta[1] = swap; pta--;
		y = cmp(pta, pta + 1); x = !y; swap = pta[y]; pta[0] = pta[x]; pta[1] = swap;
	}
}

template <class T, class CMP>
void parity_swap_eight(T *data, CMP cmp)
{
	T swap[8], *ptl, *ptr, *pts;
	unsigned char x, y;

	ptl = data;
	y = cmp(ptl, ptl + 1); x = !y; swap[0] = ptl[y]; ptl[0] = ptl[x]; ptl[1] = swap[0]; ptl += 2;
	y = cmp(ptl, ptl + 1); x = !y; swap[0] = ptl[y]; ptl[0] = ptl[x]; ptl[1] = swap[0]; ptl += 2;
	y = cmp(ptl, ptl + 1); x = !y; swap[0] = ptl[y]; ptl[0] = ptl[x]; ptl[1] = swap[0]; ptl += 2;
	y = cmp(ptl, ptl + 1); x = !y; swap[0] = ptl[y]; ptl[0] = ptl[x]; ptl[1] = swap[0];

	if (cmp(data + 1, data + 2) && cmp(data + 3, data + 4) && cmp(data + 5, data + 6))
	{
		return;
	}
	parity_merge_two(data + 0, swap + 0, x, y, ptl, ptr, pts, cmp);
	parity_merge_two(data + 4, swap + 4, x, y, ptl, ptr, pts, cmp);

	parity_merge_four(swap, data, x, y, ptl, ptr, pts, cmp);
}

template <class T, class CMP>
void parity_merge(T *dest, T *from, size_t block, size_t nmemb, CMP cmp)
{
	T *ptl, *ptr, *tpl, *tpr, *tpd, *ptd;
	unsigned char x, y;

	ptl = from;
	ptr = from + block;
	ptd = dest;
	tpl = from + block - 1;
	tpr = from + nmemb - 1;
	tpd = dest + nmemb - 1;

	for (block-- ; block ; block--)
	{
		x = cmp(ptl, ptr); y = !x; ptd[x] = *ptr; ptr += y; ptd[y] = *ptl; ptl += x; ptd++;
		x = cmp(tpl, tpr); y = !x; tpd--; tpd[x] = *tpr; tpr -= x; tpd[y] = *tpl; tpl -= y;
	}
	*ptd = cmp(ptl, ptr) ? *ptl : *ptr;
	*tpd = !cmp(tpl, tpr) ? *tpl : *tpr;
}

template <class T, class CMP>
void parity_swap_sixteen(T *data, CMP cmp)
{
	T swap[16], *ptl, *ptr, *pts;
	unsigned char x, y;

	quad_swap_four(data +  0, cmp);
	quad_swap_four(data +  4, cmp);
	quad_swap_four(data +  8, cmp);
	quad_swap_four(data + 12, cmp);

	if (cmp(data + 3, data + 4) && cmp(data + 7, data + 8) && cmp(data + 11, data + 12))
	{
		return;
	}
	parity_merge_four(data + 0, swap + 0, x, y, ptl, ptr, pts, cmp);
	parity_merge_four(data + 8, swap + 8, x, y, ptl, ptr, pts, cmp);

	parity_merge(data, swap, 8, 16, cmp);
}

template <class T, class CMP>
void tail_swap(T *data, size_t nmemb, CMP cmp)
{
	if (nmemb < 4)
	{
		bubble_sort(data, nmemb, cmp);
		return;
	}
	if (nmemb < 8)
	{
		quad_swap_four(data, cmp);
		unguarded_insert(data, 4, nmemb, cmp);
		return;
	}
	if (nmemb < 16)
	{
		parity_swap_eight(data, cmp);
		unguarded_insert(data, 8, nmemb, cmp);
		return;
	}
	parity_swap_sixteen(data, cmp);
	unguarded_insert(data, 16, nmemb, cmp);
}

// the next three functions create sorted blocks of 32 elements

template <class T, class CMP>
void parity_tail_swap_eight(T *data, CMP cmp)
{
	T swap[8], *ptl, *ptr, *pts;
	unsigned char x, y;

	if (!cmp(data + 4, data + 5)) { swap[5] = data[4]; data[4] = data[5]; data[5] = swap[5]; }
	if (!cmp(data + 6, data + 7)) { swap[7] = data[6]; data[6] = data[7]; data[7] = swap[7]; } else

		if (cmp(data + 3, data + 4) && cmp(data + 5, data + 6))
		{
			return;
		}
	swap[0] = data[0]; swap[1] = data[1]; swap[2] = data[2]; swap[3] = data[3];

	parity_merge_two(data + 4, swap + 4, x, y, ptl, ptr, pts, cmp);

	parity_merge_four(swap, data, x, y, ptl, ptr, pts, cmp);
}

template <class T, class CMP>
void parity_tail_flip_eight(T *data, CMP cmp)
{
	T swap[8], *ptl, *ptr, *pts;
	unsigned char x, y;

	if (cmp(data + 3, data + 4))
	{
		return;
	}
	swap[0] = data[0]; swap[1] = data[1]; swap[2] = data[2]; swap[3] = data[3];
	swap[4] = data[4]; swap[5] = data[5]; swap[6] = data[6]; swap[7] = data[7];

	parity_merge_four(swap, data, x, y, ptl, ptr, pts, cmp);
}

template <class T, class CMP>
void tail_merge(T *data, T *swap, size_t swap_size, size_t nmemb, size_t block, CMP cmp);

template <class T, class CMP>
size_t quad_swap(T *data, size_t nmemb, CMP cmp)
{
	T swap[32];
	size_t count, reverse, x, y;
	T *pta, *pts, *pte, tmp;

	pta = data;

	count = nmemb / 8 * 2;

	while (count--)
	{
		switch ((!cmp(pta, pta + 1)) | (!cmp(pta + 1, pta + 2)) * 2 | (!cmp(pta + 2, pta + 3)) * 4)
		{
		case 0:
			break;
		case 1:
			tmp = pta[0]; pta[0] = pta[1]; pta[1] = tmp;
			pta += 1; y = cmp(pta, pta + 1); x = !y; tmp = pta[y]; pta[0] = pta[x]; pta[1] = tmp;
			pta += 1; y = cmp(pta, pta + 1); x = !y; tmp = pta[y]; pta[0] = pta[x]; pta[1] = tmp;
			pta -= 2;
			break;
		case 2:
			tmp = pta[1]; pta[1] = pta[2]; pta[2] = tmp;
			y = cmp(pta, pta + 1); x = !y; tmp = pta[y]; pta[0] = pta[x]; pta[1] = tmp;
			pta += 2; y = cmp(pta, pta + 1); x = !y; tmp = pta[y]; pta[0] = pta[x]; pta[1] = tmp;
			pta -= 1; y = cmp(pta, pta + 1); x = !y; tmp = pta[y]; pta[0] = pta[x]; pta[1] = tmp;
			pta -= 1;
			break;
		case 3:
			tmp = pta[0]; pta[0] = pta[2]; pta[2] = tmp;
			pta += 2; y = cmp(pta, pta + 1); x = !y; tmp = pta[y]; pta[0] = pta[x]; pta[1] = tmp;
			pta -= 1; y = cmp(pta, pta + 1); x = !y; tmp = pta[y]; pta[0] = pta[x]; pta[1] = tmp;
			pta -= 1;
			break;
		case 4:
			tmp = pta[2]; pta[2] = pta[3]; pta[3] = tmp;
			pta += 1; y = cmp(pta, pta + 1); x = !y; tmp = pta[y]; pta[0] = pta[x]; pta[1] = tmp;
			pta -= 1; y = cmp(pta, pta + 1); x = !y; tmp = pta[y]; pta[0] = pta[x]; pta[1] = tmp;
			break;
		case 5:
			tmp = pta[0]; pta[0] = pta[1]; pta[1] = tmp;
			tmp = pta[2]; pta[2] = pta[3]; pta[3] = tmp;
			pta += 1; y = cmp(pta, pta + 1); x = !y; tmp = pta[y]; pta[0] = pta[x]; pta[1] = tmp;
			pta += 1; y = cmp(pta, pta + 1); x = !y; tmp = pta[y]; pta[0] = pta[x]; pta[1] = tmp;
			pta -= 2; y = cmp(pta, pta + 1); x = !y; tmp = pta[y]; pta[0] = pta[x]; pta[1] = tmp;
			break;
		case 6:
			tmp = pta[1]; pta[1] = pta[3]; pta[3] = tmp;
			y = cmp(pta, pta + 1); x = !y; tmp = pta[y]; pta[0] = pta[x]; pta[1] = tmp;
			pta += 1; y = cmp(pta, pta + 1); x = !y; tmp = pta[y]; pta[0] = pta[x]; pta[1] = tmp;
			pta -= 1;
			break;
		case 7:
			pts = pta;
			goto swapper;
		}
		count--;

		parity_tail_swap_eight(pta, cmp);

		pta += 8;

		continue;

	swapper:

		pta += 4;

		if (count--)
		{
			if (!cmp(&pta[0], &pta[1]))
			{
				if (!cmp(&pta[2], &pta[3]))
				{
					if (!cmp(&pta[1], &pta[2]))
					{
						if (!cmp(&pta[-1], &pta[0]))
						{
							goto swapper;
						}
					}
					tmp = pta[2]; pta[2] = pta[3]; pta[3] = tmp;
				}
				tmp = pta[0]; pta[0] = pta[1]; pta[1] = tmp;
			}
			else if (!cmp(&pta[2], &pta[3]))
			{
				tmp = pta[2]; pta[2] = pta[3]; pta[3] = tmp;
			}

			if (!cmp(&pta[1], &pta[2]))
			{
				tmp = pta[1]; pta[1] = pta[2]; pta[2] = tmp;

				y = cmp(pta, pta + 1); x = !y; tmp = pta[y]; pta[0] = pta[x]; pta[1] = tmp; pta += 2;
				y = cmp(pta, pta + 1); x = !y; tmp = pta[y]; pta[0] = pta[x]; pta[1] = tmp; pta -= 1;
				y = cmp(pta, pta + 1); x = !y; tmp = pta[y]; pta[0] = pta[x]; pta[1] = tmp; pta -= 1;
			}
			pte = pta - 1;

			reverse = (pte - pts) / 2;

			do
			{
				tmp = *pts; *pts++ = *pte; *pte-- = tmp;
			}
			while (reverse--);

			if (count % 2 == 0)
			{
				pta -= 4;

				parity_tail_flip_eight(pta, cmp);
			}
			else
			{
				count--;

				parity_tail_swap_eight(pta, cmp);
			}
			pta += 8;

			continue;
		}

		if (pts == data)
		{
			switch (nmemb % 8)
			{
			case 7: if (cmp(pta + 5, pta + 6)) break;
			case 6: if (cmp(pta + 4, pta + 5)) break;
			case 5: if (cmp(pta + 3, pta + 4)) break;
			case 4: if (cmp(pta + 2, pta + 3)) break;
			case 3: if (cmp(pta + 1, pta + 2)) break;
			case 2: if (cmp(pta + 0, pta + 1)) break;
			case 1: if (cmp(pta - 1, pta + 0)) break;
			case 0:
				pte = pts + nmemb - 1;

				reverse = (pte - pts) / 2;

				do
				{
					tmp = *pts; *pts++ = *pte; *pte-- = tmp;
				}
				while (reverse--);

				return 1;
			}
		}
		pte = pta - 1;

		reverse = (pte - pts) / 2;

		do
		{
			tmp = *pts; *pts++ = *pte; *pte-- = tmp;
		}
		while (reverse--);

		break;
	}

	tail_swap(pta, nmemb % 8, cmp);

	pta = data;

	for (count = nmemb / 32 ; count-- ; pta += 32)
	{
		if (cmp(pta + 7, pta + 8) && cmp(pta + 15, pta + 16) && cmp(pta + 23, pta + 24))
		{
			continue;
		}
		parity_merge(swap, pta, 8, 16, cmp);
		parity_merge(swap + 16, pta + 16, 8, 16, cmp);
		parity_merge(pta, swap, 16, 32, cmp);
	}

	if (nmemb % 32 > 8)
	{
		tail_merge(pta, swap, 32, nmemb % 32, 8, cmp);
	}
	return 0;
}

// quad merge support routines
template <class T, class CMP>
void forward_merge(T *dest, T *from, size_t block, CMP cmp)
{
	T *ptl, *ptr, *m, *e; // left, right, middle, end
	size_t x, y;

	ptl = from;
	ptr = from + block;
	m = ptr;
	e = ptr + block;

	if (cmp(m - 1, e - block / 4))
	{
		m -= 2;

		while (ptl < m)
		{
			if (cmp(ptl + 1, ptr))
			{
				*dest++ = *ptl++; *dest++ = *ptl++;
			}
			else if (!cmp(ptl, ptr + 1))
			{
				*dest++ = *ptr++; *dest++ = *ptr++;
			}
			else 
			{
				x = cmp(ptl, ptr); y = !x; dest[x] = *ptr; ptr += 1; dest[y] = *ptl; ptl += 1; dest += 2;
				x = cmp(ptl, ptr); y = !x; dest[x] = *ptr; ptr += y; dest[y] = *ptl; ptl += x; dest++;
			}
		}
		m += 2;

		while (ptl < m)
		{
			*dest++ = cmp(ptl, ptr) ? *ptl++ : *ptr++;
		}

		do *dest++ = *ptr++; while (ptr < e);
	}
	else if (!cmp(m - block / 4, e - 1))
	{
		e -= 2;

		while (ptr < e)
		{
			if (!cmp(ptl, ptr + 1))
			{
				*dest++ = *ptr++; *dest++ = *ptr++;
			}
			else if (cmp(ptl + 1, ptr))
			{
				*dest++ = *ptl++; *dest++ = *ptl++;
			}
			else 
			{
				x = cmp(ptl, ptr); y = !x; dest[x] = *ptr; ptr += 1; dest[y] = *ptl; ptl += 1; dest += 2;
				x = cmp(ptl, ptr); y = !x; dest[x] = *ptr; ptr += y; dest[y] = *ptl; ptl += x; dest++;
			}
		}
		e += 2;

		while (ptr < e)
		{
			*dest++ = !cmp(ptl, ptr) ? *ptr++ : *ptl++;
		}

		do *dest++ = *ptl++; while (ptl < m);
	}
	else
	{
		parity_merge(dest, from, block, block * 2, cmp);
	}
}

template <class T, class CMP>
void quad_merge_block(T *data, T *swap, size_t block, CMP cmp)
{
	T *pts, *c, *c_max;
	size_t block_x_2 = block * 2;

	c_max = data + block;

	if (cmp(c_max - 1, c_max))
	{
		c_max += block_x_2;

		if (cmp(c_max - 1, c_max))
		{
			c_max -= block;

			if (cmp(c_max - 1, c_max))
			{
				return;
			}
			pts = swap;

			c = data;

			do *pts++ = *c++; while (c < c_max); // step 1

			c_max = c + block_x_2;

			do *pts++ = *c++; while (c < c_max); // step 2

			return forward_merge(data, swap, block_x_2, cmp); // step 3
		}
		pts = swap;

		c = data;
		c_max = data + block_x_2;

		do *pts++ = *c++; while (c < c_max); // step 1
	}
	else
	{
		forward_merge(swap, data, block, cmp); // step 1
	}
	forward_merge(swap + block_x_2, data + block_x_2, block, cmp); // step 2

	forward_merge(data, swap, block_x_2, cmp); // step 3
}

template <class T, class CMP>
size_t quad_merge(T *data, T *swap, size_t swap_size, size_t nmemb, size_t block, CMP cmp)
{
	T *pta, *pte;

	pte = data + nmemb;

	block *= 4;

	while (block <= nmemb && block <= swap_size)
	{
		pta = data;

		do
		{
			quad_merge_block(pta, swap, block / 4, cmp);

			pta += block;
		}
		while (pta + block <= pte);

		tail_merge(pta, swap, swap_size, pte - pta, block / 4, cmp);

		block *= 4;
	}

	tail_merge(data, swap, swap_size, nmemb, block / 4, cmp);

	return block / 2;
}

template <class T, class CMP>
void partial_forward_merge(T *data, T *swap, size_t nmemb, size_t block, CMP cmp)
{
	T *r, *m, *e, *s; // right, middle, end, swap
	size_t x, y;

	r = data + block;
	e = data + nmemb - 1;

	memcpy(swap, data, block * sizeof(T));

	s = swap;
	m = swap + block - 1;

	while (s <= m - 2 && r <= e - 2)
	{
		if (!cmp(s, r + 1))
		{
			*data++ = *r++; *data++ = *r++;
		}
		else if (cmp(s + 1, r))
		{
			*data++ = *s++; *data++ = *s++;
		}
		else 
		{
			x = cmp(s, r); y = !x; data[x] = *r; r += 1; data[y] = *s; s += 1; data += 2;
			x = cmp(s, r); y = !x; data[x] = *r; r += y; data[y] = *s; s += x; data++;
		}
	}

	while (s <= m && r <= e)
	{
		*data++ = cmp(s, r) ? *s++ : *r++;
	}

	while (s <= m)
	{
		*data++ = *s++;
	}
}

template <class T, class CMP>
void partial_backward_merge(T *data, T *swap, size_t nmemb, size_t block, CMP cmp)
{
	T *m, *e, *s; // middle, end, swap
	size_t x, y;

	m = data + block - 1;
	e = data + nmemb - 1;

	if (cmp(m, m + 1))
	{
		return;
	}

	memcpy(swap, data + block, (nmemb - block) * sizeof(T));

	s = swap + nmemb - block - 1;

	while (s >= swap + 2 && m >= data + 2)
	{
		if (!cmp(m - 1, s))
		{
			*e-- = *m--;
			*e-- = *m--;
		}
		else if (cmp(m, s - 1))
		{
			*e-- = *s--;
			*e-- = *s--;
		}
		else
		{
			x = cmp(m, s); y = !x; e--; e[x] = *s; s -= 1; e[y] = *m; m -= 1; e--;
			x = cmp(m, s); y = !x; e--; e[x] = *s; s -= x; e[y] = *m; m -= y;
		}
	}

	while (s >= swap && m >= data)
	{
		*e-- = !cmp(m, s) ? *m-- : *s--;
	}

	while (s >= swap)
	{
		*e-- = *s--;
	}
}

template <class T, class CMP>
void tail_merge(T *data, T *swap, size_t swap_size, size_t nmemb, size_t block, CMP cmp)
{
	T *pta, *pte;

	pte = data + nmemb;

	while (block < nmemb && block <= swap_size)
	{
		pta = data;

		for (pta = data ; pta + block < pte ; pta += block * 2)
		{
			if (pta + block * 2 < pte)
			{
				partial_backward_merge(pta, swap, block * 2, block, cmp);

				continue;
			}
			partial_backward_merge(pta, swap, pte - pta, block, cmp);

			break;
		}
		block *= 2;
	}
}

// the next four functions provide in-place rotate merge support
template <class T>
void trinity_rotation(T *data, T *swap, size_t swap_size, size_t nmemb, size_t left)
{
	size_t bridge, right = nmemb - left;

	if (left < right)
	{
		if (left <= swap_size)
		{
			memcpy(swap, data, left * sizeof(T));
			memmove(data, data + left, right * sizeof(T));
			memcpy(data + right, swap, left * sizeof(T));
		}
		else
		{
			T *pta, *ptb, *ptc, *ptd;

			pta = data;
			ptb = pta + left;

			bridge = right - left;

			if (bridge <= swap_size && bridge > 3)
			{
				ptc = pta + right;
				ptd = ptc + left;

				memcpy(swap, ptb, bridge * sizeof(T));

				while (left--)
				{
					*--ptc = *--ptd; *ptd = *--ptb;
				}
				memcpy(pta, swap, bridge * sizeof(T));
			}
			else
			{
				ptc = ptb;
				ptd = ptc + right;

				bridge = left / 2;

				while (bridge--)
				{
					*swap = *--ptb; *ptb = *pta; *pta++ = *ptc; *ptc++ = *--ptd; *ptd = *swap;
				}

				bridge = (ptd - ptc) / 2;

				while (bridge--)
				{
					*swap = *ptc; *ptc++ = *--ptd; *ptd = *pta; *pta++ = *swap;
				}

				bridge = (ptd - pta) / 2;

				while (bridge--)
				{
					*swap = *pta; *pta++ = *--ptd; *ptd = *swap;
				}
			}
		}
	}
	else if (right < left)
	{
		if (right <= swap_size)
		{
			memcpy(swap, data + left, right * sizeof(T));
			memmove(data + right, data, left * sizeof(T));
			memcpy(data, swap, right * sizeof(T));
		}
		else
		{
			T *pta, *ptb, *ptc, *ptd;

			pta = data;
			ptb = pta + left;

			bridge = left - right;

			if (bridge <= swap_size && bridge > 3)
			{
				ptc = pta + right;
				ptd = ptc + left;

				memcpy(swap, ptc, bridge * sizeof(T));

				while (right--)
				{
					*ptc++ = *pta; *pta++ = *ptb++;
				}
				memcpy(ptd - bridge, swap, bridge * sizeof(T));
			}
			else
			{
				ptc = ptb;
				ptd = ptc + right;

				bridge = right / 2;

				while (bridge--)
				{
					*swap = *--ptb; *ptb = *pta; *pta++ = *ptc; *ptc++ = *--ptd; *ptd = *swap;
				}

				bridge = (ptb - pta) / 2;

				while (bridge--)
				{
					*swap = *--ptb; *ptb = *pta; *pta++ = *--ptd; *ptd = *swap;
				}

				bridge = (ptd - pta) / 2;

				while (bridge--)
				{
					*swap = *pta; *pta++ = *--ptd; *ptd = *swap;
				}
			}
		}
	}
	else
	{
		T *pta, *ptb;

		pta = data;
		ptb = pta + left;

		while (left--)
		{
			*swap = *pta; *pta++ = *ptb; *ptb++ = *swap;
		}
	}
}

template <class T, class CMP>
size_t monobound_binary_first(T *data, T *value, size_t top, CMP cmp)
{
	T *end;
	size_t mid;

	end = data + top;

	while (top > 1)
	{
		mid = top / 2;

		if (cmp(value, end - mid))
		{
			end -= mid;
		}
		top -= mid;
	}

	if (cmp(value, end - 1))
	{
		end--;
	}
	return (end - data);
}

template <class T, class CMP>
void blit_merge_block(T *data, T *swap, size_t swap_size, size_t block, size_t right, CMP cmp)
{
	size_t left;

	if (cmp(data + block - 1, data + block))
	{
		return;
	}

	left = monobound_binary_first(data + block, data + block / 2, right, cmp);

	right -= left;

	block /= 2;

	if (left)
	{
		trinity_rotation(data + block, swap, swap_size, block + left, block);

		if (left <= swap_size)
		{
			partial_backward_merge(data, swap, block + left, block, cmp);
		}
		else if (block <= swap_size)
		{
			partial_forward_merge(data, swap, block + left, block, cmp);
		}
		else
		{
			blit_merge_block(data, swap, swap_size, block, left, cmp);
		}
	}

	if (right)
	{
		if (right <= swap_size)
		{
			partial_backward_merge(data + block + left, swap, block + right, block, cmp);
		}
		else if (block <= swap_size)
		{
			partial_forward_merge(data + block + left, swap, block + right, block, cmp);
		}
		else
		{
			blit_merge_block(data + block + left, swap, swap_size, block, right, cmp);
		}
	}
}

template <class T, class CMP>
void blit_merge(T *data, T *swap, size_t swap_size, size_t nmemb, size_t block, CMP cmp)
{
	T *pta, *pte;

	pte = data + nmemb;

	while (block < nmemb)
	{
		pta = data;

		for (pta = data ; pta + block < pte ; pta += block * 2)
		{
			if (pta + block * 2 < pte)
			{
				blit_merge_block(pta, swap, swap_size, block, block, cmp);

				continue;
			}
			blit_merge_block(pta, swap, swap_size, block, pte - pta - block, cmp);

			break;
		}
		block *= 2;
	}
}

template <class T, class CMP>
void quadsort(T *data, size_t nmemb, CMP cmp)
{
	if (nmemb < 32)
	{
		tail_swap(data, nmemb, cmp);
	}
	else if (quad_swap(data, nmemb, cmp) == 0)
	{
		T *swap;
		size_t swap_size = 32;

		while (swap_size * 4 <= nmemb)
		{
			swap_size *= 4;
		}
		swap = (T *)malloc(swap_size * sizeof(T));

		if (swap == NULL)
		{
			swap = (T *)malloc((swap_size = 512) * sizeof(T));

			if (swap == NULL)
			{
				T stack[32];

				tail_merge(data, stack, 32, nmemb, 32, cmp);

				blit_merge(data, stack, 32, nmemb, 64, cmp);

				return;
			}
		}
		quad_merge(data, swap, swap_size, nmemb, 32, cmp);

		blit_merge(data, swap, swap_size, nmemb, swap_size * 2, cmp);

		free(swap);
	}
}

template <class T, class CMP>
void quadsort_swap(T *data, T *swap, size_t swap_size, size_t nmemb, CMP cmp)
{
	if (nmemb < 32)
	{
		tail_swap(data, nmemb, cmp);
	}
	else if (quad_swap(data, nmemb, cmp) == 0)
	{
		size_t block;

		block = quad_merge(data, swap, swap_size, nmemb, 32, cmp);

		blit_merge(data, swap, swap_size, nmemb, block, cmp);
	}
}

#endif
