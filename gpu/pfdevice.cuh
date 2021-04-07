/***********************************************************************[pfdevice.h]
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

#ifndef __SIGMA_DEVICE_
#define __SIGMA_DEVICE_

#include "pfgrid.cuh"
#include "pfatomics.cuh"
#include "pfsimp.cuh"

namespace pFROST {

	namespace SIGmA {

		// [0] = hse_limit, [1] = bce_limit,
		// [2] = ere_limit, [3] = xor_max_arity
		__constant__ uint32 dc_limits[NLIMITS];
		// number of buckets to store 'SCLAUSE'
		__constant__ int dc_nbuckets = sizeof(SCLAUSE) / sizeof(uint32);

		template<class T>
		struct DEFAULT_CMP {
			_PFROST_D_ bool operator () (T& x, T& y) {
				return x < y;
			}
		};

		_PFROST_H_D_ void pLit(const uint32& l) { printf("%d", SIGN(l) ? -int(ABS(l)) : ABS(l)); }

		_PFROST_H_D_ void pSharedClause(uint32* c, const int& sz)
		{
			printf("Shared Clause (size = %d)->[", sz);
			for (int k = 0; k < sz; k++) pLit(c[k]), printf("  ");
			printf("]\n");
		}

		_PFROST_H_D_ void pClauseSet(CNF& cnf, OL& poss, OL& negs)
		{
			for (uint32 i = 0; i < poss.size(); i++) {
				printf("c | "), cnf[poss[i]].print();
			}
			for (uint32 i = 0; i < negs.size(); i++) {
				printf("c | "), cnf[negs[i]].print();
			}
		}

		_PFROST_H_D_ void pClauseSet(CNF& cnf, OL& ol)
		{
			for (uint32 i = 0; i < ol.size(); i++) {
				SCLAUSE& c = cnf[ol[i]];
				printf("c | c(%lld)->", uint64(ol[i])), c.print();
			}
		}

		_PFROST_H_D_ bool checkMolten(CNF& cnf, OL& poss, OL& negs)
		{
			for (uint32 i = 0; i < poss.size(); i++)
				if (cnf[poss[i]].molten()) return false;
			for (uint32 i = 0; i < negs.size(); i++)
				if (cnf[negs[i]].molten()) return false;
			return true;
		}

		_PFROST_H_D_ bool checkDeleted(CNF& cnf, OL& poss, OL& negs)
		{
			for (uint32 i = 0; i < poss.size(); i++)
				if (cnf[poss[i]].deleted()) return false;
			for (uint32 i = 0; i < negs.size(); i++)
				if (cnf[negs[i]].deleted()) return false;
			return true;
		}

		_PFROST_H_D_ void printGate(CNF& cnf, OL& poss, OL& negs)
		{
			for (uint32 i = 0; i < poss.size(); i++) {
				if (cnf[poss[i]].molten()) {
					printf("c | ");
					cnf[poss[i]].print();
				}
			}
			for (uint32 i = 0; i < negs.size(); i++) {
				if (cnf[negs[i]].molten()) {
					printf("c | ");
					cnf[negs[i]].print();
				}
			}
		}

		template<class T, class CMP>
		_PFROST_H_D_ bool devIsSorted(T* d, const int& sz, CMP cmp) {
			for (int i = 1; i < sz; i++)
				if (cmp(d[i], d[i - 1])) return false;
			return true;
		}

		template<typename T, typename CMP>
		_PFROST_D_ void cswap(T& x, T& y, CMP cmp)
		{
			bool which = cmp(x, y);
			T ta = which ? x : y;
			T tb = which ? y : x;
			x = ta, y = tb;
		}

		template<typename T, typename CMP>
		_PFROST_D_ void sort2(T& x, T& y, CMP cmp)
		{
			cswap(x, y, cmp);
		}

		template<typename T, typename CMP>
		_PFROST_D_ void sort3(T& x, T& y, T& z, CMP cmp)
		{
			cswap(y, z, cmp);
			cswap(x, z, cmp);
			cswap(x, y, cmp);
		}

		template<typename T, typename CMP>
		_PFROST_D_ void sort4(T* d, CMP cmp)
		{
			cswap(d[0], d[1], cmp);
			cswap(d[2], d[3], cmp);
			cswap(d[0], d[2], cmp);
			cswap(d[1], d[3], cmp);
			cswap(d[1], d[2], cmp);
		}

		template<typename T, typename CMP>
		_PFROST_D_ void sort5(T* d, CMP cmp)
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
		_PFROST_D_ void sort6(T* d, CMP cmp)
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
			if (size == 2) sort2(data[0], data[1], cmp);
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
			devSort(data, size, DEFAULT_CMP<T>());
			assert(devIsSorted(data, size, DEFAULT_CMP<T>()));
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
#define MIN(x, y) (x < y ? x : y)
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
		_PFROST_D_ void devSwap(T& a, T& b) { T c = a; a = b, b = c; }

		template<typename T>
		_PFROST_D_ void warpReduce(T* smem, T& val) {
			if (threadIdx.x < 32) {
				if (blockDim.x >= 64) val += smem[threadIdx.x + 32];
				val += __shfl_down_sync(FULLWARP, val, 16);
				val += __shfl_down_sync(FULLWARP, val, 8);
				val += __shfl_down_sync(FULLWARP, val, 4);
				val += __shfl_down_sync(FULLWARP, val, 2);
				val += __shfl_down_sync(FULLWARP, val, 1);
			}
		}

		template<typename T1, typename T2>
		_PFROST_D_ void warpReduce(T1* smem1, T1& val1, T2* smem2, T2& val2) {
			if (threadIdx.x < 32) {
				if (blockDim.x >= 64) {
					val1 += smem1[threadIdx.x + 32];
					val2 += smem2[threadIdx.x + 32];
				}
				val1 += __shfl_down_sync(FULLWARP, val1, 16);
				val2 += __shfl_down_sync(FULLWARP, val2, 16);
				val1 += __shfl_down_sync(FULLWARP, val1, 8);
				val2 += __shfl_down_sync(FULLWARP, val2, 8);
				val1 += __shfl_down_sync(FULLWARP, val1, 4);
				val2 += __shfl_down_sync(FULLWARP, val2, 4);
				val1 += __shfl_down_sync(FULLWARP, val1, 2);
				val2 += __shfl_down_sync(FULLWARP, val2, 2);
				val1 += __shfl_down_sync(FULLWARP, val1, 1);
				val2 += __shfl_down_sync(FULLWARP, val2, 1);
			}
		}

		template<typename T, typename S>
		_PFROST_D_ void loadShared(T* smem, const T& val, const S& size) {
			smem[threadIdx.x] = (threadIdx.x < size) ? val : 0;
			__syncthreads();
		}

		template<typename T1, typename T2, typename S>
		_PFROST_D_ void loadShared(T1* smem1, const T1& val1, T2* smem2, const T2& val2, const S& size) {
			if (threadIdx.x < size) { smem1[threadIdx.x] = val1, smem2[threadIdx.x] = val2; }
			else { smem1[threadIdx.x] = 0, smem2[threadIdx.x] = 0; }
			__syncthreads();
		}

		template<typename T>
		_PFROST_D_ void sharedReduce(T* smem, T& val) {
			if (blockDim.x >= 512) {
				if (threadIdx.x < 256) smem[threadIdx.x] = val = val + smem[threadIdx.x + 256];
				__syncthreads();
			}
			if (blockDim.x >= 256) {
				if (threadIdx.x < 128) smem[threadIdx.x] = val = val + smem[threadIdx.x + 128];
				__syncthreads();
			}
			if (blockDim.x >= 128) {
				if (threadIdx.x < 64) smem[threadIdx.x] = val = val + smem[threadIdx.x + 64];
				__syncthreads();
			}
		}

		template<typename T1, typename T2>
		_PFROST_D_ void sharedReduce(T1* smem1, T1& val1, T2* smem2, T2& val2) {
			if (blockDim.x >= 512) {
				if (threadIdx.x < 256) {
					smem1[threadIdx.x] = val1 = val1 + smem1[threadIdx.x + 256];
					smem2[threadIdx.x] = val2 = val2 + smem2[threadIdx.x + 256];
				}
				__syncthreads();
			}
			if (blockDim.x >= 256) {
				if (threadIdx.x < 128) {
					smem1[threadIdx.x] = val1 = val1 + smem1[threadIdx.x + 128];
					smem2[threadIdx.x] = val2 = val2 + smem2[threadIdx.x + 128];
				}
				__syncthreads();
			}
			if (blockDim.x >= 128) {
				if (threadIdx.x < 64) {
					smem1[threadIdx.x] = val1 = val1 + smem1[threadIdx.x + 64];
					smem2[threadIdx.x] = val2 = val2 + smem2[threadIdx.x + 64];
				}
				__syncthreads();
			}
		}

		template<typename T>
		_PFROST_D_ void cuVec<T>::insert(const T& val)
		{
			uint32 idx = atomicInc(&sz, cap);
			assert(idx < cap);
			_mem[idx] = val;
		}

		template<typename T>
		_PFROST_D_ void cuVec<T>::push(const T& val) {
			uint32 idx = atomicAggInc(&sz);
			assert(idx < cap);
			_mem[idx] = val;
		}

		template<typename T>
		_PFROST_D_ T* cuVec<T>::jump(const uint32& n) {
			uint32 idx = atomicAdd(&sz, n);
			assert(idx < cap);
			return _mem + idx;
		}

		_PFROST_D_ S_REF* CNF::jump(S_REF& ref, const uint32& nCls, const uint32& nLits) {
			assert(nLits >= nCls);
			S_REF regionSize = (nLits - nCls) + dc_nbuckets * nCls;
			ref = atomicAdd(&_data.size, regionSize);
			assert(ref < _data.cap);
			return cs.jump(nCls);
		}

		_PFROST_D_ uint32 rscore(const uint32& ps, const uint32& ns) {
			return (!ps || !ns) ? ps | ns : ps * ns;
		}

		_PFROST_D_ void calcSig(SCLAUSE& c)
		{
			if (c.size() <= 1) return;
			uint32 sig = 0;
			uint32* lit = c, *cend = c.end();
#pragma unroll
			while (lit != cend) sig |= MAPHASH(*lit++);
			c.set_sig(sig);
		}

		_PFROST_D_ void calcSig(uint32* data, const int& size, uint32& _sig)
		{
			if (size <= 1) return;
			assert(_sig == 0);
			uint32* end = data + size;
#pragma unroll
			while (data != end) _sig |= MAPHASH(*data++);
		}

		_PFROST_D_ void countOrgs(CNF& cnf, OL& list, uint32& orgs)
		{
			assert(!orgs);
#pragma unroll
			for (S_REF* i = list; i != list.end(); i++)
				if (cnf[*i].original()) orgs++;
		}

		_PFROST_D_ bool isTautology(const uint32& a, const uint32& b) { return (a ^ b) == NEG_SIGN; }

		_PFROST_D_ bool isTautology(const uint32& x, SCLAUSE& c1, SCLAUSE& c2)
		{
			assert(x);
			assert(!c1.deleted());
			assert(!c2.deleted());
			assert(c1.size() > 1);
			assert(c2.size() > 1);
			int n1 = c1.size(), n2 = c2.size();
			int it1 = 0, it2 = 0;
			uint32 lit1, lit2, v1, v2;
			while (it1 < n1 && it2 < n2) {
				lit1 = c1[it1], lit2 = c2[it2];
				v1 = ABS(lit1), v2 = ABS(lit2);
				if (v1 == x) it1++;
				else if (v2 == x) it2++;
				else if (isTautology(lit1, lit2)) return true;
				else if (v1 < v2) it1++;
				else if (v2 < v1) it2++;
				else { it1++; it2++; }
			}
			return false;
		}

		_PFROST_D_ bool isTautology(const uint32& x, SCLAUSE& c1, uint32* c2, const int& n2)
		{
			assert(x);
			assert(!c1.deleted());
			assert(c1.size() > 1);
			assert(n2 > 1);
			int n1 = c1.size();
			int it1 = 0, it2 = 0;
			uint32 lit1, lit2, v1, v2;
			while (it1 < n1 && it2 < n2) {
				lit1 = c1[it1], lit2 = c2[it2];
				v1 = ABS(lit1), v2 = ABS(lit2);
				if (v1 == x) it1++;
				else if (v2 == x) it2++;
				else if (isTautology(lit1, lit2)) return true; 
				else if (v1 < v2) it1++;
				else if (v2 < v1) it2++;
				else { it1++; it2++; }
			}
			return false;
		}

		_PFROST_D_ void	saveResolved(uint32*& saved, const uint32& lit)
		{
			*saved++ = lit, *saved++ = 1;
		}

		_PFROST_D_ void	saveResolved(uint32*& saved, SCLAUSE& c, const uint32& x)
		{
			uint32* first = saved, * witness = NULL;
			assert(c.original());
			uint32* k = c, * cend = c.end();
			while (k != cend) {
				const uint32 lit = *k++;
				if (lit == x) {
					witness = saved;
				}
				*saved++ = lit;
			}
			assert(witness >= first);
			if (witness != first)
				devSwap(*first, *witness);
			else
				assert(*witness == *first);
			*saved++ = c.size();
		}

		_PFROST_D_ void countLitsBefore(CNF& cnf, OL& list, uint32& nLitsBefore)
		{
#pragma unroll
			for (S_REF* i = list; i != list.end(); i++) {
				SCLAUSE& c = cnf[*i];
				if (c.original()) nLitsBefore += c.size();
			}
		}

		_PFROST_D_ void merge(const uint32& x, SCLAUSE& c1, SCLAUSE& c2, SCLAUSE* out_c)
		{
			assert(x);
			assert(c1.original());
			assert(c2.original());
			assert(c1.size() > 1);
			assert(c2.size() > 1);
			assert(out_c->empty());
			int n1 = c1.size(), n2 = c2.size();
			int it1 = 0, it2 = 0;
			uint32 lit1, lit2, v1, v2;
			while (it1 < n1 && it2 < n2) {
				lit1 = c1[it1], lit2 = c2[it2];
				v1 = ABS(lit1), v2 = ABS(lit2);
				if (v1 == x) it1++;
				else if (v2 == x) it2++;
				else if (v1 < v2) { it1++; out_c->push(lit1); }
				else if (v2 < v1) { it2++; out_c->push(lit2); }
				else { // repeated literal
					it1++, it2++;
					out_c->push(lit1);
				}
			}
			while (it1 < n1) {
				lit1 = c1[it1];
				if (ABS(lit1) == x) it1++;
				else { it1++; out_c->push(lit1); }
			}
			while (it2 < n2) {
				lit2 = c2[it2];
				if (ABS(lit2) == x) it2++;
				else { it2++; out_c->push(lit2); }
			}
			calcSig(*out_c);
			out_c->set_status(ORIGINAL);
			assert(out_c->isSorted());
			assert(out_c->hasZero() < 0);
		}

		_PFROST_D_ int merge(const uint32& x, SCLAUSE& c1, SCLAUSE& c2, uint32* out_c)
		{
			assert(x);
			assert(c1.original());
			assert(c2.original());
			assert(c1.size() > 1);
			assert(c2.size() > 1);
			int n1 = c1.size(), n2 = c2.size();
			int it1 = 0, it2 = 0;
			uint32 lit1, lit2, v1, v2;
			int len = 0;
			while (it1 < n1 && it2 < n2) {
				lit1 = c1[it1], lit2 = c2[it2];
				v1 = ABS(lit1), v2 = ABS(lit2);
				if (v1 == x) it1++;
				else if (v2 == x) it2++;
				else if (isTautology(lit1, lit2)) return 0;
				else if (v1 < v2) { it1++; out_c[len++] = lit1; }
				else if (v2 < v1) { it2++; out_c[len++] = lit2; }
				else { // repeated literal
					it1++, it2++;
					out_c[len++] = lit1;
				}
			}
			while (it1 < n1) {
				lit1 = c1[it1];
				if (ABS(lit1) == x) it1++;
				else { it1++; out_c[len++] = lit1; }
			}
			while (it2 < n2) {
				lit2 = c2[it2];
				if (ABS(lit2) == x) it2++;
				else { it2++; out_c[len++] = lit2; }
			}
			return len;
		}

		_PFROST_D_ int merge_ere(const uint32& x, SCLAUSE& c1, SCLAUSE& c2, uint32* out_c)
		{
			assert(x);
			assert(!c1.deleted());
			assert(!c2.deleted());
			assert(c1.size() > 1);
			assert(c2.size() > 1);
			int n1 = c1.size(), n2 = c2.size();
			int it1 = 0, it2 = 0;
			uint32 lit1, lit2, v1, v2;
			int len = 0;
			while (it1 < n1 && it2 < n2) {
				lit1 = c1[it1], lit2 = c2[it2];
				v1 = ABS(lit1), v2 = ABS(lit2);
				if (v1 == x) it1++;
				else if (v2 == x) it2++;
				else if (isTautology(lit1, lit2)) return 0;  
				else if (v1 < v2) { it1++; out_c[len++] = lit1; }
				else if (v2 < v1) { it2++; out_c[len++] = lit2; }
				else { // repeated literal
					it1++, it2++;
					out_c[len++] = lit1;
				}
			}
			while (it1 < n1) {
				lit1 = c1[it1];
				if (ABS(lit1) == x) it1++;
				else { it1++; out_c[len++] = lit1; }
			}
			while (it2 < n2) {
				lit2 = c2[it2];
				if (ABS(lit2) == x) it2++;
				else { it2++; out_c[len++] = lit2; }
			}
			assert(len <= SH_MAX_ERE_OUT);
			return len;
		}

		_PFROST_D_ int merge(const uint32& x, uint32* c1, const int& n1, SCLAUSE& c2, uint32* out_c)
		{
			assert(x);
			assert(n1 > 1);
			assert(c2.original());
			assert(c2.size() > 1);
			int n2 = c2.size();
			int it1 = 0, it2 = 0;
			uint32 lit1, lit2, v1, v2;
			int len = 0;
			while (it1 < n1 && it2 < n2) {
				lit1 = c1[it1], lit2 = c2[it2];
				v1 = ABS(lit1), v2 = ABS(lit2);
				if (v1 == x) it1++;
				else if (v2 == x) it2++;
				else if (isTautology(lit1, lit2)) return 0;  
				else if (v1 < v2) { it1++; out_c[len++] = lit1; }
				else if (v2 < v1) { it2++; out_c[len++] = lit2; }
				else { // repeated literal
					it1++, it2++;
					out_c[len++] = lit1;
				}
			}
			while (it1 < n1) {
				lit1 = c1[it1];
				if (ABS(lit1) == x) it1++;
				else { it1++; out_c[len++] = lit1; }
			}
			while (it2 < n2) {
				lit2 = c2[it2];
				if (ABS(lit2) == x) it2++;
				else { it2++; out_c[len++] = lit2; }
			}
			return len;
		}

		_PFROST_D_ bool isEqual(SCLAUSE& c1, uint32* c2, const int& size)
		{
			assert(!c1.deleted());
			assert(c1.size() > 1);
			assert(size > 1);
#pragma unroll
			for (int it = 0; it < size; it++) if (c1[it] != c2[it]) return false;
			return true;
		}

		_PFROST_D_ bool sub(const uint32& A, const uint32& B) { return !(A & ~B); }

		_PFROST_D_ bool selfsub(const uint32& A, const uint32& B)
		{
			uint32 B_tmp = B | ((B & 0xAAAAAAAAUL) >> 1) | ((B & 0x55555555UL) << 1);
			return !(A & ~B_tmp);
		}

		_PFROST_D_ void reduceOL(CNF& cnf, OL& ol)
		{
			if (ol.size() == 0) return;
			S_REF* i, * j, * rend = ol.end();
			for (i = ol, j = ol; i != rend; i++)
				if (!cnf[*i].deleted()) *j++ = *i;
			ol._shrink(i - j);
		}

		_PFROST_D_ void forward_equ(CNF& cnf, OT& ot, uint32* m_c, const int& m_len, const CL_ST& type)
		{
			assert(m_len > 1);
			assert(type != DELETED);
			uint32 best = *m_c, m_sig = MAPHASH(best);
			assert(best > 1);
			int minsize = ot[best].size();
			for (int k = 1; k < m_len; k++) {
				int lsize = ot[m_c[k]].size();
				if (lsize < minsize) minsize = lsize, best = m_c[k];
				m_sig |= MAPHASH(m_c[k]);
			}
			OL& minList = ot[best];
			for (S_REF* i = threadIdx.x + minList; i < minList.end(); i += blockDim.x) {
				SCLAUSE& c = cnf[*i];
				CL_ST st = c.status();
				if (m_len == c.size() && ((st & LEARNT) || (st & type)) &&
					sub(m_sig, c.sig()) && isEqual(c, m_c, m_len)) {
					c.markDeleted();
					break;
				}
			}
		}

		_PFROST_D_ void blocked_x(const uint32& x, CNF& cnf, OL& poss, OL& negs, uint32* sh_c)
		{
#pragma unroll
			for (S_REF* i = negs; i != negs.end(); i++) { // start with negs
				SCLAUSE& c_me = cnf[*i];
				if (c_me.deleted() || c_me.learnt()) continue;
				bool allTautology = true;
				int c_size = c_me.size();
				if (c_size <= SH_MAX_BCE_IN) { // use shared memory 
					c_me.shareTo(sh_c);
#pragma unroll
					for (S_REF* j = poss; j != poss.end(); j++) { // block with poss
						SCLAUSE& c_other = cnf[*j];
						if (c_other.deleted() || c_other.learnt()) continue;
						if (!isTautology(x, c_other, sh_c, c_size)) { allTautology = false; break; }
					}
				}
				else { // use global memory
#pragma unroll
					for (S_REF* j = poss; j != poss.end(); j++) { // block with poss
						SCLAUSE& c_other = cnf[*j];
						if (c_other.deleted() || c_other.learnt()) continue;
						if (!isTautology(x, c_me, c_other)) { allTautology = false; break; }
					}
				}
				if (allTautology) assert(c_me.size() > 1), c_me.markDeleted(); // all clauses but a unit can be blocked
			}
		}

	} // namespace sigma
} // namespace pfrost

#endif