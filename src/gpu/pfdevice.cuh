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
#include "pfsimptypes.h"

namespace pFROST {

	namespace SIGmA {

		#define VE_DBG 0 // set to serialize BVE
        #define IS_TAUTOLOGY(x,y)   (((x) ^ (y)) == NEG_SIGN)

		// [0] = hse_limit, [1] = bce_limit,
		// [2] = ere_limit, [3] = xor_max_arity
		// [4] = ve_clause_limit, [5] = opts.ve_lbound_en
		__constant__ uint32 dc_limits[NLIMITS];
		// number of buckets to store 'SCLAUSE'
		__constant__ int dc_nbuckets = sizeof(SCLAUSE) / sizeof(uint32);

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
				printf("c  "), cnf[poss[i]].print();
			}
			for (uint32 i = 0; i < negs.size(); i++) {
				printf("c  "), cnf[negs[i]].print();
			}
		}

		_PFROST_H_D_ void pClauseSet(CNF& cnf, OL& ol)
		{
			for (uint32 i = 0; i < ol.size(); i++) {
				SCLAUSE& c = cnf[ol[i]];
				printf("c  c(%lld)->", uint64(ol[i])), c.print();
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
					printf("c  ");
					cnf[poss[i]].print();
				}
			}
			for (uint32 i = 0; i < negs.size(); i++) {
				if (cnf[negs[i]].molten()) {
					printf("c  ");
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
			const bool which = cmp(x, y);
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
			const uint32 idx = atomicInc(&sz, cap);
			assert(idx < cap);
			_mem[idx] = val;
		}

		template<typename T>
		_PFROST_D_ void cuVec<T>::push(const T& val) {
			const uint32 idx = atomicAggInc(&sz);
			assert(idx < cap);
			_mem[idx] = val;
		}

		template<typename T>
		_PFROST_D_ T* cuVec<T>::jump(const uint32& n) {
			const uint32 idx = atomicAdd(&sz, n);
			assert(idx < cap);
			return _mem + idx;
		}

		_PFROST_D_ S_REF* CNF::jump(S_REF& ref, const uint32& nCls, const uint32& nLits) {
			assert(nLits >= nCls);
			const S_REF regionSize = nLits + dc_nbuckets * nCls;
			ref = atomicAdd(&_data.size, regionSize);
			assert(ref < _data.cap);
			return _refs.jump(nCls);
		}

		_PFROST_D_ void calcSig(SCLAUSE& c)
		{
			if (c.size() <= 1) return;
			uint32 sig = 0;
#pragma unroll
			forall_clause(c, lit) { sig |= MAPHASH(*lit); }
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

		_PFROST_D_ void freeze_binaries(CNF& cnf, OL& list)
		{
#pragma unroll
			forall_occurs(list, i) {
				SCLAUSE& c = cnf[*i];
				if (c.size() == 2) c.freeze();
			}
		}

		_PFROST_D_ bool isTautology(const uint32& x, SCLAUSE& c1, SCLAUSE& c2)
		{
			assert(x);
			assert(!c1.deleted());
			assert(!c2.deleted());
			assert(c1.size() > 1);
			assert(c2.size() > 1);
			uint32 *n1 = c1.end(), *n2 = c2.end();
			uint32 *lit1 = c1, *lit2 = c2, v1, v2;
			while (lit1 != n1 && lit2 != n2) {
				v1 = ABS(*lit1), v2 = ABS(*lit2);
				if (v1 == x) lit1++;
				else if (v2 == x) lit2++;
				else if (IS_TAUTOLOGY(*lit1, *lit2)) return true;
				else if (v1 < v2) lit1++;
				else if (v2 < v1) lit2++;
				else { lit1++, lit2++; }
			}
			return false;
		}

		_PFROST_D_ bool isTautology(const uint32& x, SCLAUSE& c1, uint32* c2, const int& n2)
		{
			assert(x);
			assert(!c1.deleted());
			assert(c1.size() > 1);
			assert(n2 > 1);
			const int n1 = c1.size();
			int it1 = 0, it2 = 0;
			uint32 lit1, lit2, v1, v2;
			while (it1 < n1 && it2 < n2) {
				lit1 = c1[it1], lit2 = c2[it2];
				v1 = ABS(lit1), v2 = ABS(lit2);
				if (v1 == x) it1++;
				else if (v2 == x) it2++;
				else if (IS_TAUTOLOGY(lit1, lit2)) return true;
				else if (v1 < v2) it1++;
				else if (v2 < v1) it2++;
				else { it1++; it2++; }
			}
			return false;
		}

		_PFROST_D_ int merge(const uint32& x, SCLAUSE& c1, SCLAUSE& c2)
		{
			assert(x);
			assert(c1.original());
			assert(c2.original());
			assert(c1.size() > 1);
			assert(c2.size() > 1);
			const int n1 = c1.size(), n2 = c2.size();
			int it1 = 0, it2 = 0;
			int len = n1 + n2 - 2;
			while (it1 < n1 && it2 < n2) {
				const uint32 lit1 = c1[it1], lit2 = c2[it2];
				const uint32 v1 = ABS(lit1), v2 = ABS(lit2);
				if (v1 == x) it1++;
				else if (v2 == x) it2++;
				else if (IS_TAUTOLOGY(lit1, lit2)) return 0;
				else if (v1 < v2) it1++;
				else if (v2 < v1) it2++;
				else { // repeated literal
					it1++, it2++;
					assert(len > 1);
					len--;
				}
			}
			assert(len > 0);
			return len;
		}

		_PFROST_D_ int merge(const uint32& x, SCLAUSE& c1, SCLAUSE& c2, uint32& unit)
		{
			assert(x);
			assert(c1.original());
			assert(c2.original());
			assert(c1.size() > 1);
			assert(c2.size() > 1);
			const int n1 = c1.size(), n2 = c2.size();
			int it1 = 0, it2 = 0;
			int len = n1 + n2 - 2;
			unit = 0;
			while (it1 < n1 && it2 < n2) {
				const uint32 lit1 = c1[it1], lit2 = c2[it2];
				const uint32 v1 = ABS(lit1), v2 = ABS(lit2);
				if (v1 == x) it1++;
				else if (v2 == x) it2++;
				else if (IS_TAUTOLOGY(lit1, lit2)) return 0;
				else if (v1 < v2) it1++;
				else if (v2 < v1) it2++;
				else { // repeated literal
					it1++, it2++;
					assert(len > 1);
					len--;
					if (len == 1) unit = lit1;
				}
			}
			assert(len > 0);
			return len;
		}

		_PFROST_D_ void merge(const uint32& x, SCLAUSE& c1, SCLAUSE& c2, SCLAUSE* out_c)
		{
			assert(x);
			assert(c1.original());
			assert(c2.original());
			assert(c1.size() > 1);
			assert(c2.size() > 1);
			assert(out_c->empty());
			const int n1 = c1.size(), n2 = c2.size();
			int it1 = 0, it2 = 0;
			while (it1 < n1 && it2 < n2) {
				const uint32 lit1 = c1[it1], lit2 = c2[it2];
				const uint32 v1 = ABS(lit1), v2 = ABS(lit2);
				if (v1 == x) it1++;
				else if (v2 == x) it2++;
				else if (v1 < v2) { it1++, out_c->push(lit1); }
				else if (v2 < v1) { it2++, out_c->push(lit2); }
				else { // repeated literal
					it1++, it2++;
					out_c->push(lit1);
				}
			}
			while (it1 < n1) {
				const uint32 lit1 = c1[it1++];
				if (NEQUAL(ABS(lit1), x)) out_c->push(lit1);
			}
			while (it2 < n2) {
				const uint32 lit2 = c2[it2++];
				if (NEQUAL(ABS(lit2), x)) out_c->push(lit2);
			}
			calcSig(*out_c);
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
			const int n1 = c1.size(), n2 = c2.size();
			int it1 = 0, it2 = 0;
			int len = 0;
			while (it1 < n1 && it2 < n2) {
				const uint32 lit1 = c1[it1], lit2 = c2[it2];
				const uint32 v1 = ABS(lit1), v2 = ABS(lit2);
				if (v1 == x) it1++;
				else if (v2 == x) it2++;
				else if (IS_TAUTOLOGY(lit1, lit2)) return 0;
				else if (v1 < v2) { it1++; out_c[len++] = lit1; }
				else if (v2 < v1) { it2++; out_c[len++] = lit2; }
				else { // repeated literal
					it1++, it2++;
					out_c[len++] = lit1;
				}
			}
			while (it1 < n1) {
				const uint32 lit1 = c1[it1++];
				if (NEQUAL(ABS(lit1), x)) out_c[len++] = lit1;
			}
			while (it2 < n2) {
				const uint32 lit2 = c2[it2++];
				if (NEQUAL(ABS(lit2), x)) out_c[len++] = lit2;
			}
			return len;
		}

		_PFROST_D_ int merge(const uint32& x, uint32* c1, const int& n1, SCLAUSE& c2, uint32* out_c)
		{
			assert(x);
			assert(n1 > 1);
			assert(c2.original());
			assert(c2.size() > 1);
			const int n2 = c2.size();
			int it1 = 0, it2 = 0;
			int len = 0;
			while (it1 < n1 && it2 < n2) {
				const uint32 lit1 = c1[it1], lit2 = c2[it2];
				const uint32 v1 = ABS(lit1), v2 = ABS(lit2);
				if (v1 == x) it1++;
				else if (v2 == x) it2++;
				else if (IS_TAUTOLOGY(lit1, lit2)) return 0;
				else if (v1 < v2) { it1++; out_c[len++] = lit1; }
				else if (v2 < v1) { it2++; out_c[len++] = lit2; }
				else { // repeated literal
					it1++, it2++;
					out_c[len++] = lit1;
				}
			}
			while (it1 < n1) {
				const uint32 lit1 = c1[it1++];
				if (NEQUAL(ABS(lit1), x)) out_c[len++] = lit1;
			}
			while (it2 < n2) {
				const uint32 lit2 = c2[it2++];
				if (NEQUAL(ABS(lit2), x)) out_c[len++] = lit2;
			}
			return len;
		}

		_PFROST_D_ void countOrgs(CNF& cnf, OL& list, uint32& orgs)
		{
			assert(!orgs);
#pragma unroll
			forall_occurs(list, i) {
				if (cnf[*i].original()) orgs++;
			}
		}

		_PFROST_D_ void countOrgs(CNF& cnf, OL& list, uint32& nClsBefore, uint32& nLitsBefore)
		{
			assert(!nClsBefore);
			assert(!nLitsBefore);
#pragma unroll
			forall_occurs(list, i) {
				const SCLAUSE& c = cnf[*i];
				if (c.original()) {
					nClsBefore++;
					nLitsBefore += c.size();
				}
			}
		}

		_PFROST_D_ void countLitsBefore(CNF& cnf, OL& list, uint32& nLitsBefore)
		{
#pragma unroll
			forall_occurs(list, i) {
				const SCLAUSE& c = cnf[*i];
				if (c.original()) nLitsBefore += c.size();
			}
		}

		_PFROST_D_ bool countResolvents(const uint32& x, const uint32& nClsBefore, CNF& cnf, OL& me, OL& other, uint32& nAddedCls, uint32& nAddedLits)
		{
			assert(x);
			assert(nClsBefore);
			assert(!nAddedLits);
			const int rlimit = dc_limits[4];
			forall_occurs(me, i) {
				SCLAUSE& ci = cnf[*i];
				if (ci.learnt()) continue;
				forall_occurs(other, j) {
					SCLAUSE& cj = cnf[*j];
					if (cj.learnt()) continue;
					const int rsize = merge(x, ci, cj);
					if (rsize > 1) {
						if (++nAddedCls > nClsBefore || (rlimit && rsize > rlimit)) return true;
						nAddedLits += rsize;
					}
				}
			}
			if (nAddedCls > ADDEDCLS_MAX || nAddedLits > ADDEDLITS_MAX) return true;
			if (dc_limits[5]) {
				uint32 nLitsBefore = 0;
				countLitsBefore(cnf, me, nLitsBefore);
				countLitsBefore(cnf, other, nLitsBefore);
				if (nAddedLits > nLitsBefore) return true;
			}
			return false;
		}

		_PFROST_D_ bool countSubstituted(const uint32& x, const uint32& nClsBefore, CNF& cnf, OL& me, OL& other, uint32& nAddedCls, uint32& nAddedLits)
		{
			assert(x);
			assert(!nAddedCls);
			assert(!nAddedLits);
			const int rlimit = dc_limits[4];
			forall_occurs(me, i) {
				SCLAUSE& ci = cnf[*i];
				if (ci.learnt()) continue;
				const bool ci_m = ci.molten();
				forall_occurs(other, j) {
					SCLAUSE& cj = cnf[*j];
					if (cj.learnt()) continue;
					if (NEQUAL(ci_m, cj.molten())) {
						const int rsize = merge(x, ci, cj);
						if (rsize > 1) {
							if (++nAddedCls > nClsBefore || (rlimit && rsize > rlimit)) return true;
							nAddedLits += rsize;
						}
					}
				}
			}
			if (nAddedCls > ADDEDCLS_MAX || nAddedLits > ADDEDLITS_MAX) return true;
			if (dc_limits[5]) {
				uint32 nLitsBefore = 0;
				countLitsBefore(cnf, me, nLitsBefore);
				countLitsBefore(cnf, other, nLitsBefore);
				if (nAddedLits > nLitsBefore) return true;
			}
			return false;
		}

		_PFROST_D_ void	saveWitness(uint32*& saved, const uint32* vorg, const uint32& witness)
		{
			assert(witness > 1);
			assert(vorg);
			const uint32 orgWitness = V2DEC(vorg[ABS(witness)], SIGN(witness));
			assert(orgWitness > 1);
			*saved++ = orgWitness;
			*saved++ = 1;
		}

		_PFROST_D_ void	saveLiteral(uint32*& saved, const uint32* vorg, const uint32& lit)
		{
			assert(lit > 1);
			assert(vorg);
			const uint32 orgLit = V2DEC(vorg[ABS(lit)], SIGN(lit));
			assert(orgLit > 1);
			*saved++ = orgLit;
		}

		_PFROST_D_ void	saveClause(uint32*& saved, SCLAUSE& c, const uint32* vorg, const uint32& witlit)
		{
			uint32* first = saved, * witness = NULL;
			assert(c.original());
			forall_clause(c, k) {
				const uint32 lit = *k;
				if (lit == witlit) {
					witness = saved;
				}
				saveLiteral(saved, vorg, lit);
			}
			assert(witness >= first);
			if (witness != first)
				devSwap(*first, *witness);
			else
				assert(*witness == *first);
			*saved++ = c.size();
		}

		_PFROST_D_ void saveResolved(const uint32& p, const uint32* vorg, CNF& cnf, OL& poss, OL& negs, cuVecU* resolved)
		{
			const uint32 n = NEG(p);
			const bool which = poss.size() > negs.size();
			if (which) {
				uint32 nsCls = 0, nsLits = 0;
				countOrgs(cnf, negs, nsCls, nsLits);
				uint32* saved = resolved->jump(nsCls + nsLits + 2);
#if VE_DBG
				printf("c  saving witness(%d) of length %d at position %d\n",
					ABS(p), nsCls + nsLits + 2, uint32(saved - resolved->data()));
#endif
#pragma unroll
				forall_occurs(negs, i) {
					SCLAUSE& c = cnf[*i];
					if (c.original()) saveClause(saved, c, vorg, n);
				}
				saveWitness(saved, vorg, p);
			}
			else {
				uint32 psCls = 0, psLits = 0;
				countOrgs(cnf, poss, psCls, psLits);
				uint32* saved = resolved->jump(psCls + psLits + 2);
#if VE_DBG
				printf("c  saving witness(%d) of length %d at position %d\n",
					ABS(p), psCls + psLits + 2, uint32(saved - resolved->data()));
#endif
#pragma unroll
				forall_occurs(poss, i) {
					SCLAUSE& c = cnf[*i];
					if (c.original()) saveClause(saved, c, vorg, p);
				}
				saveWitness(saved, vorg, n);
			}
		}

		_PFROST_D_ void saveResolved(const uint32& p, const uint32* vorg,
			const uint32& pOrgs, const uint32& nOrgs,
			CNF& cnf, OL& poss, OL& negs, cuVecU* resolved)
		{
			const uint32 n = NEG(p);
			const bool which = pOrgs > nOrgs;
			if (which) {
				uint32 nsLits = 0;
				countLitsBefore(cnf, negs, nsLits);
				uint32* saved = resolved->jump(nOrgs + nsLits + 2);
#if VE_DBG
				printf("c  saving witness(%d) of length %d at position %d\n",
					ABS(p), nOrgs + nsLits + 2, uint32(saved - resolved->data()));
#endif
#pragma unroll
				forall_occurs(negs, i) {
					SCLAUSE& c = cnf[*i];
					if (c.original()) saveClause(saved, c, vorg, n);
				}
				saveWitness(saved, vorg, p);
			}
			else {
				uint32 psLits = 0;
				countLitsBefore(cnf, poss, psLits);
				uint32* saved = resolved->jump(pOrgs + psLits + 2);
#if VE_DBG
				printf("c  saving witness(%d) of length %d at position %d\n",
					ABS(p), pOrgs + psLits + 2, uint32(saved - resolved->data()));
#endif
#pragma unroll
				forall_occurs(poss, i) {
					SCLAUSE& c = cnf[*i];
					if (c.original()) saveClause(saved, c, vorg, p);
				}
				saveWitness(saved, vorg, n);
			}
		}

		_PFROST_D_ void toblivion(const uint32& p, const uint32* vorg,
			const uint32& pOrgs, const uint32& nOrgs, 
			CNF& cnf, OL& poss, OL& negs, cuVecU* resolved)
		{
			const uint32 n = NEG(p);
			const bool which = pOrgs > nOrgs;
			if (which) {
				uint32 nsLits = 0;
				countLitsBefore(cnf, negs, nsLits);
				uint32* saved = resolved->jump(nOrgs + nsLits + 2);
#if VE_DBG
				printf("c  saving witness(%d) of length %d at position %d\n",
					ABS(p), nOrgs + nsLits + 2, uint32(saved - resolved->data()));
#endif
#pragma unroll
				forall_occurs(negs, i) {
					SCLAUSE& c = cnf[*i];
					if (c.original()) saveClause(saved, c, vorg, n);
					c.markDeleted();
				}
				saveWitness(saved, vorg, p);
			}
			else {
				uint32 psLits = 0;
				countLitsBefore(cnf, poss, psLits);
				uint32* saved = resolved->jump(pOrgs + psLits + 2);
#if VE_DBG
				printf("c  saving witness(%d) of length %d at position %d\n",
					ABS(p), pOrgs + psLits + 2, uint32(saved - resolved->data()));
#endif
#pragma unroll
				forall_occurs(poss, i) {
					SCLAUSE& c = cnf[*i];
					if (c.original()) saveClause(saved, c, vorg, p);
					c.markDeleted();
				}
				saveWitness(saved, vorg, n);
			}
			OL& other = which ? poss : negs;
#pragma unroll
			forall_occurs(other, i) cnf[*i].markDeleted();
			poss.clear(true), negs.clear(true);
		}

		_PFROST_D_ void toblivion(CNF& cnf, OL& poss, OL& negs)
		{
#pragma unroll
			forall_occurs(poss, i) cnf[*i].markDeleted();
#pragma unroll
			forall_occurs(negs, i) cnf[*i].markDeleted();
			poss.clear(true), negs.clear(true);
		}

		_PFROST_D_ bool isEqual(SCLAUSE& c1, uint32* c2, const int& size)
		{
			assert(!c1.deleted());
			assert(c1.size() > 1);
			assert(size > 1);
#pragma unroll
			for (int it = 0; it < size; it++)
				if (NEQUAL(c1[it], c2[it]))
					return false;
			return true;
		}

		_PFROST_D_ bool sub(const uint32& A, const uint32& B) { return !(A & ~B); }

		_PFROST_D_ bool selfsub(const uint32& A, const uint32& B)
		{
			uint32 B_tmp = B | ((B & 0xAAAAAAAAUL) >> 1) | ((B & 0x55555555UL) << 1);
			return !(A & ~B_tmp);
		}

		_PFROST_D_ void reduceOL(const CNF& cnf, OL& ol)
		{
			if (ol.empty()) return;
			S_REF *j = ol;
			forall_occurs(ol, i) {
				const S_REF ref = *i;
				if (!cnf[ref].deleted())
					*j++ = ref;
			}
			ol.resize(j - ol);
		}

		_PFROST_D_ void blocked_x(const uint32& x, const uint32* vorg, CNF& cnf, OL& poss, OL& negs, cuVecU* resolved, uint32* sh_c)
		{
#pragma unroll
			forall_occurs(negs, i) { // start with negs
				SCLAUSE& ci = cnf[*i];
				if (ci.deleted() || ci.learnt()) continue;
				bool allTautology = true;
				int c_size = ci.size();
				if (c_size <= SH_MAX_BCE_IN) { // use shared memory 
					ci.shareTo(sh_c);
#pragma unroll
					forall_occurs(poss, j) { // block with poss
						SCLAUSE& cj = cnf[*j];
						if (cj.deleted() || cj.learnt()) continue;
						if (!isTautology(x, cj, sh_c, c_size)) { allTautology = false; break; }
					}
				}
				else { // use global memory
#pragma unroll
					forall_occurs(poss, j) { // block with poss
						SCLAUSE& cj = cnf[*j];
						if (cj.deleted() || cj.learnt()) continue;
						if (!isTautology(x, ci, cj)) { allTautology = false; break; }
					}
				}
				if (allTautology) {
					assert(ci.size() > 1);
					uint32* saved = resolved->jump(ci.size() + 1);
					saveClause(saved, ci, vorg, NEG(V2L(x)));
					ci.markDeleted();
				}
			}
		}

	} // namespace sigma
} // namespace pfrost

#endif