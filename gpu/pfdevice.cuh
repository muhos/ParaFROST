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
#include "pfsimp.h"

namespace pFROST {

	namespace SIGmA {

#define SS_DBG 0
#define VE_DBG 0

		_PFROST_H_D_ void pLit(const uint32& l) { printf("%d", ISNEG(l) ? -int(ABS(l)) : ABS(l)); }

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
				printf("c | "), c.print();
			}
		}

#if SS_DBG
		_PFROST_D_ void DBG_updateOL(const uint32& x, CNF& cnf, OL& ol)
		{
			if (ol.size() == 0) return;
			int idx = 0, ol_sz = ol.size();
			while (idx < ol_sz) {
				SCLAUSE& c = cnf[ol[idx]];
				if (c.deleted()) ol[idx] = ol[--ol_sz];
				else if (c.molten()) { c.freeze(); ol[idx] = ol[--ol_sz]; }
				else idx++;
			}
			if (ol_sz != ol.size()) { printf("c | Updated list:\n"); pClauseSet(cnf, ol); printf("c | == End of list ==\n"); }
			ol.resize(ol_sz);
	}
#endif

		_PFROST_D_ bool isCoherent(const uint32& x, SCLAUSE& g_c, uint32* sh_c, const int& size)
		{
			uint32 it = 0;
			if (g_c.size() != size) {
				printf("c | WARNING - memory incoherency for var(%d) -> clause size not equal\n", ABS(x));
				printf("c | global clause: "), g_c.print();
				printf("c | shared clause: "), pSharedClause(sh_c, size);
				return false;
			}
			while (it < size) {
				if (sh_c[it] != g_c[it]) {
					printf("c | WARNING - memory incoherency for var(%d) -> clauses not equal\n", ABS(x));
					printf("c | global clause: "), g_c.print();
					printf("c | shared clause: "), pSharedClause(sh_c, size);
					return false;
				}
				else it++;
			}
			return true;
		}

		template<typename T>
		_PFROST_D_ void warpReduce(T* smem, T& val) {
			if (threadIdx.x < warpSize) {
				if (blockDim.x >= 64) val += smem[threadIdx.x + warpSize];
				val += __shfl_xor_sync(FULLWARP, val, 16, warpSize);
				val += __shfl_xor_sync(FULLWARP, val, 8, warpSize);
				val += __shfl_xor_sync(FULLWARP, val, 4, warpSize);
				val += __shfl_xor_sync(FULLWARP, val, 2, warpSize);
				val += __shfl_xor_sync(FULLWARP, val, 1, warpSize);
			}
		}

		template<typename T1, typename T2>
		_PFROST_D_ void warpReduce(T1* smem1, T1& val1, T2* smem2, T2& val2) {
			if (threadIdx.x < warpSize) {
				if (blockDim.x >= 64) {
					val1 += smem1[threadIdx.x + warpSize];
					val2 += smem2[threadIdx.x + warpSize];
				}
				val1 += __shfl_xor_sync(FULLWARP, val1, 16, warpSize);
				val2 += __shfl_xor_sync(FULLWARP, val2, 16, warpSize);
				val1 += __shfl_xor_sync(FULLWARP, val1, 8, warpSize);
				val2 += __shfl_xor_sync(FULLWARP, val2, 8, warpSize);
				val1 += __shfl_xor_sync(FULLWARP, val1, 4, warpSize);
				val2 += __shfl_xor_sync(FULLWARP, val2, 4, warpSize);
				val1 += __shfl_xor_sync(FULLWARP, val1, 2, warpSize);
				val2 += __shfl_xor_sync(FULLWARP, val2, 2, warpSize);
				val1 += __shfl_xor_sync(FULLWARP, val1, 1, warpSize);
				val2 += __shfl_xor_sync(FULLWARP, val2, 1, warpSize);
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
		_PFROST_D_ void cuVec<T>::pop(void) { uint32 idx = atomicSub(&sz, 1); assert(idx < UINT32_MAX); }

		template<typename T>
		_PFROST_D_ void cuVec<T>::shrink(const uint32& n) { uint32 idx = atomicSub(&sz, n); assert(idx < UINT32_MAX); }

		template<typename T>
		_PFROST_D_ T* cuVec<T>::jump(const uint32& n) {
			uint32 idx = atomicAdd(&sz, n);
			assert(idx < cap);
			return _mem + idx;
		}

		_PFROST_D_ S_REF* CNF::jump(S_REF& ref, const uint32& nCls, const uint32& nLits) {
			assert(nLits >= nCls);
			S_REF regionSize = (nLits - nCls) + sizeof(S_REF) * nCls;
			ref = atomicAdd(&_data.size, regionSize);
			assert(ref < _data.cap);
			return cs.jump(nCls);
		}

		_PFROST_D_ void devSort(uint32* data, const int& size)
		{
			if (size == 2) {
				uint32 val0 = *data, val1 = *(data + 1);
				if (val0 > val1) {
					*data = val1;
					*(data + 1) = val0;
				}
			}
			else {
				int i, j;
#pragma unroll
				for (i = 1; i < size; i++) {
					uint32 tmp = data[i];
					for (j = i; j > 0 && data[j - 1] > tmp; j--) data[j] = data[j - 1];
					data[j] = tmp;
				}
			}
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
				else if ((lit1 ^ lit2) == NEG_SIGN) return true; // tautology detected ==> abort
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
				else if ((lit1 ^ lit2) == NEG_SIGN) return true; // tautology detected ==> abort
				else if (v1 < v2) it1++;
				else if (v2 < v1) it2++;
				else { it1++; it2++; }
			}
			return false;
		}

		_PFROST_D_ bool merge(const uint32& x, SCLAUSE& c1, SCLAUSE& c2, SCLAUSE& out_c)
		{
			assert(x);
			assert(!c1.deleted());
			assert(!c2.deleted());
			assert(c1.size() > 1);
			assert(c2.size() > 1);
			out_c.clear();
			int n1 = c1.size(), n2 = c2.size();
			int it1 = 0, it2 = 0;
			uint32 lit1, lit2, v1, v2;
			while (it1 < n1 && it2 < n2) {
				lit1 = c1[it1], lit2 = c2[it2];
				v1 = ABS(lit1), v2 = ABS(lit2);
				if (v1 == x) it1++;
				else if (v2 == x) it2++;
				else if ((lit1 ^ lit2) == NEG_SIGN) { out_c.clear(); return false; }  // tautology 
				else if (v1 < v2) { it1++; out_c.push(lit1); }
				else if (v2 < v1) { it2++; out_c.push(lit2); }
				else { // repeated literal
					it1++, it2++;
					out_c.push(lit1);
				}
			}
			while (it1 < n1) {
				lit1 = c1[it1];
				if (ABS(lit1) == x) it1++;
				else { it1++; out_c.push(lit1); }
			}
			while (it2 < n2) {
				lit2 = c2[it2];
				if (ABS(lit2) == x) it2++;
				else { it2++; out_c.push(lit2); }
			}
			calcSig(out_c);
			out_c.set_status(ORIGINAL);
			assert(out_c.isSorted());
			assert(out_c.hasZero() < 0);
			return true;
		}

		_PFROST_D_ int merge(const uint32& x, SCLAUSE& c1, SCLAUSE& c2, uint32* out_c)
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
				else if ((lit1 ^ lit2) == NEG_SIGN) return 0;  // tautology 
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

		_PFROST_D_ int merge(const uint32& x, uint32* c1, const int& n1, SCLAUSE& c2, uint32* out_c)
		{
			assert(x);
			assert(n1 > 1);
			assert(!c2.deleted());
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
				else if ((lit1 ^ lit2) == NEG_SIGN) return 0;  // tautology 
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

		_PFROST_D_ bool isTautology_and(SCLAUSE& org, uint32* defs, const int& nDefs)
		{
			assert(!org.deleted());
			assert(org.size() > 1);
			assert(nDefs > 1);
#pragma unroll
			for (int i = 0; i < nDefs; i++) {
				if (org.has(defs[i])) return true;
			}
			return false;
		}

		_PFROST_D_ bool isTautology_or(SCLAUSE& org, uint32* defs, const int& nDefs)
		{
			assert(!org.deleted());
			assert(org.size() > 1);
			assert(nDefs > 1);
#pragma unroll
			for (int i = 0; i < nDefs; i++) {
				if (org.has(FLIP(defs[i]))) return true;
			}
			return false;
		}

		_PFROST_D_ bool clause_extend_and(const uint32& neg_x, SCLAUSE& org, uint32* defs, const int& nDefs, SCLAUSE& out_c)
		{
			assert(neg_x);
			assert(!org.deleted());
			assert(org.size() > 1);
			assert(nDefs > 1);
			out_c.clear();
#pragma unroll
			for (int i = 0; i < nDefs; i++)
				if (org.has(defs[i])) return false; // tautology
			// extend 
#pragma unroll
			for (int i = 0; i < org.size(); i++) {
				uint32 lit = org[i];
				if (lit == neg_x) {
#pragma unroll
					for (int k = 0; k < nDefs; k++) out_c.push(FLIP(defs[k]));
				}
				else out_c.push(lit);
			}
			// attach 
			devSort(out_c, out_c.size());
			out_c.filter();
			out_c.set_status(ORIGINAL);
			assert(out_c.hasZero() < 0);
			assert(out_c.isSorted());
			return true;
		}

		_PFROST_D_ int clause_extend_and(const uint32& neg_x, SCLAUSE& org, uint32* defs, const int& nDefs, uint32* out_c)
		{
			assert(neg_x);
			assert(!org.deleted());
			assert(org.size() > 1);
			assert(nDefs > 1);
#pragma unroll
			for (int i = 0; i < nDefs; i++)
				if (org.has(defs[i])) return 0; // tautology
			int len = 0;
#pragma unroll
			for (int i = 0; i < org.size(); i++) {
				uint32 lit = org[i];
				if (lit == neg_x) {
#pragma unroll
					for (int k = 0; k < nDefs; k++) out_c[len++] = FLIP(defs[k]);
				}
				else out_c[len++] = lit;
			}
			devSort(out_c, len);
			return len;
		}

		_PFROST_D_ bool clause_extend_or(const uint32& x, SCLAUSE& org, uint32* defs, const int& nDefs, SCLAUSE& out_c)
		{
			assert(x);
			assert(!org.deleted());
			assert(org.size() > 1);
			assert(nDefs > 1);
			out_c.clear();
#pragma unroll
			for (int i = 0; i < nDefs; i++)
				if (org.has(FLIP(defs[i]))) return false; // tautology
			// extend
#pragma unroll
			for (int i = 0; i < org.size(); i++) {
				uint32 lit = org[i];
				if (lit == x) {
#pragma unroll
					for (int k = 0; k < nDefs; k++) out_c.push(defs[k]);
				}
				else out_c.push(lit);
			}
			// attach 
			devSort(out_c, out_c.size());
			out_c.filter();
			out_c.set_status(ORIGINAL);
			assert(out_c.hasZero() < 0);
			assert(out_c.isSorted());
			return true;
		}

		_PFROST_D_ int clause_extend_or(const uint32& x, SCLAUSE& org, uint32* defs, const int& nDefs, uint32* out_c)
		{
			assert(x);
			assert(!org.deleted());
			assert(org.size() > 1);
			assert(nDefs > 1);
#pragma unroll
			for (int i = 0; i < nDefs; i++) {
				if (org.has(FLIP(defs[i]))) return 0; // tautology
			}
			int len = 0;
#pragma unroll
			for (int i = 0; i < org.size(); i++) {
				uint32 lit = org[i];
				if (lit == x) {
#pragma unroll
					for (uint32 k = 0; k < nDefs; k++) out_c[len++] = defs[k];
				}
				else out_c[len++] = lit;
			}
			devSort(out_c, len);
			return len;
		}

		_PFROST_D_ bool clause_split_and(const uint32& x, const uint32& def, SCLAUSE& org, SCLAUSE& out_c)
		{
			assert(x);
			assert(def);
			assert(!org.deleted());
			assert(org.size() > 1);
			out_c.clear();
			if (org.has(FLIP(def))) return false;  // tautology
			// split
			uint32 sig = 0;
#pragma unroll
			for (int k = 0; k < org.size(); k++) {
				uint32 lit = org[k];
				if (lit == def) continue; // repeated literal
				if (lit == x) { out_c.push(def); sig |= MAPHASH(def); }
				else { out_c.push(lit); sig |= MAPHASH(lit); }
			}
			// attach
			devSort(out_c, out_c.size());
			out_c.set_sig(sig);
			out_c.set_status(ORIGINAL);
			assert(out_c.isSorted());
			assert(out_c.hasZero() < 0);
			return true;
		}

		_PFROST_D_ int clause_split_and(const uint32& x, const uint32& def, SCLAUSE& org, uint32* out_c)
		{
			assert(x);
			assert(def);
			assert(!org.deleted());
			assert(org.size() > 1);
			if (org.has(FLIP(def))) return 0;
			int len = 0;
#pragma unroll
			for (int k = 0; k < org.size(); k++) {
				uint32 lit = org[k];
				if (lit == def) continue; // repeated literal
				if (lit == x) out_c[len++] = def;
				else out_c[len++] = lit;
			}
			devSort(out_c, len);
			return len;
		}

		_PFROST_D_ bool clause_split_or(const uint32& neg_x, const uint32& def, SCLAUSE& org, SCLAUSE& out_c)
		{
			assert(neg_x);
			assert(def);
			assert(!org.deleted());
			assert(org.size() > 1);
			out_c.clear();
			if (org.has(def)) return false;  // tautology
			// split
			uint32 sig = 0;
#pragma unroll
			for (int k = 0; k < org.size(); k++) {
				uint32 lit = org[k];
				if (lit == FLIP(def)) continue; // repeated literal
				if (lit == neg_x) { out_c.push(FLIP(def)); sig |= MAPHASH(FLIP(def)); }
				else { out_c.push(lit); sig |= MAPHASH(lit); }
			}
			// attach
			devSort(out_c, out_c.size());
			out_c.set_sig(sig);
			out_c.set_status(ORIGINAL);
			assert(out_c.isSorted());
			assert(out_c.hasZero() < 0);
			return true;
		}

		_PFROST_D_ int clause_split_or(const uint32& neg_x, const uint32& def, SCLAUSE& org, uint32* out_c)
		{
			assert(neg_x);
			assert(def);
			assert(!org.deleted());
			assert(org.size() > 1);
			if (org.has(def)) return 0;  // tautology
			uint32 len = 0;
#pragma unroll
			for (int k = 0; k < org.size(); k++) {
				uint32 lit = org[k];
				if (lit == FLIP(def)) continue; // repeated literal
				if (lit == neg_x) out_c[len++] = FLIP(def);
				else out_c[len++] = lit;
			}
			devSort(out_c, len);
			return len;
		}

		_PFROST_D_ void calcResolvents(const uint32& x, CNF& cnf, OL& poss, OL& negs, uint32& addedCls, uint64& addedLits)
		{
			assert(x);
			uint32 ps = poss.size(), ns = negs.size(), nTs = 0;
#pragma unroll
			for (uint32 i = 0; i < ps; i++) {
#pragma unroll
				for (uint32 j = 0; j < ns; j++) {
					if (isTautology(x, cnf[poss[i]], cnf[negs[j]])) nTs++;
					else addedLits += cnf[poss[i]].size() + cnf[negs[j]].size() - 2;
				}
			}
			addedCls += ps * ns - nTs;
		}

		_PFROST_D_ void deleteClauses(CNF& cnf, OL& poss, OL& negs)
		{
#pragma unroll
			for (uint32 i = 0; i < poss.size(); i++) cnf[poss[i]].markDeleted();
#pragma unroll
			for (uint32 i = 0; i < negs.size(); i++) cnf[negs[i]].markDeleted();
		}

		_PFROST_D_ void countLitsBefore(CNF& cnf, OL& poss, OL& negs, uint64& nLitsBefore)
		{
#pragma unroll
			for (uint32 i = 0; i < poss.size(); i++) nLitsBefore += cnf[poss[i]].size();
#pragma unroll
			for (uint32 i = 0; i < negs.size(); i++) nLitsBefore += cnf[negs[i]].size();
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

		_PFROST_D_ bool sub(const uint32& A, const uint32& B) { return (A & ~B) == 0; }

		_PFROST_D_ bool selfSub(const uint32& A, const uint32& B)
		{
			uint32 B_tmp = B | ((B & 0xAAAAAAAAUL) >> 1) | ((B & 0x55555555UL) << 1);
			return (A & ~B_tmp) == 0;
		}

		_PFROST_D_ bool sub(SCLAUSE& sm, SCLAUSE& lr)
		{
			assert(!sm.deleted());
			assert(!lr.deleted());
			int it1 = 0, it2 = 0;
			int sm_sz = sm.size(), lr_sz = lr.size(), sub = 0;
			while (it1 < sm_sz && it2 < lr_sz) {
				if (sm[it1] < lr[it2]) it1++;
				else if (lr[it2] < sm[it1]) it2++;
				else { sub++; it1++; it2++; }
			}
			if (sub == sm_sz) return true;
			return false;
		}

		_PFROST_D_ bool sub(SCLAUSE& sm, uint32* lr, const int& lr_sz)
		{
			assert(!sm.deleted());
			int it1 = 0, it2 = 0, sub = 0;
			int sm_sz = sm.size();
			while (it1 < sm_sz && it2 < lr_sz) {
				if (sm[it1] < lr[it2]) it1++;
				else if (lr[it2] < sm[it1]) it2++;
				else { sub++; it1++; it2++; }
			}
			if (sub == sm_sz) return true;
			return false;
		}

		_PFROST_D_ bool selfSub(const uint32& x, SCLAUSE& sm, SCLAUSE& lr)
		{
			assert(!sm.deleted());
			assert(!lr.deleted());
			assert(sm.size() > 1);
			assert(lr.size() > 1);
			int it1 = 0, it2 = 0;
			int sm_sz = sm.size(), lr_sz = lr.size(), sub = 0;
			bool self = false;
			while (it1 < sm_sz && it2 < lr_sz) {
				if (sm[it1] == FLIP(x)) it1++;
				else if (lr[it2] == x) { self = true; it2++; }
				else if (sm[it1] < lr[it2]) it1++;
				else if (lr[it2] < sm[it1]) it2++;
				else { sub++; it1++; it2++; }
			}
			if ((sub + 1) == sm_sz) {
				if (self) return true;
				else {
					while (it2 < lr_sz) {
						if (lr[it2] == x) return true;
						it2++;
					}
				}
			}
			return false;
		}

		_PFROST_D_ bool selfSub(const uint32& x, SCLAUSE& sm, uint32* lr, const int& lr_sz)
		{
			assert(!sm.deleted());
			assert(sm.size() > 1);
			assert(lr_sz > 1);
			int it1 = 0, it2 = 0;
			int sm_sz = sm.size(), sub = 0;
			bool self = false;
			while (it1 < sm_sz && it2 < lr_sz) {
				if (sm[it1] == FLIP(x)) it1++;
				else if (lr[it2] == x) { self = true; it2++; }
				else if (sm[it1] < lr[it2]) it1++;
				else if (lr[it2] < sm[it1]) it2++;
				else { sub++; it1++; it2++; }
			}
			if ((sub + 1) == sm_sz) {
				if (self) return true;
				else {
					while (it2 < lr_sz) {
						if (lr[it2] == x) return true;
						it2++;
					}
				}
			}
			return false;
		}

		_PFROST_D_ void strengthen(cuVecU* units, SCLAUSE& c, uint32* sh_c, uint32 self_lit)
		{
			assert(!c.deleted());
			assert(c.size() > 1);
			assert(self_lit);
			int n = 0;
			uint32 sig = 0;
			bool check = false;
#pragma unroll
			for (int k = 0; k < c.size(); k++) {
				uint32 lit = sh_c[k];
				if (lit != self_lit) {
					c[n] = lit, sh_c[n] = lit;
					sig |= MAPHASH(lit);
					n++;
				}
				else check = true;
			}
			assert(check);
			assert(n == c.size() - 1);
			assert(n <= SH_MAX_HSE_IN);
			assert(c.hasZero() < 0);
			assert(c.isSorted());
			c.set_sig(sig);
			c.pop();
			if (c.size() == 1) units->push(*c);
#if SS_DBG
			if (c.size() == 1) printf("c | SS(%d): ", ABS(self_lit)), c.print();
#endif
		}

		_PFROST_D_ void strengthen(cuVecU* units, SCLAUSE& c, uint32 self_lit)
		{
			assert(!c.deleted());
			assert(c.size() > 1);
			assert(self_lit);
			int n = 0;
			uint32 sig = 0;
			bool check = false;
#pragma unroll
			for (int k = 0; k < c.size(); k++) {
				uint32 lit = c[k];
				if (lit != self_lit) { c[n++] = lit, sig |= MAPHASH(lit); }
				else check = true;
			}
			assert(check);
			assert(n == c.size() - 1);
			assert(n <= SH_MAX_HSE_IN);
			assert(c.hasZero() < 0);
			assert(c.isSorted());
			c.set_sig(sig);
			c.pop();
			if (c.size() == 1) units->push(*c);
#if SS_DBG
			if (c.size() == 1) printf("c | SS(%d): ", ABS(self_lit)), c.print();
#endif
		}

		_PFROST_D_ void reduceOL(CNF& cnf, OL& ol)
		{
			if (ol.size() == 0) return;
			int idx = 0, ol_sz = ol.size();
			while (idx < ol_sz) {
				uint32& r = ol[idx];
				if (cnf[r].deleted()) r = ol[--ol_sz];
				else idx++;
			}
			ol.resize(ol_sz);
		}

		_PFROST_D_ void updateOL(CNF& cnf, OL& ol)
		{
			if (ol.size() == 0) return;
			int idx = 0, ol_sz = ol.size();
			while (idx < ol_sz) {
				uint32& r = ol[idx];
				SCLAUSE& c = cnf[r];
				if (c.deleted()) r = ol[--ol_sz];
				else if (c.molten()) c.freeze(), r = ol[--ol_sz];
				else idx++;
			}
			ol.resize(ol_sz);
		}

		_PFROST_D_ void self_sub_x(const uint32& x, CNF& cnf, OL& poss, OL& negs, cuVecU* units, uint32* sh_c)
		{
			// positives vs negatives
#pragma unroll
			for (uint32 i = 0; i < poss.size(); i++) {
				SCLAUSE& pos = cnf[poss[i]];
				if (pos.deleted()) continue;
				// self-subsumption check
				if (pos.size() > SH_MAX_HSE_IN) { // use global memory 
#pragma unroll
					for (uint32 j = 0; j < negs.size(); j++) {
						SCLAUSE& neg = cnf[negs[j]];
						if (neg.deleted()) continue;
						if (neg.size() > 1 && neg.size() < pos.size() &&
							selfSub(neg.sig(), pos.sig()) && selfSub(x, neg, pos)) {
							strengthen(units, pos, x);
							pos.melt(); // mark for fast recongnition in ot update 
						}
						else if (neg.size() > 1 && neg.size() == pos.size() &&
							selfSub(neg.sig(), pos.sig()) && selfSub(x, neg, pos)) {
							strengthen(units, pos, x);
							pos.melt();
							neg.markDeleted();
						}
					}
					// subsumption check
					for (uint32 j = 0; j < poss.size(); j++) {
						SCLAUSE& sm_c = cnf[poss[j]];
						if (!sm_c.deleted() && sm_c.size() < pos.size() &&
							sub(sm_c.sig(), pos.sig()) && sub(sm_c, pos)) {
							pos.markDeleted();
							break;
						}
					}
				}
				else { // use shared memory
					uint32* sh_pos = sh_c;
					pos.shareTo(sh_pos);
#pragma unroll
					for (uint32 j = 0; j < negs.size(); j++) {
						SCLAUSE& neg = cnf[negs[j]];
						if (neg.deleted()) continue;
						if (neg.size() > 1 && neg.size() < pos.size() &&
							selfSub(neg.sig(), pos.sig()) && selfSub(x, neg, sh_pos, pos.size())) {
							strengthen(units, pos, sh_pos, x);
							pos.melt();
						}
						else if (neg.size() > 1 && neg.size() == pos.size() &&
							selfSub(neg.sig(), pos.sig()) && selfSub(x, neg, sh_pos, pos.size())) {
							strengthen(units, pos, sh_pos, x);
							pos.melt();
							neg.markDeleted();
						}
					}
					for (uint32 j = 0; j < poss.size(); j++) {
						SCLAUSE& sm_c = cnf[poss[j]];
						assert(isCoherent(x, pos, sh_pos, pos.size()));
						if (!sm_c.deleted() && sm_c.size() < pos.size() &&
							sub(sm_c.sig(), pos.sig()) && sub(sm_c, sh_pos, pos.size())) {
							pos.markDeleted();
							break;
						}
					}
				}
			}
#if SS_DBG
			DBG_updateOL(x, cnf, poss); // discard (self)-subsumed clauses
#else
			updateOL(cnf, poss); // discard deleted or (self)-subsumed clauses
#endif
	// negatives vs positives
#pragma unroll
			for (uint32 i = 0; i < negs.size(); i++) {
				SCLAUSE& neg = cnf[negs[i]];
				if (neg.deleted()) continue;
				if (neg.size() > SH_MAX_HSE_IN) { // use global memory
					// self-subsumption check
#pragma unroll
					for (uint32 j = 0; j < poss.size(); j++) {
						SCLAUSE& pos = cnf[poss[j]];
						if (pos.deleted()) continue;
						if (pos.size() > 1 && pos.size() < neg.size() &&
							selfSub(pos.sig(), neg.sig()) && selfSub(NEG(x), pos, neg)) {
							strengthen(units, neg, NEG(x));
							neg.melt();
						}
					}
					for (uint32 j = 0; j < negs.size(); j++) {
						SCLAUSE& sm_c = cnf[negs[j]];
						if (!sm_c.deleted() && sm_c.size() < neg.size() &&
							sub(sm_c.sig(), neg.sig()) && sub(sm_c, neg)) {
							neg.markDeleted();
							break;
						}
					}
				}
				else { // use shared memory
					uint32* sh_neg = sh_c;
					neg.shareTo(sh_neg);
#pragma unroll
					for (uint32 j = 0; j < poss.size(); j++) {
						SCLAUSE& pos = cnf[poss[j]];
						if (pos.deleted()) continue;
						if (pos.size() > 1 && pos.size() < neg.size() &&
							selfSub(pos.sig(), neg.sig()) &&
							selfSub(NEG(x), pos, sh_neg, neg.size())) {
							strengthen(units, neg, sh_neg, NEG(x));
							neg.melt();
						}
					}
					// subsumption check
					for (uint32 j = 0; j < negs.size(); j++) {
						SCLAUSE& sm_c = cnf[negs[j]];
						assert(isCoherent(x, neg, sh_neg, neg.size()));
						if (!sm_c.deleted() && sm_c.size() < neg.size() &&
							sub(sm_c.sig(), neg.sig()) && sub(sm_c, sh_neg, neg.size())) {
							neg.markDeleted();
							break;
						}
					}
				}
			}
#if SS_DBG
			DBG_updateOL(x, cnf, negs); // discard (self)-subsumed clauses
#else
			updateOL(cnf, negs); // discard  deleted or (self)-subsumed clauses
#endif
		}

		_PFROST_D_ void forward_equ(CNF& cnf, OT& ot, uint32* m_c, const int64& m_sig, const uint32& m_len)
		{
			for (int k = 0; k < m_len; k++) {
				assert(m_c[k]);
				OL& ol = ot[m_c[k]];
				for (uint32 i = threadIdx.x; i < ol.size(); i += blockDim.x) {
					SCLAUSE& org = cnf[ol[i]];
					if (!org.deleted() && m_len == org.size() &&
						sub(m_sig, org.sig()) && isEqual(org, m_c, m_len)) {
						org.markDeleted();
						return;
					}
				}
			}
		}

		_PFROST_D_ void blocked_x(const uint32& x, CNF& cnf, OL& poss, OL& negs, uint32* sh_c)
		{
			uint32 ps = poss.size(), ns = negs.size(), nTs = 0;
			if (ps <= ns) {  // start with positives
#pragma unroll
				for (uint32 i = 0; i < ps; i++) {
					SCLAUSE& pos = cnf[poss[i]];
					if (!pos.deleted()) {
						int pos_size = pos.size();
						nTs = 0;
						if (pos_size <= SH_MAX_BCE_IN) { // use shared memory 
							uint32* sh_pos = sh_c;
							pos.shareTo(sh_pos);
#pragma unroll
							for (uint32 j = 0; j < ns; j++) {
								SCLAUSE& neg = cnf[negs[j]];
								if (!neg.deleted() && isTautology(x, neg, sh_pos, pos_size)) nTs++;
							}
						}
						else { // use global memory
#pragma unroll
							for (uint32 j = 0; j < ns; j++) {
								SCLAUSE& neg = cnf[negs[j]];
								if (!neg.deleted() && isTautology(x, pos, neg)) nTs++;
							}
						}
						if (nTs == ns) pos.markDeleted();
					}
				}
			}
			else { // start with negatives
#pragma unroll
				for (uint32 i = 0; i < ns; i++) {
					SCLAUSE& neg = cnf[negs[i]];
					if (!neg.deleted()) {
						int neg_size = neg.size();
						nTs = 0;
						if (neg_size <= SH_MAX_BCE_IN) { // use shared memory
							uint32* sh_neg = sh_c;
							neg.shareTo(sh_neg);
#pragma unroll
							for (uint32 j = 0; j < ps; j++) {
								SCLAUSE& pos = cnf[poss[j]];
								if (!pos.deleted() && isTautology(x, pos, sh_neg, neg_size)) nTs++;
							}
						}
						else { // use global memory
#pragma unroll
							for (uint32 j = 0; j < ps; j++) {
								SCLAUSE& pos = cnf[poss[j]];
								if (!pos.deleted() && isTautology(x, neg, pos)) nTs++;
							}
						}
						if (nTs == ps) neg.markDeleted();
					}
				}
			}
		}

		_PFROST_D_ bool substitute_AND(const uint32& x, CNF& cnf, OL& poss, OL& negs, cuVecU* units, uint32* defs, const int& nDefs, uint32* out_c)
		{
			uint32 ps = poss.size(), ns = negs.size();
			uint32 nAddedCls = 0;
			uint64 nAddedLits = 0;
			// count number of added cls & literals for negatives
#pragma unroll
			for (uint32 i = 0; i < ns; i++) {
				SCLAUSE& neg = cnf[negs[i]];
				if (!isTautology_and(neg, defs, nDefs)) {
					nAddedCls++;
					nAddedLits += neg.size() + nDefs - 1;
				}
			}
			// count number of added cls & literals for positives
#pragma unroll
			for (int d = 0; d < nDefs; d++) {
#pragma unroll
				for (uint32 i = 0; i < ps; i++) {
					SCLAUSE& pos = cnf[poss[i]];
					if (!pos.has(FLIP(defs[d]))) {
						nAddedCls++;
						nAddedLits += pos.size();
					}
				}
			}
			if (nAddedCls == 0) { // No substitutions to add
				deleteClauses(cnf, poss, negs);
				return true;
			}
			if (nAddedCls > ps + ns) return false;
#if VE_DBG
			printf("c | AND(%d): added = %d, deleted = %d\n", ABS(x), nAddedCls, ps + ns), pClauseSet(cnf, poss, negs);
#endif
			S_REF ref, checksum = 0;
			S_REF* cs = cnf.jump(ref, nAddedCls, nAddedLits);
			// substitute negatives 
#pragma unroll
			for (uint32 i = 0; i < ns && checksum < nAddedCls; i++) {
				SCLAUSE& neg = cnf[negs[i]];
				int max_sz = neg.size() + nDefs - 1;
				if (max_sz > SH_MAX_BVE_OUT) { // use global memory
					if (clause_extend_and(NEG(x), neg, defs, nDefs, cnf[ref])) {
						SCLAUSE& added = cnf[ref];
						if (added.size() == 1) units->push(*added);
						cs[checksum++] = ref, ref += added.blockSize();
					}
				}
				else { // use shared memory
					int ext_sz = 0;
					if ((ext_sz = clause_extend_and(NEG(x), neg, defs, nDefs, out_c))) {
						assert(ext_sz <= SH_MAX_BVE_OUT);
						SCLAUSE* added = cnf.cref(ref);
						added = new (added) SCLAUSE();
						added->filterCopy(out_c, ext_sz);
						if (ext_sz == 1) units->push(*out_c);
						cs[checksum++] = ref, ref += (ext_sz - 1) + sizeof(S_REF);
					}
				}
#if VE_DBG
				if (added.size() == 1) printf("== AND(%d)-> unit %d\n", ABS(x), ISNEG(*added) ? -int(ABS(*added)) : ABS(*added));
				if (added.size()) printf("c | "), added.print();
#endif
			}
			// substitute positives
#pragma unroll
			for (int d = 0; d < nDefs; d++) {
#pragma unroll
				for (uint32 i = 0; i < ps && checksum < nAddedCls; i++) {
					SCLAUSE& pos = cnf[poss[i]];
					int max_sz = pos.size();
					if (max_sz > SH_MAX_BVE_OUT) { // use global memory
						if (clause_split_and(x, defs[d], pos, cnf[ref])) {
							SCLAUSE& added = cnf[ref];
							if (added.size() == 1) units->push(*added);
							cs[checksum++] = ref, ref += added.blockSize();
						}
					}
					else { // use shared memory
						int split_sz = 0;
						if (split_sz = clause_split_and(x, defs[d], pos, out_c)) {
							assert(split_sz <= SH_MAX_BVE_OUT);
							uint32 sig = 0;
							calcSig(out_c, split_sz, sig);
							SCLAUSE* added = cnf.cref(ref);
							added = new (added) SCLAUSE(out_c, split_sz);
							added->set_sig(sig);
							if (split_sz == 1) units->push(*out_c);
							cs[checksum++] = ref, ref += (split_sz - 1) + sizeof(S_REF);
							assert(added->isSorted());
							assert(added->hasZero() < 0);
						}
					}
#if VE_DBG
					if (added.size() == 1) printf("== AND(%d)-> unit %d\n", ABS(x), ISNEG(*added) ? -int(ABS(*added)) : ABS(*added));
					if (added.size()) printf("c | "), added.print();
#endif
				}
			}
			assert(checksum == nAddedCls);
			// delete org occurs
			deleteClauses(cnf, poss, negs);
			return true; // AND-substitution successful
		}

		_PFROST_D_ bool substitute_OR(const uint32& x, CNF& cnf, OL& poss, OL& negs, cuVecU* units, uint32* defs, const int& nDefs, uint32* out_c)
		{
			uint32 ps = poss.size(), ns = negs.size();
			uint32 nAddedCls = 0;
			uint64 nAddedLits = 0;
			// count number of added cls & literals for positives
#pragma unroll
			for (uint32 i = 0; i < ps; i++) {
				SCLAUSE& pos = cnf[poss[i]];
				if (!isTautology_or(pos, defs, nDefs)) {
					nAddedCls++;
					nAddedLits += pos.size() + nDefs - 1;
				}
			}
			// count number of added cls & literals for negatives
#pragma unroll
			for (int d = 0; d < nDefs; d++) {
#pragma unroll
				for (uint32 i = 0; i < ns; i++) {
					SCLAUSE& neg = cnf[negs[i]];
					if (!neg.has(defs[d])) {
						nAddedCls++;
						nAddedLits += neg.size();
					}
				}
			}
			if (nAddedCls == 0) { // No substitutions to add
				deleteClauses(cnf, poss, negs);
				return true;
			}
			if (nAddedCls > ps + ns) return false;
#if VE_DBG
			printf("c | OR(%d): added = %d, deleted = %d\n", ABS(x), nAddedCls, ps + ns), pClauseSet(cnf, poss, negs);
#endif
			S_REF ref, checksum = 0;
			S_REF* cs = cnf.jump(ref, nAddedCls, nAddedLits);
			// substitute positives
#pragma unroll
			for (uint32 i = 0; i < ps && checksum < nAddedCls; i++) {
				SCLAUSE& pos = cnf[poss[i]];
				int max_sz = pos.size() + nDefs - 1;
				if (max_sz > SH_MAX_BVE_OUT) { // use global memory
					if (clause_extend_or(x, pos, defs, nDefs, cnf[ref])) {
						SCLAUSE& added = cnf[ref];
						if (added.size() == 1) units->push(*added);
						cs[checksum++] = ref, ref += added.blockSize();
					}
				}
				else { // use shared memory
					int ext_sz = 0;
					if ((ext_sz = clause_extend_or(x, pos, defs, nDefs, out_c))) {
						assert(ext_sz <= SH_MAX_BVE_OUT);
						SCLAUSE* added = cnf.cref(ref);
						added = new (added) SCLAUSE();
						added->filterCopy(out_c, ext_sz);
						if (ext_sz == 1) units->push(*out_c);
						cs[checksum++] = ref, ref += (ext_sz - 1) + sizeof(S_REF);
					}
				}
#if VE_DBG
				if (added.size() == 1) printf("== OR(%d)-> unit %d\n", ABS(x), ISNEG(*added) ? -int(ABS(*added)) : ABS(*added));
				if (added.size()) printf("c | "), added.print();
#endif
			}
			// substitute negatives
#pragma unroll
			for (int d = 0; d < nDefs; d++) {
#pragma unroll
				for (uint32 i = 0; i < ns && checksum < nAddedCls; i++) {
					SCLAUSE& neg = cnf[negs[i]];
					int max_sz = neg.size();
					if (max_sz > SH_MAX_BVE_OUT) { // use global memory
						if (clause_split_or(NEG(x), defs[d], neg, cnf[ref])) {
							SCLAUSE& added = cnf[ref];
							if (added.size() == 1) units->push(*added);
							cs[checksum++] = ref, ref += added.blockSize();
						}
					}
					else { // use shared memory
						int split_sz = 0;
						if (split_sz = clause_split_or(NEG(x), defs[d], neg, out_c)) {
							assert(split_sz <= SH_MAX_BVE_OUT);
							uint32 sig = 0;
							calcSig(out_c, split_sz, sig);
							SCLAUSE* added = cnf.cref(ref);
							added = new (added) SCLAUSE(out_c, split_sz);
							added->set_sig(sig);
							if (split_sz == 1) units->push(*out_c);
							cs[checksum++] = ref, ref += (split_sz - 1) + sizeof(S_REF);
							assert(added->isSorted());
							assert(added->hasZero() < 0);
						}
					}
#if VE_DBG
					if (added.size() == 1) printf("== OR(%d)-> unit %d\n", ABS(x), ISNEG(*added) ? -int(ABS(*added)) : ABS(*added));
					if (added.size()) printf("c | "), added.print();
#endif
				}
			}
			assert(checksum == nAddedCls);
			// delete org occurs
			deleteClauses(cnf, poss, negs);
			return true; // OR-substitution successful
		}

		_PFROST_D_ bool gateReasoning_x(const uint32& x, CNF& cnf, OL& poss, OL& negs, cuVecU* units, uint32* defs, uint32* out_c)
		{
			uint32 ps = poss.size(), ns = negs.size();
			uint32 imp = 0, sig = 0;
			int nImps = 0;
			// AND-gate Reasoning
			if (ns > (FAN_LMT - 1)) return false; // shared memory GUARD
#pragma unroll
			for (uint32 i = 0; i < ns; i++) {
				SCLAUSE& neg = cnf[negs[i]];
				if (neg.size() == 2) {
					if ((neg[0] ^ x) == NEG_SIGN) { // found x with negative sign
						imp = FLIP(neg[1]); // toggle implied literal sign
						out_c[nImps++] = imp;
						sig |= MAPHASH(imp);
					}
					else if ((neg[1] ^ x) == NEG_SIGN) {
						imp = FLIP(neg[0]); // toggle implied literal sign
						out_c[nImps++] = imp;
						sig |= MAPHASH(imp);
					}
				}
			}
			if (nImps > 1) {
				assert(nImps <= FAN_LMT - 1);
				out_c[nImps++] = x;
				sig |= MAPHASH(x);
				devSort(out_c, nImps);
#pragma unroll
				for (uint32 i = 0; i < ps; i++) {
					SCLAUSE& pos = cnf[poss[i]];
					int pos_size = pos.size();
					if (pos_size == nImps && sub(pos.sig(), sig) && isEqual(pos, out_c, nImps)) {
						int nDefs = 0;
#pragma unroll
						for (int j = 0; j < nImps; j++) {
							uint32 lit = out_c[j];
							if (lit != x) defs[nDefs++] = FLIP(lit);
						}
						assert(nDefs < FAN_LMT);
						assert(nDefs == pos_size - 1);
						if (substitute_AND(x, cnf, poss, negs, units, defs, nDefs, out_c)) return true;
					}
				}
			}
			// OR-gate Reasoning
			if (ps > (FAN_LMT - 1)) return false; // shared memory GUARD
			imp = 0, nImps = 0, sig = 0;
#pragma unroll
			for (uint32 i = 0; i < ps; i++) {
				SCLAUSE& pos = cnf[poss[i]];
				if (pos.size() == 2) {
					if (pos[0] == x) { // found x with positive sign
						imp = FLIP(pos[1]); // toggle implied literal sign
						out_c[nImps++] = imp;
						sig |= MAPHASH(imp);
					}
					else if (pos[1] == x) {
						imp = FLIP(pos[0]); // toggle implied literal sign
						out_c[nImps++] = imp;
						sig |= MAPHASH(imp);
					}
				}
			}
			if (nImps > 1) {
				assert(nImps <= FAN_LMT - 1);
				out_c[nImps++] = FLIP(x);
				sig |= MAPHASH(FLIP(x));
				devSort(out_c, nImps);
#pragma unroll
				for (uint32 i = 0; i < ns; i++) {
					SCLAUSE& neg = cnf[negs[i]];
					int neg_size = neg.size();
					if (neg_size == nImps && sub(neg.sig(), sig) && isEqual(neg, out_c, nImps)) {
						int nDefs = 0;
#pragma unroll
						for (int j = 0; j < nImps; j++) {
							uint32 lit = out_c[j];
							if (lit != FLIP(x)) defs[nDefs++] = lit;
						}
						assert(nDefs < FAN_LMT);
						assert(nDefs == neg_size - 1);
						if (substitute_OR(x, cnf, poss, negs, units, defs, nDefs, out_c)) return true;
					}
				}
			}
			return false;
		}

		_PFROST_D_ bool resolve_x(const uint32& x, CNF& cnf, OL& poss, OL& negs, cuVecU* units, uint32* out_c)
		{
			uint32 ps = poss.size(), ns = negs.size();
			uint32 nAddedCls = 0;
			uint64 nAddedLits = 0;
			calcResolvents(x, cnf, poss, negs, nAddedCls, nAddedLits);
			if (nAddedCls == 0) { // No resolvents to add
				deleteClauses(cnf, poss, negs);
				return true;
			}
			if (nAddedCls > ps + ns) return false;
			uint64 nLitsBefore = 0;
			countLitsBefore(cnf, poss, negs, nLitsBefore);
			if (nAddedLits > nLitsBefore) return false;
#if VE_DBG
			printf("c | Resolving %d -> adding %d, deleting %d\n", x, nAddedCls, ps + ns), pClauseSet(cnf, poss, negs);
#endif
			// resolve x
			S_REF ref, checksum = 0;
			S_REF* cs = cnf.jump(ref, nAddedCls, nAddedLits);
#pragma unroll
			for (uint32 i = 0; i < ps; i++) {
				SCLAUSE& pos = cnf[poss[i]];
#pragma unroll
				for (uint32 j = 0; j < ns && checksum < nAddedCls; j++) {
					SCLAUSE& neg = cnf[negs[j]];
					int max_sz = pos.size() + neg.size() - 2;
					if (max_sz > SH_MAX_BVE_OUT) { // use global memory
						if (merge(x, pos, neg, cnf[ref])) {
							SCLAUSE& added = cnf[ref];
							if (added.size() == 1) units->push(*added);
							cs[checksum++] = ref, ref += added.blockSize();
						}
					}
					else { // use shared memory
						int merged_sz = 0;
						if ((merged_sz = merge(x, pos, neg, out_c))) {
							assert(merged_sz <= SH_MAX_BVE_OUT);
							uint32 sig = 0;
							calcSig(out_c, merged_sz, sig);
							SCLAUSE* added = cnf.cref(ref);
							added = new (added) SCLAUSE(out_c, merged_sz);
							added->set_sig(sig);
							if (merged_sz == 1) units->push(*out_c);
							cs[checksum++] = ref, ref += (merged_sz - 1) + sizeof(S_REF);
							assert(added->isSorted());
							assert(added->hasZero() < 0);
						}
					}
				}
			}
#if VE_DBG
			printf("c | Added %d resolvents\n", checksum);
			printf("c | Units"), units->print();
#endif
			assert(checksum == nAddedCls);
			// delete org occurs
			deleteClauses(cnf, poss, negs);
			return true; // resolution successful
		}

	} // namespace sigma
} // namespace pfrost

#endif