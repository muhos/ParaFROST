/***********************************************************************[ere.cuh]
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

#ifndef __GPU_ERE_
#define __GPU_ERE_

#include "elimination.cuh"
#include "mutex.cuh"

namespace ParaFROST {

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
			else if (IS_TAUTOLOGY(lit1, lit2)) return 0;
			else if (v1 < v2) { it1++; out_c[len++] = lit1; }
			else if (v2 < v1) { it2++; out_c[len++] = lit2; }
			else { // repeated literal
				it1++, it2++;
				out_c[len++] = lit1;
			}
		}
		while (it1 < n1) {
			lit1 = c1[it1++];
			if (NEQUAL(ABS(lit1), x)) out_c[len++] = lit1;
		}
		while (it2 < n2) {
			lit2 = c2[it2++];
			if (NEQUAL(ABS(lit2), x)) out_c[len++] = lit2;
		}
		assert(len <= kOpts->ere_clause_max);
		return len;
	}

	_PFROST_D_ bool equal_ere( const SCLAUSE& c1,  const uint32* c2, const int& size)
	{
		if (c1.deleted()) return false;
		assert(c1.size() > 1);
		assert(size > 1);
		for (int it = 0; it < size; it++) 
			if (NEQUAL(c1[it], c2[it])) 
				return false;
		return true;
	}

	_PFROST_D_ bool equal_ere_vec(const SCLAUSE& c1, const uint32* c2, int size) {
		// c1 is not deleted; sizes already matched
		const uint32* __restrict__ a = c1.data();
		const uint32* __restrict__ b = c2;

		// Try 16-byte aligned path
		bool aligned = (((uintptr_t)a | (uintptr_t)b) & 0xF) == 0;
		int i = 0;

		if (aligned) {
			const uint4* A = reinterpret_cast<const uint4*>(a);
			const uint4* B = reinterpret_cast<const uint4*>(b);
			int n4 = size >> 2; // size / 4
			for (int k = 0; k < n4; ++k) {
				uint4 x = A[k], y = B[k];
				if (x.x != y.x || x.y != y.y || x.z != y.z || x.w != y.w) return false;
			}
			i = n4 << 2;
		}

		// tail
		for (; i < size; ++i) if (a[i] != b[i]) return false;
		return true;
	}

	_PFROST_D_ void forward_equ(CNF& cnf, OT& ot, cuVecB* proof, uint32* m_c, const int& m_len, const CL_ST& type)
	{
		assert(m_len > 1);
		assert(type != DELETED);
		uint32 best = *m_c, m_sig = MAPHASH(best);
		assert(best > 1);
		int minsize = ot[best].size();
		for (int k = 1; k < m_len; k++) {
			const uint32 lit = m_c[k];
			int lsize = ot[lit].size();
			if (lsize < minsize) minsize = lsize, best = lit;
			m_sig |= MAPHASH(lit);
		}
		const OL& minList = ot[best];
		assert(minsize == minList.size());
		
		// search for the clause in the smaller occurrence listcd
		for (int i = threadIdx.x; i < minsize; i += blockDim.x) {
			SCLAUSE& c = cnf[minList[i]];
			if (m_len == c.size() && (c.learnt() || (c.status() == type)) &&
				SUBSIG(m_sig, c.sig()) && equal_ere(c, m_c, m_len)) {
				assert(c.size() > 1);
				c.markDeleted();
				if (proof) {
					uint32 bytes = 0;
					countProofBytes(c, bytes);
					addr_t pstart = proof->jump(bytes);
					addr_t pdata = pstart;
					saveProofClause(pdata, c, PROOF_DELETED);
					assert(pdata == pstart + bytes);
				}
				break;
			}
		}
	}
	
	// kernel
	__global__ void ere_k(
		CNF* __restrict__ cnf,
		OT* __restrict__ ot,
		cuVecB* __restrict__ proof,
		const cuVecU* __restrict__ pVars,
		const Byte* __restrict__ eliminated)
	{
		uint32* smem = SharedMemory<uint32>();
		for_parallel_y(gid, pVars->size()) {
			const uint32 v = pVars->at(gid);
			assert(v);
			assert(!ELIMINATED(eliminated[v]));
			const uint32 p = V2L(v), n = NEG(p);
			OL& poss = (*ot)[p], & negs = (*ot)[n];
			const int ds = poss.size(), fs = negs.size();
			// do merging and apply forward equality check (on-the-fly) over resolvents
			const int clause_max = kOpts->ere_clause_max;
			if (ds && fs && ds <= kOpts->ere_max_occurs && fs <= kOpts->ere_max_occurs &&
				// first clause in sorted lists must be withing the max.
				((*cnf)[*poss].size() <= clause_max && (*cnf)[*negs].size() <= clause_max)
				) {
				forall_occurs(poss, i) {
					SCLAUSE& pos = (*cnf)[*i];
					if (pos.deleted()) continue;
					forall_occurs(negs, j) {
						SCLAUSE& neg = (*cnf)[*j];
						if (neg.deleted() || (pos.size() + neg.size() - 2) > clause_max) continue;
						uint32* m_c = smem + threadIdx.y * clause_max; // shared memory for resolvent
						int m_len = 0;
						if ((m_len = merge_ere(v, pos, neg, m_c)) > 1) {
							CL_ST type = (pos.learnt() || neg.learnt()) ? LEARNT : ORIGINAL;
							forward_equ(*cnf, *ot, proof, m_c, m_len, type);
						}
					}
				}
			}
		}
	}

} // parafrost namespace


#endif