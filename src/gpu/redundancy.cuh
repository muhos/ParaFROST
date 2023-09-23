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

	#define ERE_PUSH_LIT(LIT) \
	{ \
		lsize = ot[LIT].size(); \
		if (lsize < minsize) \
			minsize = lsize, best = LIT; \
		sig |= MAPHASH(LIT); \
		if (len < clause_max) \
			out_c[len++] = LIT; \
		else \
			return 0; \
	}

	_PFROST_D_ int _merge_ere(
		const uint32& x,
		const int clause_max,
		const OT& ot,
		SCLAUSE& c1,
		SCLAUSE& c2,
		uint32* out_c,
		uint32& sig,
		uint32& best)
	{
		assert(x);
		assert(!c1.deleted());
		assert(!c2.deleted());
		assert(c1.size() > 1);
		assert(c2.size() > 1);

		sig = 0;
		best = 0;

		int lsize = INT_MAX, minsize = INT_MAX;
		uint32* d1 = c1.data();
		uint32* d2 = c2.data();
		uint32* e1 = c1.end();
		uint32* e2 = c2.end();
		uint32 lit1, lit2, v1, v2;
		int len = 0;
		while (d1 != e1 && d2 != e2) {
			lit1 = *d1;
			lit2 = *d2;
			v1 = ABS(lit1);
			v2 = ABS(lit2);
			if (v1 == x) d1++;
			else if (v2 == x) d2++;
			else if (IS_TAUTOLOGY(lit1, lit2)) return 0;
			else if (v1 < v2) {
				ERE_PUSH_LIT(lit1);
				d1++;
			}
			else if (v2 < v1) {
				ERE_PUSH_LIT(lit2);
				d2++;
			}
			else {
				assert(lit1 == lit2);
				ERE_PUSH_LIT(lit1);
				d1++, d2++;
			}
		}

		while (d1 != e1) {
			lit1 = *d1;
			if (NEQUAL(ABS(lit1), x)) {
				ERE_PUSH_LIT(lit1);
			}
			d1++;
		}

		while (d2 != e2) {
			lit2 = *d2;
			if (NEQUAL(ABS(lit2), x)) {
				ERE_PUSH_LIT(lit2);
			}
			d2++;
		}

		return len;
	}

	// kernel
	__global__ void ere_k(
		CNF* __restrict__ cnf,
		OT* __restrict__ ot,
		cuVecB* __restrict__ proof,
		const cuVecU* __restrict__ pVars,
		const Byte* __restrict__ eliminated)
	{
		grid_t gid = global_ty;
		uint32* smem = SharedMemory<uint32>();
		while (gid < pVars->size()) {
			const uint32 v = pVars->at(gid);
			assert(v);
			assert(!ELIMINATED(eliminated[v]));
			const uint32 p = V2L(v), n = NEG(p);
			OL& poss = (*ot)[p], & negs = (*ot)[n];
			const int ds = poss.size(), fs = negs.size();

			// do merging and apply forward equality check (on-the-fly) over resolvents
			if (ds && fs && ds <= kOpts->ere_max_occurs && fs <= kOpts->ere_max_occurs) {

				const int clause_max = kOpts->ere_clause_max;

				// first clause in sorted lists is already larger than the max.
				if ((*cnf)[*poss].size() > clause_max ||
					(*cnf)[*negs].size() > clause_max)
					continue;

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
			gid += stride_y;
		}
	}

} // parafrost namespace


#endif