/***********************************************************************[subsume.cuh]
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

#ifndef __GPU_SUB_
#define __GPU_SUB_

#include "elimination.cuh"

namespace ParaFROST {

#define SS_DBG 0 // set to 1 to serialize SUB

	_PFROST_D_ bool isCoherent(SCLAUSE& g_c, uint32* sh_c, const int& size)
	{
		uint32 it = 0;
		if (g_c.size() != size) {
			printf("c  WARNING - memory incoherency for -> clause size not equal\n");
			printf("c  global clause: "), g_c.print();
			printf("c  shared clause: "), print_shared_clause(sh_c, size);
			return false;
		}
		while (it < size) {
			if (sh_c[it] != g_c[it]) {
				printf("c  WARNING - memory incoherency for -> clauses not equal\n");
				printf("c  global clause: "), g_c.print();
				printf("c  shared clause: "), print_shared_clause(sh_c, size);
				return false;
			}
			else it++;
		}
		return true;
	}

	_PFROST_D_ bool sub(SCLAUSE& subsuming, SCLAUSE& subsumed)
	{
		assert(!subsuming.deleted());
		assert(!subsumed.deleted());
		assert(subsuming.size() > 1);
		assert(subsumed.size() > 1);
		assert(subsuming.size() <= subsumed.size());
		const int size = subsuming.size();
		uint32* d1 = subsuming.data(), * d2 = subsumed.data();
		uint32* e1 = d1 + size;
		uint32* e2 = subsumed.end();
		int sub = 0;
		while (d1 != e1 && d2 != e2) {
			const uint32 lit1 = *d1, lit2 = *d2;
			if (lit1 < lit2)
				d1++;
			else if (lit2 < lit1)
				d2++;
			else {
				sub++;
				d1++, d2++;
			}
		}
		if (sub == size)
			return true;
		return false;
	}

	_PFROST_D_ bool sub(SCLAUSE& subsuming, uint32* subsumed, const int& subsumedsize)
	{
		assert(!subsuming.deleted());
		assert(subsuming.size() > 1);
		assert(subsumedsize > 1);
		assert(subsuming.size() <= subsumedsize);
		const int size = subsuming.size();
		uint32* d1 = subsuming.data(), * d2 = subsumed;
		uint32* e1 = d1 + size;
		uint32* e2 = d2 + subsumedsize;
		int sub = 0;
		while (d1 != e1 && d2 != e2) {
			const uint32 lit1 = *d1, lit2 = *d2;
			if (lit1 < lit2)
				d1++;
			else if (lit2 < lit1)
				d2++;
			else {
				sub++;
				d1++, d2++;
			}
		}
		if (sub == size)
			return true;
		return false;
	}

	_PFROST_D_ bool sub(uint32* subsuming, const int& size, SCLAUSE& subsumed)
	{
		assert(!subsumed.deleted());
		assert(size > 1);
		assert(subsumed.size() > 1);
		assert(size <= subsumed.size());
		uint32* d1 = subsuming, *d2 = subsumed.data();
		uint32* e1 = d1 + size;
		uint32* e2 = subsumed.end();
		int sub = 0;
		while (d1 != e1 && d2 != e2) {
			const uint32 lit1 = *d1, lit2 = *d2;
			if (lit1 < lit2)
				d1++;
			else if (lit2 < lit1)
				d2++;
			else {
				sub++;
				d1++, d2++;
			}
		}
		if (sub == size)
			return true;
		return false;
	}

	_PFROST_D_ bool selfsub(const uint32& x, const uint32& fx, SCLAUSE& subsuming, SCLAUSE& subsumed)
	{
		assert(!subsuming.deleted());
		assert(!subsumed.deleted());
		assert(subsuming.size() > 1);
		assert(subsumed.size() > 1);
		assert(subsuming.size() <= subsumed.size());
		const int size = subsuming.size();
		uint32* d1 = subsuming.data();
		uint32* d2 = subsumed.data();
		uint32* e1 = d1 + size;
		uint32* e2 = subsumed.end();
		int sub = 0;
		bool self = false;
		while (d1 != e1 && d2 != e2) {
			const uint32 lit1 = *d1;
			const uint32 lit2 = *d2;
			if (lit1 == fx)
				d1++;
			else if (lit2 == x) {
				self = true;
				d2++;
			}
			else if (lit1 < lit2)
				d1++;
			else if (lit2 < lit1)
				d2++;
			else {
				sub++;
				d1++, d2++;
			}
		}
		if ((sub + 1) == size) {
			if (self)
				return true;
			else {
				while (d2 != e2) {
					if (*d2 == x)
						return true;
					d2++;
				}
			}
		}
		return false;
	}

	_PFROST_D_ bool selfsub(const uint32& x, const uint32& fx, SCLAUSE& subsuming, uint32* subsumed, const int& subsumedsize)
	{
		assert(!subsuming.deleted());
		assert(subsuming.size() > 1);
		assert(subsumedsize > 1);
		assert(subsuming.size() <= subsumedsize);
		const int size = subsuming.size();
		uint32* d1 = subsuming.data();
		uint32* d2 = subsumed;
		uint32* e1 = d1 + size;
		uint32* e2 = d2 + subsumedsize;
		int sub = 0;
		bool self = false;
		while (d1 != e1 && d2 != e2) {
			const uint32 lit1 = *d1;
			const uint32 lit2 = *d2;
			if (lit1 == fx)
				d1++;
			else if (lit2 == x) {
				self = true;
				d2++;
			}
			else if (lit1 < lit2)
				d1++;
			else if (lit2 < lit1)
				d2++;
			else {
				sub++;
				d1++, d2++;
			}
		}
		if ((sub + 1) == size) {
			if (self)
				return true;
			else {
				while (d2 != e2) {
					if (*d2 == x)
						return true;
					d2++;
				}
			}
		}
		return false;
	}

	_PFROST_D_ void bumpShrunken(SCLAUSE& c)
	{
		assert(c.learnt());
		assert(c.size() > 1);
		const int old_lbd = c.lbd();
		if (old_lbd <= LBD_TIER1) return; // always keep Tier1 value
		const int new_lbd = MIN(c.size() - 1, old_lbd);
		if (new_lbd >= old_lbd) return;
		c.set_lbd(new_lbd);
		c.set_usage(USAGET3);
#if SS_DBG
		printf("c  Bumping shrunken clause with (lbd:%d, usage:%d):\n", new_lbd, c.usage());
		c.print();
#endif
	}

	_PFROST_D_ void strengthen(SCLAUSE& c, uint32* sh_c, const uint32& self)
	{
		assert(!c.deleted());
		assert(c.size() > 1);
		assert(self > 1);

		const int size = c.size();
		int n = 0;
		for (int k = 0; k < size; k++) {
			const uint32 lit = sh_c[k];
			if (NEQUAL(lit, self)) {
				c[n] = lit;
				sh_c[n] = lit;
				n++;
			}
		}

		assert(n == size - 1);
		assert(n <= SH_MAX_SUB_IN);
		assert(c.hasZero() < 0);
		assert(c.isSorted());

		c.pop();

		if (n > 1) {
			calcSig(c);
			if (c.learnt())
				bumpShrunken(c);
		}
	}

	_PFROST_D_ void strengthen(SCLAUSE& c, const uint32& self)
	{
		assert(!c.deleted());
		assert(c.size() > 1);
		assert(self > 1);

		uint32* j = c;
		forall_clause(c, k) {
			const uint32 lit = *k;
			if (NEQUAL(lit, self)) 
				*j++ = lit;
		}

		assert(c.hasZero() < 0);
		assert(c.isSorted());

		c.pop();

		if (c.size() > 1) {
			calcSig(c);
			if (c.learnt())
				bumpShrunken(c);
		}
	}

	_PFROST_D_ void updateOL(CNF& cnf, OL& ol)
	{
		if (ol.empty()) return;
		S_REF* j = ol;
		forall_occurs(ol, i) {
			SCLAUSE& c = cnf[*i];
			if (c.molten()) c.freeze();
			else if (!c.deleted()) *j++ = *i;
		}
		ol.resize(j - ol);
	}

	_PFROST_D_ void subsume(CNF& cnf, OL& list, S_REF* end, SCLAUSE& cand)
	{
		const int candsz = cand.size();
		for (S_REF* j = list; j != end; j++) {
			SCLAUSE& subsuming = cnf[*j];
			if (subsuming.deleted()) continue;
			if (cand.molten() && subsuming.size() > candsz) continue;
			if (subsuming.size() > 1 && SUBSIG(subsuming.sig(), cand.sig()) && sub(subsuming, cand)) {
				if (subsuming.learnt() && cand.original()) subsuming.set_status(ORIGINAL);
				cand.markDeleted();
				#if SS_DBG
				printf("c  subsumed"), cand.print();
				printf("c        by"), subsuming.print();
				#endif
				break;
			}
		}
	}

	_PFROST_D_ void subsume(CNF& cnf, OL& list, S_REF* end, SCLAUSE& cand, uint32* sh_cand)
	{
		const int candsz = cand.size();
		for (S_REF* j = list; j != end; j++) {
			SCLAUSE& subsuming = cnf[*j];
			if (subsuming.deleted()) continue;
			if (cand.molten() && subsuming.size() > candsz) continue;
			assert(isCoherent(cand, sh_cand, candsz));
			if (subsuming.size() > 1 && SUBSIG(subsuming.sig(), cand.sig()) && sub(subsuming, sh_cand, candsz)) {
				if (subsuming.learnt() && cand.original()) subsuming.set_status(ORIGINAL);
				cand.markDeleted();
				#if SS_DBG
				printf("c  subsumed"), cand.print();
				printf("c        by"), subsuming.print();
				#endif
				break;
			}
		}
	}

	_PFROST_D_ void selfsubsume(
		const uint32& x, 
		const uint32& fx, 
		CNF& cnf, 
		OL& list,
		SCLAUSE& cand,
		uint32& nUnits)
	{
		// try to strengthen 'cand' by removing 'x'
		const int candsz = cand.size();
		const uint32 candsig = cand.sig();
		forall_occurs(list, j) {
			SCLAUSE& subsuming = cnf[*j];
			const int subsize = subsuming.size();
			if (subsize > candsz) break;
			if (subsuming.deleted() || subsuming.molten()) continue;
			if (subsize > 1 && selfsub(subsuming.sig(), candsig) && selfsub(x, fx, subsuming, cand)) {
				strengthen(cand, x);
				cand.melt(); // mark for fast recongnition in ot update 
				if (cand.size() == 1) nUnits++;
				#if SS_DBG
				printf("c  %d strengthend ", ABS(x)), cand.print();
				#endif
				break; // cannot strengthen "cand" anymore, 'x' already removed
			}
		}
	}

	_PFROST_D_ void selfsubsume(
		const uint32& x,
		const uint32& fx, 
		CNF& cnf, 
		OL& list,
		SCLAUSE& cand,
		uint32* sh_cand,
		uint32& nUnits)
	{
		// try to strengthen 'cand' by removing 'x' using shared memory 
		const int candsz = cand.size();
		const uint32 candsig = cand.sig();
		forall_occurs(list, j) {
			SCLAUSE& subsuming = cnf[*j];
			const int subsize = subsuming.size();
			if (subsize > candsz) break;
			if (subsuming.deleted() || subsuming.molten()) continue;
			if (subsize > 1 && selfsub(subsuming.sig(), candsig) && selfsub(x, fx, subsuming, sh_cand, candsz)) {
				strengthen(cand, sh_cand, x);
				cand.melt();
				if (cand.size() == 1) nUnits++;
				#if SS_DBG
				printf("c  %d strengthend ", ABS(x)), cand.print();
				#endif
				break;
			}
		}
	}

	// kernels
	__global__ void sub_k(
		CNF* __restrict__ cnfptr,
		OT* __restrict__ ot,
		cuVecB* __restrict__ proof,
		cuVecU* __restrict__ units,
		const cuVecU* __restrict__ elected,
		const Byte* __restrict__ eliminated)
	{
		uint32* sh_cls = SharedMemory<uint32>();
		for_parallel_x (tid, elected->size()) {
			const uint32 x = (*elected)[tid];
			assert(x);
			assert(!ELIMINATED(eliminated[x]));
			const uint32 p = V2L(x), n = NEG(p);
			OL& poss = (*ot)[p];
			OL& negs = (*ot)[n];
			if (poss.size() <= kOpts->sub_max_occurs && negs.size() <= kOpts->sub_max_occurs) {
				CNF cnf = *cnfptr;
				assert(checkMolten(cnf, poss, negs));
				const uint32 n = NEG(p);
				uint32 nPosUnits = 0;
				// positives vs negatives
				forall_occurs(poss, i) {
					SCLAUSE& pos = cnf[*i];
					if (pos.size() > SUB_MAX_CL_SIZE) break;
					if (pos.deleted()) continue;
					if (pos.size() > SH_MAX_SUB_IN) { // use global memory 
						selfsubsume(p, n, cnf, negs, pos, nPosUnits);
						subsume(cnf, poss, i, pos);
					}
					else { // use shared memory
						uint32* sh_pos = sh_cls + threadIdx.x * SH_MAX_SUB_IN;
						pos.shareTo(sh_pos);
						selfsubsume(p, n, cnf, negs, pos, sh_pos, nPosUnits);
						subsume(cnf, poss, i, pos, sh_pos);
					}
				}
				// negatives vs positives
				uint32 nNegUnits = 0;
				forall_occurs(negs, i) {
					SCLAUSE& neg = cnf[*i];
					if (neg.size() > SUB_MAX_CL_SIZE) break;
					if (neg.deleted()) continue;
					if (neg.size() > SH_MAX_SUB_IN) {
						selfsubsume(n, p, cnf, poss, neg, nNegUnits);
						subsume(cnf, negs, i, neg);
					}
					else { // use shared memory
						uint32* sh_neg = sh_cls + threadIdx.x * SH_MAX_SUB_IN;
						neg.shareTo(sh_neg);
						selfsubsume(n, p, cnf, poss, neg, sh_neg, nNegUnits);
						subsume(cnf, negs, i, neg, sh_neg);
					}
				}
				// add units
				if (nPosUnits || nNegUnits) {
					uint32* ustart = units->jump(nPosUnits + nNegUnits);
					uint32* udata = ustart;
					if (nPosUnits) appendUnits(cnf, poss, udata);
					if (nNegUnits) appendUnits(cnf, negs, udata);
					assert(udata == ustart + nPosUnits + nNegUnits);
				}
				// add proof and filter the lists
				if (proof) {
					uint32 bytes = 0;
					countProofSub(cnf, poss, bytes);
					countProofSub(cnf, negs, bytes);
					addr_t pstart = proof->jump(bytes);
					addr_t pdata = pstart;
					assert(pdata);
					saveProof(cnf, poss, pdata);
					saveProof(cnf, negs, pdata);
					assert(pdata == pstart + bytes);
				}
				else {
					updateOL(cnf, poss);
					updateOL(cnf, negs);
				}
				assert(checkMolten(cnf, poss, negs));
				assert(checkDeleted(cnf, poss, negs));
			}
		}
	}

} 

#endif