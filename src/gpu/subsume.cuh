/***********************************************************************[subsume.cuh]
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

#ifndef __SIGMA_HSE_
#define __SIGMA_HSE_

#include "device.cuh"

namespace pFROST {

	namespace SIGmA {

#define SS_DBG 0 // set to 1 to serialize HSE
#define LBD_TIER1 2

#if SS_DBG
		_PFROST_D_ void DBG_updateOL(const uint32& x, CNF& cnf, OL& ol)
		{
			if (ol.size() == 0) return;
			uint32 ol_sz = ol.size();
			S_REF* i, * j, * rend = ol.end();
			for (i = ol, j = ol; i != rend; i++) {
				SCLAUSE& c = cnf[*i];
				if (c.molten() || c.deleted()) printf("c  (self)-subsumed clause "), c.print();
				if (c.molten()) c.freeze();
				else if (!c.deleted()) *j++ = *i;
			}
			ol._shrink(rend - j);
			if (ol_sz != ol.size()) { printf("c  Updated list:\n"); pClauseSet(cnf, ol); printf("c  == End of list ==\n"); }
		}
#endif

		_PFROST_D_ bool isCoherent(SCLAUSE& g_c, uint32* sh_c, const int& size)
		{
			uint32 it = 0;
			if (g_c.size() != size) {
				printf("c  WARNING - memory incoherency for -> clause size not equal\n");
				printf("c  global clause: "), g_c.print();
				printf("c  shared clause: "), pSharedClause(sh_c, size);
				return false;
			}
			while (it < size) {
				if (sh_c[it] != g_c[it]) {
					printf("c  WARNING - memory incoherency for -> clauses not equal\n");
					printf("c  global clause: "), g_c.print();
					printf("c  shared clause: "), pSharedClause(sh_c, size);
					return false;
				}
				else it++;
			}
			return true;
		}

		_PFROST_D_ bool sub(SCLAUSE& sm, SCLAUSE& lr)
		{
			assert(!sm.deleted());
			assert(!lr.deleted());
			assert(sm.size() > 1);
			assert(lr.size() > 1);
			assert(sm.size() <= lr.size());
			int it1 = 0, it2 = 0, sub = 0;
			const int sm_sz = sm.size(), lr_sz = lr.size();
			while (it1 < sm_sz && it2 < lr_sz) {
				const uint32 lit1 = sm[it1], lit2 = lr[it2];
				if (lit1 < lit2) it1++;
				else if (lit2 < lit1) it2++;
				else { sub++; it1++; it2++; }
			}
			if (sub == sm_sz) return true;
			return false;
		}

		_PFROST_D_ bool sub(SCLAUSE& sm, uint32* lr, const int& lr_sz)
		{
			assert(!sm.deleted());
			assert(sm.size() > 1);
			assert(lr_sz > 1);
			assert(sm.size() <= lr_sz);
			int it1 = 0, it2 = 0, sub = 0;
			const int sm_sz = sm.size();
			while (it1 < sm_sz && it2 < lr_sz) {
				const uint32 lit1 = sm[it1], lit2 = lr[it2];
				if (lit1 < lit2) it1++;
				else if (lit2 < lit1) it2++;
				else { sub++; it1++; it2++; }
			}
			if (sub == sm_sz) return true;
			return false;
		}

		_PFROST_D_ bool selfsub(const uint32& x, const uint32& fx, SCLAUSE& sm, SCLAUSE& lr)
		{
			assert(!sm.deleted());
			assert(!lr.deleted());
			assert(sm.size() > 1);
			assert(lr.size() > 1);
			assert(sm.size() <= lr.size());
			int it1 = 0, it2 = 0, sub = 0;
			const int sm_sz = sm.size(), lr_sz = lr.size();
			bool self = false;
			while (it1 < sm_sz && it2 < lr_sz) {
				const uint32 lit1 = sm[it1], lit2 = lr[it2];
				if (lit1 == fx) it1++;
				else if (lit2 == x) { self = true; it2++; }
				else if (lit1 < lit2) it1++;
				else if (lit2 < lit1) it2++;
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

		_PFROST_D_ bool selfsub(const uint32& x, const uint32& fx, SCLAUSE& sm, uint32* lr, const int& lr_sz)
		{
			assert(!sm.deleted());
			assert(sm.size() > 1);
			assert(lr_sz > 1);
			assert(sm.size() <= lr_sz);
			int it1 = 0, it2 = 0, sub = 0;
			const int sm_sz = sm.size();
			bool self = false;
			while (it1 < sm_sz && it2 < lr_sz) {
				const uint32 lit1 = sm[it1], lit2 = lr[it2];
				if (lit1 == fx) it1++;
				else if (lit2 == x) { self = true; it2++; }
				else if (lit1 < lit2) it1++;
				else if (lit2 < lit1) it2++;
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

		_PFROST_D_ void strengthen(cuVecU* units, SCLAUSE& c, uint32* sh_c, const uint32& self_lit)
		{
			assert(!c.deleted());
			assert(c.size() > 1);
			assert(self_lit > 1);
			int n = 0;
			for (int k = 0; k < c.size(); k++) {
				const uint32 lit = sh_c[k];
				if (lit != self_lit) {
					c[n] = lit, sh_c[n] = lit;
					n++;
				}
			}
			assert(n == c.size() - 1);
			assert(n <= SH_MAX_HSE_IN);
			assert(c.hasZero() < 0);
			assert(c.isSorted());
			c.pop();
			calcSig(c);
			if (c.size() == 1) units->push(*c);
			else if (n > 1 && c.learnt()) bumpShrunken(c);
#if SS_DBG
			if (c.size() == 1) printf("c  SS(%d): ", ABS(self_lit)), c.print();
#endif
		}

		_PFROST_D_ void strengthen(cuVecU* units, SCLAUSE& c, const uint32& self_lit)
		{
			assert(!c.deleted());
			assert(c.size() > 1);
			assert(self_lit > 1);
			int n = 0;
			for (int k = 0; k < c.size(); k++) {
				const uint32 lit = c[k];
				if (lit != self_lit) c[n++] = lit;
			}
			assert(n == c.size() - 1);
			assert(c.hasZero() < 0);
			assert(c.isSorted());
			c.pop();
			calcSig(c);
			if (n == 1) units->push(*c);
			else if (n > 1 && c.learnt()) bumpShrunken(c);
#if SS_DBG
			if (c.size() == 1) printf("c  SS(%d): ", ABS(self_lit)), c.print();
#endif
		}

		_PFROST_D_ void updateOL(CNF& cnf, OL& ol)
		{
			if (ol.empty()) return;
			S_REF *j = ol;
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
				if (subsuming.size() > 1 && sub(subsuming.sig(), cand.sig()) && sub(subsuming, cand)) {
					if (subsuming.learnt() && cand.original()) subsuming.set_status(ORIGINAL);
					cand.markDeleted();
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
				if (subsuming.size() > 1 && sub(subsuming.sig(), cand.sig()) && sub(subsuming, sh_cand, candsz)) {
					if (subsuming.learnt() && cand.original()) subsuming.set_status(ORIGINAL);
					cand.markDeleted();
					break;
				}
			}
		}

		_PFROST_D_ void selfsubsume(const uint32& x, const uint32& fx, CNF& cnf, OL& list, cuVecU* units, SCLAUSE& cand)
		{
			// try to strengthen 'cand' by removing 'x'
			const int candsz = cand.size();
			for (S_REF* j = list; j != list.end(); j++) {
				SCLAUSE& subsuming = cnf[*j];
				if (subsuming.size() > candsz) break;
				if (subsuming.deleted()) continue;
				if (subsuming.size() > 1 && selfsub(subsuming.sig(), cand.sig()) && selfsub(x, fx, subsuming, cand)) {
					strengthen(units, cand, x);
					cand.melt(); // mark for fast recongnition in ot update 
					break; // cannot strengthen "cand" anymore, 'x' already removed
				}
			}
		}

		_PFROST_D_ void selfsubsume(const uint32& x, const uint32& fx, CNF& cnf, OL& list, cuVecU* units, SCLAUSE& cand, uint32* sh_cand)
		{
			// try to strengthen 'cand' by removing 'x' using shared memory 
			const int candsz = cand.size();
			for (S_REF* j = list; j != list.end(); j++) {
				SCLAUSE& subsuming = cnf[*j];
				if (subsuming.size() > candsz) break;
				if (subsuming.deleted()) continue;
				if (subsuming.size() > 1 && selfsub(subsuming.sig(), cand.sig()) && selfsub(x, fx, subsuming, sh_cand, candsz)) {
					strengthen(units, cand, sh_cand, x);
					cand.melt();
					break;
				}
			}
		}

		_PFROST_D_ void subsume_x(const uint32& p, CNF& cnf, OL& poss, OL& negs, cuVecU* units, uint32* sh_c)
		{
			assert(checkMolten(cnf, poss, negs));
			const uint32 n = NEG(p);
			// positives vs negatives
			forall_occurs(poss, i) {
				SCLAUSE& pos = cnf[*i];
				if (pos.size() > HSE_MAX_CL_SIZE) break;
				if (pos.deleted()) continue;
				if (pos.size() > SH_MAX_HSE_IN) { // use global memory 
					selfsubsume(p, n, cnf, negs, units, pos);
					subsume(cnf, poss, i, pos);
				}
				else { // use shared memory
					uint32* sh_pos = sh_c;
					pos.shareTo(sh_pos);
					selfsubsume(p, n, cnf, negs, units, pos, sh_pos);
					subsume(cnf, poss, i, pos, sh_pos);
				}
			}
#if SS_DBG
			DBG_updateOL(p, cnf, poss); // discard (self)-subsumed clauses
#else
			updateOL(cnf, poss); // discard deleted or strengthened clauses from 'p' list
#endif
			// negatives vs positives
			forall_occurs(negs, i) {
				SCLAUSE& neg = cnf[*i];
				if (neg.size() > HSE_MAX_CL_SIZE) break;
				if (neg.deleted()) continue;
				if (neg.size() > SH_MAX_HSE_IN) {
					selfsubsume(n, p, cnf, poss, units, neg);
					subsume(cnf, negs, i, neg);
				}
				else { // use shared memory
					uint32* sh_neg = sh_c;
					neg.shareTo(sh_neg);
					selfsubsume(n, p, cnf, poss, units, neg, sh_neg);
					subsume(cnf, negs, i, neg, sh_neg);
				}
			}
#if SS_DBG
			DBG_updateOL(p, cnf, negs); 
#else
			updateOL(cnf, negs); 
#endif
			assert(checkMolten(cnf, poss, negs));
			assert(checkDeleted(cnf, poss, negs));
		}
		
	} // sigma namespace
} // parafrost namespace

#endif