/***********************************************************************[pfhse.h]
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

#include "pfdevice.cuh"

namespace pFROST {

	namespace SIGmA {

#define SS_DBG 0

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
			assert(c.hasZero() < 0);
			assert(c.isSorted());
			c.set_sig(sig);
			c.pop();
			if (c.size() == 1) units->push(*c);
#if SS_DBG
			if (c.size() == 1) printf("c | SS(%d): ", ABS(self_lit)), c.print();
#endif
		}

		_PFROST_D_ void updateOL(CNF& cnf, OL& ol)
		{
			if (ol.size() == 0) return;
			uint32* r = ol, *rend = ol.end();
			while (r != rend) {
				SCLAUSE& c = cnf[*r];
				if (c.deleted()) *r = *--rend;
				else if (c.molten()) c.freeze(), *r = *--rend;
				else r++;
			}
			assert(r == rend);
			assert(rend >= ol);
			ol.resize(rend - ol);
		}

		_PFROST_D_ void self_sub_x(const uint32& x, CNF& cnf, OL& poss, OL& negs, cuVecU* units, uint32* sh_c)
		{
			// positives vs negatives
#pragma unroll
			for (uint32* i = poss; i != poss.end(); i++) {
				SCLAUSE& pos = cnf[*i];
				if (pos.deleted()) continue;
				// self-subsumption check
				if (pos.size() > SH_MAX_HSE_IN) { // use global memory 
#pragma unroll
					for (uint32* j = negs; j != negs.end(); j++) {
						SCLAUSE& neg = cnf[*j];
						if (neg.deleted() || neg.learnt()) continue;
						if (neg.size() > 1 && neg.size() < pos.size() &&
							selfSub(neg.sig(), pos.sig()) && selfSub(x, neg, pos)) {
							strengthen(units, pos, x);
							pos.melt(); // mark for fast recongnition in ot update 
						}
						else if (neg.size() > 1 && neg.size() == pos.size() && pos.original() &&
							selfSub(neg.sig(), pos.sig()) && selfSub(x, neg, pos)) {
							strengthen(units, pos, x);
							pos.melt();
							neg.markDeleted();
						}
					}
					// subsumption check
					for (uint32* j = poss; j != poss.end(); j++) {
						SCLAUSE& sm_c = cnf[*j];
						if (sm_c.original() && sm_c.size() < pos.size() &&
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
					for (uint32* j = negs; j != negs.end(); j++) {
						SCLAUSE& neg = cnf[*j];
						if (neg.deleted() || neg.learnt()) continue;
						if (neg.size() > 1 && neg.size() < pos.size() &&
							selfSub(neg.sig(), pos.sig()) && selfSub(x, neg, sh_pos, pos.size())) {
							strengthen(units, pos, sh_pos, x);
							pos.melt();
						}
						else if (neg.size() > 1 && neg.size() == pos.size() && pos.original() &&
							selfSub(neg.sig(), pos.sig()) && selfSub(x, neg, sh_pos, pos.size())) {
							strengthen(units, pos, sh_pos, x);
							pos.melt();
							neg.markDeleted();
						}
					}
					for (uint32* j = poss; j != poss.end(); j++) {
						SCLAUSE& sm_c = cnf[*j];
						assert(isCoherent(x, pos, sh_pos, pos.size()));
						if (sm_c.original() && sm_c.size() < pos.size() &&
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
			for (uint32* i = negs; i != negs.end(); i++) {
				SCLAUSE& neg = cnf[*i];
				if (neg.deleted()) continue;
				if (neg.size() > SH_MAX_HSE_IN) {
					// self-subsumption check
#pragma unroll
					for (uint32* j = poss; j != poss.end(); j++) {
						SCLAUSE& pos = cnf[*j];
						if (pos.deleted() || pos.learnt()) continue;
						if (pos.size() > 1 && pos.size() < neg.size() &&
							selfSub(pos.sig(), neg.sig()) && selfSub(NEG(x), pos, neg)) {
							strengthen(units, neg, NEG(x));
							neg.melt();
						}
					}
					// subsumption check
					for (uint32* j = negs; j < negs.end(); j++) {
						SCLAUSE& sm_c = cnf[*j];
						if (sm_c.original() && sm_c.size() < neg.size() &&
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
					for (uint32* j = poss; j != poss.end(); j++) {
						SCLAUSE& pos = cnf[*j];
						if (pos.deleted() || pos.learnt()) continue;
						if (pos.size() > 1 && pos.size() < neg.size() &&
							selfSub(pos.sig(), neg.sig()) &&
							selfSub(NEG(x), pos, sh_neg, neg.size())) {
							strengthen(units, neg, sh_neg, NEG(x));
							neg.melt();
						}
					}
					// subsumption check
					for (uint32* j = negs; j < negs.end(); j++) {
						SCLAUSE& sm_c = cnf[*j];
						assert(isCoherent(x, neg, sh_neg, neg.size()));
						if (sm_c.original() && sm_c.size() < neg.size() &&
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
	}
}

#endif