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

#ifndef __HSE_
#define __HSE_

#include "pfsolve.h" 
using namespace pFROST;

namespace SIGmA {

#define HSE_DBG 0
#define SUB_DBG 0
#define HSE_MAX_CL_SIZE 1000

	struct CNF_CMP_KEY {
		bool operator () (S_REF x, S_REF y) {
			if (x->size() != y->size()) return x->size() < y->size();
			else if (x->lit(0) != y->lit(0)) return x->lit(0) < y->lit(0);
			else if (x->lit(1) != y->lit(1)) return x->lit(1) < y->lit(1);
			else if (x->size() > 2 && x->back() != y->back()) return x->back() < y->back();
			else if (x->sig() != y->sig()) return x->sig() < y->sig();
			return x < y;
		}
	};
	struct CNF_CMP_SZ {
		bool operator () (S_REF x, S_REF y) {
			return x->size() < y->size();
		}
	};

	inline bool checkDeleted(OL& poss, OL& negs)
	{
		for (int i = 0; i < poss.size(); i++)
			if (poss[i]->deleted()) return false;
		for (int i = 0; i < negs.size(); i++)
			if (negs[i]->deleted()) return false;
		return true;
	}

	inline bool isEqual(const SCLAUSE& c1, const Lits_t& c2)
	{
		assert(!c1.deleted());
		assert(c1.size() == c2.size());
		int it = 0;
		while (it < c2.size()) {
			if (c1[it] != c2[it]) return false;
			else it++;
		}
		return true;
	}

	inline void updateOL(OL& ol)
	{
		if (ol.empty()) return;
		S_REF* i, * j, * rend = ol.end();
		for (i = ol, j = i; i != rend; i++) {
			SCLAUSE& c = **i;
			if (c.molten()) c.freeze();
			else if (!c.deleted()) *j++ = *i;
		}
		ol.shrink(int(rend - j));
	}

	inline bool sub(const uint32& A, const uint32& B) { return !(A & ~B); }

	inline bool selfsub(const uint32& A, const uint32& B)
	{
		uint32 B_tmp = B | ((B & 0xAAAAAAAAUL) >> 1) | ((B & 0x55555555UL) << 1);
		return !(A & ~B_tmp);
	}

	inline bool sub(const S_REF sm, const S_REF lr)
	{
		assert(!sm->deleted());
		assert(!lr->deleted());
		assert(sm->size() > 1);
		assert(lr->size() > 1);
		assert(sm->size() <= lr->size());
		int it1 = 0, it2 = 0, sub = 0;
		while (it1 < sm->size() && it2 < lr->size()) {
			if (sm->lit(it1) < lr->lit(it2)) it1++;
			else if (lr->lit(it2) < sm->lit(it1)) it2++;
			else { sub++; it1++; it2++; }
		}
		if (sub == sm->size())
			return true;
		return false;
	}

	inline bool sub(const Lits_t& sm, const SCLAUSE& lr)
	{
		assert(!lr.deleted());
		assert(sm.size() > 1);
		assert(lr.size() > 1);
		assert(sm.size() <= lr.size());
		int it1 = 0, it2 = 0, sub = 0;
		while (it1 < sm.size() && it2 < lr.size()) {
			if (sm[it1] < lr[it2]) it1++;
			else if (lr[it2] < sm[it1]) it2++;
			else { sub++; it1++; it2++; }
		}
		if (sub == sm.size())
			return true;
		return false;
	}

	inline bool selfsub(const uint32& x, const S_REF sm, const S_REF lr)
	{
		assert(!sm->deleted());
		assert(!lr->deleted());
		assert(sm->size() > 1);
		assert(lr->size() > 1);
		assert(sm->size() <= lr->size());
		int it1 = 0, it2 = 0, sub = 0;
		bool self = false;
		while (it1 < sm->size() && it2 < lr->size()) {
			if (sm->lit(it1) == FLIP(x)) it1++;
			else if (lr->lit(it2) == x) { self = true; it2++; }
			else if (sm->lit(it1) < lr->lit(it2)) it1++;
			else if (lr->lit(it2) < sm->lit(it1)) it2++;
			else { sub++; it1++; it2++; }
		}
		if ((sub + 1) == sm->size()) {
			if (self) return true;
			else {
				while (it2 < lr->size()) {
					if (lr->lit(it2) == x) return true;
					it2++;
				}
			}
		}
		return false;
	}

	inline void self_sub_x(const uint32& x, OL& poss, OL& negs)
	{
		// positives vs negatives
		for (int i = 0; i < poss.size(); i++) {
			S_REF pos = poss[i];
			if (pos->size() > HSE_MAX_CL_SIZE) break;
			if (pos->deleted()) continue;
			// self-subsumption check
			for (int j = 0; j < negs.size(); j++) {
				S_REF neg = negs[j];
				if (neg->size() > pos->size()) break;
				if (neg->deleted()) continue;
				if (neg->size() > 1 && selfsub(neg->sig(), pos->sig()) && selfsub(x, neg, pos)) {
#if HSE_DBG
					PFLCLAUSE(1, (*pos), " Clause ");
					PFLCLAUSE(1, (*neg), " Strengthened by ");
#endif 
					pfrost->strengthen(pos, x);
					pos->melt(); // mark for fast recongnition in ot update 
					break; // cannot strengthen "pos" anymore, 'x' already removed
				}
			}
			// subsumption check
			for (int j = 0; j < i; j++) {
				S_REF sm_c = poss[j];
				if (sm_c->deleted()) continue;
				if (pos->molten() && sm_c->size() > pos->size()) continue;
				if (sm_c->size() > 1 && sub(sm_c->sig(), pos->sig()) && sub(sm_c, pos)) {
					if (sm_c->learnt() && pos->original()) sm_c->set_status(ORIGINAL);
#if HSE_DBG
					PFLCLAUSE(1, (*pos), " Clause ");
					PFLCLAUSE(1, (*sm_c), " Subsumed by ");
#endif 
					pfrost->removeClause(pos);
					break;
				}
			}
		}
		updateOL(poss);
		// negatives vs positives
		for (int i = 0; i < negs.size(); i++) {
			S_REF neg = negs[i];
			if (neg->size() > HSE_MAX_CL_SIZE) break;
			if (neg->deleted()) continue;
			// self-subsumption check
			for (int j = 0; j < poss.size(); j++) {
				S_REF pos = poss[j];
				if (pos->size() >= neg->size()) break;
				if (pos->deleted()) continue;
				if (pos->size() > 1 && selfsub(pos->sig(), neg->sig()) && selfsub(NEG(x), pos, neg)) {
#if HSE_DBG
					PFLCLAUSE(1, (*neg), " Clause ");
					PFLCLAUSE(1, (*pos), " Strengthened by ");
#endif 
					pfrost->strengthen(neg, NEG(x));
					neg->melt();
					break;
				}
			}
			// subsumption check
			for (int j = 0; j < i; j++) {
				S_REF sm_c = negs[j];
				if (sm_c->deleted()) continue;
				if (neg->molten() && sm_c->size() > neg->size()) continue;
				if (sm_c->size() > 1 && sub(sm_c->sig(), neg->sig()) && sub(sm_c, neg)) {
					if (sm_c->learnt() && neg->original()) sm_c->set_status(ORIGINAL);
#if HSE_DBG
					PFLCLAUSE(1, (*neg), " Clause ");
					PFLCLAUSE(1, (*sm_c), " Subsumed by ");
#endif 
					pfrost->removeClause(neg);
					break;
				}
			}
		}
		updateOL(negs);
		assert(checkDeleted(poss, negs));
	}


}

#endif