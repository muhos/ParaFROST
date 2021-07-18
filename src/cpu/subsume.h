/***********************************************************************[subsume.h]
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

#ifndef __SUB_
#define __SUB_

#include "simplify.h" 
using namespace pFROST;

#define HSE_MAX_CL_SIZE 1000

inline void updateOL(OL& ol)
{
	if (ol.empty()) return;
	S_REF* j = ol;
	forall_occurs(ol, i) {
		SCLAUSE& c = **i;
		if (c.molten()) c.freeze();
		else if (!c.deleted()) *j++ = *i;
	}
	ol.resize(int(j - ol));
}

inline bool sub(SCLAUSE& sm, SCLAUSE& lr)
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

inline bool selfsub(const uint32& x, const uint32& fx, SCLAUSE& sm, SCLAUSE& lr)
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

inline bool subsume(OL& list, S_REF* end, SCLAUSE& cand)
{
	const int candsz = cand.size();
	for (S_REF* j = list; j != end; j++) {
		SCLAUSE& subsuming = **j;
		if (subsuming.deleted()) continue;
		if (cand.molten() && subsuming.size() > candsz) continue;
		if (subsuming.size() > 1 && sub(subsuming.sig(), cand.sig()) && sub(subsuming, cand)) {
			if (subsuming.learnt() && cand.original()) subsuming.set_status(ORIGINAL);
			pfrost->removeClause(cand);
			PFLCLAUSE(4, cand, " Clause ");
			PFLCLAUSE(4, subsuming, " Subsumed by ");
			return true;
		}
	}
	return false;
}

inline bool selfsubsume(const uint32& x, const uint32& fx, OL& list, SCLAUSE& cand)
{
	// try to strengthen 'cand' by removing 'x'
	const int candsz = cand.size();
	const uint32 candsig = cand.sig();
	forall_occurs(list, j) {
		SCLAUSE& subsuming = **j;
		const int subsize = subsuming.size();
		if (subsize > candsz) break;
		if (subsuming.deleted()) continue;
		if (subsize > 1 && selfsub(subsuming.sig(), candsig) && selfsub(x, fx, subsuming, cand)) {
			PFLCLAUSE(4, cand, " Clause ");
			pfrost->strengthen(cand, x);
			cand.melt(); // mark for fast recongnition in ot update 
			PFLCLAUSE(4, subsuming, " Strengthened by ");
			return true; // cannot strengthen "cand" anymore, 'x' already removed
		}
	}
	return false;
}

inline void self_sub_x(const uint32& p, OL& poss, OL& negs, SUBSTATS& substats)
{
	CHECKLIT(p);
	assert(checkMolten(poss, negs));
	const uint32 n = NEG(p);
	// positives vs negatives
	forall_occurs(poss, i) {
		SCLAUSE& pos = **i;
		if (pos.size() > HSE_MAX_CL_SIZE) break;
		if (pos.deleted()) continue;
#ifdef STATISTICS
		if (selfsubsume(p, n, negs, pos)) substats.strengthened++;
		if (subsume(poss, i, pos)) substats.subsumed++;
#else 
		selfsubsume(p, n, negs, pos);
		subsume(poss, i, pos);
#endif
	}
	updateOL(poss);
	// negatives vs positives
	forall_occurs(negs, i) {
		SCLAUSE& neg = **i;
		if (neg.size() > HSE_MAX_CL_SIZE) break;
		if (neg.deleted()) continue;
#ifdef STATISTICS
		if (selfsubsume(n, p, poss, neg)) substats.strengthened++;
		if (subsume(negs, i, neg)) substats.subsumed++;
#else
		selfsubsume(n, p, poss, neg);
		subsume(negs, i, neg);
#endif
	}
	updateOL(negs);
	assert(checkMolten(poss, negs));
	assert(checkDeleted(poss, negs));
}

#endif