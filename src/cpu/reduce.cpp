/***********************************************************************[reduce.cpp]
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

//#include <execution>
#include "solve.h"
using namespace pFROST;

struct LEARNT_CMP {
	const CMM& cm;
	LEARNT_CMP(const CMM& _cm) : cm(_cm) {}
	bool operator () (const C_REF& a, const C_REF& b) const {
		const CLAUSE& x = cm[a], & y = cm[b];
		if (x.lbd() > y.lbd()) return true;
		if (x.lbd() < y.lbd()) return false;
		return x.size() > y.size();
	}
};

bool ParaFROST::chronoHasRoot() {
	if (!opts.chrono_en || !DL()) return true;
	for (uint32 i = dlevels[1]; i < trail.size(); i++) {
		uint32 lit = trail[i];
		assert(!unassigned(lit));
		if (!l2dl(lit)) {
			PFLOG2(2, " Found root unit(%d) due to chronological backtracking", l2i(lit));
			backtrack();
			if (BCP()) { learnEmpty(); return false; }
			return true;
		}
	}
	return true;	
}

void ParaFROST::reduce()
{
	assert(sp->propagated == trail.size());
	assert(conflict == NOREF);
	assert(cnfstate == UNSOLVED);
	assert(learnts.size());
	stats.reduces++;
	if (!chronoHasRoot()) return;
	if (canSubsume()) subsume();
	bool shrunken = shrink();
	protectReasons();
	reduceLearnts();
	recycle();
	unprotectReasons();
	INCREASE_LIMIT(this, reduce, stats.reduces, nbylogn, false);
	if (shrunken && canMap()) map(); // "recycle" must be called beforehand
}

void ParaFROST::reduceLearnts()
{
	assert(reduced.empty());
	assert(learnts.size());
	reduced.reserve(learnts.size());
	C_REF* end = learnts.end();
	for (C_REF* i = learnts; i != end; i++) {
		const C_REF r = *i;
		if (cm.deleted(r)) continue;
		CLAUSE& c = cm[r];
		assert(!c.deleted());
		assert(c.learnt());
		if (c.reason()) continue;
		if (c.hyper()) {
			assert(c.size() <= 3);
			if (c.usage()) c.warm();
			else {
				removeClause(c, r);
				if (c.binary()) stats.binary.reduced++;
				else stats.ternary.reduced++;
			}
			continue;
		}
		if (c.keep()) continue;
		if (c.usage()) {
			c.warm();
			if (c.lbd() <= opts.lbd_tier2) continue;
		}
		assert(c.size() > 2);
		reduced.push(r);
	}
	if (reduced.size()) {
		C_REF pivot = opts.reduce_perc * reduced.size();
		PFLOGN2(2, " Reducing learnt database up to (%zd clauses)..", pivot);
		end = reduced.end();
		C_REF* head = reduced.data();
		//std::stable_sort(std::execution::par_unseq, reduced.data(), reduced.end(), LEARNT_CMP(cm));
		std::stable_sort(head, end, LEARNT_CMP(cm));
		// remove unlucky learnts from database
		C_REF *tail = reduced + pivot;
		for (; head != tail; head++) {
			const C_REF r = *head;
			assert(cm[r].learnt());
			assert(!cm[r].reason());
			assert(!cm[r].keep());
			assert(cm[r].lbd() > opts.lbd_tier1);
			assert(cm[r].size() > 2);
			removeClause(cm[r], r);
		}
		limit.keptsize = 0, limit.keptlbd = 0;
		for (head = tail; head != end; head++) {
			CLAUSE& c = cm[*head];
			if (c.lbd() > limit.keptlbd) limit.keptlbd = c.lbd();
			if (c.size() > limit.keptsize) limit.keptsize = c.size();
		}
		reduced.clear();
		PFLENDING(2, 5, "(kept lbd: %d, size: %d)", limit.keptlbd, limit.keptsize);
	}
}