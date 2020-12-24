/***********************************************************************[pfreduce.cpp]
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

#include "pfsolve.h"
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

struct LEARNT_SZ_CMP {
	const CMM& cm;
	LEARNT_SZ_CMP(const CMM& _cm) : cm(_cm) {}
	bool operator () (const C_REF& a, const C_REF& b) const {
		return cm[a].size() > cm[b].size();
	}
};

inline void	ParaFROST::reduceWeight(double& val) {
	double orgSize = orgs.size();
	if (orgSize > 1e5) {
		val *= log(orgSize / 1e4) / log(10);
		if (val < 1.0) val = 1.0;
	}
}

bool ParaFROST::BCPChronoRoot() {
	if (!opts.chrono_en || !DL()) return true;
	for (uint32 i = dlevels[1]; i < trail.size(); i++) {
		uint32 lit = trail[i];
		assert(!unassigned(lit));
		if (!l2dl(lit)) {
			PFLOG2(2, " Found root unit(%d) due to chronological backtracking", l2i(lit));
			backtrack();
			if (BCP()) { cnfstate = UNSAT; return false; }
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
	if (!BCPChronoRoot()) return;
	bool shrunken = shrink();
	protectReasons();
	reduceLearnts();
	recycle();
	unprotectReasons();
	double current_inc = opts.reduce_inc * (stats.reduces + 1);
	reduceWeight(current_inc);
	lrn.lastreduce = nConflicts;
	lrn.reduce_conf_max = nConflicts + current_inc;
	PFLOG2(2, " reduce limit increased to %lld conflicts by a weight %.2f", lrn.reduce_conf_max, current_inc);
	if (shrunken && canMap()) map(); // "recycle" must be called beforehand
}

void ParaFROST::reduceLearnts(const bool& sizeonly)
{
	assert(reduced.empty());
	reduced.reserve(learnts.size());
	// keep recently used, tier1, and reason learnts otherwise reduce
	for (C_REF* r = learnts; r != learnts.end(); r++) {
		CLAUSE& c = cm[*r];
		if (c.deleted()) continue;
		assert(c.learnt());
		if (c.reason()) continue;
		if (c.usage()) { c.warm(); continue; }
		if (c.keep()) continue;
		assert(c.size() > 2);
		reduced.push(*r);
	}
	if (reduced.size()) {
		uint32 pivot = opts.reduce_perc * reduced.size();
		PFLOGN2(2, " Reducing learnt database up to (%d clauses)..", pivot);
		if (sizeonly) std::stable_sort(reduced.data(), reduced.end(), LEARNT_SZ_CMP(cm));
		else std::stable_sort(reduced.data(), reduced.end(), LEARNT_CMP(cm));
		// remove unlucky learnts from database
		C_REF* h, * t = reduced + pivot;
		for (h = reduced; h != t; h++) {
			assert(cm[*h].learnt());
			assert(!cm[*h].reason());
			assert(!cm[*h].keep());
			assert(cm[*h].lbd() > opts.lbd_tier1);
			assert(cm[*h].size() > 2);
			removeClause(*h);
		}
		lrn.keptsize = 0, lrn.keptlbd = 0;
		C_REF* e = reduced.end();
		for (h = t; h != e; h++) {
			CLAUSE& c = cm[*h];
			if (c.lbd() > lrn.keptlbd) lrn.keptlbd = c.lbd();
			if (c.size() > lrn.keptsize) lrn.keptsize = c.size();
		}
		reduced.clear();
		PFLENDING(2, 5, "(kept lbd: %d, size: %d)", lrn.keptlbd, lrn.keptsize);
	}
}

void ParaFROST::reduceTop(const bool& sizeonly) {
	if (learnts.empty()) return;
	assert(sp->propagated == trail.size());
	assert(conflict == NOREF);
	assert(cnfstate == UNSOLVED);
	stats.reduces++;
	reduceLearnts(sizeonly);
	filter(learnts);
	double current_inc = opts.reduce_inc * (stats.reduces + 1);
	reduceWeight(current_inc);
	lrn.lastreduce = nConflicts;
	lrn.reduce_conf_max = nConflicts + current_inc;
	PFLOG2(2, " reduce limit increased to %lld conflicts by a weight %.2f", lrn.reduce_conf_max, current_inc);
}