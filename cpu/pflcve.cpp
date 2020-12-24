/***********************************************************************[pflcve.cpp]
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

#include "pfsimp.h"

using namespace pFROST;
using namespace SIGmA;

void ParaFROST::varReorder()
{
	PFLOGN2(2, " Finding eligible variables for LCVE..");
	if (opts.profile_simp) timer.pstart();
	eligible.resize(inf.maxVar);
	occurs.resize(inf.maxVar + 1);
	assert(!scnf.empty());
	histSimp(scnf, true);
	uint32* scores = sp->tmp_stack;
	for (uint32 v = 1; v <= inf.maxVar; v++) eligible[v - 1] = v, scores[v] = rscore(v);
	rSort(eligible, LCV_CMP(scores), LCV_RANK(scores));
	if (opts.profile_simp) timer.pstop(), timer.vo += timer.pcpuTime();
	PFLDONE(2, 5);
	if (verbose >= 3) {
		PFLOG0(" Eligible variables:");
		for (uint32 i = 0; i < eligible.size(); i++) {
			uint32 v = eligible[i];
			PFLOG1("  e[%d]->(v: %d, p: %d, n: %d, s: %d)", i, v, occurs[v].ps, occurs[v].ns, scores[v]);
		}
	}
}

bool ParaFROST::LCVE()
{
	// reorder variables
	varReorder();
	// extended LCVE
	PFLOGN2(2, " Electing variables in phase-%d..", phase);
	PVs.clear();
	for (uint32 i = 0; i < eligible.size(); i++) {
		uint32 cand = eligible[i];
		assert(cand && cand <= inf.maxVar);
		if (sp->vstate[cand] == FROZEN || sp->vstate[cand] == MELTED) continue;
		if (sp->frozen[cand]) continue;
		uint32 p = V2L(cand), n = NEG(p);
		// occur. lists may have deleted occurs
		uint32 poss_sz = (uint32)ot[p].size(), negs_sz = (uint32)ot[n].size();
		assert(poss_sz >= occurs[cand].ps);
		assert(negs_sz >= occurs[cand].ns);
		if (occurs[cand].ps == 0 && occurs[cand].ns == 0) continue;
		uint32 pos_temp = opts.mu_pos << mu_inc, neg_temp = opts.mu_neg << mu_inc;
		if (poss_sz >= pos_temp && negs_sz >= neg_temp) break;
		assert(sp->vstate[cand] == ACTIVE);
		PVs.push(cand);
		depFreeze(ot[p], cand, pos_temp, neg_temp);
		depFreeze(ot[n], cand, pos_temp, neg_temp);
	}
	assert(verifyLCVE());
	PFLENDING(2, 5, "(%d elected)", PVs.size());
	if (verbose >= 3) { PFLOGN0(" PLCVs "); printVars(PVs, PVs.size(), 'v'); }
	memset(sp->frozen, 0, inf.maxVar + 1ULL);
	if (PVs.size() < opts.lcve_min) {
		if (verbose > 1) PFLOGW("parallel variables not enough -> skip BVE");
		return false;
	}
	return true;
}

inline void ParaFROST::depFreeze(const OL& ol, const uint32& cand, const uint32& p_temp, const uint32& n_temp)
{
	for (int i = 0; i < ol.size(); i++) {
		SCLAUSE& c = *ol[i];
		if (c.deleted()) continue;
		for (int k = 0; k < c.size(); k++) {
			uint32 v = ABS(c[k]);
			if (v != cand && (occurs[v].ps < p_temp || occurs[v].ns < n_temp)) sp->frozen[v] = 1;
		}
	}
}