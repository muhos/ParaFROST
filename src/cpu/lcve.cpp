/***********************************************************************[lcve.cpp]
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

#include "simplify.h"

using namespace pFROST;

void ParaFROST::varReorder()
{
	PFLOGN2(2, " Finding eligible variables for LCVE..");
	if (opts.profile_simp) timer.pstart();
	eligible.resize(inf.maxVar);
	occurs.resize(inf.maxVar + 1);
	assert(!scnf.empty());
	histSimp(scnf, true);
	uint32* scores = sp->tmpstack;
	forall_variables(v) {
		eligible[v - 1] = v, scores[v] = prescore(v);
	}
	rSort(eligible, LCV_CMP(scores), LCV_RANK(scores));
	if (opts.profile_simp) timer.pstop(), timer.vo += timer.pcpuTime();
	PFLDONE(2, 5);
	if (verbose > 3) {
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
	sp->stacktail = sp->tmpstack;
	for (uint32 i = 0; i < eligible.size(); i++) {
		const uint32 cand = eligible[i];
		CHECKVAR(cand);
		if (iassumed(cand)) continue;
		if (sp->vstate[cand].state) continue;
		if (sp->frozen[cand]) continue;
		const uint32 p = V2L(cand), n = NEG(p);
		const uint32 poss_sz = (uint32)ot[p].size(), negs_sz = (uint32)ot[n].size();
		assert(poss_sz >= occurs[cand].ps);
		assert(negs_sz >= occurs[cand].ns);
		if (occurs[cand].ps == 0 && occurs[cand].ns == 0) continue;
		const uint32 pos_temp = opts.mu_pos << mu_inc, neg_temp = opts.mu_neg << mu_inc;
		if (poss_sz >= pos_temp && negs_sz >= neg_temp) break;
		assert(!sp->vstate[cand].state);
		PVs.push(cand);
		depFreeze(ot[p], cand, pos_temp, neg_temp);
		depFreeze(ot[n], cand, pos_temp, neg_temp);
	}
	assert(verifyLCVE());
	clearFrozen();
	PFLENDING(2, 5, "(%d elected)", PVs.size());
	if (verbose > 3) { PFLOGN0(" PLCVs "); printVars(PVs, PVs.size(), 'v'); }
	if (PVs.size() < opts.lcve_min) {
		if (verbose > 1) PFLOGW("parallel variables not enough -> skip BVE");
		return false;
	}
	return true;
}

inline void ParaFROST::depFreeze(OL& ol, const uint32& cand, const uint32& p_temp, const uint32& n_temp)
{
	LIT_ST* frozen = sp->frozen;
	uint32*& frozen_stack = sp->stacktail;
	forall_occurs(ol, i) {
		S_REF c = *i;
		if (c->deleted()) continue;
		forall_clause((*c), k) {
			const uint32 v = ABS(*k);
			CHECKVAR(v);
			if (!frozen[v] && NEQUAL(v, cand) && (occurs[v].ps < p_temp || occurs[v].ns < n_temp)) {
				frozen[v] = 1;
				assert(frozen_stack < sp->tmpstack + inf.maxVar);
				*frozen_stack++ = v;
			}
		}
	}
}