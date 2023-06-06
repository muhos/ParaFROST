/***********************************************************************[mdm.cpp]
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

#include "solve.h" 
#include "mdmrank.h"
#include "mdmlimit.h"
#include "mdmassign.h"

using namespace ParaFROST;

inline bool	Solver::verifyMDM()
{
	for (uint32 i = sp->propagated; i < trail.size(); i++) {
		uint32 v = ABS(trail[i]);
		if (sp->frozen[v]) {
			PFLOG0("");
			PFLOGEN("decision(%d) is elected and frozen", v);
			printWatched(v);
			return false;
		}
	}
	return true;
}

inline bool	Solver::verifySeen()
{
	for (uint32 v = 0; v <= inf.maxVar; v++) {
		if (sp->seen[v]) {
			PFLOG0("");
			PFLOGEN("seen(%d) is not unseen", v);
			printWatched(v);
			return false;
		}
	}
	return true;
}

inline bool Solver::valid(const LIT_ST* values, WL& ws)
{
	forall_watches(ws, i) {
		const WATCH w = *i;
		// clause satisfied
		if (values[w.imp] > 0) continue;
		// if 'w.imp' not satisfied then it's an implication of 'cand'
		if (w.binary())
			return false;
		// there cannot be falsified literal as watched,
		// so validating starts from 'c + 2'
		CLAUSE& c = cm[w.ref];
		assert(c.size() > 2);
		bool satisfied = false, unAssigned = false;
		uint32* k = c + 2, * cend = c.end();
		while (k != cend && !satisfied && !unAssigned) {
			CHECKLIT(*k);
			const LIT_ST val = values[*k];
			if (UNASSIGNED(val)) unAssigned = true;
			else if (val) satisfied = true;
			k++;
		}
		if (!satisfied && !unAssigned)
			return false;
	}
	return true;
}

inline bool Solver::depFreeze(const uint32& cand, const LIT_ST* values, LIT_ST* frozen, uint32*& stack, WL& ws)
{
	forall_watches(ws, i) {
		const WATCH w = *i;
		if (values[w.imp] > 0) continue;
		CLAUSE& c = cm[w.ref];
		uint32* lits = c.data();
		uint32 othervar = ABS(lits[0]) ^ ABS(lits[1]) ^ cand;
		if (sp->seen[othervar]) return false;
		if (!frozen[othervar]) {
			frozen[othervar] = 1;
			assert(stack < sp->tmpstack + inf.maxVar);
			*stack++ = othervar;
		}
		uint32* k = lits + 2, * cend = c.end();
		while (k != cend) {
			othervar = ABS(*k++);
			if (sp->seen[othervar]) return false;
			if (!frozen[othervar]) {
				frozen[othervar] = 1;
				assert(stack < sp->tmpstack + inf.maxVar);
				*stack++ = othervar;
			}
		}
	}
	return true;
}

inline void Solver::MDMAssume(const LIT_ST* values, LIT_ST* frozen, uint32*& tail, bool& skip)
{
	assert(sp->stacktail == sp->tmpstack);
	int level = DL();
	while (level < assumptions.size()) {
		const uint32 a = assumptions[level];
		CHECKLIT(a);
		assert(ifrozen[ABS(a)]);
		const uint32 cand = ABS(a);
		const LIT_ST val = values[a];
		if (UNASSIGNED(val)) {
			level++;
			if (!depFreeze(cand, values, frozen, tail, wt[a]))
				skip = true;
			enqueueDecision(a);
			sp->seen[cand] = 1;
		}
		else if (val) incDL(), level = DL();
		else {
			ianalyze(FLIP(a));
			cnfstate = UNSAT;
			clearMDM();
			return;
		}
	}
	stats.decisions.massumed += trail.size() - sp->propagated;
}

void Solver::MDMInit()
{
	if (!last.mdm.rounds) return;

	assert(inf.unassigned);
	assert(sp->propagated == trail.size());
	assert(conflict == NOREF);

	if (opts.mdm_walk_en) {
		stats.mdm.walks++;
		walk();
	}

	// check if formula is solved by walk strategy
	if (!UNSOLVED(cnfstate)) return;

	stats.mdm.calls++;

	PFLOG2(2, " MDM %d: electing decisions at decaying round %d..", stats.mdm.calls, last.mdm.rounds);

	eligible_initial;

	mdm_prefetch(values, states, frozen, tail);

	bool skip = false;
	if (assumptions.size())
		MDMAssume(values, frozen, tail, skip);

	if (!skip) {
		forall_vector(uint32, eligible, evar) {
			const uint32 cand = *evar;
			CHECKVAR(cand);
			if (frozen[cand] || states[cand].state || iassumed(cand))
				continue;
			const uint32 dec = V2DEC(cand, sp->psaved[cand]);
			mdm_assign(cand, dec);
		}
	}

	mdm_update;

	if (opts.mdm_vsids_pumps || opts.mdm_vmtf_pumps)
		pumpFrozen();

	clearMDM(), eligible.clear(true), occurs.clear(true);

	printStats(1, 'm', CMDM);
}

void Solver::MDM()
{
	assert(inf.unassigned);
	assert(sp->propagated == trail.size());
	assert(conflict == NOREF);
	assert(UNSOLVED(cnfstate));

	stats.mdm.calls++;
	PFLOG2(2, " MDM %d: electing decisions at decaying round %d..", stats.mdm.calls, last.mdm.rounds);

	if (vsidsEnabled()) {
		eligible_heap(vsids);
	}
	else {
		eligible_queue(vmtf);
	}

	assert(eligible.size() >= 1);

	const bool firstround = last.mdm.rounds == opts.mdm_rounds;
	if (opts.mdm_walk_en && !DL() && firstround) {
		stats.mdm.walks++;
		walk();
	}

	// check if formula is solved by walk strategy
	if (!UNSOLVED(cnfstate)) {
		eligible.clear();
		return;
	}

	mdm_prefetch(values, states, frozen, tail);

	bool skip = false;
	if (assumptions.size())
		MDMAssume(values, frozen, tail, skip);

	if (!skip) {
		const bool targeting = useTarget();
		forall_vector(uint32, eligible, evar) {
			const uint32 cand = *evar;
			CHECKVAR(cand);
			if (frozen[cand] || states[cand].state || iassumed(cand)) continue;
			uint32 dec = makeAssign(cand, targeting);
			mdm_assign(cand, dec);
		}
	}

	mdm_update;

	if (opts.mdm_vsids_pumps || opts.mdm_vmtf_pumps)
		pumpFrozen();

	clearMDM();

	printStats(firstround, 'm', CMDM);
}