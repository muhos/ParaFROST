/***********************************************************************[pfmdm.cpp]
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

#include "pfsort.h"
#include "pfsolve.h" 
using namespace pFROST;

void ParaFROST::varOrder()
{
	PFLOGN2(2, " Finding eligible decisions at initial round..");
	eligible.resize(inf.maxVar);
	occurs.resize(inf.maxVar + 1);
	histBins(true);
	hist(orgs);
	uint32* scores = sp->tmp_stack;
	for (uint32 v = 1; v <= inf.maxVar; v++) eligible[v - 1] = v, scores[v] = rscore(v);
	if (mcv_en) Sort(eligible, MCV_CMP(scores));
	else Sort(eligible, LCV_CMP(scores));
	PFLDONE(2, 5);
	if (verbose >= 3) {
		PFLOG0(" Eligible decisions:");
		for (uint32 i = 0; i < eligible.size(); i++) {
			uint32 v = eligible[i];
			PFLOG1("  e[%d]->(v: %d, p: %d, n: %d, s: %d)", i, v, occurs[v].ps, occurs[v].ns, scores[v]);
		}
	}
}

void ParaFROST::histBins(const bool& rst)
{
	if (wtBin.size() == 0) return;
	if (rst) { for (uint32 i = 0; i < occurs.size(); i++) occurs[i] = { 0 , 0 }; }
	if (sp->simplified) wtBin.recycle();
	for (uint32 v = 1; v <= inf.maxVar; v++) {
		uint32 p = v2l(v), n = neg(p);
		occurs[v].ps += wtBin[n].size(), occurs[v].ns += wtBin[p].size();
	}
}

void ParaFROST::pumpFrozenHeap(const uint32& lit)
{
	WL& ws = sp->simplified ? wt.getClean(lit) : wt[lit];
	if (ws.empty()) return;
	uint32 v = l2a(lit);
	assert(!sp->frozen[v]);
	double norm_act = (double)sp->level[v] / lrn.numMDs;
	for (WATCH* w = ws; w != ws.end(); w++) {
		CLAUSE& c = cm[w->ref];
		assert(!c.deleted());
		uint32 v0 = l2a(*c), frozen_v = (v0 == v) ? l2a(c[1]) : v0;
		if (varAct[frozen_v] == 0.0) varBumpHeap(frozen_v, norm_act);
	}
}

void ParaFROST::pumpFrozenQue(const uint32& lit)
{
	WL& ws = sp->simplified ? wt.getClean(lit) : wt[lit];
	if (ws.empty()) return;
	uint32 v = l2a(lit);
	assert(!sp->frozen[v]);
	for (WATCH* w = ws; w != ws.end(); w++) {
		CLAUSE& c = cm[w->ref];
		assert(!c.deleted());
		uint32 v0 = l2a(*c), frozen_v = (v0 == v) ? l2a(c[1]) : v0;
		if (sp->frozen[frozen_v]) {
			analyzed.push(frozen_v);
			sp->frozen[frozen_v] = 0;
		}
	}
}

void ParaFROST::pumpFrozen()
{
	if (lrn.numMDs == 0) return;
	assert(trail.size() > sp->propagated);
	if (verbose == 4) PFLOG1(" Pumping frozen variables..");
	else PFLOGN2(2, " Pumping frozen variables..");
	uint32* start = trail + sp->propagated;
	if (vsids() && mdm_vsids_pumps) {
		for (uint32* assign = trail.end() - 1; assign >= start; assign--)
			pumpFrozenHeap(*assign);
		mdm_vsids_pumps--;
	}
	else if (mdm_vmfq_pumps) { // VMFQ (requires special handling)
		assert(analyzed.empty());
		assert(inf.maxVar >= lrn.numMDs);
		analyzed.reserve(inf.maxVar - lrn.numMDs);
		for (uint32* assign = trail.end() - 1; assign >= start; assign--)
			pumpFrozenQue(*assign);
		for (uint32* v = analyzed.end() - 1; v >= analyzed; v--)
			varBumpQueueNU(*v);
		if (!sp->locked[*analyzed]) vQueue.update(*analyzed, varBumps[*analyzed]);
		analyzed.clear();
		mdm_vmfq_pumps--;
	}
	PFLDONE(2, 4);
	if (verbose >= 3 && vsids()) printHeap();
}

void ParaFROST::MDMInit()
{
	assert(!satisfied());
	assert(sp->propagated == trail.size());
	assert(conflict == NOREF);
	assert(cnfstate == UNSOLVED);
	if (!lrn.rounds) return;
	stats.mdm_calls++;
	varOrder(); // initial variable ordering
	if (verbose == 4) PFLOG1(" Electing decisions (rounds=%d)..", lrn.rounds);
	else PFLOGN2(2, " Electing decisions (rounds=%d)..", lrn.rounds);
	for (uint32 i = 0; i < inf.maxVar; i++) {
		uint32 cand = eligible[i];
		assert(cand && cand <= inf.maxVar);
		if (sp->frozen[cand] || sp->locked[cand]) continue;
		uint32 p = v2l(cand), n = neg(p);
		uint32 dec = wt[p].size() >= wt[n].size() ? p : n;
		if (wtBin[dec].size()) continue;
		if (depFreeze(wt[dec], cand)) {
			incDL();
			enqueue(dec, DL());
			sp->seen[cand] = 1;
		}
	}
	assert(verifyMDM());
	for (uint32 i = sp->propagated; i < trail.size(); i++) sp->seen[l2a(trail[i])] = 0;
	assert(verifySeen());
	lrn.numMDs = trail.size() - sp->propagated;
	lrn.nRefVars = inf.maxVar - lrn.numMDs, lrn.rounds--;
	stats.n_mds += lrn.numMDs;
	PFLENDING(2, 4, "(%d elected)", lrn.numMDs);
	if (mdm_vsids_pumps || mdm_vmfq_pumps) pumpFrozen();
	memset(sp->frozen, 0, inf.maxVar + 1LL);
	printStats(1, 'm');
}

void ParaFROST::MDM()
{
	assert(!satisfied());
	assert(sp->propagated == trail.size());
	assert(conflict == NOREF);
	assert(cnfstate == UNSOLVED);
	if (vsids()) eligibleVSIDS();
	else if (!mdmvsidsonly_en) eligibleVMFQ();
	assert(eligible.size() >= 1);
	if (verbose == 4) PFLOG1(" Electing decisions (rounds=%d)...", lrn.rounds);
	else PFLOGN2(2, " Electing decisions (rounds=%d)...", lrn.rounds);
	stats.mdm_calls++;
	for (uint32 i = 0; i < eligible.size(); i++) {
		uint32 cand = eligible[i];
		assert(cand && cand <= inf.maxVar);
		if (sp->frozen[cand] || sp->locked[cand]) continue;
		uint32 dec = what(cand, useTarget());
		if (valid(wtBin[dec], wt[dec], cand) && depFreeze(wt[dec], cand)) {
			incDL();
			enqueue(dec, DL());
			sp->seen[cand] = 1;
		}
	}
	assert(verifyMDM());
	for (uint32 i = sp->propagated; i < trail.size(); i++) sp->seen[l2a(trail[i])] = 0;
	assert(verifySeen());
	lrn.numMDs = trail.size() - sp->propagated;
	lrn.nRefVars = inf.maxVar - lrn.numMDs, lrn.rounds--;
	stats.n_mds += lrn.numMDs;
	PFLENDING(2, 4, "(%d elected)", lrn.numMDs);
	if (mdm_vmfq_pumps) pumpFrozen();
	memset(sp->frozen, 0, inf.maxVar + 1LL);
	printStats(lrn.rounds == 0, 'm');
}

bool ParaFROST::valid(WL& wsBins, WL& ws, const uint32& cand)
{
	// check binary clauses
	for (WATCH* w = wsBins; w != wsBins.end(); w++) {
		CLAUSE& c = cm[w->ref];
		assert(c.size() == 2);
		if (value(c[0]) <= 0 || value(c[1]) <= 0) return false; // binary not satisfied
	}
	// check k-size clauses
	for (WATCH* w = ws; w != ws.end(); w++) {
		CLAUSE& c = cm[w->ref];
		assert(c.size() > 2);
		assert(w->imp);
		if (isTrue(w->imp)) continue; // clause satisfied
		bool satisfied = false, unAssigned = false;
		uint32* k = c + 2, * cend = c.end();
		while (k != cend && !satisfied && !unAssigned) {
			LIT_ST val = value(*k);
			if (val > 0) satisfied = true; // clause satisfied
			else if (val < 0) unAssigned = true;
			k++;
		}
		if (!satisfied && !unAssigned) return false;
	}
	return true;
}

bool ParaFROST::depFreeze(WL& ws, const uint32& cand)
{
	for (WATCH* w = ws; w != ws.end(); w++) {
		CLAUSE& c = cm[w->ref];
		assert(c.size() > 2);
		if (isTrue(w->imp)) continue;
		uint32 v0 = l2a(c[0]), frozen_v = (v0 == cand) ? l2a(c[1]) : v0;
		if (sp->seen[frozen_v]) return false;
		sp->frozen[frozen_v] = 1;
		uint32* k = c + 2, * cend = c.end();
		while (k != cend) {
			uint32 v = l2a(*k++);
			if (sp->seen[v]) return false;
			sp->frozen[v] = 1;
		}
	}
	return true;
}

void ParaFROST::eligibleVSIDS()
{
	PFLOGN2(2, " Finding eligible decisions at MDM round %d..", lrn.rounds);
	hist(orgs, true);
	hist(learnts);
	uint32* scores = sp->tmp_stack;
	eligible.clear();
	for (uint32 i = 0; i < vHeap.size(); i++) {
		uint32 v = vHeap[i];
		if (sp->locked[v]) continue;
		eligible.push(v);
		scores[v] = rscore(v);
	}
	assert(eligible.size() >= 1);
	Sort(eligible, KEY_CMP_ACTIVITY<double, uint32>(varAct, scores));
	PFLDONE(2, 5);
	if (verbose >= 3) {
		PFLOG0(" Eligible decisions:");
		for (uint32 i = 0; i < eligible.size(); i++)
			PFLOG1("  e[%d]->(s: %d, a: %g)", eligible[i], scores[eligible[i]], varAct[eligible[i]]);
	}
}

void ParaFROST::eligibleVMFQ()
{
	PFLOGN2(2, " Finding eligible decisions at MDM round %d..", lrn.rounds);
	eligible.clear();
	uint32 free = vQueue.free();
	assert(free);
	while (free) {
		if (!sp->locked[free]) eligible.push(free);
		free = vQueue.previous(free);
	}
	assert(eligible.size() >= 1);
	Sort(eligible, KEY_CMP_BUMP<int64>(varBumps));
	PFLDONE(2, 5);
	if (verbose >= 3) {
		PFLOG0(" Eligible decisions:");
		for (uint32 i = 0; i < eligible.size(); i++)
			PFLOG1("  e[%d]->(b: %lld)", eligible[i], varBumps[eligible[i]]);
	}
}
