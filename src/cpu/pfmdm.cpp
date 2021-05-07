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

inline bool	ParaFROST::verifyMDM() {
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

inline bool	ParaFROST::verifySeen() {
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

inline void ParaFROST::pumpFrozenHeap(const uint32& lit)
{
	CHECKLIT(lit);
	WL& ws = wt[lit];
	if (ws.empty()) return;
	uint32 v = ABS(lit);
	assert(!sp->frozen[v]);
	double norm_act = (double)sp->level[v] / last.mdm.decisions;
	forall_watches(ws, w) {
		if (cm.deleted(w->ref)) continue;
		uint32 frozen_v;
		if (w->binary()) frozen_v = ABS(w->imp);
		else {
			CLAUSE& c = cm[w->ref];
			frozen_v = ABS(c[0]) ^ ABS(c[1]) ^ v;
			assert(frozen_v != v);	
		}
		CHECKVAR(frozen_v);
		if (activity[frozen_v] == 0) varBumpHeap(frozen_v, norm_act);
	}
}

inline void ParaFROST::pumpFrozenQue(const uint32& lit)
{
	CHECKLIT(lit);
	WL& ws = wt[lit];
	if (ws.empty()) return;
	uint32 v = ABS(lit);
	assert(!sp->frozen[v]);
	forall_watches(ws, w) {
		if (cm.deleted(w->ref)) continue;
		uint32 frozen_v;
		if (w->binary()) frozen_v = ABS(w->imp);
		else {
			CLAUSE& c = cm[w->ref];
			frozen_v = ABS(c[0]) ^ ABS(c[1]) ^ v;
			assert(frozen_v != v);
		}
		CHECKVAR(frozen_v);
		if (sp->frozen[frozen_v]) {
			analyzed.push(frozen_v);
			sp->frozen[frozen_v] = 0;
		}
	}
}

void ParaFROST::pumpFrozen()
{
	if (!last.mdm.decisions) return;
	assert(trail.size() > sp->propagated);
	if (verbose == 4) PFLOG1(" Pumping frozen variables..");
	else PFLOGN2(2, " Pumping frozen variables..");
	uint32* start = trail + sp->propagated;
	if (vsidsEnabled() && opts.mdm_vsids_pumps) {
		for (uint32* assign = trail.end() - 1; assign >= start; assign--)
			pumpFrozenHeap(*assign);
		opts.mdm_vsids_pumps--;
	}
	else if (opts.mdm_vmfq_pumps) { // VMFQ (requires special handling)
		assert(analyzed.empty());
		assert(inf.maxVar >= last.mdm.decisions);
		analyzed.reserve(inf.maxVar - last.mdm.decisions);
		for (uint32* assign = trail.end() - 1; assign >= start; assign--)
			pumpFrozenQue(*assign);
		if (analyzed.size()) {
			for (uint32* v = analyzed.end() - 1; v >= analyzed; v--)
				varBumpQueueNU(*v);
			uint32 first = *analyzed;
			if (UNASSIGNED(sp->value[V2L(first)])) vmtf.update(first, bumps[first]);
			analyzed.clear();
			opts.mdm_vmfq_pumps--;
		}
	}
	PFLDONE(2, 4);
	if (verbose >= 3 && vsidsEnabled()) printHeap();
}

void ParaFROST::histCNF(BCNF& cnf, const bool& reset) {
	if (cnf.empty()) return;
	if (reset) { 
		for (uint32 i = 0; i < occurs.size(); i++) 
			occurs[i] = { 0 , 0 }; 
	}
	forall_cnf(cnf, i) {
		countoccurs(cm[*i]);
	}
	assert(occurs[0].ps == 0 && occurs[0].ns == 0);
}

void ParaFROST::varOrder()
{
	PFLOGN2(2, " Finding eligible decisions at initial round..");
	assert(!learnts.size());
	histCNF(orgs, true);
	uint32* scores = sp->tmp_stack;
	forall_variables(v) { 
		eligible[v - 1] = v, scores[v] = prescore(v);
	}
	if (opts.mdm_mcv_en) rSort(eligible, MCV_CMP(scores), MCV_RANK(scores));
	else rSort(eligible, LCV_CMP(scores), LCV_RANK(scores));
	PFLDONE(2, 5);
	if (verbose >= 3) {
		PFLOG0(" Eligible decisions:");
		for (uint32 i = 0; i < eligible.size(); i++) {
			uint32 v = eligible[i];
			PFLOG1("  e[%d]->(v: %d, p: %d, n: %d, s: %d)", i, v, occurs[v].ps, occurs[v].ns, scores[v]);
		}
	}
}

void ParaFROST::MDMInit()
{
	assert(inf.unassigned);
	assert(sp->propagated == trail.size());
	assert(conflict == NOREF);
	assert(cnfstate == UNSOLVED);
	if (!last.mdm.rounds) return;
	stats.mdmcalls++;
	eligible.resize(inf.maxVar);
	occurs.resize(inf.maxVar + 1);
	varOrder(); // initial variable ordering
	if (verbose == 4) PFLOG1(" Electing decisions (rounds=%d)..", last.mdm.rounds);
	else PFLOGN2(2, " Electing decisions (rounds=%d)..", last.mdm.rounds);
	for (uint32 i = 0; i < inf.maxVar; i++) {
		uint32 cand = eligible[i];
		CHECKVAR(cand);
		if (iassumed(cand)) continue;
		if (sp->frozen[cand] || sp->vstate[cand].state) continue;
		uint32 p = V2L(cand), n = NEG(p);
		if (UNASSIGNED(sp->value[p])) {
			uint32 dec = wt[p].size() >= wt[n].size() ? p : n;
			if (valid(wt[dec]) && depFreeze(wt[dec], cand)) {
				enqueueDecision(dec);
				sp->seen[cand] = 1;
			}
		}
	}
	assert(verifyMDM());
	for (uint32 i = sp->propagated; i < trail.size(); i++) sp->seen[ABS(trail[i])] = 0;
	assert(verifySeen());
	last.mdm.decisions = trail.size() - sp->propagated;
	last.mdm.unassigned = inf.maxVar - last.mdm.decisions;
	last.mdm.rounds--;
	stats.decisions.multiple += last.mdm.decisions;
	PFLENDING(2, 4, "(%d elected)", last.mdm.decisions);
	if (opts.mdm_vsids_pumps || opts.mdm_vmfq_pumps) pumpFrozen();
	memset(sp->frozen, 0, inf.maxVar + 1LL);
	printStats(1, 'm');
	eligible.clear(true);
	occurs.clear(true);
}

void ParaFROST::MDM()
{
	assert(inf.unassigned);
	assert(sp->propagated == trail.size());
	assert(conflict == NOREF);
	assert(cnfstate == UNSOLVED);
	eligible.clear(true);
	occurs.resize(inf.maxVar + 1);
	// TODO: use assumptions in MDM as well
	// by introducing eligibleAssumptions();
	if (vsidsEnabled()) eligibleVSIDS();
	else if (!opts.mdmvsidsonly_en) eligibleVMFQ();
	assert(eligible.size() >= 1);
	if (verbose == 4) PFLOG1(" Electing decisions (rounds=%d)...", last.mdm.rounds);
	else PFLOGN2(2, " Electing decisions (rounds=%d)...", last.mdm.rounds);
	stats.mdmcalls++;
	for (uint32 i = 0; i < eligible.size(); i++) {
		uint32 cand = eligible[i];
		assert(cand && cand <= inf.maxVar);
		if (sp->frozen[cand] || sp->vstate[cand].state) continue;
		uint32 dec = makeAssign(cand, useTarget());
		if (valid(wt[dec]) && depFreeze(wt[dec], cand)) {
			enqueueDecision(dec);
			sp->seen[cand] = 1;
		}
	}
	assert(verifyMDM());
	for (uint32 i = sp->propagated; i < trail.size(); i++) sp->seen[ABS(trail[i])] = 0;
	assert(verifySeen());
	last.mdm.decisions = trail.size() - sp->propagated;
	last.mdm.unassigned = inf.maxVar - last.mdm.decisions;
	last.mdm.rounds--;
	stats.decisions.multiple += last.mdm.decisions;
	PFLENDING(2, 4, "(%d elected)", last.mdm.decisions);
	if (opts.mdm_vmfq_pumps) pumpFrozen();
	memset(sp->frozen, 0, inf.maxVar + 1LL);
	if (opts.mdmfusem_en && !last.mdm.rounds) {
		printStats(1, 'm');
		MDMFuseMaster();
	}
}

inline bool ParaFROST::valid(WL& ws)
{
	forall_watches(ws, i) {
		const WATCH w = *i;
		assert(w.imp);
		if (isTrue(w.imp)) continue; // clause satisfied
		if (w.binary()) return false; // if 'w.imp' not satisfied then it's an implication of 'cand'
		else {
			// there cannot be falsified literal as watched,
			// so validating starts from 'c + 2'
			CLAUSE& c = cm[w.ref];
			assert(c.size() > 2);
			bool satisfied = false, unAssigned = false;
			uint32* k = c + 2, *cend = c.end();
			while (k != cend && !satisfied && !unAssigned) {
				CHECKLIT(*k);
				const LIT_ST val = sp->value[*k];
				if (UNASSIGNED(val)) unAssigned = true;
				else if (val) satisfied = true;
				k++;
			}
			if (!satisfied && !unAssigned) return false;
		}
	}
	return true;
}

inline bool ParaFROST::depFreeze(WL& ws, const uint32& cand)
{
	forall_watches(ws, i) {
		const WATCH w = *i;
		CLAUSE& c = cm[w.ref];
		assert(c.size() > 1);
		if (isTrue(w.imp)) continue; 
		assert(!w.binary());
		const uint32 other_w = ABS(c[0]) ^ ABS(c[1]) ^ cand;
		if (sp->seen[other_w]) return false;
		sp->frozen[other_w] = 1;
		uint32* k = c + 2, *cend = c.end();
		while (k != cend) {
			const uint32 v = ABS(*k++);
			if (sp->seen[v]) return false;
			sp->frozen[v] = 1;
		}
	}
	return true;
}

void ParaFROST::eligibleVSIDS()
{
	PFLOGN2(2, " Finding VSIDS eligible decisions at MDM round %d..", last.mdm.rounds);
	histCNF(orgs, true);
	histCNF(learnts);
	uint32 *scores = sp->tmp_stack;
	for (uint32 i = 0; i < vsids.size(); i++) {
		const uint32 v = vsids[i];
		if (sp->vstate[v].state) continue;
		if (UNASSIGNED(sp->value[V2L(v)])) {
			eligible.push(v);
			scores[v] = prescore(v);
		}
	}
	assert(eligible.size() >= 1);
	Sort(eligible, KEY_CMP_ACTIVITY(activity, scores));
	PFLDONE(2, 5);
	if (verbose >= 3) {
		PFLOG0(" Eligible decisions:");
		for (uint32 i = 0; i < eligible.size(); i++)
			PFLOG1("  e[%d]->(s: %d, a: %g)", eligible[i], scores[eligible[i]], activity[eligible[i]]);
	}
}

void ParaFROST::eligibleVMFQ()
{
	PFLOGN2(2, " Finding VMFQ eligible decisions at MDM round %d..", last.mdm.rounds);
	uint32 free = vmtf.free();
	assert(free);
	while (free) {
		if (!sp->vstate[free].state && UNASSIGNED(sp->value[V2L(free)])) eligible.push(free);
		free = vmtf.previous(free);
	}
	assert(eligible.size() >= 1);
	rSort(eligible, KEY_CMP_BUMP(bumps), KEY_RANK_BUMP(bumps));
	PFLDONE(2, 5);
	if (verbose >= 3) {
		PFLOG0(" Eligible decisions:");
		for (uint32 i = 0; i < eligible.size(); i++)
			PFLOG1("  e[%d]->(b: %lld)", eligible[i], bumps[eligible[i]]);
	}
}

void ParaFROST::MDMFuseMaster() {
	if (opts.mdm_rounds && stats.conflicts >= limit.mdm) {
		last.mdm.rounds = opts.mdm_rounds;
		INCREASE_LIMIT(this, mdm, stats.mdmcalls, nlognlogn, true);
	}
}

void ParaFROST::MDMFuseSlave() {
	if (opts.mdm_rounds && opts.mdm_div) {
		int q = int(stats.conflicts / opts.mdm_div) + 1;
		if (stats.restart.all % q == 0) {
			last.mdm.unassigned = 0, last.mdm.rounds = opts.mdm_rounds;
			PFLOG2(2, " Starts: %lld, conflicts: %lld, dividor: %d, q: %d, rounds: %d", 
				stats.restart.all, stats.conflicts, opts.mdm_div, q, last.mdm.rounds);
		}
		if (stats.conflicts % opts.mdm_freq == 0) opts.mdm_div += opts.mdm_sinc;
	}
}