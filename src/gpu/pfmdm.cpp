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

inline bool	ParaFROST::verifyMDM() 
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

inline bool	ParaFROST::verifySeen()
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

inline void	ParaFROST::clearMDM() 
{
	assert(verifyMDM());
	uint32* start = trail + sp->propagated, *end = trail.end();
	while (start != end)
		sp->seen[ABS(*start++)] = 0;
	assert(verifySeen());
	assert((sp->stacktail - sp->tmp_stack) <= (inf.maxVar - last.mdm.decisions));
	clearFrozen();
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

void ParaFROST::varOrder()
{
	PFLOGN2(2, " Finding eligible decisions at initial round..");
	histCNF(orgs, true);
	histCNF(learnts);
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
	if (!last.mdm.rounds) return;
	assert(inf.unassigned);
	assert(sp->propagated == trail.size());
	assert(conflict == NOREF);
	assert(cnfstate == UNSOLVED);
	stats.mdm.calls++;
	PFLOG2(2, " MDM %d: electing decisions at decaying round %d..", stats.mdm.calls, last.mdm.rounds);
	eligible.resize(inf.maxVar);
	occurs.resize(inf.maxVar + 1);
	varOrder(); // initial variable ordering
	sp->stacktail = sp->tmp_stack;
	if (opts.mdmassume_en && assumptions.size()) {
		assert(sp->stacktail == sp->tmp_stack);
		int level = DL();
		while (level < assumptions.size()) {
			const uint32 a = assumptions[level];
			CHECKLIT(a);
			assert(ifrozen[ABS(a)]);
			const uint32 cand = ABS(a);
			const LIT_ST val = sp->value[a];
			if (UNASSIGNED(val)) {
				level++;
				uint32 dec = a;
				if (valid(wt[dec]) && depFreeze(wt[dec], cand)) {
					enqueueDecision(dec);
					sp->seen[cand] = 1;
					stats.decisions.massumed++;
				}
			}
			else if (!val) {
				ianalyze(FLIP(a));
				cnfstate = UNSAT;
				clearMDM(), eligible.clear(true), occurs.clear(true);
				return;
			}
		}
	}
	else 
		assert(sp->stacktail == sp->tmp_stack);
	for (uint32 i = 0; i < inf.maxVar; i++) {
		uint32 cand = eligible[i];
		CHECKVAR(cand);
		if (sp->frozen[cand] || sp->vstate[cand].state || iassumed(cand)) continue;
		const uint32 p = V2L(cand), n = NEG(p);
		if (UNASSIGNED(sp->value[p])) {
			const uint32 dec = wt[p].size() >= wt[n].size() ? p : n;
			if (valid(wt[dec]) && depFreeze(wt[dec], cand)) {
				enqueueDecision(dec);
				sp->seen[cand] = 1;
			}
		}
	}
	last.mdm.decisions = trail.size() - sp->propagated;
	last.mdm.unassigned = inf.maxVar - last.mdm.decisions;
	last.mdm.rounds--;
	stats.decisions.multiple += last.mdm.decisions;
	PFLOG2(2, " MDM %d: %d decisions are elected (%.2f%%)",
		stats.mdm.calls, last.mdm.decisions, percent(last.mdm.decisions, maxActive()));
	if (opts.mdm_vsids_pumps || opts.mdm_vmfq_pumps) pumpFrozen();
	clearMDM(), eligible.clear(true), occurs.clear(true);
	printStats(1, 'm', CMDM);
}

void ParaFROST::MDM()
{
	const bool vsidsActive = vsidsEnabled();
	if (opts.mdmvsidsonly_en && !vsidsActive) {
		last.mdm.rounds--;
		return;
	}
	assert(inf.unassigned);
	assert(sp->propagated == trail.size());
	assert(conflict == NOREF);
	assert(cnfstate == UNSOLVED);
	stats.mdm.calls++;
	PFLOG2(2, " MDM %d: electing decisions at decaying round %d..", stats.mdm.calls, last.mdm.rounds);
	eligible.clear(true);
	if (vsidsActive) eligibleVSIDS();
	else eligibleVMFQ();
	assert(eligible.size() >= 1);
	sp->stacktail = sp->tmp_stack;
	if (opts.mdmassume_en && assumptions.size()) {
		assert(sp->stacktail == sp->tmp_stack);
		int level = DL();
		while (level < assumptions.size()) {
			const uint32 a = assumptions[level];
			CHECKLIT(a);
			assert(ifrozen[ABS(a)]);
			const uint32 cand = ABS(a);
			const LIT_ST val = sp->value[a];
			if (UNASSIGNED(val)) {
				level++;
				uint32 dec = a;
				if (valid(wt[dec]) && depFreeze(wt[dec], cand)) {
					enqueueDecision(dec);
					sp->seen[cand] = 1;
					stats.decisions.massumed++;
				}
			}
			else if (!val) {
				ianalyze(FLIP(a));
				cnfstate = UNSAT;
				clearMDM();
				return;
			}
		}
	}
	else 
		assert(sp->stacktail == sp->tmp_stack);
	for (uint32 i = 0; i < eligible.size(); i++) {
		uint32 cand = eligible[i];
		CHECKVAR(cand);
		if (sp->frozen[cand] || sp->vstate[cand].state || iassumed(cand)) continue;
		uint32 dec = makeAssign(cand, useTarget());
		if (valid(wt[dec]) && depFreeze(wt[dec], cand)) {
			enqueueDecision(dec);
			sp->seen[cand] = 1;
		}
	}
	last.mdm.decisions = trail.size() - sp->propagated;
	last.mdm.unassigned = inf.maxVar - last.mdm.decisions;
	last.mdm.rounds--;
	stats.decisions.multiple += last.mdm.decisions;
	PFLOG2(2, " MDM %d: %d decisions are elected (%.2f%%)", 
		stats.mdm.calls, last.mdm.decisions, percent(last.mdm.decisions, maxActive()));
	if (opts.mdm_vsids_pumps || opts.mdm_vmfq_pumps) pumpFrozen();
	clearMDM();
	printStats(last.mdm.rounds == opts.mdm_rounds, 'm', CMDM);
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
			uint32* k = c + 2, * cend = c.end();
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
	LIT_ST* frozen = sp->frozen;
	uint32*& frozen_stack = sp->stacktail;
	forall_watches(ws, i) {
		const WATCH w = *i;
		if (isTrue(w.imp)) continue;
		assert(!w.binary());
		CLAUSE& c = cm[w.ref];
		uint32 othervar = ABS(c[0]) ^ ABS(c[1]) ^ cand;
		if (sp->seen[othervar]) return false;
		if (!frozen[othervar]) {
			frozen[othervar] = 1;
			assert(frozen_stack < sp->tmp_stack + inf.maxVar);
			*frozen_stack++ = othervar;
		}
		uint32* k = c + 2, * cend = c.end();
		while (k != cend) {
			othervar = ABS(*k++);
			if (sp->seen[othervar]) return false;
			if (!frozen[othervar]) {
				frozen[othervar] = 1;
				assert(frozen_stack < sp->tmp_stack + inf.maxVar);
				*frozen_stack++ = othervar;
			}
		}
	}
	return true;
}

void ParaFROST::eligibleVSIDS()
{
	PFLOGN2(2, "  finding VSIDS eligible decisions..");
	stats.mdm.vsids++;
	occurs.resize(inf.maxVar + 1);
	histCNF(orgs, true);
	histCNF(learnts);
	uint32* scores = sp->tmp_stack;
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
		PFLOG0("  eligible decisions:");
		for (uint32 i = 0; i < eligible.size(); i++)
			PFLOG1("  e[%d]->(s: %d, a: %g)", eligible[i], scores[eligible[i]], activity[eligible[i]]);
	}
}

void ParaFROST::eligibleVMFQ()
{
	PFLOGN2(2, "  finding VMFQ eligible decisions..");
	stats.mdm.vmtf++;
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
		PFLOG0("  eligible decisions:");
		for (uint32 i = 0; i < eligible.size(); i++)
			PFLOG1("  e[%d]->(b: %lld)", eligible[i], bumps[eligible[i]]);
	}
}

bool ParaFROST::canMMD()
{
	if (!opts.mdm_rounds) return false;
	const bool enough = varsEnough();
	const bool rounds = last.mdm.rounds;
	if (enough && !rounds && stats.conflicts >= limit.mdm) {
		last.mdm.rounds = opts.mdm_rounds;
		INCREASE_LIMIT(this, mdm, stats.mdm.calls, linear, true);
	}
	return enough && rounds;
}