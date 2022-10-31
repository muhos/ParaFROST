/***********************************************************************[mdmrank.cpp]
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
using namespace ParaFROST;


void Solver::varOrder()
{
	PFLOGN2(2, " Finding eligible decisions at initial round..");
	histCNF(orgs, true);
	histCNF(learnts);
	uint32* scores = sp->tmpstack;
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

void Solver::eligibleVSIDS()
{
	PFLOGN2(2, "  finding VSIDS eligible decisions..");
	stats.mdm.vsids++;
	occurs.resize(inf.maxVar + 1);
	histCNF(orgs, true);
	histCNF(learnts);
	uint32* scores = sp->tmpstack;
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

void Solver::eligibleVMFQ()
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