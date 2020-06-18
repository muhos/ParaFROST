/***********************************************************************[pfopts.h]
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

void ParaFROST::initQueue() {
	if (verbose == 4) PFLOG2(2, "  Initializing VMFQ Queue with %d variables..", inf.maxVar);
	else PFLOGN2(2, "  Initializing VMFQ Queue with %d variables..", inf.maxVar);
	for (uint32 v = 1; v <= inf.maxVar; v++)
		vQueue.init(v), vQueue.update(v, (varBumps[v] = ++lrn.bumped));
	PFLDONE(2, 4);
}

void ParaFROST::initHeap() {
	PFLOGN2(2, "  Initializing VSIDS Heap with %d variables..", inf.maxVar);
	for (uint32 v = 1; v <= inf.maxVar; v++)
		vHeap.insert(v);
	PFLDONE(2, 5);
}

void ParaFROST::savePhases(LIT_ST* to, const uint32& off, const uint32& size)
{
	assert(size);
	assert(off < NOREF);
	assert(off < size);
	uint32* lit = trail + off, *end = trail + size;
	for (; lit != end; lit++) 
		to[l2a(*lit)] = sign(*lit);
}

void ParaFROST::savePhases(LIT_ST* to)
{
	for (uint32 v = 1; v <= inf.maxVar; v++)
		to[v] = sp->value[v2l(v)];
}

void ParaFROST::savePhases(const int& bt_level) {

	LIT_ST reset = (lrn.lastrephased && nConflicts > lrn.rephase_last_max);
	if (reset) {
		lrn.target = 0;
		if (lrn.lastrephased == BESTPHASE) lrn.best = 0;
	}
	uint32 pivot = trail_lens[bt_level];
	if (pivot > lrn.target) {
		savePhases(sp->ptarget, pivot, trail.size());
		lrn.target = pivot;
	}
	if (pivot > lrn.best) {
		savePhases(sp->pbest, 0, pivot);
		lrn.best = pivot;
	}
	if (reset) lrn.lastrephased = 0;
}

void ParaFROST::rephase() {

	stats.n_rephs++;
	cancelAssigns();
	memset(sp->ptarget, UNDEFINED, inf.maxVar + 1ULL);
	lrn.lastrephased = 0, lrn.target = 0;
	int64 count = lrn.rephased[lrn.stable]++;
	LIT_ST which = UNDEFINED;
	if ((stable_en && vsidsonly_en) || !stable_en) {
		which = count % 5;
		if (which == 0 || which == 2 || which == 4) varBestPhase();
		else if (which == 1) varFlipPhase();
		else if (which == 3) varOrgPhase();
	}
	else if (lrn.stable) {
		which = count % 3;
		if (which == 0 || which == 2) varBestPhase();
		else if (which == 1) varOrgPhase();
	}
	else {
		assert(!lrn.stable);
		which = count % 4;
		if (which == 0 || which == 2) varFlipPhase();
		else if (which == 1 || which == 3) varBestPhase();
	}
	assert(lrn.lastrephased);
	lrn.rephase_last_max = nConflicts;
	lrn.rephase_conf_max = lrn.rephase_last_max + rephase_inc * (stats.n_rephs + 1);
}

uint32 ParaFROST::nextVSIDS()
{
	uint32 cand = 0;
	while (!vHeap.empty()) {
		cand = vHeap.top();
		assert(cand && cand <= inf.maxVar);
		if (!sp->locked[cand]) break;
		vHeap.pop();
	}
	assert(cand);
	PFLOG2(4, " Next heap choice %d, activity %g", cand, varAct[cand]);
	return cand;
}

uint32 ParaFROST::nextVMFQ()
{
	LIT_ST assigned = UNDEFINED;
	uint32 free = vQueue.free();
	assert(free);
	assigned = sp->locked[free];
	while (sp->locked[free]) free = vQueue.previous(free);
	assert(!assigned || assigned == 1);
	if (assigned) vQueue.update(free, varBumps[free]);
	PFLOG2(4, " Next queue choice %d, bumped %lld", free, varBumps[free]);
	return free;
}

void ParaFROST::decide()
{
	assert(trail.size() < inf.maxVar - inf.maxMelted);
	assert(sp->propagated == trail.size());
	assert(conflict == NOREF);
	assert(cnfstate == UNSOLVED);
	uint32 cand = vsids() ? nextVSIDS() : nextVMFQ();
	uint32 dec = what(cand, useTarget());
	incDL();
	enqueue(dec, DL());
	stats.n_fuds++;
}