/***********************************************************************[pfvmap.cpp]
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

void ParaFROST::map(BCNF& cnf) {
	assert(!vmap.empty());
	for (uint32 i = 0; i < cnf.size(); i++) {
		CLAUSE& c = cm[cnf[i]];
		assert(!c.deleted());
		vmap.mapClause(c);
	}
}

void ParaFROST::map(WL& ws) {
	if (ws.empty()) return;
	for (WATCH* w = ws; w != ws.end(); w++)
		w->imp = vmap.mapLit(w->imp);
}

void ParaFROST::map(WT& wt) {
	if (wt.empty()) return;
	assert(!vmap.empty());
	for (uint32 v = 1; v <= inf.maxVar; v++) {
		uint32 mVar = vmap.mapped(v);
		if (mVar) {
			uint32 p = V2L(v), n = NEG(p);
			uint32 mpos = V2L(mVar), mneg = NEG(mpos);
			if (mVar != v) { // map watch lists
				wt[mpos].copyFrom(wt[p]);
				wt[mneg].copyFrom(wt[n]);
			}
			map(wt[mpos]), map(wt[mneg]); // then map watch imps
		}
	}
	wt.resize(V2L(vmap.size()));
	wt.shrinkCap();
}

void ParaFROST::map(const bool& sigmified)
{
	assert(!satisfied());
	assert(conflict == NOREF);
	assert(!DL());
	assert(trail.size() == sp->propagated);
	stats.mappings++;
	int64 memBefore = sysMemUsed();
	vmap.initiate(sp);
	// map original literals with current values
	vmap.mapOrgs(model.lits);
	// map clauses and watch tables
	if (!sigmified) map(orgs), map(learnts), map(wt);
	else mapped = true, newBeginning(), mapped = false;
	// map trail, queue and heap
	vmap.mapShrinkLits(trail);
	sp->propagated = trail.size();
	const uint32 firstFrozen = vmap.firstL0();
	vmfq.map(*vmap, firstFrozen);
	vmap.mapShrinkVars(vmfq.data());
	vmap.mapShrinkVars(bumps);
	uVec1D tmp;
	while (vsids.size()) {
		uint32 x = vsids.pop();
		if (x == firstFrozen) continue;
		uint32 mx = vmap.mapped(x);
		if (mx) tmp.push(mx);
	}
	vmap.mapShrinkVars(activity);
	vsids.rebuild(tmp);
	// map search space
	SP* newSP = new SP(vmap.size());
	vmap.mapSP(newSP);
	delete sp;
	sp = newSP;
	// update phase-saving counters and
	sp->trailpivot = 0, lrn.best = lrn.target = 0;
	for (uint32 v = 1; v <= vmap.numVars(); v++) {
		if (sp->vstate[v] != ACTIVE) continue;
		if (sp->pbest[v] != UNDEFINED) lrn.best++;
		if (sp->ptarget[v] != UNDEFINED) lrn.target++;
	}
	// reset markers
	stats.marker = 0;
	memset(sp->marks, UNDEFINED, vmap.size());
	PFLOG2(2, " Variable mapping compressed %d to %d, saving %.2f KB of memory",
		inf.maxVar, vmap.numVars(), double(abs(memBefore - sysMemUsed())) / KBYTE);
	inf.maxVar = vmap.numVars();
	inf.nDualVars = V2L(inf.maxVar + 1);
	inf.maxFrozen = inf.maxMelted = 0;
	int64 current_inc = opts.map_inc * (stats.mappings + 1LL);
	lrn.map_conf_max = nConflicts + current_inc;
	PFLOG2(2, " map limit increased to %lld conflicts by a value %lld", lrn.map_conf_max, current_inc);
	vmap.destroy();
}