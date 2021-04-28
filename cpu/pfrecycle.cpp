/***********************************************************************[pfrecycle.cpp]
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

inline void ParaFROST::moveClause(C_REF& r, CMM& newBlock)
{
	assert(r < cm.size());
	CLAUSE& c = cm[r];
	assert(!c.deleted());
	if (c.moved()) { r = c.ref(); return; }
	r = newBlock.alloc(c);
	c.set_ref(r);
}

inline void	ParaFROST::moveWatches(WL& ws, CMM& new_cm)
{
	forall_watches(ws, w) {
		moveClause(w->ref, new_cm);
	}
}

inline void	ParaFROST::recycleWL(const uint32& lit, const bool& priorbins)
{
	CHECKLIT(lit);
	WL& ws = wt[lit], hypers, saved;
	if (ws.empty()) { ws.clear(true); return; }
	const uint32 flit = FLIP(lit);
	WATCH *j = ws;
	forall_watches(ws, i) {
		WATCH w = *i;
		const C_REF r = w.ref;
		assert(r != NOREF);
		if (cm.deleted(r)) continue;
		CLAUSE& c = cm[r];
		w.imp = c[0] ^ c[1] ^ flit;
		w.size = c.size();
		if (c.binary()) {
			if (c.hyper()) hypers.push(w);
			else *j++ = w;
		}
		else if (c.original()) {
			if (priorbins) saved.push(w);
			else *j++ = w;
		}
	}
	ws.resize(int(j - ws));
	forall_watches(hypers, i) ws.push(*i);
	if (priorbins) {
		forall_watches(saved, i) ws.push(*i);
		saved.clear(true);
	}
	ws.shrinkCap();
	hypers.clear(true);
}

void ParaFROST::protectReasons() 
{
	VSTATE* states = sp->vstate;
	C_REF* sources = sp->source;
	forall_vector(uint32, trail, t) {
		const uint32 lit = *t, v = ABS(lit);
		if (states[v].state) continue;
		assert(!unassigned(lit));
		assert(sp->level[v]);
		const C_REF r = sources[v];
		if (!REASON(r)) continue;
		assert(!cm[r].reason());
		cm[r].markReason();
	}
}

void ParaFROST::unprotectReasons() 
{
	VSTATE* states = sp->vstate;
	C_REF* sources = sp->source;
	forall_vector(uint32, trail, t) {
		const uint32 lit = *t, v = ABS(lit);
		if (states[v].state) continue;
		assert(!unassigned(lit));
		assert(sp->level[v]);
		const C_REF r = sources[v];
		if (!REASON(r)) continue;
		assert(cm[r].reason());
		cm[r].initReason();
	}
}

void ParaFROST::recycleWT() 
{
	forall_variables(v) {
		uint32 p = V2L(v), n = NEG(p);
		recycleWL(p, false), recycleWL(n, false);
	}
	forall_cnf(learnts, i) {
		const C_REF r = *i;
		assert(r < cm.size());
		if (cm.deleted(r)) continue;
		CLAUSE& c = cm[r];
		if (c.binary()) continue;
		sortClause(c);
		attachWatch(r, c);
	}
}

void ParaFROST::recycle(CMM& new_cm)
{
	reduced.clear(true);
	analyzed.clear(true);
	recycleWT();
	for (uint32 q = vmtf.last(); q; q = vmtf.previous(q)) {
		const uint32 lit = makeAssign(q), flit = FLIP(lit);
		moveWatches(wt[lit], new_cm);
		moveWatches(wt[flit], new_cm);
	}
	C_REF* sources = sp->source;
	int* levels = sp->level;
	uint32* tend = trail.end();
	for (uint32* t = trail; t != tend; t++) {
		const uint32 lit = *t, v = ABS(lit);
		C_REF& r = sources[v];
		if (r == NOREF) continue;
		if (!levels[v]) { r = NOREF; continue; }
		assert(r < cm.size());
		if (cm.deleted(r)) r = NOREF;
		else moveClause(r, new_cm);
	}
	filter(orgs, new_cm);
	filter(learnts, new_cm);
	orgs.shrinkCap();
}

void ParaFROST::recycle() 
{
	assert(sp->propagated == trail.size());
	assert(conflict == NOREF);
	assert(cnfstate == UNSOLVED);
	shrink();
	if (canCollect()) {
		stats.recycle.hard++;
		PFLOGN2(2, " Recycling garbage..");
		CMM new_cm(cm.size() - cm.garbage());
		recycle(new_cm);
		PFLGCMEM(2, cm, new_cm);
		new_cm.migrateTo(cm);
		PFLDONE(2, 5);
	}
	else {
		stats.recycle.soft++;
		recycleWT();
		filter(learnts);
	}
}

void ParaFROST::filter(BCNF& cnf) 
{
	if (cnf.empty()) return;
	C_REF* j = cnf;
	forall_cnf(cnf, i) {
		const C_REF r = *i;
		if (cm.deleted(r)) continue;
		*j++ = r;
	}
	assert(j >= cnf);
	cnf.resize(uint32(j - cnf));
}

void ParaFROST::filter(BCNF& cnf, CMM& new_cm)
{
	if (cnf.empty()) return;
	C_REF* j = cnf;
	forall_cnf(cnf, i) {
		C_REF r = *i;
		if (cm.deleted(r)) continue;
		moveClause(r, new_cm);
		*j++ = r; // must follow moveClause
	}
	assert(j >= cnf);
	cnf.resize(uint32(j - cnf));
}

