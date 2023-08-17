/***********************************************************************[recycle.cpp]
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

inline void Solver::moveClause(C_REF& r, CMM& newBlock)
{
	assert(r < cm.size());
	CLAUSE& c = cm[r];
	assert(!c.deleted());
	if (c.moved()) { r = c.ref(); return; }
	r = newBlock.alloc(c);
	c.set_ref(r);
}

inline void	Solver::moveWatches(WL& ws, CMM& new_cm)
{
	forall_watches(ws, w) {
		moveClause(w->ref, new_cm);
	}
	ws.shrinkCap();
}

inline void	Solver::recycleWL(const uint32& lit)
{
	CHECKLIT(lit);
	WL& ws = wt[lit], hypers;
	if (ws.empty()) return;
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
		else if (c.original()) *j++ = w;
	}
	ws.resize(int(j - ws));
	forall_watches(hypers, i) ws.push(*i);
	hypers.clear(true);
}

void Solver::markReasons() 
{
	const VSTATE* states = sp->vstate;
	const C_REF* sources = sp->source;
	forall_vector(uint32, trail, t) {
		const uint32 lit = *t, v = ABS(lit);
		if (states[v].state) continue;
		assert(!unassigned(lit));
		assert(sp->level[v]);
		const C_REF r = sources[v];
		if (REASON(r)) {
			assert(!cm[r].reason());
			cm[r].markReason();
		}
	}
}

void Solver::unmarkReasons() 
{
	const VSTATE* states = sp->vstate;
	const C_REF* sources = sp->source;
	forall_vector(uint32, trail, t) {
		const uint32 lit = *t, v = ABS(lit);
		if (states[v].state) continue;
		assert(!unassigned(lit));
		assert(sp->level[v]);
		const C_REF r = sources[v];
		if (REASON(r)) {
			assert(cm[r].reason());
			cm[r].initReason();
		}
	}
}

void Solver::recycleWT() 
{
	forall_variables(v) {
		uint32 p = V2L(v), n = NEG(p);
		recycleWL(p), recycleWL(n);
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

void Solver::recycle(CMM& new_cm)
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
	forall_vector(uint32, trail, t) {
		const uint32 lit = *t, v = ABS(lit);
		C_REF& r = sources[v];
		if (REASON(r)) {
			if (levels[v]) {
				assert(r < cm.size());
				if (cm.deleted(r)) r = NOREF;
				else moveClause(r, new_cm);
			}
			else r = NOREF;
		}
	}
	filter(orgs, new_cm);
	filter(learnts, new_cm);
	orgs.shrinkCap();
}

void Solver::recycle() 
{
	assert(sp->propagated == trail.size());
	assert(conflict == NOREF);
	assert(UNSOLVED(cnfstate));
	shrink();
	if (canCollect()) {
		PFLOGN2(2, " Recycling garbage..");
		stats.recycle.hard++;
		assert(cm.size() >= cm.garbage());
		const size_t bytes = cm.size() - cm.garbage();
		CMM new_cm(bytes);
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

void Solver::filter(BCNF& cnf) 
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

void Solver::filter(BCNF& cnf, CMM& new_cm)
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

