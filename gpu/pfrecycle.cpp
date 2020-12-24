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

inline void ParaFROST::moveClause(C_REF& r, CMM& newBlock) {
	assert(r < cm.size());
	CLAUSE& c = cm[r];
	assert(!c.deleted());
	if (c.moved()) { r = c.ref(); return; }
	r = newBlock.alloc(c);
	c.set_ref(r);
}

inline void	ParaFROST::recycleWL(WL& ws, CMM& new_cm) {
	for (WATCH* w = ws; w != ws.end(); w++)
		moveClause(w->ref, new_cm);
}

inline void	ParaFROST::recycleWL(const uint32& lit) {
	WL& ws = wt[lit], saved;
	if (ws.empty()) { ws.clear(true); return; }
	WATCH* i, *j = ws, *end = ws.end();
	for (i = j; i != end; i++) {
		WATCH w = *i;
		assert(w.ref != NOREF);
		const CLAUSE& c = cm[w.ref];
		if (c.deleted()) continue;
		w.size = c.size();
		int litpos = (c[0] == FLIP(lit));
		assert(c[!litpos] == FLIP(lit));
		w.imp = c[litpos];
		if (w.binary()) *j++ = w;
		else saved.push(w);
	}
	ws.resize(int(j - ws));
	for (WATCH* s = saved; s != saved.end(); s++) ws.push(*s);
	if (ws.empty()) ws.clear(true);
	saved.clear(true);
}

void ParaFROST::protectReasons() {
	for (uint32 i = 0; i < trail.size(); i++) {
		uint32 lit = trail[i], v = ABS(lit);
		if (sp->vstate[v] != ACTIVE) continue;
		assert(!unassigned(lit));
		assert(sp->level[v]);
		C_REF r = sp->source[v];
		if (!REASON(r)) continue;
		assert(!cm[r].reason());
		cm[r].markReason();
	}
}

void ParaFROST::unprotectReasons() {
	for (uint32 i = 0; i < trail.size(); i++) {
		uint32 lit = trail[i], v = ABS(lit);
		if (sp->vstate[v] != ACTIVE) continue;
		assert(!unassigned(lit));
		assert(sp->level[v]);
		C_REF r = sp->source[v];
		if (!REASON(r)) continue;
		assert(cm[r].reason());
		cm[r].initReason();
	}
}

void ParaFROST::recycleWT() {
	for (uint32 v = 1; v <= inf.maxVar; v++) {
		uint32 p = V2L(v), n = NEG(p);
		recycleWL(p), recycleWL(n);
	}
}

void ParaFROST::recycle(CMM& new_cm)
{
	recycleWT();
	for (uint32 q = vmfq.last(); q; q = vmfq.previous(q)) {
		uint32 lit = makeAssign(q), flit = FLIP(lit);
		recycleWL(wt[lit], new_cm);
		recycleWL(wt[flit], new_cm);
	}
	uint32 count = 0;
	for (uint32* t = trail; t != trail.end(); t++) {
		uint32 lit = *t, v = ABS(lit);
		C_REF& r = sp->source[v];
		if (r == NOREF) continue;
		if (!sp->level[v]) { r = NOREF; continue; }
		assert(r < cm.size());
		if (cm[r].deleted()) { r = NOREF; continue; }
		assert(cm[r].reason());
		moveClause(r, new_cm);
		count++;
	}
	PFPRINT(2, 5, "(updated %d sources)", count);
	filter(orgs, new_cm);
	filter(learnts, new_cm);
}

void ParaFROST::recycle() {
	assert(sp->propagated == trail.size());
	assert(conflict == NOREF);
	assert(cnfstate == UNSOLVED);
	stats.recyclings++;
	shrink();
	if (canCollect()) {
		PFLOGN2(2, " Recycling garbage..");
		CMM new_cm(cm.size() - cm.garbage());
		recycle(new_cm);
		PFLGCMEM(2, cm, new_cm);
		new_cm.migrate(cm);
		PFLDONE(2, 5);
	}
	else recycleWT(), filter(learnts);
}

//=========================//
// shrinking CNF routines  //
//=========================//
bool ParaFROST::shrink()
{
	if (sp->simplified >= inf.maxFrozen) return false;
	sp->simplified = inf.maxFrozen;
	stats.shrinkages++;
	PFLOGN2(2, " Shrinking CNF..");
	assert(trail.size());
	assert(conflict == NOREF);
	assert(cnfstate == UNSOLVED);
	assert(sp->propagated == trail.size());
	assert(!unassigned(trail.back()));
	int64 beforeCls = maxClauses(), beforeLits = maxLiterals();
	shrink(orgs);
	shrink(learnts);
	PFLENDING(2, 5, "(-%lld clauses, -%lld literals)", 
		beforeCls - maxClauses(), beforeLits - maxLiterals());
	return true;
}

void ParaFROST::shrink(BCNF& cnf) {
	if (cnf.empty()) return;
	C_REF* i, * j, * end = cnf.end();
	for (i = cnf, j = i; i != end; i++) {
		C_REF r = *i;
		CLAUSE& c = cm[r];
		assert(!c.moved());
		if (c.deleted()) continue;
		CL_ST st = rooted(c);
		if (st > 0) removeClause(r);
		else if (!st) {
			shrinkClause(r);
			*j++ = r;
		}
		else *j++ = r;
	}
	assert(j >= cnf);
	cnf.resize(uint32(j - cnf));
}

void ParaFROST::filter(BCNF& cnf) {
	if (cnf.empty()) return;
	C_REF* i, * j, * end = cnf.end();
	for (i = cnf, j = i; i != end; i++) {
		C_REF r = *i;
		if (cm[r].deleted()) continue;
		*j++ = r;
	}
	assert(j >= cnf);
	cnf.resize(uint32(j - cnf));
}

void ParaFROST::filter(BCNF& cnf, BCNF& dest, const CL_ST& t) {
	if (cnf.empty()) return;
	C_REF* i, * j, * end = cnf.end();
	for (i = cnf, j = i; i != end; i++) {
		C_REF r = *i;
		CL_ST st = cm[r].status();
		if (st & DELETED) continue;
		if (st & t) {
			dest.push(r);
			continue;
		}
		*j++ = r;
	}
	assert(j >= cnf);
	cnf.resize(uint32(j - cnf));
}

void ParaFROST::filter(BCNF& cnf, CMM& new_cm) {
	if (cnf.empty()) return;
	C_REF* i, * j, * end = cnf.end();
	for (i = cnf, j = i; i != end; i++) {
		C_REF r = *i;
		if (cm[r].deleted()) continue;
		moveClause(r, new_cm);
		*j++ = r; // must follow moveClause
	}
	assert(j >= cnf);
	cnf.resize(uint32(j - cnf));
}

