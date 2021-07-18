/***********************************************************************[uip.cpp]
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
using namespace pFROST;

inline int ParaFROST::calcLBD(CLAUSE& c) 
{
	int lbd = 0;
	if (c.binary()) lbd = (l2dl(c[0]) != l2dl(c[1])) + 1;
	else {
		const uint64 marker = ++stats.marker;
		forall_clause(c, k) {
			const int litLevel = l2dl(*k);
			if (NEQUAL(sp->board[litLevel], marker)) { 
				sp->board[litLevel] = marker; 
				lbd++; 
			}
		}
	}
	return lbd;
}

inline void	ParaFROST::bumpClause(CLAUSE& c)
{
	assert(c.learnt());
	assert(c.size() > 2);
	const bool hyper = c.hyper();
	if (!hyper && c.keep()) return;
	const CL_ST used = c.usage();
	c.initTier3();
	if (hyper) return;
	const int old_lbd = c.lbd();
	const int new_lbd = calcLBD(c);
	if (new_lbd < old_lbd) { // update old LBD
		if (new_lbd <= opts.lbd_tier1) c.set_keep(true);
		else if (old_lbd > opts.lbd_tier2 && new_lbd <= opts.lbd_tier2) c.initTier2();
		c.set_lbd(new_lbd);
		PFLCLAUSE(4, c, "  bumping clause with lbd %d ", new_lbd);
	}
	else if (used && old_lbd <= opts.lbd_tier2) c.initTier2();
}

inline void	ParaFROST::analyzeLit(const uint32& lit, int& track, int& size)
{
	CHECKLIT(lit);
	assert(isFalse(lit));
	const uint32 v = ABS(lit);
	const int litlevel = sp->level[v];
	if (litlevel) {
		size++;
		if (!sp->seen[v]) {
			sp->resolventsize++;
			const int level = DL();
			sp->seen[v] = ANALYZED_M;
			analyzed.push(v);
			assert(litlevel <= level);
			if (litlevel == level) track++;
			else if (litlevel < level) {
				sp->seen[v] = REMOVABLE_M;
				learntC.push(lit);
				VSTATE& vstate = sp->vstate[litlevel];
				assert(vstate.dlcount < MAX_DLC);
				if (vstate.dlcount == 2 || vstate.dlcount++) return;
				lbdlevels.push(litlevel);
			}
		}
	}
}

inline bool ParaFROST::analyzeReason(const C_REF& ref, const uint32& parent, int& track) 
{
	CHECKLIT(parent);
	CLAUSE& reason = cm[ref];
	PFLCLAUSE(4, reason, "  analyzing %d reason", l2i(parent));
	sp->reasonsize = 1;
	sp->conflictdepth++;
	if (reason.binary()) 
		analyzeLit((reason[0] ^ reason[1] ^ parent), track, sp->reasonsize);
	else {
		if (reason.learnt()) bumpClause(reason);
		forall_clause(reason, k) {
			const uint32 lit = *k;
			if (NEQUAL(lit, parent))
				analyzeLit(lit, track, sp->reasonsize);
		}
	}
	assert(sp->resolventsize > 0);
	sp->resolventsize--;
	if (sp->reasonsize > 2 && sp->resolventsize < sp->reasonsize) {
		assert(sp->resolventsize >= 0);
		assert(!reason.binary());
		assert(!reason.deleted());
		strengthenOTF(reason, ref, parent);
		if (sp->conflictdepth == 1 && sp->resolventsize < sp->conflictsize) {
			assert(sp->conflictsize > 2);
			assert(conflict != NOREF);
			assert(ref != conflict);
			CLAUSE& subsumed = cm[conflict];
			assert(reason.size() <= subsumed.size());
			if (reason.original() || subsumed.learnt()) {
				PFLCLAUSE(4, subsumed, "  found subsumed conflict");
				removeClause(subsumed, conflict);
				stats.subsume.subsumedfly++;
			}
		}
		conflict = ref;
		return true;
	}
	return false;
}

inline void ParaFROST::strengthenOTF(CLAUSE& c, const C_REF& ref, const uint32& self)
{
	CHECKLIT(self);
	assert(c.size() > 2);
	assert(c[0] == self || c[1] == self);
	PFLOG2(4, "  parent(%d) is strengthening last reason clause", l2i(self));
	stats.subsume.strengthenedfly++;
	const uint32 other = c[0] ^ c[1] ^ self;
	c[0] = other, c[1] = self;
	assert(c[0] != c[1]);
	detachWatch(FLIP(self), ref);
	if (opts.proof_en) proof.shrinkClause(c, self);
	uint32* j = c + 1, * end = c.end();
	for (uint32* i = c + 2; i != end; i++) {
		const uint32 lit = *i;
		assert(isFalse(lit));
		if (l2dl(lit))
			*j++ = *i;
	}
	const int removed = int(end - j);
	shrinkClause(c, removed);
	int maxPos = 1;
	int maxLevel = l2dl(c[1]);
	const int size = c.size();
	for (int i = 2; i < size; i++) {
		const uint32 other = c[i];
		const int otherLevel = l2dl(other);
		if (otherLevel > maxLevel) {
			maxPos = i;
			maxLevel = otherLevel;
		}
	}
	if (NEQUAL(maxPos, 1)) swap(c[1], c[maxPos]);
	const uint32 first = c[0];
	const uint32 second = c[1];
	attachWatch(second, first, ref, size);
	WL& ws = wt[FLIP(first)];
	forall_watches(ws, i) {
		if (i->ref == ref) {
			i->imp = second;
			i->size = size;
			break;
		}
	}
}

bool ParaFROST::finduip()
{
	assert(UNSOLVED(cnfstate));
	assert(conflict != NOREF);
	const int level = DL();
	sp->learntLBD = UNDEFINED;
	sp->reasonsize = 0;
	sp->conflictsize = 0;
	sp->resolventsize = 0;
	sp->conflictdepth = 0;
	learntC.push(0);
	int track = 0;
	// analyze conflict
	CLAUSE& c = cm[conflict];
	PFLCLAUSE(4, c, "  analyzing conflict");
	if (c.binary()) {
		analyzeLit(c[0], track, sp->conflictsize);
		analyzeLit(c[1], track, sp->conflictsize);
	}
	else {
		if (c.learnt()) bumpClause(c);
		forall_clause(c, k) {
			analyzeLit(*k, track, sp->conflictsize);
		}
	}
	// analyze reasons
	assert(sp->conflictsize == sp->resolventsize);
	assert(!sp->reasonsize);
	assert(!sp->conflictdepth);
	uint32 parent = 0;
	uint32* tail = trail.end();
	while (true) {
		while (!parent) { // find next implication 
			assert(tail > trail.data());
			const uint32 lit = *--tail, v = ABS(lit);
			if (sp->seen[v] && sp->level[v] == level) 
				parent = lit;
		}
		assert(track);
		if (!--track) break;
		CHECKLIT(parent);
		const C_REF reason = l2r(parent);
		assert(REASON(reason));
		if (analyzeReason(reason, parent, track)) {
			learntC.clear();
			clearLevels();
			clearAnalyzed();
			return true;
		}
		parent = 0;
	}
	assert(learntC[0] == 0);
	learntC[0] = FLIP(parent);
	PFLLEARNT(this, 3);
	sp->learntLBD = (int)lbdlevels.size();
	assert(sp->learntLBD >= 0);
	assert(sp->learntLBD <= learntC.size());
	PFLOG2(4, "  LBD of learnt clause = %d", sp->learntLBD);
	return false;
}