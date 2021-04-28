/***********************************************************************[pfuip.cpp]
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

inline int ParaFROST::calcLBD(CLAUSE& c) {
	int lbd = 0;
	if (c.binary()) lbd = (l2dl(c[0]) != l2dl(c[1])) + 1;
	else {
		uint64 marker = ++stats.marker;
		forall_clause(c, k) {
			const int litLevel = l2dl(*k);
			if (sp->board[litLevel] != marker) { sp->board[litLevel] = marker; lbd++; }
		}
	}
	return lbd;
}

inline void	ParaFROST::bumpClause(CLAUSE& c) {
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

inline void	ParaFROST::analyzeLit(const uint32& lit, int& track) {
	CHECKLIT(lit);
	assert(isFalse(lit));
	const uint32 v = ABS(lit);
	const int litlevel = sp->level[v];
	if (litlevel && !sp->seen[v]) {
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

inline void ParaFROST::analyzeReason(const C_REF& r, const uint32& parent, int& track) {
	CLAUSE& c = cm[r];
	PFLCLAUSE(4, c, "  analyzing %d %s", parent ? l2i(parent) : l2i(*c), parent ? "reason" : "conflict");
	if (c.binary()) {
		if (parent) analyzeLit((c[0] ^ c[1] ^ parent), track);
		else { analyzeLit(c[0], track), analyzeLit(c[1], track); }
	}
	else {
		if (c.learnt()) bumpClause(c);
		forall_clause(c, k) {
			const uint32 lit = *k;
			if (NEQUAL(lit, parent))
				analyzeLit(lit, track);
		}
	}
}

void ParaFROST::finduip()
{
	assert(cnfstate);
	assert(conflict != NOREF);
	const int level = DL();
	sp->learntLBD = UNDEFINED;
	learntC.push(0);
	uint32 parent = 0;
	int track = 0;
	uint32* tail = trail.end();
	C_REF reason = conflict;
	while (true) {
		assert(reason != NOREF);
		analyzeReason(reason, parent, track);
		parent = 0;
		while (!parent) { // find next implication 
			assert(tail > trail.data());
			const uint32 lit = *--tail, v = ABS(lit);
			if (sp->seen[v] && sp->level[v] == level) parent = lit;
		}
		if (!--track) break;
		CHECKLIT(parent);
		reason = l2r(parent);
	}
	assert(learntC[0] == 0);
	learntC[0] = FLIP(parent);
	PFLLEARNT(this, 3);
	sp->learntLBD = (int)lbdlevels.size();
	assert(sp->learntLBD >= 0);
	assert(sp->learntLBD <= learntC.size());
	PFLOG2(4, "  LBD of learnt clause = %d", sp->learntLBD);
}