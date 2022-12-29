/***********************************************************************[minimize.cpp]
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

#include "solve.hpp"

using namespace ParaFROST;

inline int Solver::calcLBD()
{
	assert(!learntC.empty());
	if (learntC.size() == 1) return 0;
	else if (learntC.size() == 2) return (l2dl(learntC[0]) != l2dl(learntC[1]));
	else {
		int lbd = 0;
		const uint64 marker = ++stats.marker;
		forall_clause(learntC, k) {
			const int litLevel = l2dl(*k);
			if (NEQUAL(sp->board[litLevel], marker)) { 
				sp->board[litLevel] = marker;
				lbd++;
			}
		}
		return lbd - 1;
	}
}

void Solver::minimizebin()
{
	const uint32 uip = FLIP(learntC[0]);
	const uint64 unmarker = stats.marker++, marker = stats.marker;
	uint32 *j, *end = learntC.end();
	for (j = learntC + 1; j != end; j++)
		sp->board[ABS(*j)] = marker;
	WL& ws = wt[uip];
	int nLitsRem = 0;
	const int size = learntC.size() - 1;
	forall_watches(ws, i) {
		const WATCH w = *i;
		if (w.binary()) {
			const uint32 other = w.imp, v = ABS(other);
			if (sp->board[v] == marker && isTrue(other)) {
				sp->board[v] = unmarker;
				sp->seen[v] = 0; // see no evil!
				nLitsRem++;
			}
		}
		if (nLitsRem == size) break; // bail out early
	}
	if (nLitsRem) {
		assert(end == learntC.end());
		j = end;
		uint32* i, *newend = j - nLitsRem;
		for (i = learntC + 1; i != newend; i++) {
			const uint32 lit = *i;
			if (NEQUAL(sp->board[ABS(lit)], marker))
				*i-- = *--j;
		}
		learntC.shrink(nLitsRem);
		sp->learntLBD = calcLBD();
		assert(sp->learntLBD <= opts.lbd_tier2);
		PFLOG2(4, " Learnt clause strengthened by %d binaries with new LBD of %d", nLitsRem, sp->learntLBD);
	}
}

bool Solver::minimize(const uint32& lit, const int& depth)
{
	CHECKLIT(lit);
	if (depth >= opts.minimize_depth) return false;
	const uint32 v = ABS(lit);
	const int litlevel = sp->level[v];
	if (!litlevel || (REMOVABLE(sp->seen[v]) && depth)) return true;
	const C_REF r = sp->source[v];
	if (DECISION(r) || POISONED(sp->seen[v])) return false;
	if (sp->vstate[litlevel].dlcount < 2) return false;
	assert(REASON(r));
	CLAUSE& c = cm[r];
	PFLCLAUSE(4, c, "  checking %d reason", -l2i(lit));
	bool gone = true;
	uint32* cend = c.end();
	for (uint32* k = c; gone && k != cend; k++) {
		const uint32 other = *k;
		if (NEQUAL(other, lit))
			gone = minimize(FLIP(other), depth + 1);
	}
	if (depth) sp->seen[v] = gone ? REMOVABLE_M : POISONED_M;
	else assert(REMOVABLE(sp->seen[v]));
	minimized.push(v);
	return gone;
}

void Solver::minimize()
{
	assert(learntC.size() > 1);
	assert(minimized.empty());
	uint32 *j = learntC + 1, *end = learntC.end();
	for (uint32* i = j; i != end; i++) {
		const uint32 lit = *i;
		if (!minimize(FLIP(lit)))
			*j++ = lit;
	}
	int shrunken = int(end - j);
	PFLOG2(4, " Learnt clause minimized by %d literals", shrunken);
	learntC.shrink(shrunken);
	clearMinimized();
}