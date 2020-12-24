/***********************************************************************[pfminimize.cpp]
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

void ParaFROST::minimizeBin()
{
	markLearnt();
	uint32 uip = FLIP(learntC[0]);
	WL& ws = wt[uip];
	int nLitsRem = 0;
	for (WATCH* w = ws; w != ws.end(); w++) {
		if (!w->binary()) continue;
		uint32 other = w->imp;
		uint32 v = ABS(other);
		if (sp->marks[v] == !SIGN(other)) {
			sp->marks[v] = UNDEFINED; // unmark
			nLitsRem++;
		}
	}
	if (nLitsRem) {
		uint32* i, * j = learntC.end(), * newend = j - nLitsRem;
		for (i = learntC + 1; i != newend; i++) {
			uint32 lit = *i;
			if (sp->marks[ABS(lit)] != SIGN(lit))
				swap(*--j, *i--);
		}
		learntC.shrink(nLitsRem);
		for (i = learntC + 1; i != learntC.end(); i++)
			sp->marks[ABS(*i)] = UNDEFINED; // unmark
		PFLOG2(4, " Learnt clause minimized by %d literals using binary strengthening", nLitsRem);
	}
	else
		unmarkLearnt();
}

bool ParaFROST::minimize(const uint32& lit, const int& depth)
{
	register uint32 v = ABS(lit);
	C_REF r = sp->source[v];
	if (!sp->level[v] || REMOVABLE(sp->seen[v]) || KEPT(sp->seen[v])) return true;
	if (DECISION(r) || POISONED(sp->seen[v]) || sp->level[v] == DL()) return false;
	if (depth > opts.minimize_depth) return false;
	assert(r != NOREF);
	CLAUSE& c = cm[r];
	PFLCLAUSE(4, c, "  checking %d reason", -l2i(lit));
	bool gone = true;
	if (c.binary()) gone = minimize(FLIP(c[0] ^ c[1] ^ lit), depth + 1);
	else {
		for (uint32* k = c; gone && k != c.end(); k++) {
			uint32 other = *k;
			if (other != lit) gone = minimize(FLIP(other), depth + 1);
		}
	}
	sp->seen[v] = gone ? REMOVABLE_M : POISONED_M;
	minimized.push(v);
	if (!depth) PFLOG2(4, "  self-subsuming %d %s", l2i(lit), gone ? "succeeded" : "failed");
	return gone;
}

void ParaFROST::minimize()
{
	assert(learntC.size() > 1);
	assert(minimized.empty());
	uint32* i = learntC, * j = i, * end = learntC.end();
	for (; i != end; i++) if (!minimize(FLIP(*i))) sp->seen[ABS(*j++ = *i)] = KEEP_M;
	PFLOG2(4, " Learnt clause minimized by %d literals", int(learntC.end() - j));
	int newSize = int(j - learntC);
	learntC.resize(newSize);
	if (newSize > 1 && newSize <= opts.minimizebin_max) minimizeBin();
}