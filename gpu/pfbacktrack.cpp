/***********************************************************************[pfbacktrack.cpp]
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

inline void ParaFROST::savePhases(const int& bt_level) {

	LIT_ST reset = (lrn.lastrephased && nConflicts > lrn.rephase_last_max);
	if (reset) {
		lrn.target = 0;
		if (lrn.lastrephased == BESTPHASE) lrn.best = 0;
	}
	if (sp->trailpivot > lrn.target) {
		lrn.target = sp->trailpivot;
		savePhases(sp->ptarget);
	}
	if (sp->trailpivot > lrn.best) {
		lrn.best = sp->trailpivot;
		savePhases(sp->pbest);
	}
	if (reset) lrn.lastrephased = 0;
}

inline void	ParaFROST::cancelAssign(const uint32& lit) {
	assert(lit > 1);
	sp->value[lit] = UNDEFINED;
	sp->value[FLIP(lit)] = UNDEFINED;
	PFLOG2(4, " Literal %d@%d cancelled", l2i(lit), l2dl(lit));
}

void ParaFROST::backtrack(const int& bt_level)
{
	if (DL() == bt_level) return;
	int pivot = bt_level + 1;
	PFLOG2(3, " Backtracking to level %d, at trail index %d", bt_level, dlevels[pivot]);
	savePhases(bt_level);
	uint32 from = dlevels[pivot], i = from, j = from;
	while (i < trail.size()) {
		uint32 lit = trail[i++], v = ABS(lit);
		if (sp->level[v] > bt_level) {
			cancelAssign(lit);
			sp->locked[v] = 0;
			if (!vsids.has(v)) vsids.insert(v);
			if (vmfq.bumped() < bumps[v]) vmfq.update(v, bumps[v]);
		}
		else assert(opts.chrono_en), trail[j++] = lit;
	}
	PFLOG2(3, " %d literals kept (%d are saved) and %d are cancelled", j, j - from, trail.size() - j);
	trail.resize(j);
	if (sp->propagated > from) sp->propagated = from;
	if (sp->trailpivot > from) sp->trailpivot = from;
	dlevels.resize(pivot);
	assert(DL() == bt_level);
}