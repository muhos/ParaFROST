/***********************************************************************[backtrack.cpp]
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

inline void Solver::savePhases() {
    const LIT_ST reset = (last.rephase.type && stats.conflicts > last.rephase.conflicts);
    if (!probed) {
        if (reset) last.rephase.target = 0;
        if (sp->trailpivot > last.rephase.target) {
            last.rephase.target = sp->trailpivot;
            memcpy(sp->ptarget, sp->psaved, inf.maxVar + 1ULL);
        }
        if (sp->trailpivot > last.rephase.best) {
            last.rephase.best = sp->trailpivot;
            memcpy(sp->pbest, sp->psaved, inf.maxVar + 1ULL);
        }
        sp->trailpivot = 0;
    }
    if (reset) last.rephase.type = 0;
}

inline void Solver::cancelAssign(const uint32& lit) {
    CHECKLIT(lit);
    assert(inf.unassigned < inf.maxVar);
    sp->value[lit] = UNDEFINED;
    sp->value[FLIP(lit)] = UNDEFINED;
    inf.unassigned++;
    PFLOG2(4, "  literal %d@%d cancelled", l2i(lit), l2dl(lit));
}

void Solver::backtrack(const int& jmplevel) {
    if (DL() == jmplevel) return;
    const uint32 pivot = jmplevel + 1;
    PFLOG2(3, " Backtracking to level %d, at trail index %d", jmplevel, dlevels[pivot]);
    savePhases();
    const uint32 from = dlevels[pivot];
    uint32 *i = trail + from, *j = i, *end = trail.end();
    if (stable) {
        while (i != end) {
            const uint32 lit = *i++, v = ABS(lit);
            if (sp->level[v] > jmplevel) {
                cancelAssign(lit);
                if (!vsids.has(v)) vsids.insert(v);
            } else
                *j++ = lit;
        }
    } else {
        while (i != end) {
            const uint32 lit = *i++, v = ABS(lit);
            if (sp->level[v] > jmplevel) {
                cancelAssign(lit);
                if (vmtf.bumped() < bumps[v]) vmtf.update(v, bumps[v]);
            } else
                *j++ = lit;
        }
    }
    const uint32 remained = uint32(j - trail);
    PFLOG2(3, "  %d literals kept (%d are saved) and %zd are cancelled", remained, remained - from, end - j);
    trail.resize(remained);
    if (sp->propagated > from) sp->propagated = from;
    dlevels.resize(pivot);
    assert(DL() == jmplevel);
}