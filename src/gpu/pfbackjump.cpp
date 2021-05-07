/***********************************************************************[pfbackjump.cpp]
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

int ParaFROST::where()
{
	const int level = DL();
	int jmplevel = UNDEFINED;
	if (learntC.size() == 1) jmplevel = 0;
	else if (learntC.size() == 2) jmplevel = l2dl(learntC[1]);
	else {
		assert(learntC.size() > 2);
		assert(level > 1);
		int chronolevel = level - 1;
		uint32* maxPos = learntC + 1, maxLit = *maxPos, * end = learntC.end();
		jmplevel = l2dl(maxLit);
		for (uint32* k = learntC + 2; k != end; k++) {
			const uint32 lit = *k;
			const int litLevel = l2dl(lit);
			if (jmplevel >= litLevel) continue;
			jmplevel = litLevel;
			maxLit = lit;
			maxPos = k;
			if (litLevel == chronolevel) break;
		}
		*maxPos = learntC[1];
		learntC[1] = maxLit;
		PFLLEARNT(this, 3);
	}
	assert(level > jmplevel);
	if (opts.chrono_en && level - jmplevel > opts.chrono_min) {
		jmplevel = level - 1;
		stats.backtrack.chrono++;
		PFLOG2(3, "  forced chronological backtracking to level %d", jmplevel);
	}
	else if (opts.chrono_en && opts.chronoreuse_en) {
		uint32 bestvar = 0, bestpos = 0;
		if (vsidsEnabled()) {
			VSIDS_CMP hcmp(activity);
			for (uint32 i = dlevels[jmplevel + 1]; i < trail.size(); i++) {
				const uint32 v = ABS(trail[i]);
				if (bestvar && !hcmp(bestvar, v)) continue;
				bestvar = v;
				bestpos = i;
			}
		}
		else {
			for (uint32 i = dlevels[jmplevel + 1]; i < trail.size(); i++) {
				const uint32 v = ABS(trail[i]);
				if (bestvar && bumps[bestvar] >= bumps[v]) continue;
				bestvar = v;
				bestpos = i;
			}
		}
		CHECKVAR(bestvar);
		PFLOG2(4, "  found best variable %d at trail position %d", bestvar, bestpos);
		int old_jmplevel = jmplevel;
		while (jmplevel < level - 1 && dlevels[jmplevel + 1] <= bestpos) jmplevel++;
		if (old_jmplevel == jmplevel) stats.backtrack.nonchrono++;
		else stats.backtrack.chrono++;
	}
	else stats.backtrack.nonchrono++;
	assert(jmplevel != UNDEFINED);
	return jmplevel;
}

C_REF ParaFROST::backjump() {
	assert(trail.size());
	const int jmplevel = where();
	backtrack(jmplevel);
	if (learntC.size() == 1) {
		enqueue(learntC[0]);
		stats.units.learnt++;
	}
	else {
		if (opts.proof_en) proof.addClause(learntC);
		C_REF r = newClause(learntC, true);
		enqueue(*learntC, jmplevel, r);
		return r;
	}
	return NOREF;
}