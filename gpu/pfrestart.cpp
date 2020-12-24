/***********************************************************************[pfrestart.cpp]
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

bool ParaFROST::vibrate() {
	if (!opts.stable_en) return false;
	if (vsidsOnly()) return true;
	if (nConflicts >= lrn.stable_conf_max) {
		lrn.stable = !lrn.stable;
		opts.stabrestart_inc *= opts.stabrestart_rate;
		if (opts.stabrestart_inc > int64(2e9))
			opts.stabrestart_inc = int64(2e9);
		lrn.stable_conf_max = nConflicts + opts.stabrestart_inc;
		if (lrn.stable_conf_max <= nConflicts)
			lrn.stable_conf_max = nConflicts + 1;
		PFLOG2(2, " stable restarts limit increased to %lld conflicts", lrn.stable_conf_max);
		lbdrest.swap();
		printStats();
	}
	return lrn.stable;
}

int ParaFROST::reuse() {
	bool stable = vsidsEnabled();
	uint32 cand = stable ? nextVSIDS() : nextVMFQ();
	assert(cand);
	int currLevel = DL(), target = 0;
	if (stable) {
		HEAP_CMP hcmp(activity);
		double candAct = activity[cand];
		while (target < currLevel) {
			uint32 pivot = dlevels[target + 1];
			if (hcmp(cand, ABS(trail[pivot])))
				target++;
			else break;
		}
	}
	else {
		int64 candBump = bumps[cand];
		while (target < currLevel) {
			uint32 pivot = dlevels[target + 1];
			if (candBump < bumps[ABS(trail[pivot])])
				target++;
			else break;
		}
	}
	if (target) stats.reuses++;
	return target;
}

void ParaFROST::restart()
{
	assert(sp->propagated == trail.size());
	assert(conflict == NOREF);
	assert(cnfstate == UNSOLVED);
	if (lrn.stable) stats.stab_restarts++;
	starts++;
	backtrack(opts.reusetrail_en ? reuse() : 0);
	if (opts.mdmfuses_en) MDMFuseSlave();
	lrn.restarts_conf_max = nConflicts + opts.restart_inc;
}

// Donald Knuth version of Luby restart as in Cadical ..
// implemented here to avoid overwhelming warnings about "-u"
void LUBYREST::update() {
	if (!period || restart) return;
	if (--countdown) return;
	if ((u & -u) == v) u++, v = 1;
	else v <<= 1;
	if (limited && v >= limit) u = v = 1;
	countdown = v * period;
	restart = true;
}