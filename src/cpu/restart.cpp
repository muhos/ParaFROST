/***********************************************************************[restart.cpp]
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

bool Solver::canRestart()
{
	if (!DL()) return false;
	if (stats.conflicts < limit.restart.conflicts) return false;
	vibrate();
	if (stable) return lubyrest;
	return lbdrest.restart(stable);
}

void Solver::restart()
{
	assert(sp->propagated == trail.size());
	assert(conflict == NOREF);
	assert(UNSOLVED(cnfstate));
	stats.restart.all++;
	backtrack(reuse());
	if (stable) stats.restart.stable++;
	else updateUnstableLimit();
}

void Solver::updateUnstableLimit()
{
	assert(!stable);
	uint64 increase = opts.restart_inc - 1;
	increase += logn(stats.restart.all);
	limit.restart.conflicts = stats.conflicts + increase;
}

int Solver::reuse() 
{
	bool stable = vsidsEnabled();
	uint32 cand = stable ? nextVSIDS() : nextVMFQ();
	assert(cand);
	int currLevel = DL(), target = 0;
	if (stable) {
		VSIDS_CMP hcmp(activity);
		while (target < currLevel) {
			uint32 pivot = dlevels[target + 1];
			if (hcmp(cand, ABS(trail[pivot])))
				target++;
			else break;
		}
	}
	else {
		uint64 candBump = bumps[cand];
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

#if defined(_WIN32)
#pragma warning(push)
#pragma warning(disable : 4146)
#endif

void LUBYREST::update()
{
	if (!period || restart) return;
	assert(countdown > 0);
	if (--countdown) return;
	if ((u & -u) == v) u++, v = 1;
	else v <<= 1;
	assert(v);
	assert((UINT64_MAX / v) >= period);
	countdown = v * period;
	restart = true;
	// reset if limit is reached
	if (limited && countdown > limit) {
		u = v = 1;
		countdown = period;
	}
}

#if defined(_WIN32)
#pragma warning(pop)
#endif