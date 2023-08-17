/***********************************************************************[mode.cpp]
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


void Solver::vibrate() {

	assert(UNSOLVED(cnfstate));
	if (last.rephase.type) return;
	if (!opts.stable_en) return;
	if (stable) {
		if (stats.searchticks < limit.mode.ticks) return;
	}
	else if (stats.conflicts < limit.mode.conflicts) return;
	if (stable) unstableMode();
	else stableMode();
}

void Solver::stableMode()
{
	assert(!stable);
	stats.stablemodes++;
	stable = true;
	updateModeLimit();
	if (opts.luby_inc) 
		lubyrest.enable(opts.luby_inc, opts.luby_max);
	printStats();
	updateHeap();
}

void Solver::unstableMode()
{
	assert(stable);
	stats.unstablemodes++;
	stable = false;
	updateModeLimit();
	updateQueue();
	updateUnstableLimit();
}

void Solver::updateModeLimit()
{
	assert(opts.stable_en);
	if (stable) {
		assert(stats.searchticks > stats.mode.ticks);
		uint64 ticks = stats.searchticks - stats.mode.ticks;
		ticks *= opts.stable_rate;
		limit.mode.ticks = stats.searchticks + ticks;
		PFLOG2(2, " Stable mode limit increased to %lld ticks", limit.mode.ticks);
	}
	else {
		const uint64 scaled = quadratic(stats.stablemodes) + 1;
		const uint64 increase = opts.mode_inc * scaled;
		limit.mode.conflicts = stats.conflicts + increase;
		PFLOG2(2, " Unstable mode limit increased to %lld conflicts", limit.mode.conflicts);
	}
	last.rephase.type = 0;
	stats.mode.ticks = stats.searchticks;
}
