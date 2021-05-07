/***********************************************************************[pfrephase.cpp]
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

inline void	ParaFROST::varOrgPhase() {
	stats.rephase.org++;
	memset(sp->psaved, opts.polarity, inf.maxVar + 1ULL);
	last.rephase.type = ORGPHASE;
}


inline void	ParaFROST::varInvPhase() {
	stats.rephase.inv++;
	memset(sp->psaved, !opts.polarity, inf.maxVar + 1ULL);
	last.rephase.type = INVPHASE;
}

inline void	ParaFROST::varFlipPhase() {
	stats.rephase.flip++;
	LIT_ST* end = sp->psaved + inf.maxVar + 1;
	for (LIT_ST* s = sp->psaved; s != end; s++)
		*s ^= NEG_SIGN;
	last.rephase.type = FLIPPHASE;
}

inline void	ParaFROST::varBestPhase() {
	stats.rephase.best++;
	forall_variables(v) {
		if (sp->pbest[v] != UNDEFINED)
			sp->psaved[v] = sp->pbest[v];
	}
	last.rephase.type = BESTPHASE;
}

inline void	ParaFROST::varRandPhase() {
	stats.rephase.random++;
	forall_variables(v) {
		sp->psaved[v] = random.genbool() ? 1 : 0;
	}
	last.rephase.type = RANDPHASE;
}

void ParaFROST::rephase()
{
	rootify();
	assert(cnfstate == UNSOLVED);
	stats.rephase.all++;
	memset(sp->ptarget, UNDEFINED, inf.maxVar + 1ULL);
	last.rephase.target = 0;
	const uint64 count = last.rephase.count++;
	if (!count) varOrgPhase();
	else if (count == 1) varInvPhase();
	else {
		switch ((count - 2) % 8)
		{
		default:
		case 0:
			varBestPhase();
			break;
		case 1:
			varOrgPhase();
			break;
		case 2:
			varBestPhase();
			break;
		case 3:
			varInvPhase();
			break;
		case 4:
			varBestPhase();
			break;
		case 5:
			varRandPhase();
			break;
		case 6:
			varBestPhase();
			break;
		case 7:
			varFlipPhase();
			break;
		}
	}
	assert(last.rephase.type);
	last.rephase.conflicts = stats.conflicts;
	INCREASE_LIMIT(this, rephase, stats.rephase.all, linear, false);
}