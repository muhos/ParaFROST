/***********************************************************************[rephase.cpp]
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

#include "solver.h"
using namespace ParaFROST;

inline void	Solver::varOrgPhase() 
{
#ifdef STATISTICS
	stats.rephase.org++;
#endif
	memset(sp->psaved, opts.polarity, inf.maxVar + 1ULL);
	last.rephase.type = ORGPHASE;
}


inline void	Solver::varInvPhase() 
{
#ifdef STATISTICS
	stats.rephase.inv++;
#endif
	memset(sp->psaved, !opts.polarity, inf.maxVar + 1ULL);
	last.rephase.type = INVPHASE;
}

inline void	Solver::varFlipPhase()
{
#ifdef STATISTICS
	stats.rephase.flip++;
#endif
	LIT_ST* end = sp->psaved + inf.maxVar + 1;
	for (LIT_ST* s = sp->psaved; s != end; s++)
		*s ^= NEG_SIGN;
	last.rephase.type = FLIPPHASE;
}

inline void	Solver::varRandPhase()
{
#ifdef STATISTICS
	stats.rephase.random++;
#endif
	forall_variables(v) {
		sp->psaved[v] = random.brand();
	}
	last.rephase.type = RANDPHASE;
}

inline void	Solver::varBestPhase() 
{
#ifdef STATISTICS
	stats.rephase.best++;
#endif
	forall_variables(v) {
		if (!UNASSIGNED(sp->pbest[v]))
			sp->psaved[v] = sp->pbest[v];
	}
	last.rephase.type = BESTPHASE;
}

inline void	Solver::varWalkPhase() 
{
	walk();
	if (last.rephase.type == WALKPHASE) autarky();
	else {
		assert(!last.rephase.type);
		varBestPhase();
	}
}

void Solver::rephase()
{
	rootify();
	assert(UNSOLVED(cnfstate));
	stats.rephase.all++;
	memset(sp->ptarget, UNDEFINED, inf.maxVar + 1ULL);
	last.rephase.target = 0;
	const uint64 count = last.rephase.count++;
	// same order as Kissat solver but with
	// 'nlognlogn' interval and short steps
	if (!count) varOrgPhase();
	else if (count == 1) varInvPhase();
	else {
		switch ((count - 2) % 12)
		{
		default:
		case 0:
			varBestPhase();
			break;
		case 1:
			varWalkPhase();
			break;
		case 2:
			varOrgPhase();
			break;
		case 3:
			varBestPhase();
			break;
		case 4:
			varWalkPhase();
			break;
		case 5:
			varInvPhase();
			break;
		case 6:
			varBestPhase();
			break;
		case 7:
			varWalkPhase();
			break;
		case 8:
			varRandPhase();
			break;
		case 9:
			varBestPhase();
			break;
		case 10:
			varWalkPhase();
			break;
		case 11:
			varFlipPhase();
			break;
		}
	}
	assert(last.rephase.type);
	last.rephase.conflicts = stats.conflicts;
	INCREASE_LIMIT(rephase, stats.rephase.all, nlognlogn, false);
}