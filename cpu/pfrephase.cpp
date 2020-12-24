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
#include "pfrandom.h"
using namespace pFROST;

inline void	ParaFROST::varOrgPhase() {
	memset(sp->psaved, opts.polarity, inf.maxVar + 1ULL);
	lrn.lastrephased = ORGPHASE;
}


inline void	ParaFROST::varInvPhase() {
	memset(sp->psaved, !opts.polarity, inf.maxVar + 1ULL);
	lrn.lastrephased = INVPHASE;
}

inline void	ParaFROST::varFlipPhase() {
	LIT_ST* end = sp->psaved + inf.maxVar + 1;
	for (LIT_ST* s = sp->psaved; s != end; s++)
		*s ^= NEG_SIGN;
	lrn.lastrephased = FLIPPHASE;
}

inline void	ParaFROST::varBestPhase() {
	for (uint32 v = 1; v <= inf.maxVar; v++)
		if (sp->pbest[v] != UNDEFINED)
			sp->psaved[v] = sp->pbest[v];
	lrn.lastrephased = BESTPHASE;
}

inline void	ParaFROST::varRandPhase() {
	stats.n_randrephs++;
	RANDOM random(opts.seed);
	random += stats.n_randrephs;
	for (uint32 v = 1; v <= inf.maxVar; v++)
		sp->psaved[v] = random.generate_bool() ? 0 : 1;
	lrn.lastrephased = RANDPHASE;
}

void ParaFROST::rephase() {

	stats.n_rephs++;
	backtrack();
	memset(sp->ptarget, UNDEFINED, inf.maxVar + 1ULL);
	lrn.target = 0;
	int64 count = lrn.rephased[lrn.stable]++;
	bool single = opts.stable_en && opts.vsidsonly_en ? true : !opts.stable_en;
	if (single) {
		switch (count % 8) {
			case 0: varInvPhase(); break;
			case 1: varBestPhase(); break;
			case 2: varFlipPhase(); break;
			case 3: varBestPhase(); break;
			case 4: varRandPhase(); break;
			case 5: varBestPhase(); break;
			case 6: varOrgPhase(); break;
			case 7: varBestPhase(); break;
			default: lrn.lastrephased = 0; break;
		}
	}
	else if (lrn.stable) {
		if (!count)	 varOrgPhase();
		else if (count == 1) varInvPhase();
		else {
			switch ((count - 2) % 4) {
				case 0:  varBestPhase(); break;
				case 1:  varOrgPhase(); break;
				case 2:  varBestPhase(); break;
				case 3:  varInvPhase(); break;
				default: lrn.lastrephased = 0; break;
			}
		}
	}
	else if (!lrn.stable) {
		if (!count) varFlipPhase();
		else {
			switch ((count - 1) % 4) {
				case 0:  varRandPhase(); break;
				case 1:  varBestPhase(); break;
				case 2:  varFlipPhase(); break;
				case 3:  varBestPhase(); break;
				default: lrn.lastrephased = 0; break;
			}
		}
	}
	assert(lrn.lastrephased);
	int64 current_inc = opts.rephase_inc * (stats.n_rephs + 1);
	lrn.rephase_last_max = nConflicts;
	lrn.rephase_conf_max = lrn.rephase_last_max + current_inc;
}