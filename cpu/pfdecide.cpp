/***********************************************************************[pfdecide.cpp]
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

uint32 ParaFROST::nextVSIDS()
{
	uint32 cand = 0;
	while (!vsids.empty()) {
		cand = vsids.top();
		assert(cand && cand <= inf.maxVar);
		if (!sp->locked[cand]) break;
		vsids.pop();
	}
	assert(cand);
	PFLOG2(4, " Next heap choice %d, activity %g", cand, activity[cand]);
	return cand;
}

uint32 ParaFROST::nextVMFQ()
{
	LIT_ST assigned = UNDEFINED;
	uint32 free = vmfq.free();
	assert(free);
	assigned = sp->locked[free];
	while (sp->locked[free]) free = vmfq.previous(free);
	assert(!assigned || assigned == 1);
	if (assigned) vmfq.update(free, bumps[free]);
	PFLOG2(4, " Next queue choice %d, bumped %lld", free, bumps[free]);
	return free;
}

uint32 ParaFROST::makeAssign(const uint32& v, const bool& tphase) {
	assert(v < NOVAR);
	LIT_ST pol = UNDEFINED;
	if (UNASSIGNED(pol) && tphase) pol = sp->ptarget[v];
	if (UNASSIGNED(pol)) pol = sp->psaved[v];
	if (UNASSIGNED(pol)) pol = opts.polarity;
	assert(pol >= 0);
	return V2DEC(v, pol);
}

void ParaFROST::decide()
{
	assert(trail.size() < inf.maxVar - inf.maxMelted);
	assert(sp->propagated == trail.size());
	assert(conflict == NOREF);
	assert(cnfstate == UNSOLVED);
	uint32 cand = vsidsEnabled() ? nextVSIDS() : nextVMFQ();
	uint32 dec = makeAssign(cand, useTarget());
	incDL();
	enqueue(dec, DL());
	stats.n_fuds++;
}