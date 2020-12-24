/***********************************************************************[pfsclause.cpp]
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

void ParaFROST::newClause(SCLAUSE& s)
{
	if (s.deleted()) return;
	assert(s.size() > 1);
	assert(!s.molten());
	int sz = s.size();
	assert(sz > 1);
	// NOTE: 's' should be used with 'sp' before any mapping is done
	if (stats.sigmifications > 1 && s.added()) markSubsume(s);
	C_REF r = cm.alloc(sz);
	CLAUSE& new_c = cm[r];
	if (mapped) vmap.mapClause(new_c, s);
	else new_c.copyLitsFrom(s);
	assert(sz == new_c.size());
	assert(new_c.keep());
	assert(new_c[0] > 1 && new_c[1] > 1);
	assert(new_c[0] <= NOVAR && new_c[1] <= NOVAR);
	new_c.set_status(s.status());
	if (s.learnt()) {
		assert(s.lbd());
		assert(s.usage() <= USAGET2);
		assert(!s.added());
		int lbd = s.lbd();
		if (sz > 2 && lbd > opts.lbd_tier1)
			new_c.set_keep(0);
		new_c.set_lbd(lbd);
		new_c.set_usage(s.usage());
		learnts.push(r);
		inf.nLearntLits += sz;
	}
	else {
		assert(new_c.original());
		orgs.push(r);
		inf.nLiterals += sz;
	}
}

void ParaFROST::markSubsume(SCLAUSE& s) {
	assert(s.added());
	for (uint32* k = s; k != s.end(); k++)
		sp->subsume[ABS(*k)] = 1;
}