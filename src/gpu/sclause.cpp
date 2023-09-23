/***********************************************************************[sclause.cpp]
Copyright(c) 2020, Muhammad Osama - Anton Wijs,
Copyright(c) 2022-present, Muhammad Osama.

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

void Solver::newClause(SCLAUSE& s)
{
	int size = s.size();
	assert(size > 1);
	assert(!s.deleted());
	assert(!s.molten());	
	// NOTE: 's' should be used with 'sp' before any mapping is done
	if (stats.sigma.calls > 1 && s.added()) markSubsume(s);
	C_REF r = cm.alloc(size);
	CLAUSE& new_c = cm[r];
	if (mapped) vmap.mapClause(new_c, s);
	else new_c.copyLitsFrom(s);
	assert(size == new_c.size());
	assert(new_c.keep());
	assert(new_c[0] > 1 && new_c[1] > 1);
	assert(new_c[0] <= NOVAR && new_c[1] <= NOVAR);
	assert(!new_c.deleted());
	assert(s.status() == ORIGINAL || s.status() == LEARNT);
	if (s.learnt()) {
		assert(s.lbd());
		assert(s.usage() <= USAGET2);
		assert(!s.added());
		new_c.markLearnt();
		int lbd = s.lbd();
		if (size > 2 && lbd > opts.lbd_tier1) new_c.set_keep(0);
		new_c.set_lbd(lbd);
		new_c.set_usage(s.usage());
		learnts.push(r);
		stats.literals.learnt += size;
	}
	else {
		assert(new_c.original());
		orgs.push(r);
		stats.literals.original += size;
	}
}

void Solver::markSubsume(SCLAUSE& s) {
	assert(s.added());
	forall_clause(s, k) {
		markSubsume(*k);
	}
}