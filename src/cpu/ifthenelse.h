/***********************************************************************[ifthenelse.h]
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

#ifndef __ITE_
#define __ITE_

#include "simplify.h"
using namespace ParaFROST;

inline S_REF fast_equality_check(OT& ot, uint32 x, uint32 y, uint32 z)
{
	if (ot[y].size() > ot[z].size()) swap(y, z);
	if (ot[x].size() > ot[y].size()) swap(x, y);
	OL& list = ot[x];
	sort3(x, y, z);
	assert(x <= y && y <= z && x <= z);
	forall_occurs(list, i) {
		SCLAUSE& c = **i;
		if (c.molten()) continue;
		assert(c.isSorted());
		if (c.original() && c.size() == 3 &&
			c[0] == x && c[1] == y && c[2] == z) return *i;
	}
	return NULL;
}

inline bool find_ITE_gate(const uint32& dx, const int& nOrgCls, OT& ot, int& nAddedCls)
{
	const uint32 fx = FLIP(dx);
	assert(checkMolten(ot[dx], ot[fx]));
	OL& itarget = ot[dx];
	S_REF* end = itarget.end();
	for (S_REF* i = itarget; i != end; i++) {
		SCLAUSE& ci = **i;
		if (ci.original() && ci.size() == 3) {
			uint32 xi = ci[0], yi = ci[1], zi = ci[2];
			if (yi == dx) swap(xi, yi);
			if (zi == dx) swap(xi, zi);
			assert(xi == dx);
			for (S_REF* j = i + 1; j != end; j++) {
				SCLAUSE& cj = **j;
				if (cj.original() && cj.size() == 3) {
					assert(cj.original());
					uint32 xj = cj[0], yj = cj[1], zj = cj[2];
					if (yj == dx) swap(xj, yj);
					if (zj == dx) swap(xj, zj);
					assert(xj == dx);
					if (ABS(yi) == ABS(zj)) swap(yj, zj);
					if (ABS(zi) == ABS(zj)) continue;
					if (yi != FLIP(yj)) continue;
					S_REF d1 = fast_equality_check(ot, fx, yi, FLIP(zi));
					if (!d1) continue;
					S_REF d2 = fast_equality_check(ot, fx, yj, FLIP(zj));
					if (!d2) continue;
					assert(d1->original());
					assert(d2->original());
					// mark gate clauses
					ci.melt(), cj.melt();
					d1->melt(), d2->melt();
					// check resolvability
					uint32 v = ABS(dx);
					OL& otarget = ot[fx];
					nAddedCls = 0;
					if (countSubstituted(v, nOrgCls, itarget, otarget, nAddedCls)) {
						ci.freeze(), cj.freeze();
						d1->freeze(), d2->freeze();
						return false;
					}
					// can be substituted
					if (verbose >= 4) {
						PFLOG1(" Gate %d = ITE(%d, %d, %d) found ==> added = %d, deleted = %d",
							l2i(dx), -int(ABS(yi)), -int(ABS(zi)), -int(ABS(zj)), 
							nAddedCls, itarget.size() + otarget.size());
						printGate(itarget, otarget);
					}
					return true;
				}
			}
		}
	}
	return false;
}

#endif