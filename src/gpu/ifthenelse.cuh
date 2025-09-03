/***********************************************************************[ite.cuh]
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

#ifndef __ITE_
#define __ITE_

#include "elimination.cuh"

namespace ParaFROST {

	#define ITE_DBG 0

	_PFROST_D_ S_REF fast_equality_check(CNF& cnf, OT& ot, uint32 x, uint32 y, uint32 z) 
	{
		if (ot[y].size() > ot[z].size()) devSwap(y, z);
		if (ot[x].size() > ot[y].size()) devSwap(x, y);

		OL& list = ot[x];

		sort3(x, y, z, GPU_DEFAULT_CMP<uint32>());

		assert(x <= y && y <= z && x <= z);

		forall_occurs(list, i) {
			SCLAUSE& c = cnf[*i];
			if (c.molten()) continue;
			assert(c.isSorted());
			const uint32* lits = c.data();
			if (c.original() && c.size() == 3 && 
				lits[0] == x && lits[1] == y && lits[2] == z)
				return *i;
		}

		return GNOREF;
	}

	_PFROST_D_ bool find_ite_gate(
		const uint32& dx, OL& dx_list, 
		const uint32& fx, OL& fx_list,
		const uint32& nOrgCls,
		CNF& cnf,
		OT& ot,
		uint32& nElements,
		uint32& nAddedCls, 
		uint32& nAddedLits)
	{
		assert(dx > 1);

		assert(dx_list.size());

		if (cnf[dx_list.back()].size() == 2)
			return false;

		assert(checkMolten(cnf, dx_list, fx_list));

		const uint32 v = ABS(dx);

		S_REF* end = dx_list.end();
		for (S_REF* i = dx_list; i != end; ++i) {
			SCLAUSE& ci = cnf[*i];
			if (ci.original() && ci.size() == 3) {
				uint32 xi = ci[0], yi = ci[1], zi = ci[2];
				if (yi == dx) devSwap(xi, yi);
				if (zi == dx) devSwap(xi, zi);
				assert(xi == dx);
				for (S_REF* j = i + 1; j != end; j++) {
					SCLAUSE& cj = cnf[*j];
					if (cj.original() && cj.size() == 3) {
						assert(cj.original());
						uint32 xj = cj[0], yj = cj[1], zj = cj[2];
						if (yj == dx) devSwap(xj, yj);
						if (zj == dx) devSwap(xj, zj);
						assert(xj == dx);
						if (ABS(yi) == ABS(zj)) devSwap(yj, zj);
						if (ABS(zi) == ABS(zj)) continue;
						if (yi != FLIP(yj)) continue;
						
						S_REF r1 = fast_equality_check(cnf, ot, fx, yi, FLIP(zi));
						if (r1 == GNOREF) continue;
						S_REF r2 = fast_equality_check(cnf, ot, fx, yj, FLIP(zj));
						if (r2 == GNOREF) continue;
						assert(cnf[r1].original());
						assert(cnf[r2].original());

						// mark gate clauses
						ci.melt(), cj.melt();
						cnf[r1].melt(), cnf[r2].melt();

						// check resolvability
						nElements = 0, nAddedCls = 0, nAddedLits = 0;
						if (countSubstituted(v, nOrgCls, cnf, dx_list, fx_list, nElements, nAddedCls, nAddedLits)) {
							ci.freeze(), cj.freeze();
							cnf[r1].freeze(), cnf[r2].freeze();
							return false;
						}

						// can be substituted
						#if ITE_DBG
						printf("c  Gate %d = ITE(%d, %d, %d) found ==> added = %d, deleted = %d\n",
							ABS(dx), -int(ABS(yi)), -int(ABS(zi)), -int(ABS(zj)),
							nAddedCls, dx_list.size() + fx_list.size());
						print_gate(cnf, dx_list, fx_list);
						#endif

						return true;
					}
				}
			}
		}
		return false;
	}

} 

#endif