/***********************************************************************[pfite.cuh]
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

#ifndef __SIGMA_ITE_
#define __SIGMA_ITE_

#include "pfdevice.cuh"

namespace pFROST {

	namespace SIGmA {

		_PFROST_D_ S_REF fast_equality_check(CNF& cnf, OT& ot, uint32 x, uint32 y, uint32 z) {
			if (ot[y].size() > ot[z].size()) devSwap(y, z);
			if (ot[x].size() > ot[y].size()) devSwap(x, y);
			OL& list = ot[x];
			sort3(x, y, z, GPU_DEFAULT_CMP<uint32>());
			assert(x <= y && y <= z && x <= z);
#pragma unroll
			forall_occurs(list, i) {
				SCLAUSE& c = cnf[*i];
				if (c.molten()) continue;
				assert(c.isSorted());
				if (c.original() && c.size() == 3 &&
					c[0] == x && c[1] == y && c[2] == z) return *i;
			}
			return GNOREF;
		}

		_PFROST_D_ bool find_ite_gate(const uint32& dx, const uint32& nOrgCls, CNF& cnf, OT& ot, uint32& nAddedCls, uint32& nAddedLits)
		{
			assert(dx > 1);
			assert(checkMolten(cnf, ot[dx], ot[FLIP(dx)]));
			OL& itarget = ot[dx];
			forall_occurs(itarget, i) {
				SCLAUSE& ci = cnf[*i];
				if (ci.learnt() || ci.size() < 3 || ci.size() > 3) continue;
				assert(ci.original());
				uint32 xi = ci[0], yi = ci[1], zi = ci[2];
				if (yi == dx) devSwap(xi, yi);
				if (zi == dx) devSwap(xi, zi);
				assert(xi == dx);
				for (S_REF* j = i + 1; j != itarget.end(); j++) {
					SCLAUSE& cj = cnf[*j];
					if (cj.learnt() || cj.size() < 3 || cj.size() > 3) continue;
					assert(cj.original());
					uint32 xj = cj[0], yj = cj[1], zj = cj[2];
					if (yj == dx) devSwap(xj, yj);
					if (zj == dx) devSwap(xj, zj);
					assert(xj == dx);
					if (ABS(yi) == ABS(zj)) devSwap(yj, zj);
					if (ABS(zi) == ABS(zj)) continue;
					if (yi != FLIP(yj)) continue;
					const uint32 f_dx = FLIP(dx);
					S_REF r1 = fast_equality_check(cnf, ot, f_dx, yi, FLIP(zi));
					if (r1 == GNOREF) continue;
					S_REF r2 = fast_equality_check(cnf, ot, f_dx, yj, FLIP(zj));
					if (r2 == GNOREF) continue;
					assert(cnf[r1].original());
					assert(cnf[r2].original());
					// mark gate clauses
					ci.melt(), cj.melt();
					cnf[r1].melt(), cnf[r2].melt();
					// check resolvability
					const uint32 v = ABS(dx);
					OL& otarget = ot[f_dx];
					nAddedCls = 0, nAddedLits = 0;
					if (countSubstituted(v, nOrgCls, cnf, itarget, otarget, nAddedCls, nAddedLits)) {
						ci.freeze(), cj.freeze();
						cnf[r1].freeze(), cnf[r2].freeze();
						return false;
					}
					// can be substituted
#if VE_DBG
					printf("c  Gate %d = ITE(%d, %d, %d) found ==> added = %d, deleted = %d", ABS(dx), ABS(yi), ABS(zi), ABS(zj), nAddedCls, itarget.size() + otarget.size());
					printGate(cnf, itarget, otarget);
#endif
					return true;
				}
			}
			return false;
		}

	} // sigma namespace
} // parafrost namespace


#endif