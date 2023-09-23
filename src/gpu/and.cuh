/***********************************************************************[and.cuh]
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

#ifndef __AND_
#define __AND_

#include "elimination.cuh"

namespace ParaFROST {

	#define AND_DBG 0

	_PFROST_D_ int find_fanin(const uint32& gate_out, CNF& cnf, OL& list, uint32* out_c, uint32& sig)
	{
		assert(gate_out > 1);
		sig = 0;
		int nImps = 0;
		forall_occurs(list, i) {
			SCLAUSE& c = cnf[*i];
			assert(!c.molten());
			if (c.original() && c.size() == 2) {
				const uint32 imp = FLIP(c[0] ^ c[1] ^ gate_out);
				out_c[nImps++] = imp;
				sig |= MAPHASH(imp);
				c.melt(); // mark as gate clause
			}
		}
		return nImps;
	}

	_PFROST_D_ bool find_ao_gate(
		const uint32& dx, OL& dx_list, 
		const uint32& fx, OL& fx_list,
		const uint32& nOrgCls,
		CNF& cnf, 
		uint32* out_c, 
		uint32& nElements,
		uint32& nAddedCls,
		uint32& nAddedLits)
	{
		assert(dx > 1);
		assert(dx_list.size());
		assert(fx_list.size());

		if (cnf[*dx_list].size() > 2 || cnf[fx_list.back()].size() < 3)
			return false;

		assert(checkMolten(cnf, dx_list, fx_list));

		uint32 sig;
		int nImps = find_fanin(dx, cnf, dx_list, out_c, sig);

		if (nImps > 1) {

			const uint32 x = ABS(dx);

			out_c[nImps++] = fx;

			sig |= MAPHASH(fx);

			devSort(out_c, nImps);

			forall_occurs(fx_list, i) {

				SCLAUSE& c = cnf[*i];

				if (c.original() && c.size() == nImps && SUBSIG(c.sig(), sig) && isEqual(c, out_c, nImps)) {

					c.melt(); // mark as fanout clause

					// check resolvability
					nElements = 0, nAddedCls = 0, nAddedLits = 0;
					if (countSubstituted(x, nOrgCls, cnf, dx_list, fx_list, nElements, nAddedCls, nAddedLits)) { 
						c.freeze();
						break; 
					}

					// can be substituted
					#if AND_DBG
					printf("c  Gate %d = %s(", ABS(dx), SIGN(dx) ? "AND" : "OR");
					for (int k = 0; k < nImps; k++) {
						if (ABS(out_c[k]) == ABS(dx)) continue;
						printf(" %d", ABS(out_c[k]));
						if (k < nImps - 1) printf(",");
					}
					printf(" ) found ==> added = %d, deleted = %d\n", nAddedCls, dx_list.size() + fx_list.size());
					printGate(cnf, dx_list, fx_list);
					#endif

					return true;
				}
			}
		}

		freezeBinaries(cnf, dx_list);

		return false;
	}

} 


#endif