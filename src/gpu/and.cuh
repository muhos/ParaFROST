/***********************************************************************[and.cuh]
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

#ifndef __GPU_AND_
#define __GPU_AND_

#include "elimination.cuh"

namespace ParaFROST {

	#define AND_DBG 0

	_PFROST_D_ int find_fanin(const uint32& gate_out, CNF& cnf, OL& list, uint32* out_c, uint32& sig)
	{
		assert(gate_out > 1);
		sig = 0;
		int nImps = 0;
#pragma unroll
		forall_occurs(list, i) {
			SCLAUSE& c = cnf[*i];
			if (c.learnt()) continue;
			assert(!c.molten());
			if (c.size() == 2) {
				const uint32 imp = FLIP(c[0] ^ c[1] ^ gate_out);
				out_c[nImps++] = imp;
				sig |= MAPHASH(imp);
				c.melt(); // mark as gate clause
			}
		}
		return nImps;
	}

	_PFROST_D_ bool find_ao_gate(
		const uint32& dx,
		const uint32& nOrgCls,
		CNF& cnf, 
		OT& ot,
		uint32* out_c, 
		uint32& nElements,
		uint32& nAddedCls,
		uint32& nAddedLits)
	{
		assert(dx > 1);
		assert(checkMolten(cnf, ot[dx], ot[FLIP(dx)]));
		const uint32 x = ABS(dx);
		uint32 sig;
		// (-) ==> look for AND , (+) ==> look for OR
#if AND_DBG
		const char* type = SIGN(dx) ? "AND" : "OR";
#endif
		OL& itarget = ot[dx];
		int nImps = find_fanin(dx, cnf, itarget, out_c, sig);
		if (nImps > 1) {
			uint32 f_dx = FLIP(dx);
			out_c[nImps++] = f_dx;
			sig |= MAPHASH(f_dx);
			devSort(out_c, nImps);
			OL& otarget = ot[f_dx];
			forall_occurs(otarget, i) {
				SCLAUSE& c = cnf[*i];
				if (c.learnt()) continue;
				if (c.size() == nImps && sub(c.sig(), sig) && isEqual(c, out_c, nImps)) {
					c.melt(); // mark as fanout clause
					// check resolvability
					nElements = 0, nAddedCls = 0, nAddedLits = 0;
					if (countSubstituted(x, nOrgCls, cnf, itarget, otarget, nElements, nAddedCls, nAddedLits)) { 
						c.freeze();
						break; 
					}

					// can be substituted
					#if AND_DBG
					printf("c  Gate %d = %s(", ABS(dx), type);
					for (int k = 0; k < nImps; k++) {
						if (ABS(out_c[k]) == ABS(dx)) continue;
						printf(" %d", ABS(out_c[k]));
						if (k < nImps - 1) printf(",");
					}
					printf(" ) found ==> added = %d, deleted = %d\n", nAddedCls, itarget.size() + otarget.size());
					printGate(cnf, itarget, otarget);
					#endif

					return true;
				}
			}
		}
		freezeBinaries(cnf, itarget);
		return false;
	}

} // parafrost namespace


#endif