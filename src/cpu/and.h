/***********************************************************************[and.h]
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

#ifndef __AND_
#define __AND_

#include "simplify.h" 
using namespace pFROST;

namespace SIGmA {

	inline void find_fanin(const uint32& gate_out, OL& list, Lits_t& out_c, uint32& sig)
	{
		assert(gate_out > 1);
		out_c.clear();
		sig = 0;
		uint32 imp = 0;
		for (S_REF* i = list; i != list.end(); i++) {
			SCLAUSE& c = **i;
			assert(!c.molten());
			if (c.original() && c.size() == 2) {
				imp = FLIP(c[0] ^ c[1] ^ gate_out);
				out_c.push(imp);
				sig |= MAPHASH(imp);
				c.melt(); // mark as gate clause
			}
		}
	}

	inline bool find_AO_gate(const uint32& dx, const int& nOrgCls, OT& ot, Lits_t& out_c, int& nAddedCls)
	{
		assert(dx > 1);
		assert(checkMolten(ot[dx], ot[FLIP(dx)]));
		out_c.clear();
		uint32 sig, x = ABS(dx);
		// (-) ==> look for AND , (+) ==> look for OR
#if VE_DBG
		const char* type = SIGN(dx) ? "AND" : "OR";
#endif
		OL& itarget = ot[dx];
		find_fanin(dx, itarget, out_c, sig);
		if (out_c.size() > 1) {
			uint32 f_dx = FLIP(dx);
			out_c.push(f_dx);
			sig |= MAPHASH(f_dx);
			Sort(out_c, LESS<uint32>());
			OL& otarget = ot[f_dx];
			for (int i = 0; i < otarget.size(); i++) {
				SCLAUSE& c = *otarget[i];
				if (c.original() && c.size() == out_c.size() && sub(c.sig(), sig) && isEqual(c, out_c)) {
					c.melt(); // mark as fanout clause
					// check resolvability
					nAddedCls = 0;
					if (countSubstituted(x, nOrgCls, itarget, otarget, nAddedCls)) { c.freeze(); break; }
					// can be substituted
#if VE_DBG
					PFLOGN1(" Gate %d = %s(", ABS(dx), type);
					for (int k = 0; k < out_c.size(); k++) {
						if (ABS(out_c[k]) == ABS(dx)) continue;
						fprintf(stdout, " %d", ABS(out_c[k]));
						if (k < out_c.size() - 1) putc(',', stdout);
					}
					fprintf(stdout, " ) found ==> added = %d, deleted = %d\n", nAddedCls, itarget.size() + otarget.size());
					printGate(itarget, otarget);
#endif
					return true;
				}
			}
		}
		freezeBinaries(itarget);
		return false;
	}

}

#endif