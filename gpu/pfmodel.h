/***********************************************************************[pfsolvertypes.h]
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

#ifndef __VMODEL_
#define __VMODEL_

#include "pfvec.h"
#include "pfdefs.h"

namespace pFROST {
	namespace SIGmA {
		struct MODEL {
			Vec<LIT_ST> value;
			uVec1D lits, resolved;
			uint32 maxVar;
			MODEL() : maxVar(0) {}
			void		init() {
				assert(inf.maxVar);
				PFLOG2(2, " Initially mapping original variables to literals..");
				maxVar = inf.maxVar;
				lits.resize(maxVar + 1), lits[0] = 0;
				for (uint32 v = 1; v <= inf.maxVar; v++) lits[v] = v2l(v);
			}
			void		print() {
				PFLMH('v');
				for (uint32 v = 1; v <= maxVar; v++)
					PFLMLIT(v, value[v]);
				putc('\n', stdout);
				if (!quiet_en) PFLOG0("");
			}
			void		extend(LIT_ST*);
			__forceinline
			bool		satisfied(const uint32& lit) { return value[l2a(lit)] == !sign(lit); }
		};
	}
}

#endif
