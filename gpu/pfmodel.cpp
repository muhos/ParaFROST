/***********************************************************************[pfvec.h]
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

#include "pfmodel.h"

using namespace pFROST;
using namespace SIGmA;

void MODEL::extend(LIT_ST* currValue)
{
	uint32 updated = 0;
	value.resize(maxVar + 1, 0);
	for (uint32 v = 1; v <= maxVar; v++) {
		uint32 mlit = lits[v];
		if (mlit) {
			value[v] = currValue[mlit];
			updated++;
		}
	}
	PFLOG2(2, " Extending model updated %d mapped values.", updated);
	if (!resolved.empty()) {
		uint32 before = updated;
		PFLOGN2(2, " Extending model with eliminated variables..");
		if (verbose == 4) putc('\n', stdout);
		uint32* x = resolved.end(), lit;
		while (x != resolved) {
			bool sat = false;
			PFLOGN2(4, " Checking clause: ");
			while (lit = *--x) {
				if (verbose == 4) printLit(lit);
				if (sat) continue;
				if (satisfied(*x)) sat = true;
			}
			if (sat) while (*--x);
			else while (lit = *--x) value[l2a(*x)] = !sign(*x), updated++;
			if (verbose == 4) putc('\n', stdout);
		}
		PFLENDING(2, 4, "(%d updated)", updated - before);
	}
}
