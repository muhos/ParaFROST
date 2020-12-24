/***********************************************************************[pfextend.cpp]
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

#ifndef __EXTEND_
#define __EXTEND_

#include "pfmodel.h"

using namespace pFROST;

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
    if (resolved.size()) {
        uint32 before = updated;
        PFLOGN2(2, " Extending model with eliminated variables..");
        uint32 *x = resolved.end() - 1, k;
        while (x > resolved) {
            bool unsat = true;
            for (k = *x--; k > 1; k--, x--) {
                if (satisfied(*x)) { unsat = false; break; }
            }
            if (unsat) value[ABS(*x)] = !SIGN(*x), updated++;
            x -= k;
        }
        PFLENDING(2, 5, "(%d updated)", updated - before);
    }
}

#endif
