/***********************************************************************[pfrandom.h]
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

#ifndef __RANDOM_
#define __RANDOM_

#include <cassert>
#include "pfdtypes.h"

namespace pFROST {

    // Cadical-style random generator
    class RANDOM {
        uint64 state;
        void    add             (uint64 val) {
            if (!(state += val)) state = 1;
            next();
        }

    public:
        RANDOM                  () : state(0) {};
        RANDOM                  (const uint64& seed) : state(seed) { }
        void operator +=        (uint64 val) { add(val); }
        uint64  seed            () const { return state; }
        uint64  next            () {
            state *= 6364136223846793005ull;
            state += 1442695040888963407ull;
            assert(state);
            return state;
        }
        uint32  generate        () { next(); return state >> 32; }
        bool    generate_bool   () { return generate() < 2147483648u; }
    };

}

#endif
