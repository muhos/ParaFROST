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

    public:
                       RANDOM          () : state(0) {};
                       RANDOM          (const uint64& seed) : state(seed) {}
        inline void    init            (const uint64& seed) { state = seed; }
        inline uint64  seed            () const { return state; }
        inline uint64  next64          () {
            state *= 6364136223846793005ul;
            state += 1442695040888963407ul;
            assert(state);
            return state;
        }
        inline uint32  next32          () {
            return next64() >> 32;
        }
        inline bool    genbool          () {
            const uint32 next = next32();
            const double fraction = next / 4294967296.0;
            assert(fraction >= 0 && fraction < 1);
            const uint32 value = uint32(2 * fraction);
            assert(value < 2);
            return value;
        }
    };

}

#endif
