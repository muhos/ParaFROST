/***********************************************************************[pfatomics.cuh]
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

#ifndef __GL_ATOMIC_
#define __GL_ATOMIC_

#include "pfcudefs.h"
#include "pfdtypes.h"
#include "pfwarp.cuh"

namespace pFROST {

    namespace SIGmA {

        _PFROST_D_ uint32 atomicAggInc(uint32* counter) {
            uint32 mask = __activemask(), total = __popc(mask), laneMask;
            laneMask_lt(laneMask);
            uint32 prefix = __popc(mask & laneMask);
            int lowest_lane = __ffs(mask) - 1;
            uint32 warpRes = 0;
            if (prefix == 0) warpRes = atomicAdd(counter, total);
            warpRes = __shfl_sync(mask, warpRes, lowest_lane);
            uint32 teread_offset = prefix + warpRes;
            return teread_offset;
        }

        _PFROST_D_ int atomicAggInc(int* counter) {
            uint32 mask = __activemask(), total = __popc(mask), laneMask;
            laneMask_lt(laneMask);
            uint32 prefix = __popc(mask & laneMask);
            int lowest_lane = __ffs(mask) - 1;
            int warpRes = 0;
            if (prefix == 0) warpRes = atomicAdd(counter, total);
            warpRes = __shfl_sync(mask, warpRes, lowest_lane);
            int teread_offset = prefix + warpRes;
            return teread_offset;
        }

    }

}

#endif