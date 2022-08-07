/***********************************************************************[atomics.cuh]
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

#include "definitions.cuh"
#include "datatypes.h"
#include "warp.cuh"

namespace pFROST {

	template<class T>
	_PFROST_D_ T atomicAggInc(T* counter) {
		const uint32 mask = __activemask(), total = __popc(mask);
		uint32 laneMask;
		laneMask_lt(laneMask);
		const T prefix = (T)__popc(mask & laneMask);
		const int lowest_lane = __ffs(mask) - 1;
		T warpRes = prefix ? 0 : atomicAdd(counter, total);
		warpRes = __shfl_sync(mask, warpRes, lowest_lane);
		return (prefix + warpRes);
	}

	template<class T, class R>
	_PFROST_D_ void atomicAggMax(T* counter, const R ref) {
		const uint32 mask = __activemask(), max_id = (32 - __clz(mask)) - 1;
		uint32 lane_id;
		laneId(lane_id);
		if (lane_id == max_id)
			atomicMax(counter, ref);
	}

}

#endif