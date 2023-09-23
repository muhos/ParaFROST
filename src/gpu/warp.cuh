/***********************************************************************[warp.cuh]
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

#ifndef __GL_WARP_
#define __GL_WARP_

#include "definitions.cuh"
#include "datatypes.hpp"

namespace ParaFROST {

#define FULLWARP 0xFFFFFFFFU

	_PFROST_IN_D_ void laneId(uint32& id) {
		asm("mov.u32 %0, %%laneid;" : "=r"(id));
	}

	_PFROST_IN_D_ void laneMask_lt(uint32& lanemask) {
		asm("mov.u32 %0, %%lanemask_lt;" : "=r"(lanemask));
	}

	_PFROST_IN_D_ uint32 blockWarpId() {
		return threadIdx.x >> 5;
	}

	_PFROST_IN_D_ uint32 blockWarps() {
		return blockDim.x >> 5;
	}

}

#endif