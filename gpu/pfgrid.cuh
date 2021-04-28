/***********************************************************************[pfgrid.cuh]
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

#ifndef __GRID_INFO_
#define __GRID_INFO_

#include "pfdefinitions.cuh"
#include "pfdatatypes.h"

namespace pFROST {

	namespace SIGmA {
		// x
		_PFROST_D_ uint32 global_bx		() { return blockDim.x * blockIdx.x; }
		_PFROST_D_ uint32 global_bx_off	() { return (blockDim.x << 1) * blockIdx.x; }
		_PFROST_D_ uint32 global_tx		() { return global_bx() + threadIdx.x; }
		_PFROST_D_ uint32 global_tx_off	() { return global_bx_off() + threadIdx.x; }
		_PFROST_D_ uint32 stride_x		() { return blockDim.x * gridDim.x; }
		_PFROST_D_ uint32 stride_x_off	() { return (blockDim.x << 1) * gridDim.x; }
		// y
		_PFROST_D_ uint32 global_by		() { return blockDim.y * blockIdx.y; }
		_PFROST_D_ uint32 global_ty		() { return global_by() + threadIdx.y; }
		_PFROST_D_ uint32 stride_y		() { return blockDim.y * gridDim.y; }
	}
}

#endif