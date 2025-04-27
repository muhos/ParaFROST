/***********************************************************************[elimination.cuh]
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

#ifndef __GPU_REDUCE_
#define __GPU_REDUCE_

#include "grid.cuh"

namespace ParaFROST {

	template<typename T>
	_PFROST_IN_D_ void warpReduce(volatile T* smem, T& val)
	{
		if (threadIdx.x < 32) {
			if (blockDim.x >= 64) val += smem[threadIdx.x + 32];
			unsigned mask = __activemask();
			val += __shfl_down_sync(mask, val, 16);
			val += __shfl_down_sync(mask, val, 8);
			val += __shfl_down_sync(mask, val, 4);
			val += __shfl_down_sync(mask, val, 2);
			val += __shfl_down_sync(mask, val, 1);
		}
	}

	template<typename T1, typename T2>
	_PFROST_IN_D_ void warpReduce(volatile T1* smem1, T1& val1, volatile T2* smem2, T2& val2) 
	{
		if (threadIdx.x < 32) {
			if (blockDim.x >= 64) {
				val1 += smem1[threadIdx.x + 32];
				val2 += smem2[threadIdx.x + 32];
			}
			unsigned mask = __activemask();
			if (blockDim.x >= 32) {
				val1 += __shfl_down_sync(mask, val1, 16);
				val2 += __shfl_down_sync(mask, val2, 16);
			}
			if (blockDim.x >= 16) {
				val1 += __shfl_down_sync(mask, val1, 8);
				val2 += __shfl_down_sync(mask, val2, 8);
			}
			if (blockDim.x >= 8) {
				val1 += __shfl_down_sync(mask, val1, 4);
				val2 += __shfl_down_sync(mask, val2, 4);
			}
			if (blockDim.x >= 4) {
				val1 += __shfl_down_sync(mask, val1, 2);
				val2 += __shfl_down_sync(mask, val2, 2);
			}
			if (blockDim.x >= 2) {
				val1 += __shfl_down_sync(mask, val1, 1);
				val2 += __shfl_down_sync(mask, val2, 1);
			}
		}
	}

	template<typename T, typename S>
	_PFROST_IN_D_ void loadShared(T* smem, const T& val, const S& size) 
	{
		smem[threadIdx.x] = (threadIdx.x < size) ? val : 0;
		__syncthreads();
	}

	template<typename T1, typename T2, typename S>
	_PFROST_IN_D_ void loadShared(T1* smem1, const T1& val1, T2* smem2, const T2& val2, const S& size) 
	{
		if (threadIdx.x < size) { smem1[threadIdx.x] = val1, smem2[threadIdx.x] = val2; }
		else { smem1[threadIdx.x] = 0, smem2[threadIdx.x] = 0; }
		__syncthreads();
	}

	template<typename T>
	_PFROST_IN_D_ void sharedReduce(T* smem, T& val) {
		if (blockDim.x >= 512) {
			if (threadIdx.x < 256) smem[threadIdx.x] = val = val + smem[threadIdx.x + 256];
			__syncthreads();
		}
		if (blockDim.x >= 256) {
			if (threadIdx.x < 128) smem[threadIdx.x] = val = val + smem[threadIdx.x + 128];
			__syncthreads();
		}
		if (blockDim.x >= 128) {
			if (threadIdx.x < 64) smem[threadIdx.x] = val = val + smem[threadIdx.x + 64];
			__syncthreads();
		}
	}

	template<typename T1, typename T2>
	_PFROST_IN_D_ void sharedReduce(T1* smem1, T1& val1, T2* smem2, T2& val2)
	{
		if (blockDim.x >= 512) {
			if (threadIdx.x < 256) {
				smem1[threadIdx.x] = val1 = val1 + smem1[threadIdx.x + 256];
				smem2[threadIdx.x] = val2 = val2 + smem2[threadIdx.x + 256];
			}
			__syncthreads();
		}
		if (blockDim.x >= 256) {
			if (threadIdx.x < 128) {
				smem1[threadIdx.x] = val1 = val1 + smem1[threadIdx.x + 128];
				smem2[threadIdx.x] = val2 = val2 + smem2[threadIdx.x + 128];
			}
			__syncthreads();
		}
		if (blockDim.x >= 128) {
			if (threadIdx.x < 64) {
				smem1[threadIdx.x] = val1 = val1 + smem1[threadIdx.x + 64];
				smem2[threadIdx.x] = val2 = val2 + smem2[threadIdx.x + 64];
			}
			__syncthreads();
		}
	}

} 

#endif