/***********************************************************************[grid.cuh]
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

#ifndef __GRID_INFO_
#define __GRID_INFO_

#include "definitions.cuh"
#include "datatypes.hpp"

namespace ParaFROST {

	// x
	typedef uint32 grid_t;

	#define global_bx		(grid_t)(blockDim.x * blockIdx.x)
	#define global_bx_off	(grid_t)((blockDim.x << 1) * blockIdx.x)
	#define global_tx		(grid_t)(global_bx + threadIdx.x)
	#define global_tx_off	(grid_t)(global_bx_off + threadIdx.x)
	#define stride_x        (grid_t)(blockDim.x * gridDim.x)
	#define stride_x_off	(grid_t)((blockDim.x << 1) * gridDim.x)
	// y
	#define global_by		(grid_t)(blockDim.y * blockIdx.y)
	#define global_ty		(grid_t)(global_by + threadIdx.y)
	#define stride_y		(grid_t)(blockDim.y * gridDim.y)

	#define ROUNDUP(DATALEN, DIVISOR) ((grid_t(DATALEN) + (DIVISOR) - 1) / (DIVISOR))

	#define for_parallel_x(IDX, SIZE) \
		for (grid_t IDX = global_tx, stride = stride_x, data_size = grid_t(SIZE); IDX < data_size; IDX += stride)

	#define for_parallel_x_double(IDX, SIZE) \
		for (grid_t IDX = global_tx_off, stride = stride_x_off, data_size = grid_t(SIZE); IDX < data_size; IDX += stride)

	#define for_parallel_x_tiled(BIDX, SIZE) \
		for (grid_t BIDX = blockIdx.x, MAX_BLOCKS_X = ROUNDUP(SIZE, blockDim.x); BIDX < MAX_BLOCKS_X; BIDX += gridDim.x)

	#define for_parallel_y(IDX, SIZE) \
		for (grid_t IDX = global_ty, stride = stride_y, data_size = grid_t(SIZE); IDX < data_size; IDX += stride)

	#define for_parallel_y_off(IDX, OFF, SIZE) \
		for (grid_t IDX = global_ty + OFF, stride = stride_y, data_size = grid_t(SIZE); IDX < data_size; IDX += stride)

	#define for_parallel_y_tiled(BIDX, SIZE) \
		for (grid_t BIDX = blockIdx.y, MAX_BLOCKS_Y = ROUNDUP(SIZE, blockDim.y); BIDX < MAX_BLOCKS_Y; BIDX += gridDim.y)


	#define WARP_ALIGN_DOWN(V,W)   ( ((V) / (W)) * (W) )
	#define WARP_ALIGN_UP(V,W)     ( ( ROUNDUP((V), (W)) * (W) ) )
	#define TPB_2D(BLOCK2D)        ( ((BLOCK2D).x) * ((BLOCK2D).y) )

	#define ROUNDUPBLOCKS(DATALEN, NTHREADS)  ROUNDUP((DATALEN), (NTHREADS))

	#define LOOP_SHRINK_YX_TO_FIT(BLOCK2D, SMEM_EXPR, BLOCKCAP, MAXSMEM, WARP)                  \
	do {                                                                                        \
		(BLOCK2D).x = MIN((BLOCK2D).x, WARP_ALIGN_DOWN((BLOCKCAP), (WARP)));                    \
		(BLOCK2D).x = MAX((BLOCK2D).x, (WARP));                                                 \
		(BLOCK2D).x = WARP_ALIGN_DOWN((BLOCK2D).x, (WARP));                                     \
		while ( ((size_t)(SMEM_EXPR) > (size_t)(MAXSMEM)) ||                                    \
				( TPB_2D(BLOCK2D) > (BLOCKCAP) ) ) {                                            \
			if ((BLOCK2D).y > 1) {                                                              \
				(BLOCK2D).y >>= 1;                                                              \
				continue;                                                                       \
			}                                                                                   \
			uint32 _maxX_ = WARP_ALIGN_DOWN((BLOCKCAP), (WARP));                                \
			(BLOCK2D).x = MIN((BLOCK2D).x, _maxX_);                                             \
			(BLOCK2D).x = MAX((BLOCK2D).x, (WARP));                                             \
			(BLOCK2D).x = WARP_ALIGN_DOWN((BLOCK2D).x, (WARP));                                 \
			if ((size_t)(SMEM_EXPR) > (size_t)(MAXSMEM)) {                                      \
				LOGERROR("insufficient shared memory");                                         \
			}                                                                                   \
		}                                                                                       \
	} while (0)

	#define LOOP_OVERSUB_Y(BLOCK2D, NVARS, MINTHREADS, TARGETBLOCKS,                             \
						SMEM_EXPR, BLOCKCAP, MAXSMEM, REALBLOCKS_LVAL)                           \
	do {                                                                                         \
		(REALBLOCKS_LVAL) = ROUNDUP((NVARS), (BLOCK2D).y);                                       \
		while ( (REALBLOCKS_LVAL) < (grid_t)(TARGETBLOCKS) && (BLOCK2D).y > 1 ) {                \
			unsigned _nextY_ = MAX((MINTHREADS), ((BLOCK2D).y >> 1));                            \
			if (_nextY_ == (BLOCK2D).y) break;                                                   \
			unsigned _prevY_ = (BLOCK2D).y;                                                      \
			(BLOCK2D).y = _nextY_;                                                               \
			if ( ((size_t)(SMEM_EXPR) > (size_t)(MAXSMEM)) ||                                    \
				( TPB_2D(BLOCK2D) > (BLOCKCAP) ) ) {                                             \
				(BLOCK2D).y = _prevY_;                                                           \
				break;                                                                           \
			}                                                                                    \
			(REALBLOCKS_LVAL) = ROUNDUP((NVARS), (BLOCK2D).y);                                   \
		}                                                                                        \
	} while (0)

	#define LOOP_SHRINK_NT_FOR_SMEM(NTHREADS, MINTHREADS, WARP, MAXSMEM, SMEM_EXPR)              \
	do {                                                                                         \
		(NTHREADS) = MAX((MINTHREADS), WARP_ALIGN_DOWN((NTHREADS), (WARP)));                     \
		size_t _smem_ = (size_t)(SMEM_EXPR);                                                     \
		while (_smem_ > (size_t)(MAXSMEM) && (NTHREADS) > (MINTHREADS)) {                        \
			uint32 _next_ = MAX((MINTHREADS), WARP_ALIGN_DOWN(((NTHREADS) >> 1), (WARP)));       \
			if (_next_ == (NTHREADS)) break;                                                     \
			uint32 _save_ = (NTHREADS); (NTHREADS) = _next_;                                     \
			size_t _probe_ = (size_t)(SMEM_EXPR); (NTHREADS) = _save_;                           \
			if (_probe_ > (size_t)(MAXSMEM)) break;                                              \
			(NTHREADS) = _next_;                                                                 \
			_smem_ = _probe_;                                                                    \
		}                                                                                        \
	} while (0)

	#define LOOP_OVERSUB_NT(NVARS, NTHREADS, MINTHREADS, WARP, TARGETBLOCKS,                      \
							MAXSMEM, SMEM_EXPR, REALBLOCKS_LVAL)                                  \
	do {                                                                                          \
		(REALBLOCKS_LVAL) = ROUNDUP((NVARS), (NTHREADS));                                         \
		while ( (REALBLOCKS_LVAL) < (grid_t)(TARGETBLOCKS) && (NTHREADS) > (MINTHREADS) ) {       \
			uint32 _next_ = MAX((MINTHREADS), WARP_ALIGN_DOWN(((NTHREADS) >> 1), (WARP)));        \
			if (_next_ == (NTHREADS)) break;                                                      \
			uint32 _save_ = (NTHREADS); (NTHREADS) = _next_;                                      \
			size_t _probe_ = (size_t)(SMEM_EXPR); (NTHREADS) = _save_;                            \
			if (_probe_ > (size_t)(MAXSMEM)) break;                                               \
			(NTHREADS) = _next_;                                                                  \
			(REALBLOCKS_LVAL) = ROUNDUP((NVARS), (NTHREADS));                                     \
		}                                                                                         \
	} while (0)

	#define OPTIMIZEBLOCKS(NVARS, NTHREADS, SMEMSIZE)                                            \
	grid_t nBlocks = 0, hardCap = 0;                                                             \
	do {                                                                                         \
		assert(NVARS);                                                                           \
		assert(NTHREADS);                                                                        \
		const uint32 warp     = (uint32)devProp.warpSize;                                        \
		const uint32 blockCap = (uint32)MIN(devProp.maxThreadsPerBlock, 1024);                   \
		const size_t maxSmem  = (devProp.sharedMemPerBlockOptin ?                                \
								devProp.sharedMemPerBlockOptin : devProp.sharedMemPerBlock);     \
		const uint32 SMs      = (uint32)devProp.multiProcessorCount;                             \
		const uint32 minThr   = (uint32)warp;                                                    \
		(NTHREADS) = MIN((uint32)(NTHREADS), blockCap);                                          \
		(NTHREADS) = MAX((uint32)(NTHREADS), minThr);                                            \
		LOOP_SHRINK_NT_FOR_SMEM((NTHREADS), minThr, warp, maxSmem, (SMEMSIZE));                  \
		grid_t realBlocks;                                                                       \
		const grid_t targetBlocks = (grid_t)SMs * 8;                                             \
		LOOP_OVERSUB_NT((NVARS), (NTHREADS), minThr, warp, targetBlocks,                         \
						maxSmem, (SMEMSIZE), realBlocks);                                        \
		const uint32 maxGrid  = (uint32)devProp.maxGridSize[0];                                  \
		hardCap = MIN(targetBlocks, maxGrid);                                                    \
		nBlocks = MIN(realBlocks, hardCap);                                                      \
	} while (0)

	#define OPTIMIZEBLOCKSELIM(NVARS, NTHREADS, SMEMSIZE, MINOPTS)                                 \
	grid_t nBlocks = 0, hardCap = 0;                                                               \
	do {                                                                                           \
		assert(NVARS);                                                                             \
		assert(NTHREADS);                                                                          \
		const grid_t MINTHREADS = MINOPTS ## _min_threads;                                         \
		const uint32 warp     = (uint32)devProp.warpSize;                                          \
		const uint32 blockCap = (uint32)MIN(devProp.maxThreadsPerBlock, 1024);                     \
		const size_t maxSmem  = (devProp.sharedMemPerBlockOptin ?                                  \
								devProp.sharedMemPerBlockOptin : devProp.sharedMemPerBlock);       \
		const uint32 SMs      = (uint32)devProp.multiProcessorCount;                               \
		const uint32 minThr   = (uint32)MAX(MINTHREADS, warp);                                     \
		(NTHREADS) = MIN((uint32)(NTHREADS), blockCap);                                            \
		(NTHREADS) = MAX((uint32)(NTHREADS), minThr);                                              \
		LOOP_SHRINK_NT_FOR_SMEM((NTHREADS), minThr, warp, maxSmem, (SMEMSIZE));                    \
		grid_t realBlocks;                                                                         \
		const grid_t targetBlocks = (grid_t)SMs * 8;                                               \
		LOOP_OVERSUB_NT((NVARS), (NTHREADS), minThr, warp, targetBlocks,                           \
						maxSmem, (SMEMSIZE), realBlocks);                                          \
		const uint32 maxGrid  = (uint32)devProp.maxGridSize[0];                                    \
		hardCap = MIN(targetBlocks, maxGrid);                                                      \
		nBlocks = MIN(realBlocks, hardCap);                                                        \
	} while (0)

	#define OPTIMIZEBLOCKSERE(NVARS, BLOCK2D, SMEM, MINOPTS)                                       \
	grid_t nBlocks = 0, hardCap = 0;                                                               \
	do {                                                                                           \
		assert((NVARS) > 0);                                                                       \
		assert((BLOCK2D).x > 0);                                                                   \
		assert((BLOCK2D).y > 0);                                                                   \
		const grid_t MINTHREADS = MINOPTS ## _min_threads;                                         \
		const uint32 warp     = (uint32)devProp.warpSize;                                          \
		const uint32 blockCap = (uint32)MIN(devProp.maxThreadsPerBlock, 1024);                     \
		const size_t maxSmem  = (devProp.sharedMemPerBlockOptin ?                                  \
								devProp.sharedMemPerBlockOptin : devProp.sharedMemPerBlock);       \
		const uint32 SMs      = (uint32)devProp.multiProcessorCount;                               \
		const grid_t minThreads = (grid_t)MAX(MINTHREADS, warp);                                   \
		grid_t realBlocks;                                                                         \
		LOOP_SHRINK_YX_TO_FIT((BLOCK2D), (SMEM), blockCap, maxSmem, warp);                         \
		const grid_t targetBlocks = (grid_t)SMs * 8;                                               \
		LOOP_OVERSUB_Y((BLOCK2D), (NVARS), minThreads, targetBlocks,                               \
					(SMEM), blockCap, maxSmem, realBlocks);                                        \
		const uint32 maxGrid  = (uint32)devProp.maxGridSize[1];                                    \
		hardCap = MIN(targetBlocks, maxGrid);                                                      \
		nBlocks = MIN(realBlocks, hardCap);                                                        \
	} while (0)

	#define OPTIMIZEBLOCKS2(DATALEN, NTHREADS)                               	\
			assert(DATALEN);                                                 	\
			assert(NTHREADS);                                                	\
			assert(maxGPUThreads);                                           	\
			grid_t nBlocks = MAXREDUCEBLOCKS + 1;								\
			grid_t multiplier = 0;												\
			uint32 blockSize = 0;                                      			\
			while (MAXREDUCEBLOCKS < nBlocks && blockSize < 1024) {				\
				blockSize = BLOCK1D << multiplier;								\
				assert(blockSize);												\
				const grid_t REALTHREADS = (blockSize) << 1;					\
				const grid_t REALBLOCKS = ROUNDUPBLOCKS(DATALEN, REALTHREADS);  \
				const grid_t MAXBLOCKS = maxGPUThreads / REALTHREADS;			\
				nBlocks = MIN(REALBLOCKS, MAXBLOCKS); 							\
				multiplier++;													\
			}  																	\

	// macros for shared memory calculation
    #define OPTIMIZESHARED(NTHREADS, MINCAP)                 	\
            assert(MINCAP);                                  	\
            assert(NTHREADS);                                	\
            assert(maxGPUSharedMem);                         	\
            const size_t smemSize = (NTHREADS) * (MINCAP);   	\

}

#endif