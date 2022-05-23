/***********************************************************************[grid.cuh]
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

#include "definitions.cuh"
#include "datatypes.h"

namespace pFROST {

	// x
	typedef uint32 gridtype;

	#define global_bx		(gridtype)(blockDim.x * blockIdx.x)
	#define global_bx_off	(gridtype)((blockDim.x << 1) * blockIdx.x)
	#define global_tx		(gridtype)(global_bx + threadIdx.x)
	#define global_tx_off	(gridtype)(global_bx_off + threadIdx.x)
	#define stride_x        (gridtype)(blockDim.x * gridDim.x)
	#define stride_x_off	(gridtype)((blockDim.x << 1) * gridDim.x)
	// y
	#define global_by		(gridtype)(blockDim.y * blockIdx.y)
	#define global_ty		(gridtype)(global_by + threadIdx.y)
	#define stride_y		(gridtype)(blockDim.y * gridDim.y)


	// macros for blocks calculation
	#define ROUNDUPBLOCKS(DATALEN, NTHREADS)							     \
			(((DATALEN) + (NTHREADS) - 1) / (NTHREADS))

	#define OPTIMIZEBLOCKS(DATALEN, NTHREADS)                                \
			assert(DATALEN);                                                 \
			assert(NTHREADS);                                                \
			assert(maxGPUThreads);                                           \
			const gridtype REALBLOCKS = ROUNDUPBLOCKS(DATALEN, NTHREADS);    \
			const gridtype MAXBLOCKS = maxGPUThreads / NTHREADS;             \
			const gridtype nBlocks = MIN(REALBLOCKS, MAXBLOCKS);             \

	#define OPTIMIZEBLOCKSELIM(NVARS, MAXTHREADS, MINOPTS)                   \
			assert(NVARS);                                                   \
			assert(MAXTHREADS);                                              \
			assert(maxGPUThreads);                                           \
			const gridtype MINTHREADS = MINOPTS ## _min_threads;		     \
			gridtype nThreads = MAXTHREADS;								     \
			gridtype realBlocks = ROUNDUPBLOCKS(NVARS, nThreads);            \
			const gridtype MAXBLOCKS = maxGPUThreads / MAXTHREADS;           \
			const gridtype MINBLOCKS =										 \
				  gridtype(MAXBLOCKS * (MINOPTS ## _min_blocks));		     \
			while (nThreads > MINTHREADS && realBlocks <= MINBLOCKS) {	     \
				nThreads >>= 1;											     \
				realBlocks = ROUNDUPBLOCKS(NVARS, nThreads);			     \
			}															     \
			const gridtype nBlocks = MIN(realBlocks, MAXBLOCKS);             \

	#define OPTIMIZEBLOCKSERE(NVARS, BLOCK2D, MINOPTS)                       \
			assert(NVARS);                                                   \
			assert(BLOCK2D.x);                                               \
			assert(BLOCK2D.y);                                               \
			assert(maxGPUThreads);                                           \
			const gridtype MINTHREADS = MINOPTS ## _min_threads;		     \
			gridtype realBlocks = ROUNDUPBLOCKS(NVARS, BLOCK2D.y);		     \
			const gridtype MAXBLOCKS = maxGPUThreads /					     \
										(BLOCK2D.x * BLOCK2D.y);             \
			const gridtype MINBLOCKS =										 \
				  gridtype(MAXBLOCKS * (MINOPTS ## _min_blocks));		     \
			while (BLOCK2D.y > MINTHREADS && realBlocks <= MINBLOCKS) {	     \
				BLOCK2D.y >>= 1;											 \
				realBlocks = ROUNDUPBLOCKS(NVARS, BLOCK2D.y);			     \
			}																 \
			const gridtype nBlocks = MIN(realBlocks, MAXBLOCKS);             \

	#define OPTIMIZEBLOCKS2(DATALEN, NTHREADS)                               \
			assert(DATALEN);                                                 \
			assert(NTHREADS);                                                \
			assert(maxGPUThreads);                                           \
			const gridtype REALTHREADS = (NTHREADS) << 1;					 \
			const gridtype REALBLOCKS = ROUNDUPBLOCKS(DATALEN, REALTHREADS); \
			const gridtype MAXBLOCKS = maxGPUThreads / REALTHREADS;          \
			const gridtype nBlocks = MIN(REALBLOCKS, MAXBLOCKS);             \

	// macros for shared memory calculation
    #define OPTIMIZESHARED(NTHREADS, MINCAP)                 \
            assert(MINCAP);                                  \
            assert(NTHREADS);                                \
            assert(maxGPUSharedMem);                         \
            const size_t smemSize = (NTHREADS) * (MINCAP);   \
            assert(maxGPUSharedMem >= smemSize);             \

}

#endif