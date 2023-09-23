/***********************************************************************[simplify.cuh]
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

#ifndef __GPU_SIMP_
#define __GPU_SIMP_

#include <cuda_runtime.h>
#include <thrust/sort.h>
#include <thrust/binary_search.h>
#include <thrust/adjacent_difference.h>
#include <thrust/iterator/counting_iterator.h>
#include "thrustalloc.cuh"
#include "variables.cuh"
#include "histogram.cuh"
#include "options.cuh"
#include "printer.cuh"
#include "memory.cuh"
#include "timer.cuh"
#include "cache.cuh"
#include "proof.cuh"
#include "table.cuh"
#include "key.cuh"
#include "cnf.cuh"
#include "vstate.hpp"

namespace ParaFROST {
	//======================================================//
	//                GPU Wrappers Declaration              //
	//======================================================//
	void printConstants();
	void initSharedMem();
	void initDevOpts();
	void initDevVorg(const cuHist&);
	void mapFrozenAsync(VARS*, const uint32&);
	void cuMemSetAsync(addr_t, const Byte&, const size_t&);
	void copyIf(uint32*, CNF*);
	void copyIfAsync(uint32*, CNF*);
	void calcScores(VARS*, uint32*);
	void calcScores(VARS*, uint32*, OT*);
	void prepareCNFAsync(CNF*, const cudaStream_t&);
	void createOTAsync(CNF*, OT*, const bool&);
	void reduceOTAsync(CNF*, OT*, const bool&);
	void veAsync(CNF*, OT*, VARS*, cudaStream_t*, cuVecB*, cuMM&, const cuHist&, const bool&);
	void veResizeCNFAsync(CNF*, const cuHist&);
	void subAsync(CNF*, OT*, VARS*, cuVecB*);
	void bceAsync(CNF*, OT*, VARS*, cuVecB*);
	void ereAsync(CNF*, OT*, VARS*, cuVecB*);
	void parcountCls(CNF*);
	void parcountAll(CNF*);
	void parcountLits(CNF*);

}

#endif 
