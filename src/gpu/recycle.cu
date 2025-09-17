/***********************************************************************[recycle.cu]
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

#include "grid.cuh"
#include "timer.cuh"
#include "memory.cuh"
#include "options.cuh"
#include "definitions.hpp"
#include <cub/device/device_scan.cuh>
#include <cub/device/device_select.cuh>

using namespace cub;

namespace ParaFROST {

	//=================================//
	//	CNF Garbage Collection on GPU  //
	//=================================//
	__global__ void scatter_k(const CNF* __restrict__ src, S_REF* __restrict__ scatter, addr_t __restrict__ stencil, const size_t max_size)
	{
		for_parallel_x(tid, src->size()) {
			const SCLAUSE& c = src->clause(tid);
			assert(tid < max_size);
			if (c.deleted()) { stencil[tid] = 0, scatter[tid] = 0; }
			else { 
				stencil[tid] = 1, scatter[tid] = c.size() + DC_NBUCKETS;
				assert(c.size() == scatter[tid] - DC_NBUCKETS);
			}
		}
	}

	__global__ void compact_k(CNF* __restrict__ src, CNF* __restrict__ dest, const S_REF* __restrict__ scatter, const addr_t __restrict__ stencil)
	{
		for_parallel_x(tid, src->size()) {
			if (stencil[tid]) {
				const S_REF new_r = scatter[tid];
				new (dest->cref(new_r)) SCLAUSE(src->clause(tid));
				assert(src->clause(tid).size() == dest->cref(new_r)->size());
				assert(src->clause(tid).capacity() == dest->cref(new_r)->capacity());
			}
			else assert(src->clause(tid).deleted());
		}
	}

	void cuMM::compactCNF(CNF* src, CNF* dest)
	{
		if (gopts.profile_gpu) cutimer->start();
		assert(src);
		assert(dest);

		const uint32 old_size = pinned_cnf->size();

		assert(old_size <= nscatters);
		assert(SCLAUSEBUCKETS == sizeof(SCLAUSE) / sizeof(uint32));

		LOG2(2, " Compacting simplified CNF (%d to %d) on GPU..", old_size, inf.numClauses);

		const S_REF data_size = REGIONBUCKETS(inf.numClauses, inf.numLiterals);
		resizeCNFAsync(dest, data_size, inf.numClauses);
		grid_t nThreads(BLOCK1D);
		OPTIMIZEBLOCKS(old_size, nThreads, 0);
		scatter_k << <nBlocks, nThreads >> > (src, d_scatter, d_stencil, nscatters);
		if (gopts.sync_always) {
			LASTERR("Scattering CNF failed");
			SYNC(0);
		}

		size_t ebytes = 0, fbytes = 0;
		DeviceScan::ExclusiveSum(NULL, ebytes, d_scatter, d_scatter, old_size);
		uint32* tmp = resizeLits(ebytes);
		if (!tmp) { throw MEMOUTEXCEPTION(); }
		DeviceScan::ExclusiveSum(tmp, ebytes, d_scatter, d_scatter, old_size);

		// *tmp will hold the new size, so tmp + 1 will be the start of the temporary array
		DeviceSelect::Flagged(NULL, fbytes, d_scatter, d_stencil, d_refs_mem, tmp, old_size);
		tmp = resizeLits(fbytes);
		if (!tmp) { throw MEMOUTEXCEPTION(); }
		DeviceSelect::Flagged(tmp + 1, fbytes, d_scatter, d_stencil, d_refs_mem, tmp, old_size);

		compact_k << <nBlocks, BLOCK1D >> > (src, dest, d_scatter, d_stencil);

		pinned_cnf->resize(data_size, inf.numClauses);

		LASTERR("CNF compact failed");
		SYNC(0);

		if (gopts.profile_gpu) cutimer->stop(), cutimer->gc += cutimer->gpuTime();
	}

}