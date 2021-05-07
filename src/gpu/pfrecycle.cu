/***********************************************************************[pfrecycle.cu]
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

#include "pfmemory.cuh"
#include "pfdevice.cuh"
#include <cub/device/device_scan.cuh>
#include <cub/device/device_select.cuh>
using namespace cub;

namespace pFROST {

	namespace SIGmA {
		//=================================//
		//	CNF Garbage Collection on GPU  //
		//=================================//
		__global__ void scatter_k(const CNF* __restrict__ src, S_REF* __restrict__ scatter, addr_t __restrict__ stencil)
		{
			uint32 tid = global_tx;
			while (tid < src->size()) {
				const SCLAUSE& c = src->clause(tid);
				if (c.deleted()) stencil[tid] = 0, scatter[tid] = 0;
				else stencil[tid] = 1, scatter[tid] = c.size() + dc_nbuckets;
				tid += stride_x;
			}
		}

		__global__ void compact_k(CNF* __restrict__ src, CNF* __restrict__ dest, const S_REF* __restrict__ scatter, const addr_t __restrict__ stencil)
		{
			uint32 tid = global_tx;
			while (tid < src->size()) {
				if (stencil[tid]) {
					const S_REF new_r = scatter[tid];
					new (dest->cref(new_r)) SCLAUSE(src->clause(tid));
					assert(src->clause(tid).size() == dest->cref(new_r)->size());
					assert(src->clause(tid).capacity() == dest->cref(new_r)->capacity());
				}
				tid += stride_x;
			}
		}

		void cuMM::scatterCNF(CNF* src, S_REF* scatter, Byte* stencil)
		{
			assert(src);
			assert(scatter);
			assert(stencil);
			assert(maxGPUThreads);
			const uint32 old_size = pinned_cnf->size();
			uint32 nThreads = BLOCK1D, maxBlocks = maxGPUThreads / nThreads;
			uint32 nBlocks = MIN((old_size + nThreads - 1) / nThreads, maxBlocks);
			scatter_k << <nBlocks, nThreads >> > (src, scatter, stencil);
			if (sync_always) {
				LOGERR("Scattering CNF failed");
				sync();
			}
		}

		void cuMM::compactCNF(CNF* src, CNF* dest)
		{
			assert(src);
			assert(dest);
			const uint32 old_size = pinned_cnf->size();
			assert(old_size <= nscatters);
			assert(hc_nbuckets == sizeof(SCLAUSE) / sizeof(uint32));
			PFLOG2(2, " Compacting simplified CNF (%d to %d) on GPU..", old_size, inf.nClauses);
			const S_REF data_size = inf.nClauses * hc_nbuckets + inf.nLiterals;
			resizeCNFAsync(dest, data_size, inf.nClauses);
			scatterCNF(src, d_scatter, d_stencil);
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
			OPTIMIZEBLOCKS(old_size, BLOCK1D);
			compact_k << <nBlocks, BLOCK1D >> > (src, dest, d_scatter, d_stencil);
			pinned_cnf->resize(data_size, inf.nClauses);
			LOGERR("CNF compact failed");
			sync();
		}

	}
}