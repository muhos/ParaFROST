/***********************************************************************[pfmemory.cu]
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
		//=============================//
		//	  CUDA memory management   //
		//=============================//
		const size_t hc_srsize = sizeof(S_REF);
		const size_t hc_scsize = sizeof(SCLAUSE);
		const size_t hc_otsize = sizeof(OT);
		const size_t hc_olsize = sizeof(OL);
		const size_t hc_cnfsize = sizeof(CNF);
		const size_t hc_varsize = sizeof(uint32);
		const size_t hc_gstsize = sizeof(GSTATS);
		const size_t hc_cuvecsize = sizeof(cuVecU);

		__global__ void resizeCNF_k(CNF* cnf, S_REF d_size, uint32 cs_size)
		{
			cnf->resize(d_size, cs_size);
		}

		__global__ void assignListPtrs(OT* __restrict__ ot, const uint32* __restrict__ hist, const S_REF* __restrict__ segs, uint32 size)
		{
			uint32 tid = global_tx;
			while (tid < size) {
				assert(segs[tid] < UINT32_MAX);
				(*ot)[tid].alloc(ot->data(segs[tid]), hist[tid]);
				tid += stride_x;
			}
		}

		void cuMM::resizeCNFAsync(CNF* dcnf, const S_REF& data_size, const uint32& cs_size)
		{
			assert(dcnf);
			assert(data_size);
			assert(cs_size);
			resizeCNF_k << <1, 1 >> > (dcnf, data_size, cs_size);
			if (sync_always) {
				LOGERR("Resizing CNF failed");
				sync();
			}
		}

		uint32* cuMM::resizeLits(const size_t& min_lits)
		{
			assert(min_lits);
			const size_t min_cap = min_lits * hc_varsize;
			if (litsPool.cap < min_cap) {
				DFREE(litsPool);
				assert(litsPool.mem == NULL);
				if (!hasDeviceMem(min_cap, "Literals")) return NULL;
				CHECK(cudaMalloc((void**)&litsPool.mem, min_cap));
				litsPool.thrust_lits = t_iptr(litsPool.mem);
				litsPool.cap = min_cap;
				litsPool.size = min_lits;
			}
			return litsPool.mem;
		}

		addr_t cuMM::allocTemp(const size_t& min_cap)
		{
			assert(min_cap);
			if (tmpPool.cap < min_cap) {
				DFREE(tmpPool);
				assert(tmpPool.mem == NULL);
				if (!hasDeviceMem(min_cap, "Temporary")) return NULL;
				CHECK(cudaMalloc((void**)&tmpPool.mem, min_cap));
				tmpPool.cap = min_cap;
			}
			return tmpPool.mem;
		}

		bool cuMM::allocHist(cuHist& cuhist)
		{
			assert(inf.nDualVars == V2L(inf.maxVar + 1ULL));
			const size_t segBytes = inf.nDualVars * hc_srsize;
			const size_t histBytes = inf.nDualVars * hc_varsize;
			const size_t varsBytes = (inf.maxVar + 1) * hc_varsize;
			const size_t min_cap = segBytes + histBytes + varsBytes;
			assert(min_cap);
			if (hhistPool.cap < histBytes) {
				if (hhistPool.cap) {
					assert(hhistPool.mem != NULL);
					assert(cuhist.h_hist != NULL);
					CHECK(cudaFreeHost(hhistPool.mem));
					hhistPool.mem = NULL;
					hhistPool.cap = 0;
				}
				assert(hhistPool.mem == NULL);
				CHECK(cudaHostAlloc((void**)&hhistPool.mem, histBytes, cudaHostAllocDefault));
				cuhist.h_hist = (uint32*)hhistPool.mem;
				hhistPool.cap = histBytes;
			}
			if (histPool.cap < min_cap) {
				DFREE(histPool);
				assert(histPool.mem == NULL);
				if (!hasDeviceMem(min_cap, "Histogram")) return false;
				CHECK(cudaMalloc((void**)&histPool.mem, min_cap));
				// NOTE: d_segs, d_hist used internally by OT allocation and externally
				//       by BVE for calculating resolvents offsets (memory reuse)
				addr_t ea = histPool.mem;
				cuhist.d_segs = d_segs = (S_REF*)ea, ea += segBytes;
				cuhist.d_hist = d_hist = (uint32*)ea, ea += histBytes;
				cuhist.d_vorg = (uint32*)ea, ea += varsBytes;
				assert(ea == histPool.mem + min_cap);
				cuhist.thrust_hist = t_iptr(d_hist);
				histPool.cap = min_cap;
			}
			return true;
		}

		bool cuMM::allocAux(const size_t& clsCap)
		{
			assert(clsCap && clsCap <= UINT32_MAX);
			const size_t scatterBytes = clsCap * hc_srsize;
			const size_t min_cap = scatterBytes + clsCap;
			assert(min_cap);
			if (pinned_cnf == NULL)
				CHECK(cudaHostAlloc((void**)&pinned_cnf, hc_cnfsize, cudaHostAllocDefault));
			if (auxPool.cap < min_cap) {
				DFREE(auxPool);
				assert(auxPool.mem == NULL);
				if (!hasDeviceMem(min_cap, "Auxiliary")) return false;
				CHECK(cudaMalloc((void**)&auxPool.mem, min_cap));
				d_scatter = (S_REF*)auxPool.mem;
				d_stencil = auxPool.mem + scatterBytes;
				auxPool.cap = min_cap;
				nscatters = clsCap;
			}
			return true;
		}

		bool cuMM::resizeCNF(CNF*& cnf, const size_t& clsCap, const size_t& litsCap)
		{
			assert(clsCap && clsCap <= UINT32_MAX);
			assert(litsCap && litsCap <= UINT32_MAX);
			assert(litsCap >= clsCap);
			const size_t csBytes = clsCap * hc_srsize;
			const size_t dataBytes = clsCap * hc_scsize + litsCap * hc_bucket;
			assert(dataBytes % hc_bucket == 0);
			const size_t min_cap = hc_cnfsize + dataBytes + csBytes;
			assert(min_cap);
			if (cnfPool.cap == 0) {
				assert(cnf == NULL);
				assert(cnfPool.mem == NULL);
				if (!hasUnifiedMem(min_cap, "CNF")) return false;
				CHECK(cudaMallocManaged((void**)&cnfPool.mem, min_cap));
				if (devProp.major > 5) {
					PFLOGN2(2, " Advising GPU driver to favor global over system memory..");
					CHECK(cudaMemAdvise(cnfPool.mem, min_cap, cudaMemAdviseSetPreferredLocation, MASTER_GPU));
					PFLDONE(2, 5);
				}
				cnf = (CNF*)cnfPool.mem;
				const S_REF data_cap = S_REF(dataBytes / hc_bucket);
				new (cnf) CNF(data_cap, uint32(clsCap));
				d_cnf_mem = cnf->data().mem, d_refs_mem = cnf->refsData();
				cnfPool.cap = min_cap;
			}
			else {
				assert(cnf != NULL);
				assert(cnfPool.mem != NULL);
				if (!hasUnifiedMem(min_cap, "CNF")) return false;
				cacheCNFPtr(cnf);
				addr_t newMem = NULL;
				CHECK(cudaMallocManaged((void**)&newMem, min_cap));
				sync();
				if (devProp.major > 5) {
					PFLOGN2(2, " Advising GPU driver to favor global over system memory..");
					CHECK(cudaMemAdvise(newMem, min_cap, cudaMemAdviseSetPreferredLocation, MASTER_GPU));
					CHECK(cudaMemPrefetchAsync(newMem, min_cap, MASTER_GPU));
					PFLDONE(2, 5);
				}
				CNF* tmp_cnf = (CNF*)newMem;
				const S_REF data_cap = S_REF(dataBytes / hc_bucket);
				new (tmp_cnf) CNF(data_cap, uint32(clsCap));
				d_cnf_mem = tmp_cnf->data().mem, d_refs_mem = tmp_cnf->refsData();
				if (profile_gpu) cutimer->start();
				if (gc_gpu) compactCNF(cnf, tmp_cnf);
				else {
					sync(), tmp_cnf->copyFrom(cnf);
					pinned_cnf->resize(tmp_cnf->data().size, tmp_cnf->size());
				}
				if (profile_gpu) cutimer->stop(), cutimer->gc += cutimer->gpuTime();
				FREE(cnfPool);
				cnfPool.mem = newMem;
				cnfPool.cap = min_cap;
				cnf = tmp_cnf;
			}
			return true;
		}

		bool cuMM::resizeOTAsync(OT*& ot, const size_t& min_lits, const cudaStream_t& _s)
		{
			assert(d_hist != NULL);
			assert(d_segs != NULL);
			assert(min_lits && min_lits <= UINT32_MAX);
			uint32* tmp = resizeLits(min_lits);
			if (!tmp) return false;
			size_t ebytes = 0;
			DeviceScan::ExclusiveSum(NULL, ebytes, d_hist, d_segs, inf.nDualVars, _s);
			DeviceScan::ExclusiveSum(tmp, ebytes, d_hist, d_segs, inf.nDualVars, _s);
			if (!otBlocks) otBlocks = MIN((inf.nDualVars + BLOCK1D - 1) / BLOCK1D, maxGPUThreads / BLOCK1D);
			const size_t min_cap = hc_otsize + inf.nDualVars * hc_olsize + min_lits * hc_srsize;
			assert(min_cap);
			if (otPool.cap < min_cap) { // realloc
				FREE(otPool);
				assert(otPool.mem == NULL);
				if (!hasUnifiedMem(min_cap, "OT")) return false;
				CHECK(cudaMallocManaged((void**)&otPool.mem, min_cap));
				if (devProp.major > 5) {
					PFLOGN2(2, " Advising GPU driver to favor global over system memory..");
					CHECK(cudaMemAdvise(otPool.mem, min_cap, cudaMemAdviseSetPreferredLocation, MASTER_GPU));
					CHECK(cudaMemPrefetchAsync(otPool.mem, min_cap, MASTER_GPU, _s));
					PFLDONE(2, 5);
				}
				ot = (OT*)otPool.mem;
				LOGERR("Exclusively scanning histogram failed");
				sync(_s); // needed for calling the next constructor on host
				new (ot) OT(inf.nDualVars);
				assignListPtrs << <otBlocks, BLOCK1D, 0, _s >> > (ot, d_hist, d_segs, inf.nDualVars);
				otPool.cap = min_cap;
			}
			else
				assignListPtrs << <otBlocks, BLOCK1D, 0, _s >> > (ot, d_hist, d_segs, inf.nDualVars);
			if (sync_always) {
				LOGERR("Occurrence lists allocation failed");
				sync(_s);
			}
			return true;
		}

		bool cuMM::allocVars(VARS*& vars, const size_t& resCap)
		{
			assert(varsPool.cap == 0);
			assert(vars == NULL);
			assert(resCap && resCap <= UINT32_MAX);
			vars = new VARS();
			const size_t uintVec_sz = inf.maxVar * hc_varsize;
			const size_t varsize = inf.maxVar + 1;
			const size_t scores_sz = varsize * hc_varsize;
			const size_t resolved_sz = resCap * hc_varsize;
			const size_t min_cap = hc_gstsize + scores_sz + resolved_sz + (hc_cuvecsize + uintVec_sz) * 3;
			assert(min_cap);
			if (!hasUnifiedMem(min_cap, "Fixed")) return false;
			CHECK(cudaMallocManaged((void**)&varsPool.mem, min_cap));
			addr_t ea = varsPool.mem, end = ea + min_cap;
			vars->gstats = (GSTATS*)ea, ea += hc_gstsize;
			vars->pVars = (cuVecU*)ea, ea += hc_cuvecsize;
			vars->units = (cuVecU*)ea, ea += hc_cuvecsize;
			vars->resolved = (cuVecU*)ea, ea += hc_cuvecsize;
			uint32* uintPtr = (uint32*)ea;
			vars->pVars->alloc(uintPtr, inf.maxVar), uintPtr += inf.maxVar; d_units = uintPtr;
			vars->units->alloc(uintPtr, inf.maxVar), uintPtr += inf.maxVar;
			vars->eligible = uintPtr, uintPtr += inf.maxVar;
			vars->scores = uintPtr, uintPtr += varsize;
			vars->resolved->alloc(uintPtr, uint32(resCap)), uintPtr += resCap;
			assert((addr_t)uintPtr == end);
			varsPool.cap = min_cap;
			if (devProp.major > 5) {
				PFLOGN2(2, " Advising GPU driver to favor global over system memory..");
				CHECK(cudaMemAdvise(vars->gstats, hc_gstsize, cudaMemAdviseSetPreferredLocation, MASTER_GPU));
				addr_t tmpPtr = ea + uintVec_sz; // skip pVars
				CHECK(cudaMemAdvise(tmpPtr, end - tmpPtr, cudaMemAdviseSetPreferredLocation, MASTER_GPU));
				CHECK(cudaMemPrefetchAsync(tmpPtr, end - tmpPtr, MASTER_GPU));
				PFLDONE(2, 5);
			}
			if (pinned_units == NULL) {
				CHECK(cudaHostAlloc((void**)&pinned_units, uintVec_sz, cudaHostAllocDefault));
				vars->cachedUnits = pinned_units;
			}
			else vars->cachedUnits = pinned_units;
			return true;
		}

		void cuMM::createMirror(CNF*& hcnf, const size_t& clsCap, const size_t& litsCap)
		{
			assert(clsCap && clsCap <= UINT32_MAX);
			assert(litsCap && litsCap <= UINT32_MAX);
			assert(litsCap >= clsCap);
			const size_t csBytes = clsCap * hc_srsize;
			const size_t dataBytes = clsCap * hc_scsize + litsCap * hc_bucket;
			assert(dataBytes % hc_bucket == 0);
			const size_t min_cap = hc_cnfsize + dataBytes + csBytes;
			assert(min_cap);
			if (hcnfPool.cap < min_cap) {
				hcnfPool.cap = min_cap;
				pfralloc(hcnfPool.mem, hcnfPool.cap);
				hcnf = (CNF*)hcnfPool.mem;
			}
			const S_REF data_cap = S_REF(dataBytes / hc_bucket);
			new (hcnf) CNF(data_cap, uint32(clsCap));
		}

		void cuMM::mirrorCNF(CNF*& hcnf)
		{
			assert(cnfPool.cap);
			assert(cnfPool.mem != NULL);
			CHECK(cudaMemcpy(hcnf, cnfPool.mem, hc_cnfsize, cudaMemcpyDeviceToHost));
			const size_t csBytes = hcnf->size() * hc_srsize;
			const size_t dataBytes = hcnf->data().size * hc_bucket;
			const size_t min_cap = hc_cnfsize + dataBytes + csBytes;
			assert(min_cap <= cnfPool.cap);
			if (hcnfPool.cap < min_cap) {
				hcnfPool.cap = min_cap;
				pfralloc(hcnfPool.mem, hcnfPool.cap);
				hcnf = (CNF*)hcnfPool.mem;
			}
			hcnf->fixPointer(); // replace device with host pointers
		}

		void cuMM::freeVars() {
			PFLOGN2(2, " Freeing up fixed unified memory..");
			FREE(varsPool);
			d_units = NULL;
			PFLENDING(2, 5, "(remaining: %lld)", cap);
		}

		void cuMM::freeCNF() {
			PFLOGN2(2, " Freeing up CNF unified memory..");
			FREE(cnfPool);
			d_cnf_mem = NULL, d_refs_mem = NULL;
			PFLENDING(2, 5, "(remaining: %lld)", cap);
		}

		void cuMM::freeOT() {
			PFLOGN2(2, " Freeing up OT unified memory..");
			FREE(otPool);
			PFLENDING(2, 5, "(remaining: %lld)", cap);
		}

		void cuMM::freeFixed() {
			PFLOGN2(2, " Freeing up device memory..");
			if (auxPool.mem != NULL) {
				CHECK(cudaFree(auxPool.mem)), auxPool.mem = NULL;
				_free += auxPool.cap, auxPool.cap = 0;
				d_scatter = NULL, d_stencil = NULL;
				nscatters = 0;
			}
			if (histPool.mem != NULL) {
				CHECK(cudaFree(histPool.mem)), histPool.mem = NULL;
				_free += histPool.cap, histPool.cap = 0;
				d_segs = NULL, d_hist = NULL;
			}
			if (litsPool.mem != NULL) {
				CHECK(cudaFree(litsPool.mem)), litsPool.mem = NULL;
				_free += litsPool.cap, litsPool.cap = litsPool.size = 0;
			}
			if (tmpPool.mem != NULL) {
				CHECK(cudaFree(tmpPool.mem)), tmpPool.mem = NULL;
				_free += tmpPool.cap, tmpPool.cap = 0;
			}
			PFLENDING(2, 5, "(released cap: %lld)", dcap);
			dcap = 0;
		}

		void cuMM::freePinned() {
			if (pinned_cnf != NULL) CHECK(cudaFreeHost(pinned_cnf)), pinned_cnf = NULL;
			if (pinned_units != NULL) CHECK(cudaFreeHost(pinned_units)), pinned_units = NULL;
			if (hhistPool.mem != NULL) {
				CHECK(cudaFreeHost(hhistPool.mem)), hhistPool.mem = NULL;
				hhistPool.cap = 0;
			}
		}

		void cuMM::breakMirror() {
			if (hcnfPool.mem != NULL) {
				std::free(hcnfPool.mem), hcnfPool.mem = NULL;
				hcnfPool.cap = 0;
			}
		}

		void TCA::destroy() {
			for (freeBlock_t::iterator i = freeBlocks.begin(); i != freeBlocks.end(); i++) {
				thrust::cuda_cub::free(thrust::cuda_cub::pointer<void>(i->second));
				i->second = NULL;
			}
			for (allocBlock_t::iterator i = allocBlocks.begin(); i != allocBlocks.end(); i++)
				thrust::cuda_cub::free(thrust::cuda_cub::pointer<void>(i->first));
			freeBlocks.clear();
			allocBlocks.clear();
			used = 0;
		}

		char* TCA::allocate(int64 new_cap) {
			char* result = NULL;
			freeBlock_t::iterator freeBlock = freeBlocks.lower_bound(new_cap);
			// found free block
			if (freeBlock != freeBlocks.end()) {
				result = freeBlock->second;
				new_cap = freeBlock->first;
				freeBlocks.erase(freeBlock);
			}
			// no free blocks, allocate new one
			else {
				try {
					result = thrust::cuda_cub::malloc<char>(new_cap).get();
					used += new_cap;
				}
				catch (std::runtime_error&) {
					PFLOGE("cannot allocate new memory block for Thrust");
					throw;
				}
			}
			assert(result);
			allocBlocks.insert(std::make_pair(result, new_cap)); // cache new block
			return result;
		}

		void TCA::deallocate(char* ptr, size_t) {
			allocBlock_t::iterator allocBlock = allocBlocks.find(ptr);
			if (allocBlock == allocBlocks.end()) throw INVALID_PTR(ptr);
			const int64 new_cap = allocBlock->second;
			allocBlocks.erase(allocBlock);
			freeBlocks.insert(std::make_pair(new_cap, ptr)); // cache free block
		}

	}
}