/***********************************************************************[pfcualloc.cu]
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

#include "pfsolve.h"
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

		__global__ void resizeCNF_k(CNF* cnf, CNF* hcnf) { cnf->resize(hcnf->data().size, hcnf->size()); }
		__global__ void resizeCNF_k(CNF* cnf, S_REF d_size, uint32 cs_size) { cnf->resize(d_size, cs_size); }
		__global__ void scatter_k(CNF* src, S_REF* scatter, addr_t stencil) {
			uint32 tid = global_tx();
			while (tid < src->size()) {
				SCLAUSE& c = src->clause(tid);
				if (c.deleted()) stencil[tid] = 0, scatter[tid] = 0;
				else stencil[tid] = 1, scatter[tid] = (c.size() - 1) + dc_nbuckets;
				tid += stride_x();
			}
		}
		__global__ void compact_k(CNF* src, CNF* dest, S_REF* scatter, addr_t stencil) {
			uint32 tid = global_tx();
			while (tid < src->size()) {
				if (stencil[tid]) {
					S_REF new_r = scatter[tid];
					new (dest->cref(new_r)) SCLAUSE(src->clause(tid));
					assert(src->clause(tid).size() == dest->cref(new_r)->size());
				}
				tid += stride_x();
			}
		}
		__global__ void assignListPtrs(OT* ot, uint32* hist, S_REF* segs, uint32 size)
		{
			uint32 tid = global_tx();
			while (tid < size) {
				assert(segs[tid] < UINT32_MAX);
				(*ot)[tid].alloc(ot->data(segs[tid]), hist[tid]);
				tid += stride_x();
			}
		}
		
		void cuMM::compactCNF(CNF* src, CNF* dest) {
			uint32 old_size = pinned_cnf->size();
			assert(old_size <= nscatters);
			assert(hc_nbuckets == sizeof(SCLAUSE) / sizeof(uint32));
			PFLOGN2(2, " Compacting simplified CNF (%d to %d) on GPU..", old_size, inf.nClauses);
			S_REF data_size = inf.nClauses * hc_nbuckets + (inf.nLiterals - inf.nClauses);
			resizeCNF_k << <1, 1 >> > (dest, data_size, inf.nClauses);
			size_t tb = 0, ftb = 0;
			uint32* ts = d_lits;
			uint32 nTereads = BLOCK1D, maxBlocks = maxGPUTereads / nTereads;
			uint32 nBlocks = std::min((old_size + nTereads - 1) / nTereads, maxBlocks);
			scatter_k << <nBlocks, nTereads >> > (src, d_scatter, d_stencil);
			DeviceScan::ExclusiveSum(NULL, tb, d_scatter, d_scatter, old_size), assert(tb <= litsbytes);
			DeviceScan::ExclusiveSum(ts, tb, d_scatter, d_scatter, old_size);
			DeviceSelect::Flagged(NULL, ftb, d_scatter, d_stencil, d_cs_mem, ts, old_size), assert(ftb <= litsbytes);
			DeviceSelect::Flagged(ts + 1, ftb, d_scatter, d_stencil, d_cs_mem, ts, old_size);
			compact_k << <nBlocks, nTereads >> > (src, dest, d_scatter, d_stencil);
			pinned_cnf->resize(data_size, inf.nClauses);
			LOGERR(" CNF compact failed");
			sync();
			PFLDONE(2, 5);
		}

		bool cuMM::allocVars(VARS*& vars, const size_t& resCap) {
			assert(vars == NULL);
			assert(resCap && resCap <= UINT32_MAX);
			vars = new VARS();
			size_t uintVec_sz = inf.maxVar * hc_varsize;
			size_t varsize = inf.maxVar + 1;
			size_t scores_sz = varsize * hc_varsize;
			size_t resolved_sz = resCap * hc_varsize;
			size_t newCap = hc_gstsize + scores_sz + resolved_sz + (hc_cuvecsize + uintVec_sz) * 3;
			varsPool.cap = newCap, cap += newCap;
			assert(varsPool.cap);
			if (!hasUnifiedMem("Fixed")) { cap -= newCap; return false; }
			CHECK(cudaMallocManaged((void**)&varsPool.mem, varsPool.cap));
			addr_t ea = varsPool.mem, end = ea + varsPool.cap;
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

		bool cuMM::allocHist(cuHist& cuhist, const size_t& litsCap)
		{
			assert(litsCap && litsCap <= UINT32_MAX);
			assert(inf.nDualVars == V2L(inf.maxVar + 1ULL));
			litsbytes = litsCap * hc_varsize;
			size_t segBytes = inf.nDualVars * hc_srsize;
			size_t histBytes = inf.nDualVars * hc_varsize;
			size_t newCap = segBytes + histBytes + litsbytes;
			assert(newCap);
			if (hhistPool.cap < histBytes) {
				if (hhistPool.cap) {
					assert(hhistPool.mem != NULL);
					assert(cuhist.h_hist != NULL);
					CHECK(cudaFreeHost(hhistPool.mem));
					hhistPool.mem = NULL;
				}
				assert(hhistPool.mem == NULL);
				CHECK(cudaHostAlloc((void**)&hhistPool.mem, histBytes, cudaHostAllocDefault));
				cuhist.h_hist = (uint32*)hhistPool.mem;
				hhistPool.cap = histBytes;
			}
			if (histPool.cap < newCap) {
				if (histPool.cap) {
					assert(histPool.mem != NULL);
					assert(d_hist != NULL);
					assert(d_segs != NULL);
					assert(cuhist.d_hist != NULL);
					assert(cuhist.d_segs != NULL);
					assert(cuhist.d_lits != NULL);
					CHECK(cudaFree(histPool.mem));
					histPool.mem = NULL;
					dcap -= histPool.cap;
					assert(dcap >= 0);
					_free += histPool.cap;
				}
				assert(histPool.mem == NULL);
				dcap += newCap;
				if (!hasDeviceMem("Histogram")) { dcap -= newCap; return false; }
				CHECK(cudaMalloc((void**)&histPool.mem, newCap));
				// NOTE: d_segs, d_hist used internally by OT allocation and externally by BVE for calculating resolvents offsets (memory reuse)
				//		 d_lits is used as temporary storage as well for CUB routines
				addr_t ea = histPool.mem;
				cuhist.d_segs = d_segs = (S_REF*)ea, ea += segBytes;
				cuhist.d_hist = d_hist = (uint32*)ea, ea += histBytes;
				cuhist.d_lits = d_lits = (uint32*)ea, ea += litsbytes;
				assert(ea == histPool.mem + newCap);
				cuhist.thrust_hist = t_iptr(d_hist);
				cuhist.thrust_lits = t_iptr(d_lits);
				histPool.cap = newCap;
			}
			return true;
		}

		bool cuMM::allocAux(const size_t& clsCap)
		{
			assert(clsCap && clsCap <= UINT32_MAX);
			nscatters = clsCap;
			size_t scatterBytes = nscatters * hc_srsize;
			size_t newCap = scatterBytes + nscatters;
			assert(newCap);
			if (pinned_cnf == NULL) 
				CHECK(cudaHostAlloc((void**)&pinned_cnf, hc_cnfsize, cudaHostAllocDefault));
			if (auxPool.cap < newCap) {
				if (auxPool.cap) {
					assert(auxPool.mem != NULL);
					assert(d_scatter != NULL);
					CHECK(cudaFree(auxPool.mem));
					auxPool.mem = NULL;
					dcap -= auxPool.cap;
					assert(dcap >= 0);
					_free += auxPool.cap;
				}
				assert(auxPool.mem == NULL);
				dcap += newCap;
				if (!hasDeviceMem("Auxiliary ")) { dcap -= newCap; return false; }
				CHECK(cudaMalloc((void**)&auxPool.mem, newCap));
				d_scatter = (S_REF*)auxPool.mem;
				d_stencil = auxPool.mem + scatterBytes;
				auxPool.cap = newCap;
			}
			return true;
		}

		bool cuMM::resizeCNF(CNF*& cnf, const size_t& clsCap, const size_t& litsCap) {
			assert(clsCap && clsCap <= UINT32_MAX);
			assert(litsCap && litsCap <= UINT32_MAX);
			assert(litsCap >= clsCap);
			size_t csBytes = clsCap * hc_srsize;
			size_t dataBytes = clsCap * hc_scsize + (litsCap - clsCap) * hc_bucket;
			assert(dataBytes % hc_bucket == 0);
			size_t newCap = hc_cnfsize + dataBytes + csBytes;
			assert(newCap);
			if (cnfPool.cap == 0) {
				assert(cnf == NULL);
				assert(cnfPool.mem == NULL);
				cnfPool.cap = newCap, cap += cnfPool.cap;
				if (!hasUnifiedMem("CNF")) { cap -= newCap; return false; }
				CHECK(cudaMallocManaged((void**)&cnfPool.mem, cnfPool.cap));
				if (devProp.major > 5) {
					PFLOGN2(2, " Advising GPU driver to favor global over system memory..");
					CHECK(cudaMemAdvise(cnfPool.mem, cnfPool.cap, cudaMemAdviseSetPreferredLocation, MASTER_GPU));
					PFLDONE(2, 5);
				}
				cnf = (CNF*)cnfPool.mem;
				S_REF data_cap = S_REF(dataBytes / hc_bucket);
				new (cnf) CNF(data_cap, uint32(clsCap));
				d_cnf_mem = cnf->data().mem, d_cs_mem = cnf->csData();
			}
			else {
				assert(cnf != NULL);
				assert(cnfPool.mem != NULL);
				cap -= cnfPool.cap;
				assert(cap >= 0);
				_free += cnfPool.cap;
				cap += newCap;
				if (!hasUnifiedMem("CNF")) { cap -= newCap; return false; }
				cacheCNFPtr(cnf);
				addr_t newMem = NULL;
				CHECK(cudaMallocManaged((void**)&newMem, newCap));
				sync();
				if (devProp.major > 5) {
					PFLOGN2(2, " Advising GPU driver to favor global over system memory..");
					CHECK(cudaMemAdvise(newMem, newCap, cudaMemAdviseSetPreferredLocation, MASTER_GPU));
					CHECK(cudaMemPrefetchAsync(newMem, newCap, MASTER_GPU));
					PFLDONE(2, 5);
				}
				CNF* tmp_cnf = (CNF*)newMem;
				S_REF data_cap = S_REF(dataBytes / hc_bucket);
				new (tmp_cnf) CNF(data_cap, uint32(clsCap));
				d_cnf_mem = tmp_cnf->data().mem, d_cs_mem = tmp_cnf->csData();
				if (profile_gpu) cutimer->start();
				if (gc_par) compactCNF(cnf, tmp_cnf);
				else {
					sync(), tmp_cnf->copyFrom(cnf);
					pinned_cnf->resize(tmp_cnf->data().size, tmp_cnf->size());
				}
				if (profile_gpu) cutimer->stop(), cutimer->gc += cutimer->gpuTime();
				CHECK(cudaFree(cnfPool.mem));
				cnfPool.mem = newMem, cnfPool.cap = newCap, cnf = tmp_cnf;
			}
			return true;
		}

		bool cuMM::resizeOTAsync(OT*& ot, const size_t& litsCap, const cudaStream_t& _s) {
			assert(d_hist != NULL);
			assert(d_segs != NULL);
			assert(litsCap && litsCap <= UINT32_MAX);
			size_t tb = 0;
			DeviceScan::ExclusiveSum(NULL, tb, d_hist, d_segs, inf.nDualVars, _s), assert(tb <= litsbytes);
			DeviceScan::ExclusiveSum(d_lits, tb, d_hist, d_segs, inf.nDualVars, _s);
			if (!otBlocks) otBlocks = std::min((inf.nDualVars + BLOCK1D - 1) / BLOCK1D, maxGPUTereads / BLOCK1D);
			size_t newCap = hc_otsize + inf.nDualVars * hc_olsize + litsCap * hc_srsize;
			assert(newCap);
			if (otPool.cap < newCap) { // realloc
				if (otPool.cap) {
					assert(otPool.mem != NULL);
					assert(ot != NULL);
					CHECK(cudaFree(otPool.mem));
					otPool.mem = NULL;
					cap -= otPool.cap;
					assert(cap >= 0);
					_free += otPool.cap;
				}
				assert(otPool.mem == NULL);
				cap += newCap;
				if (!hasUnifiedMem("OT")) { cap -= newCap; return false; }
				CHECK(cudaMallocManaged((void**)&otPool.mem, newCap));
				if (devProp.major > 5) {
					PFLOGN2(2, " Advising GPU driver to favor global over system memory..");
					CHECK(cudaMemAdvise(otPool.mem, newCap, cudaMemAdviseSetPreferredLocation, MASTER_GPU));
					CHECK(cudaMemPrefetchAsync(otPool.mem, newCap, MASTER_GPU, _s));
					PFLDONE(2, 5);
				}
				ot = (OT*)otPool.mem;
				LOGERR("Summing histogram failed");
				sync(_s); // needed for calling the next constructor on host
				new (ot) OT(inf.nDualVars);
				assignListPtrs << <otBlocks, BLOCK1D, 0, _s >> > (ot, d_hist, d_segs, inf.nDualVars);
				otPool.cap = newCap;
			}
			else 
				assignListPtrs << <otBlocks, BLOCK1D, 0, _s >> > (ot, d_hist, d_segs, inf.nDualVars);
			if (sync_always) {
				LOGERR("Occurrence lists allocation failed");
				sync(_s);
			}
			return true;
		}

		void cuMM::createMirror(CNF*& hcnf, const size_t& clsCap, const size_t& litsCap)
		{
			assert(clsCap && clsCap <= UINT32_MAX);
			assert(litsCap && litsCap <= UINT32_MAX);
			assert(litsCap >= clsCap);
			size_t csBytes = clsCap * hc_srsize;
			size_t dataBytes = clsCap * hc_scsize + (litsCap - clsCap) * hc_bucket;
			assert(dataBytes % hc_bucket == 0);
			size_t newCap = hc_cnfsize + dataBytes + csBytes;
			assert(newCap);
			if (hcnfPool.cap < newCap) {
				hcnfPool.cap = newCap;
				pfalloc(hcnfPool.mem, hcnfPool.cap);
				hcnf = (CNF*)hcnfPool.mem;
			}
			S_REF data_cap = S_REF(dataBytes / hc_bucket);
			new (hcnf) CNF(data_cap, uint32(clsCap));
		}

		void cuMM::mirrorCNF(CNF*& hcnf)
		{
			assert(cnfPool.cap);
			assert(cnfPool.mem != NULL);
			CHECK(cudaMemcpy(hcnf, cnfPool.mem, hc_cnfsize, cudaMemcpyDeviceToHost));
			size_t csBytes = hcnf->size() * hc_srsize;
			size_t dataBytes = hcnf->data().size * hc_bucket;
			size_t newCap = hc_cnfsize + dataBytes + csBytes;
			assert(newCap <= cnfPool.cap);
			if (hcnfPool.cap < newCap) {
				hcnfPool.cap = newCap;
				pfalloc(hcnfPool.mem, hcnfPool.cap);
				hcnf = (CNF*)hcnfPool.mem;
			}
			hcnf->fixPointer(); // replace device with host pointers
		}

		void cuMM::resizeCNFAsync(CNF* dcnf, const S_REF& data_size, const uint32& cs_size)
		{
			assert(dcnf != NULL);
			assert(data_size);
			assert(cs_size);
			resizeCNF_k << <1, 1 >> > (dcnf, data_size, cs_size);
			if (sync_always) {
				LOGERR("Resizing CNF failed");
				sync();
			}
		}

		void cuMM::freeVars() {
			if (varsPool.mem != NULL) {
				d_units = NULL;
				CHECK(cudaFree(varsPool.mem)), varsPool.mem = NULL;
				cap -= varsPool.cap;
				assert(cap >= 0);
				_free += varsPool.cap;
				varsPool.cap = 0;
			}
		}

		void cuMM::freeCNF() {
			if (cnfPool.mem != NULL) {
				d_cnf_mem = NULL, d_cs_mem = NULL;
				CHECK(cudaFree(cnfPool.mem)), cnfPool.mem = NULL;
				cap -= cnfPool.cap;
				assert(cap >= 0);
				_free += cnfPool.cap;
				cnfPool.cap = 0;
			}
		}

		void cuMM::freeOT() {
			if (otPool.mem != NULL) {
				CHECK(cudaFree(otPool.mem)), otPool.mem = NULL;
				cap -= otPool.cap;
				assert(cap >= 0);
				_free += otPool.cap;
				otPool.cap = 0;
			}
		}

		void cuMM::freeFixed() {
			if (auxPool.mem != NULL) {
				CHECK(cudaFree(auxPool.mem)), auxPool.mem = NULL;
				_free += auxPool.cap, auxPool.cap = 0;
			}
			if (histPool.mem != NULL) {
				CHECK(cudaFree(histPool.mem)), histPool.mem = NULL;
				_free += histPool.cap, histPool.cap = 0;
			}
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
					PFLOGE("cannot allocate new memory block for Thrust.");
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
			int64 new_cap = allocBlock->second;
			allocBlocks.erase(allocBlock);
			freeBlocks.insert(std::make_pair(new_cap, ptr)); // cache free block
		}
	}
}