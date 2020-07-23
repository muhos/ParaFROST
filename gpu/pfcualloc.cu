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
#include <thrust/scan.h>
#include <cub/device/device_scan.cuh>
#include <cub/device/device_select.cuh>
using namespace cub;

namespace pFROST {

	namespace SIGmA {
		//=============================//
		//	  CUDA memory management   //
		//=============================//
		__global__ void resizeCNF_k(CNF* cnf, CNF* hcnf) { cnf->resize(hcnf->data().size, hcnf->size()); }
		__global__ void resizeCNF_k(CNF* cnf, uint32 d_size, uint32 cs_size) { cnf->resize(d_size, cs_size); }
		__global__ void scatter_k(CNF* src, uint32* scatter, addr_t stencil) {
			uint32 tid = global_tx();
			while (tid < src->size()) {
				SCLAUSE& c = src->clause(tid);
				if (c.deleted()) stencil[tid] = 0, scatter[tid] = 0;
				else stencil[tid] = 1, scatter[tid] = c.blockSize();
				tid += stride_x();
			}
		}
		__global__ void compact_k(CNF* src, CNF* dest, uint32* scatter, addr_t stencil) {
			uint32 tid = global_tx();
			while (tid < src->size()) {
				if (stencil[tid]) {
					uint32 new_r = scatter[tid];
					new (dest->cref(new_r)) SCLAUSE(src->clause(tid));
					assert(src->clause(tid).size() == dest->cref(new_r)->size());
				}
				tid += stride_x();
			}
		}
		__global__ void assignListPtrs_k0(uint32* refs, uint32* hist, uint32 size)
		{
			uint32 v = global_tx();
			while (v < size) {
				uint32 p = V2D(v + 1), n = NEG(p);
				refs[v] = hist[p] + hist[n];
				v += stride_x();
			}
		}
		__global__ void assignListPtrs_k1(OT* ot, uint32* hist, uint32* refs, uint32 size)
		{
			uint32 v = global_tx();
			while (v < size) {
				uint32 p = V2D(v + 1), n = NEG(p);
				uint32 ps = hist[p], ns = hist[n];
				uint32* head = ot->data(refs[v]);
				(*ot)[p].alloc(head, ps);
				(*ot)[n].alloc(head + ps, ns);
				v += stride_x();
			}
		}
		
		void cuMM::compactCNF(CNF* src, CNF* dest) {
			uint32 old_size = pinned_cnf->size();
			assert(old_size <= nscatters);
			PFLOGN2(2, " Compacting simplified CNF (%d to %d) on GPU..", old_size, inf.nClauses);
			S_REF data_size = inf.nClauses * sizeof(S_REF) + (inf.nLiterals - inf.nClauses);
			resizeCNF_k << <1, 1 >> > (dest, data_size, inf.nClauses);
			size_t tb = 0, ftb = 0;
			uint32* ts = d_lits;
			uint32 nThreads = BLOCK1D, maxBlocks = maxGPUThreads / nThreads;
			uint32 nBlocks = std::min((old_size + nThreads - 1) / nThreads, maxBlocks);
			scatter_k << <nBlocks, nThreads >> > (src, d_scatter, d_stencil);
			DeviceScan::ExclusiveSum(NULL, tb, d_scatter, d_scatter, old_size), assert(tb <= litsbytes);
			DeviceScan::ExclusiveSum(ts, tb, d_scatter, d_scatter, old_size);
			DeviceSelect::Flagged(NULL, ftb, d_scatter, d_stencil, d_cs_mem, ts, old_size), assert(ftb <= litsbytes);
			DeviceSelect::Flagged(ts + 1, ftb, d_scatter, d_stencil, d_cs_mem, ts, old_size);
			compact_k << <nBlocks, nThreads >> > (src, dest, d_scatter, d_stencil);
			pinned_cnf->resize(data_size, inf.nClauses);
			LOGERR(" CNF compact failed");
			sync();
			PFLDONE(2, 5);
		}

		bool cuMM::allocVars(VARS*& vars, const uint32& resCap) {
			assert(vars == NULL);
			vars = new VARS();
			size_t gsts_sz = sizeof(GSTATS);
			size_t cuVec_sz = sizeof(cuVecU);
			size_t uintVec_sz = inf.maxVar * sizeof(uint32);
			size_t varsize = inf.maxVar + 1;
			size_t scores_sz = varsize * sizeof(uint32);
			size_t resolved_sz = resCap * sizeof(uint32);
			size_t new_cap = gsts_sz + scores_sz + resolved_sz + (cuVec_sz + uintVec_sz) * 3;
			varsPool.cap = new_cap, cap += new_cap;
			assert(varsPool.cap);
			if (!hasUnifiedMem("Fixed")) return false;
			CHECK(cudaMallocManaged((void**)&varsPool.mem, varsPool.cap));
			addr_t ea = varsPool.mem, end = ea + varsPool.cap;
			vars->gstats = (GSTATS*)ea, ea += gsts_sz;
			vars->pVars = (cuVecU*)ea, ea += cuVec_sz;
			vars->units = (cuVecU*)ea, ea += cuVec_sz;
			vars->resolved = (cuVecU*)ea, ea += cuVec_sz;
			uint32* uintPtr = (uint32*)ea;
			vars->pVars->alloc(uintPtr, inf.maxVar), uintPtr += inf.maxVar; d_units = uintPtr;
			vars->units->alloc(uintPtr, inf.maxVar), uintPtr += inf.maxVar;
			vars->eligible = uintPtr, uintPtr += inf.maxVar;
			vars->scores = uintPtr, uintPtr += varsize;
			vars->resolved->alloc(uintPtr, resCap), uintPtr += resCap;
			assert((addr_t)uintPtr == end);
			if (devProp.major > 5) {
				PFLOGN2(2, " Advising GPU driver to favor global over system memory..");
				CHECK(cudaMemAdvise(vars->gstats, gsts_sz, cudaMemAdviseSetPreferredLocation, MASTER_GPU));
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

		bool cuMM::allocHist(cuHist& cuhist, const uint32& litsCap)
		{
			assert(litsCap);
			assert(inf.nDualVars == v2l(inf.maxVar + 1ULL));
			cuhist.litsbytes = litsbytes = litsCap * sizeof(uint32);
			size_t hsumBytes = inf.maxVar * sizeof(uint32);
			size_t histBytes = inf.nDualVars * sizeof(uint32);
			size_t newCap = hsumBytes + histBytes + litsbytes;
			assert(newCap);
			if (histPool.cap < newCap) {
				if (histPool.cap) {
					assert(histPool.mem != NULL);
					assert(d_hist != NULL);
					assert(d_hsum != NULL);
					assert(cuhist.h_hist != NULL);
					assert(cuhist.d_hist != NULL);
					assert(cuhist.d_hsum != NULL);
					assert(cuhist.d_lits != NULL);
					CHECK(cudaFree(histPool.mem));
					CHECK(cudaFreeHost(hhistPool.mem));
					hhistPool.mem = NULL;
					histPool.mem = NULL;
				}
				assert(histPool.mem == NULL);
				dcap -= histPool.cap, _free += histPool.cap, dcap += newCap;
				if (!hasDeviceMem("Histogram")) return false;
				CHECK(cudaHostAlloc((void**)&hhistPool.mem, histBytes, cudaHostAllocDefault));
				CHECK(cudaMalloc((void**)&histPool.mem, newCap));
				// NOTE: d_hsum, d_hist used internally by OT allocation and externally by BVE for different purpose (memory reuse)
				//		 d_lits is used as temporary storage as well for CUB routines
				cuhist.h_hist = hhistPool.mem;
				cuhist.d_hsum = d_hsum = histPool.mem; 
				cuhist.d_hist = d_hist = d_hsum + inf.maxVar;
				cuhist.d_lits = d_lits = d_hist + inf.nDualVars;
				cuhist.thrust_hist = t_iptr(d_hist);
				cuhist.thrust_lits = t_iptr(cuhist.d_lits);
				histPool.cap = newCap;
				hhistPool.cap = histBytes;
			}
			return true;
		}

		bool cuMM::allocAux(const uint32& clsCap)
		{
			assert(clsCap);
			nscatters = clsCap;
			uint32 maxSize = nscatters;
			size_t newCap = maxSize * sizeof(uint32) + maxSize;
			assert(newCap);
			if (pinned_cnf == NULL) 
				CHECK(cudaHostAlloc((void**)&pinned_cnf, sizeof(CNF), cudaHostAllocDefault));
			if (auxPool.cap < newCap) {
				if (auxPool.cap) {
					assert(auxPool.mem != NULL);
					assert(d_scatter != NULL);
					CHECK(cudaFree(auxPool.mem));
					auxPool.mem = NULL;
				}
				assert(auxPool.mem == NULL);
				dcap -= auxPool.cap, _free += auxPool.cap, dcap += newCap;
				if (!hasDeviceMem("Auxiliary ")) return false;
				CHECK(cudaMalloc((void**)&auxPool.mem, newCap));
				d_scatter = auxPool.mem;
				d_stencil = addr_t(d_scatter + maxSize);
				auxPool.cap = newCap;
			}
			return true;
		}

		bool cuMM::resizeCNF(CNF*& cnf, const uint32& clsCap, const uint32& litsCap) {
			assert(clsCap);
			assert(litsCap);
			assert(litsCap >= clsCap);
			size_t bucket = sizeof(S_REF);
			size_t cnfBytes = sizeof(CNF);
			size_t csBytes = clsCap * bucket;
			size_t dataBytes = clsCap * sizeof(SCLAUSE) + (litsCap - clsCap) * sizeof(uint32);
			assert(dataBytes % bucket == 0);
			size_t newCap = cnfBytes + dataBytes + csBytes;
			assert(newCap);
			if (cnfPool.cap == 0) {
				assert(cnf == NULL);
				assert(cnfPool.mem == NULL);
				assert(sizeof(SCLAUSE) / sizeof(S_REF) == sizeof(uint32));
				cnfPool.cap = newCap, cap += cnfPool.cap;
				if (!hasUnifiedMem("CNF")) return false;
				CHECK(cudaMallocManaged((void**)&cnfPool.mem, cnfPool.cap));
				if (devProp.major > 5) {
					PFLOGN2(2, " Advising GPU driver to favor global over system memory..");
					CHECK(cudaMemAdvise(cnfPool.mem, cnfPool.cap, cudaMemAdviseSetPreferredLocation, MASTER_GPU));
					PFLDONE(2, 5);
				}
				cnf = (CNF*)cnfPool.mem;
				S_REF data_cap = S_REF(dataBytes / bucket);
				new (cnf) CNF(data_cap, clsCap);
				d_cnf_mem = cnf->data().mem, d_cs_mem = cnf->csData();
			}
			else {
				assert(cnf != NULL);
				assert(cnfPool.mem != NULL);
				cap -= cnfPool.cap, _free += cnfPool.cap, cap += newCap;
				if (!hasUnifiedMem(__func__)) return false;
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
				S_REF data_cap = S_REF(dataBytes / bucket);
				new (tmp_cnf) CNF(data_cap, clsCap);
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

		bool cuMM::resizeOTAsync(OT*& ot, const uint32& litsCap, const cudaStream_t& _s) {
			assert(d_hist != NULL);
			assert(d_hsum != NULL);
			assert(litsCap);
			if (!otBlocks) otBlocks = std::min((inf.maxVar + BLOCK1D - 1) / BLOCK1D, maxGPUThreads / BLOCK1D);
			assignListPtrs_k0 << <otBlocks, BLOCK1D, 0, _s >> > (d_hsum, d_hist, inf.maxVar);
			size_t tb = 0;
			DeviceScan::ExclusiveSum(NULL, tb, d_hsum, d_hsum, inf.maxVar, _s), assert(tb <= litsbytes);
			size_t newCap = sizeof(OT) + inf.nDualVars * sizeof(OL) + litsCap * sizeof(uint32);
			assert(newCap);
			if (otPool.cap < newCap) { // realloc
				if (otPool.cap) {
					assert(otPool.mem != NULL);
					assert(ot != NULL);
					CHECK(cudaFree(otPool.mem));
					otPool.mem = NULL;
				}
				assert(otPool.mem == NULL);
				cap -= otPool.cap, _free += otPool.cap, cap += newCap;
				if (!hasUnifiedMem("OT")) return false;
				CHECK(cudaMallocManaged((void**)&otPool.mem, newCap));
				if (devProp.major > 5) {
					PFLOGN2(2, " Advising GPU driver to favor global over system memory..");
					CHECK(cudaMemAdvise(otPool.mem, newCap, cudaMemAdviseSetPreferredLocation, MASTER_GPU));
					CHECK(cudaMemPrefetchAsync(otPool.mem, newCap, MASTER_GPU, _s));
					PFLDONE(2, 5);
				}
				ot = (OT*)otPool.mem;
				sync(_s); // needed for calling the next constructor on host
				new (ot) OT(inf.nDualVars);
				DeviceScan::ExclusiveSum(d_lits, tb, d_hsum, d_hsum, inf.maxVar, _s);
				assignListPtrs_k1 << <otBlocks, BLOCK1D, 0, _s >> > (ot, d_hist, d_hsum, inf.maxVar);
				otPool.cap = newCap;
			}
			else {
				DeviceScan::ExclusiveSum(d_lits, tb, d_hsum, d_hsum, inf.maxVar, _s);
				assignListPtrs_k1 << <otBlocks, BLOCK1D, 0, _s >> > (ot, d_hist, d_hsum, inf.maxVar);
			}
			if (sync_always) sync(_s);
			return true;
		}

		void cuMM::createMirror(CNF*& hcnf, const uint32& clsCap, const uint32& litsCap)
		{
			assert(clsCap);
			assert(litsCap);
			assert(litsCap >= clsCap);
			size_t bucket = sizeof(S_REF);
			size_t cnfBytes = sizeof(CNF);
			size_t csBytes = clsCap * bucket;
			size_t dataBytes = clsCap * sizeof(SCLAUSE) + (litsCap - clsCap) * sizeof(uint32);
			assert(dataBytes % bucket == 0);
			size_t newCap = cnfBytes + dataBytes + csBytes;
			assert(newCap);
			if (hcnfPool.cap < newCap) {
				hcnfPool.cap = newCap;
				pfalloc(hcnfPool.mem, hcnfPool.cap);
				hcnf = (CNF*)hcnfPool.mem;
			}
			S_REF data_cap = S_REF(dataBytes / bucket);
			new (hcnf) CNF(data_cap, clsCap);
		}

		void cuMM::mirrorCNF(CNF*& hcnf)
		{
			assert(cnfPool.cap);
			assert(cnfPool.mem != NULL);
			CHECK(cudaMemcpy(hcnf, cnfPool.mem, sizeof(CNF), cudaMemcpyDeviceToHost));
			size_t bucket = sizeof(S_REF);
			size_t cnfBytes = sizeof(CNF);
			size_t csBytes = hcnf->size() * bucket;
			size_t dataBytes = hcnf->data().size * bucket;
			size_t newCap = cnfBytes + dataBytes + csBytes;
			assert(newCap <= cnfPool.cap);
			if (hcnfPool.cap < newCap) {
				hcnfPool.cap = newCap;
				pfalloc(hcnfPool.mem, hcnfPool.cap);
				hcnf = (CNF*)hcnfPool.mem;
			}
			hcnf->fixPointer(); // replace device with host pointers
		}

		void cuMM::resizeCNFAsync(CNF* dcnf, CNF* hcnf)
		{
			assert(dcnf != NULL);
			assert(hcnf != NULL);
			assert(hcnf->data().size);
			assert(hcnf->size());
			resizeCNF_k << <1, 1 >> > (dcnf, hcnf->data().size, hcnf->size());
			if (sync_always) sync();
		}

		void cuMM::freeDynamic() {
			if (!cap) return;
			d_units = NULL, d_cnf_mem = NULL, d_cs_mem = NULL;
			if (otPool.mem != NULL) CHECK(cudaFree(otPool.mem)), otPool.mem = NULL;
			if (varsPool.mem != NULL) CHECK(cudaFree(varsPool.mem)), varsPool.mem = NULL;
			if (cnfPool.mem != NULL) CHECK(cudaFree(cnfPool.mem)), cnfPool.mem = NULL;
			varsPool.cap = 0, cnfPool.cap = 0, otPool.cap = 0, cap = 0, _free = _tot - dcap;
		}

		void cuMM::freeFixed() {
			if (!dcap) return;
			if (pinned_cnf != NULL) CHECK(cudaFreeHost(pinned_cnf)), pinned_cnf = NULL;
			if (pinned_units != NULL) CHECK(cudaFreeHost(pinned_units)), pinned_units = NULL;
			if (hhistPool.mem != NULL) CHECK(cudaFreeHost(hhistPool.mem)), hhistPool.mem = NULL;
			if (auxPool.mem != NULL) CHECK(cudaFree(auxPool.mem)), auxPool.mem = NULL;
			if (histPool.mem != NULL) CHECK(cudaFree(histPool.mem)), histPool.mem = NULL;
			hhistPool.cap = 0, histPool.cap = 0, auxPool.cap = 0;
			dcap = 0, _free = _tot - cap;
		}

		void cuMM::breakMirror() {
			if (hcnfPool.mem != NULL) std::free(hcnfPool.mem), hcnfPool.mem = NULL;
			hcnfPool.cap = 0;
		}

	}
}