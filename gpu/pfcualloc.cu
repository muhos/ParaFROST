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

namespace pFROST {

	namespace SIGmA {
		//=============================//
		//	   GPU memory management   //
		//=============================//
		__global__ void resizeCNF_k(CNF* cnf, uint32 data_size, uint32 cs_size) { cnf->resize(data_size, cs_size); }
		__global__ void resetOTCap_k(OT* ot) { ot->resetCap(); }
		__global__ void assignListPtrs_k(OT* ot, uint32* hist, uint32 size)
		{
			uint32 v = global_tx();
			while (v < size) {
				uint32 p = V2D(v + 1), n = NEG(p);
				uint32 ps = hist[p], ns = hist[n];
				uint32* head = ot->data(ps + ns);
				(*ot)[p].alloc(head, ps);
				(*ot)[n].alloc(head + ps, ns);
				v += stride_x();
			}
		}

		bool cuMM::allocFixed(VARS*& vars, const uint32& resCap) {
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
			if (!hasFreeMem("Fixed")) return false;
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
			vars->cachedUnits = new uint32[inf.maxVar];
			return true;
		}

		void cuMM::resizeHostCNF(CNF*& cnf, const uint32& clsCap, const uint32& litsCap)
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
				cnf = (CNF*)hcnfPool.mem;
				S_REF data_cap = S_REF(dataBytes / bucket);
				cnf = new (cnf) CNF(data_cap, clsCap);
			}
		}

		void cuMM::mirrorCNF(CNF*& cnf)
		{
			assert(cnfPool.cap);
			assert(cnfPool.mem != NULL);
			hcnfPool.cap = cnfPool.cap;
			pfalloc(hcnfPool.mem, hcnfPool.cap);
			cnf = (CNF*)hcnfPool.mem;
			cnf = new (cnf) CNF();
			CHECK(cudaMemcpy(cnf, cnfPool.mem, sizeof(CNF), cudaMemcpyDeviceToHost));
			cnf->fixPointer();
		}

		void cuMM::resizeCNFAsync(CNF* dcnf, CNF* hcnf) 
		{
			assert(dcnf != NULL);
			assert(hcnf != NULL);
			assert(hcnf->data().size);
			assert(hcnf->size());
			resizeCNF_k << <1, 1 >> > (dcnf, hcnf->data().size, hcnf->size());
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
				if (!hasFreeMem("CNF")) return false;
				CHECK(cudaMallocManaged((void**)&cnfPool.mem, cnfPool.cap));
				if (devProp.major > 5) {
					PFLOGN2(2, " Advising GPU driver to favor global over system memory..");
					CHECK(cudaMemAdvise(cnfPool.mem, cnfPool.cap, cudaMemAdviseSetPreferredLocation, MASTER_GPU));
					PFLDONE(2, 5);
				}
				cnf = (CNF*)cnfPool.mem;
				S_REF data_cap = S_REF(dataBytes / bucket);
				cnf = new (cnf) CNF(data_cap, clsCap);
				d_cnf_mem = cnf->data().mem, d_cs_mem = cnf->csData();
			}
			else {
				assert(cnf != NULL);
				assert(cnfPool.mem != NULL);
				cap -= cnfPool.cap, _free += cnfPool.cap, cap += newCap;
				if (!hasFreeMem(__func__)) return false;
				addr_t newMem = NULL;
				CHECK(cudaMallocManaged((void**)&newMem, newCap));
				CNF* tmp_cnf = (CNF*)newMem;
				S_REF data_cap = S_REF(dataBytes / bucket);
				tmp_cnf = new (tmp_cnf) CNF(data_cap, clsCap);
				d_cnf_mem = tmp_cnf->data().mem, d_cs_mem = tmp_cnf->csData();
				tmp_cnf->copyFrom(cnf);
				CHECK(cudaFree(cnfPool.mem));
				cnfPool.mem = newMem, cnfPool.cap = newCap, cnf = tmp_cnf;
				prefetchCNF();
			}
			return true;
		}

		bool cuMM::resizeOTAsync(OT*& ot, uint32* d_hist, const uint32& litsCap, const cudaStream_t& _s) {
			assert(d_hist != NULL);
			assert(litsCap);
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
				if (!hasFreeMem("OT")) return false;
				CHECK(cudaMallocManaged((void**)&otPool.mem, newCap));
				if (devProp.major > 5) {
					PFLOGN2(2, " Advising GPU driver to favor global over system memory..");
					CHECK(cudaMemAdvise(otPool.mem, newCap, cudaMemAdviseSetPreferredLocation, MASTER_GPU));
					CHECK(cudaMemPrefetchAsync(otPool.mem, newCap, MASTER_GPU, _s));
					PFLDONE(2, 5);
				}
				ot = (OT*)otPool.mem;
				ot = new (ot) OT(inf.nDualVars);
				if (!otBlocks) otBlocks = MIN((inf.maxVar + BLOCK1D - 1) / BLOCK1D, maxGPUThreads / BLOCK1D);
				assignListPtrs_k << <otBlocks, BLOCK1D, 0, _s >> > (ot, d_hist, inf.maxVar);
				otPool.cap = newCap;
			}
			else assignListPtrs_k << <otBlocks, BLOCK1D, 0, _s >> > (ot, d_hist, inf.maxVar);
			return true;
		}

		void cuMM::resetOTCapAsync(OT* ot, const cudaStream_t& _s) { 
			assert(ot != NULL);
			resetOTCap_k << <1, 1, 0, _s >> > (ot);
		}

		void cuMM::destroy() {
			if (!cap) return;
			d_units = NULL, d_cnf_mem = NULL, d_cs_mem = NULL;
			if (otPool.mem != NULL) CHECK(cudaFree(otPool.mem)), otPool.mem = NULL;
			if (varsPool.mem != NULL) CHECK(cudaFree(varsPool.mem)), varsPool.mem = NULL;
			if (cnfPool.mem != NULL) CHECK(cudaFree(cnfPool.mem)), cnfPool.mem = NULL;
			breakMirror();
			varsPool.cap = 0, cnfPool.cap = 0, otPool.cap = 0, cap = 0, _free = _tot;
		}

		void cuMM::breakMirror() {
			if (hcnfPool.mem != NULL) std::free(hcnfPool.mem), hcnfPool.mem = NULL;
			hcnfPool.cap = 0;
		}
	}
}