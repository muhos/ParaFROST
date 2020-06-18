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
#include "pfgrid.cuh"
#include "pfatomics.cuh"

namespace pFROST {

	namespace SIGmA {
		//=============================//
		//	   GPU memory management   //
		//=============================//
		_PFROST_D_ uint32* OT::data(const uint32& offset) { return occurs + atomicAdd(&maxEntries, offset); }
		//_PFROST_D_ void CNF::newClause(SCLAUSE& src) {
		//	assert(src.size() > 0);
		//	S_REF newC = cls + atomicAggInc(&n_cls);
		//	newC->set_ptr(lits + atomicAdd(&n_lits, src.size()));
		//	newC->copyFrom(src);
		//}
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
		//__global__ void shrink_k(CNF* dest, CNF* src)
		//{
		//	uint32 i = global_tx();
		//	while (i < src->size()) {
		//		SCLAUSE& s = (*src)[i];
		//		if (s.status() == LEARNT || s.status() == ORIGINAL) dest->newClause(s);
		//		i += stride_x();
		//	}
		//}

		bool cuMM::allocPV(VARS*& vars) {
			assert(varsPool.mem == NULL);
			vars = new VARS();
			assert(vars->pVars == NULL);
			assert(vars->units == NULL);
			assert(vars->scores == NULL);
			assert(vars->gstats == NULL);
			size_t gsts_sz = sizeof(GSTATS);
			size_t cuVec_sz = sizeof(cuVecU);
			size_t uintVec_sz = inf.maxVar * sizeof(uint32);
			size_t varsize = inf.maxVar + 1;
			size_t scores_sz = varsize * sizeof(uint32);
			varsPool.cap = gsts_sz + cuVec_sz * 2 + uintVec_sz * 3 + scores_sz + varsize;
			assert(varsPool.cap);
			cap += varsPool.cap;
			if (!hasFreeMem("Variables")) return false;
			CHECK(cudaMallocManaged((void**)&varsPool.mem, varsPool.cap));
			addr_t ea = varsPool.mem, end = ea + varsPool.cap;
			vars->gstats = (GSTATS*)ea, ea += gsts_sz;
			vars->pVars = (cuVecU*)ea, ea += cuVec_sz;
			vars->units = (cuVecU*)ea, ea += cuVec_sz;
			uint32* uintPtr = (uint32*)ea;
			vars->pVars->alloc(uintPtr, inf.maxVar), uintPtr += inf.maxVar; d_units = uintPtr;
			vars->units->alloc(uintPtr, inf.maxVar), uintPtr += inf.maxVar;
			vars->eligible = uintPtr, uintPtr += inf.maxVar;
			vars->scores = uintPtr, uintPtr += varsize;
			LIT_ST* litPtr = (LIT_ST*)uintPtr;
			vars->vstate = litPtr, litPtr += varsize;
			assert((addr_t)litPtr == end);
			if (devProp.major > 5) {
				PFLOGN2(2, " Advising GPU driver to favor global over system memory..");
				CHECK(cudaMemAdvise(vars->gstats, gsts_sz, cudaMemAdviseSetPreferredLocation, MASTER_GPU));
				addr_t tmpPtr = ea + uintVec_sz; // skip pVars too
				CHECK(cudaMemAdvise(tmpPtr, uintVec_sz * 3 + varsize, cudaMemAdviseSetPreferredLocation, MASTER_GPU));
				CHECK(cudaMemPrefetchAsync(tmpPtr, uintVec_sz * 3 + varsize, MASTER_GPU));
				PFLDONE(2, 5);
			}
			vars->cachedUnits = new uint32[inf.maxVar];
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
				if (!hasFreeMem("CNF")) return false;
				CHECK(cudaMallocManaged((void**)&cnfPool.mem, cnfPool.cap));
				CNF* scnf = (CNF*)cnfPool.mem;
				S_REF data_cap = S_REF(dataBytes / bucket);
				cnf = new (scnf) CNF(data_cap, clsCap);
				d_cnf_mem = cnf->data().mem, d_cs_mem = cnf->csData();
			}
			else {
				assert(cnf != NULL);
				assert(cnfPool.mem != NULL);
				if (cnfPool.cap >= newCap)
					cnf->shrink();
				else {
					cap -= cnfPool.cap, cap += newCap;
					if (!hasFreeMem(__func__)) return false;
					addr_t newMem = NULL;
					CHECK(cudaMallocManaged((void**)&newMem, newCap));
					CNF* tmp_cnf = (CNF*)newMem;
					S_REF data_cap = S_REF(dataBytes / bucket);
					tmp_cnf = new (tmp_cnf) CNF(data_cap, clsCap);
					tmp_cnf->copyFrom(cnf);
					CHECK(cudaFree(cnfPool.mem));
					cnfPool.mem = newMem, cnfPool.cap = newCap, cnf = tmp_cnf;
					d_cnf_mem = cnf->data().mem, d_cs_mem = cnf->csData();
					prefetchCNF();
				}
			}
			return true;
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
				cap -= otPool.cap, cap += newCap;
				if (!hasFreeMem("OT")) return false;
				CHECK(cudaMallocManaged((void**)&otPool.mem, newCap));
				if (devProp.major > 5) {
					PFLOGN2(2, " Advising GPU driver to favor global over system memory..");
					CHECK(cudaMemAdvise(otPool.mem, newCap, cudaMemAdviseSetPreferredLocation, MASTER_GPU));
					CHECK(cudaMemPrefetchAsync(otPool.mem, newCap, MASTER_GPU));
					PFLDONE(2, 5);
				}
				ot = (OT*)otPool.mem;
				ot = new (ot) OT(inf.nDualVars);
				if (!otBlocks) otBlocks = MIN((inf.maxVar + BLOCK1D - 1) / BLOCK1D, maxGPUThreads / BLOCK1D);
				assignListPtrs_k << <otBlocks, BLOCK1D >> > (ot, d_hist, inf.maxVar);
				otPool.cap = newCap;
			}
			else assignListPtrs_k << <otBlocks, BLOCK1D, 0, _s >> > (ot, d_hist, inf.maxVar);
			return true;
		}

		void cuMM::resetOTCapAsync(OT* ot, const cudaStream_t& _s) { 
			if (ot == NULL) return;
			resetOTCap_k << <1, 1, 0, _s >> > (ot);
		}

		void cuMM::destroy() {
			if (!cap) return;
			d_units = NULL, d_cnf_mem = NULL, d_cs_mem = NULL;
			if (otPool.mem != NULL) CHECK(cudaFree(otPool.mem)), otPool.mem = NULL;
			if (varsPool.mem != NULL) CHECK(cudaFree(varsPool.mem)), varsPool.mem = NULL;
			if (cnfPool.mem != NULL) CHECK(cudaFree(cnfPool.mem)), cnfPool.mem = NULL;
			cap = 0;
		}

		void cuMM::breakMirror() {
			if (hcnfPool.mem != NULL) std::free(hcnfPool.mem), hcnfPool.mem = NULL;
		}
	}
}