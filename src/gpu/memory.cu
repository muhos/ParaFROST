/***********************************************************************[memory.cu]
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

#include "primitives.cuh"
#include "memory.cuh"
#include "options.cuh"
#include <cub/device/device_scan.cuh>
#include <cub/device/device_select.cuh>

using namespace cub;
using namespace ParaFROST;

//=============================//
//    CUDA memory management   //
//=============================//
const size_t hc_srsize = sizeof(S_REF);
const size_t hc_scsize = sizeof(SCLAUSE);
const size_t hc_otsize = sizeof(OT);
const size_t hc_olsize = sizeof(OL);
const size_t hc_cnfsize = sizeof(CNF);
const size_t hc_varsize = sizeof(uint32);
const size_t hc_cuvecsize = sizeof(cuVecU);

__global__ void resizeCNF_k(CNF* cnf, const S_REF d_size, const uint32 cs_size)
{
	cnf->resize(d_size, cs_size);
}

__global__ void assignListPtrs(OT* __restrict__ ot, const uint32* __restrict__ hist, const S_REF* __restrict__ segs, const uint32 size)
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
	if (gopts.sync_always) {
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

bool cuMM::allocHist(cuHist& cuhist, const bool& proofEnabled)
{
	assert(inf.nDualVars == V2L(inf.maxVar + 1ULL));
	const size_t segBytes = inf.nDualVars * hc_srsize;
	const size_t histBytes = inf.nDualVars * hc_varsize;
	const size_t varsBytes = (inf.maxVar + 1) * hc_varsize;
	size_t min_cap = segBytes + histBytes + varsBytes;
	if (proofEnabled) 
		min_cap += inf.nDualVars;
	assert(min_cap);
	if (histPool.cap < min_cap) {
		DFREE(histPool);
		assert(histPool.mem == NULL);
		if (!hasDeviceMem(min_cap, "Histogram")) return false;
		CHECK(cudaMalloc((void**)&histPool.mem, min_cap));
		// NOTE: d_segs, d_hist used internally by OT allocation and externally
		//       by BVE for calculating resolvents offsets (memory reuse)
		//		 lbyte is used for proof byte counting
		addr_t ea = histPool.mem;
		cuhist.d_segs = d_segs = (S_REF*)ea, ea += segBytes;
		cuhist.d_hist = d_hist = (uint32*)ea, ea += histBytes;
		cuhist.d_vorg = (uint32*)ea, ea += varsBytes;
		if (proofEnabled) {
			cuhist.d_lbyte = ea;
			ea += inf.nDualVars;
		}
		assert(ea == histPool.mem + min_cap);
		cuhist.thrust_hist = t_iptr(d_hist);
		histPool.cap = min_cap;
	}
	return true;
}

bool cuMM::allocVars(VARS*& vars, const size_t& resolvedCap)
{
	assert(varsPool.cap == 0);
	assert(vars == NULL);
	assert(resolvedCap && resolvedCap < UINT32_MAX);
	vars = new VARS();
	const size_t uintVec_sz = inf.maxVar * hc_varsize;
	const size_t varsize = inf.maxVar + 1;
	const size_t scores_sz = varsize * hc_varsize;
	const size_t resolved_sz = resolvedCap * hc_varsize;
	size_t min_cap = hc_cuvecsize * 3;                             // headers: (pVars + units + resolved) 
	min_cap += uintVec_sz * 3 + scores_sz + resolved_sz + varsize; // data:    (pVars + units + eligible) + scores + resolved + eliminated
	assert(min_cap);
	if (!hasUnifiedMem(min_cap, "Fixed")) return false;
	CHECK(cudaMallocManaged((void**)&varsPool.mem, min_cap));
	addr_t ea = varsPool.mem, end = ea + min_cap;
	vars->pVars = (cuVecU*)ea, ea += hc_cuvecsize;
	vars->units = (cuVecU*)ea, ea += hc_cuvecsize;
	vars->resolved = (cuVecU*)ea, ea += hc_cuvecsize;
	uint32* uintPtr = (uint32*)ea;
	vars->pVarsData = uintPtr;
	vars->pVarsSize = (uint32*)(vars->pVars + sizeof(uint32*));
	vars->pVars->alloc(uintPtr, inf.maxVar), uintPtr += inf.maxVar, d_units = uintPtr;
	vars->units->alloc(uintPtr, inf.maxVar), uintPtr += inf.maxVar;
	vars->eligible = uintPtr, uintPtr += inf.maxVar;
	vars->scores = uintPtr, uintPtr += varsize;
	vars->resolved->alloc(uintPtr, uint32(resolvedCap)), uintPtr += resolvedCap;
	Byte* bytePtr = (Byte*)uintPtr;
	vars->eliminated = bytePtr, bytePtr += varsize;
	assert(bytePtr == end);
	varsPool.cap = min_cap;
	#if !defined(_WIN32)
	if (devProp.major > 5) {
		PFLOGN2(2, " Advising GPU driver to favor global over system memory..");
		addr_t tmpPtr = ea + uintVec_sz; // skip pVars
		CHECK(cudaMemAdvise(tmpPtr, end - tmpPtr, cudaMemAdviseSetPreferredLocation, MASTER_GPU));
		CHECK(cudaMemPrefetchAsync(tmpPtr, end - tmpPtr, MASTER_GPU));
		PFLDONE(2, 5);
	}
	#endif
	return true;
}

bool cuMM::allocPinned(VARS* vars, cuHist& cuhist)
{
	assert(vars);
	assert(inf.nDualVars == V2L(inf.maxVar + 1ULL));
	const size_t elimBytes = inf.maxVar + 1;
	const size_t unitBytes = inf.maxVar * hc_varsize;
	const size_t histBytes = inf.nDualVars * hc_varsize;
	size_t min_cap = hc_cnfsize + elimBytes + unitBytes + histBytes;
	assert(min_cap);
	if (pinnedPool.cap) {
		assert(pinnedPool.mem);
		CHECK(cudaFreeHost(pinnedPool.mem));
		pinnedPool.mem = NULL;
		pinnedPool.cap = 0;
	}
	assert(pinnedPool.mem == NULL);
	cudaError_t retVal = cudaHostAlloc((void**)&pinnedPool.mem, min_cap, cudaHostAllocDefault);
	if (retVal != cudaSuccess || retVal == cudaErrorMemoryAllocation) {
		PFLOGW("CUDA runtime failure due to %s", cudaGetErrorString(retVal));
		return false;
	}
	addr_t ea = pinnedPool.mem;
	pinned_cnf = (CNF*)ea, ea += hc_cnfsize;
	cuhist.h_hist = (uint32*)ea, ea += histBytes;
	vars->cachedUnits = (uint32*)ea, ea += unitBytes;
	vars->cachedEliminated = ea, ea += elimBytes;
	assert(ea == (pinnedPool.mem + min_cap));
	pinnedPool.cap = min_cap;
	return true;
}

bool cuMM::allocAux(const size_t& clsCap)
{
	assert(clsCap && clsCap <= UINT32_MAX);
	const size_t scatterBytes = clsCap * hc_srsize;
	const size_t min_cap = scatterBytes + clsCap;
	assert(min_cap);
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
		#if !defined(_WIN32)
		if (devProp.major > 5) {
			PFLOGN2(2, " Advising GPU driver to favor global over system memory..");
			CHECK(cudaMemAdvise(cnfPool.mem, min_cap, cudaMemAdviseSetPreferredLocation, MASTER_GPU));
			PFLDONE(2, 5);
		}
		#endif
		cnf = (CNF*)cnfPool.mem;
		const S_REF data_cap = S_REF(dataBytes / hc_bucket);
		new (cnf) CNF(data_cap, uint32(clsCap));
		d_cnf_mem = cnf->data().mem, d_refs_mem = cnf->refsData();
		cnfPool.cap = min_cap;
	}
	else {
		assert(cnf);
		assert(cnfPool.mem);
		if (!hasUnifiedMem(min_cap, "CNF")) return false;
		cacheCNFPtr(cnf);
		addr_t newMem = NULL;
		CHECK(cudaMallocManaged((void**)&newMem, min_cap));
		sync();
		#if !defined(_WIN32)
		if (devProp.major > 5) {
			PFLOGN2(2, " Advising GPU driver to favor global over system memory..");
			CHECK(cudaMemAdvise(newMem, min_cap, cudaMemAdviseSetPreferredLocation, MASTER_GPU));
			CHECK(cudaMemPrefetchAsync(newMem, min_cap, MASTER_GPU));
			PFLDONE(2, 5);
		}
		#endif
		CNF* tmp_cnf = (CNF*)newMem;
		const S_REF data_cap = S_REF(dataBytes / hc_bucket);
		new (tmp_cnf) CNF(data_cap, uint32(clsCap));
		d_cnf_mem = tmp_cnf->data().mem, d_refs_mem = tmp_cnf->refsData();
		if (gopts.profile_gpu) cutimer->start();
		if (gopts.gc_gpu) compactCNF(cnf, tmp_cnf);
		else {
			sync(), tmp_cnf->copyFrom(cnf);
			pinned_cnf->resize(tmp_cnf->data().size, tmp_cnf->size());
		}
		if (gopts.profile_gpu) cutimer->stop(), cutimer->gc += cutimer->gpuTime();
		FREE(cnfPool);
		cnfPool.mem = newMem;
		cnfPool.cap = min_cap;
		cnf = tmp_cnf;
	}
	return true;
}

bool cuMM::resizeOTAsync(OT*& ot, const size_t& min_lits, const cudaStream_t& _s)
{
	assert(d_hist);
	assert(d_segs);
	assert(min_lits && min_lits <= UINT32_MAX);
	uint32* tmp = resizeLits(min_lits);
	if (!tmp) return false;
	size_t ebytes = 0;
	DeviceScan::ExclusiveSum(NULL, ebytes, d_hist, d_segs, inf.nDualVars, _s);
	DeviceScan::ExclusiveSum(tmp, ebytes, d_hist, d_segs, inf.nDualVars, _s);
	OPTIMIZEBLOCKS(inf.nDualVars, BLOCK1D);
	const size_t min_cap = hc_otsize + inf.nDualVars * hc_olsize + min_lits * hc_srsize;
	assert(min_cap);
	if (otPool.cap < min_cap) { // realloc
		FREE(otPool);
		assert(otPool.mem == NULL);
		if (!hasUnifiedMem(min_cap, "OT")) return false;
		CHECK(cudaMallocManaged((void**)&otPool.mem, min_cap));
		#if !defined(_WIN32)
		if (devProp.major > 5) {
			PFLOGN2(2, " Advising GPU driver to favor global over system memory..");
			CHECK(cudaMemAdvise(otPool.mem, min_cap, cudaMemAdviseSetPreferredLocation, MASTER_GPU));
			CHECK(cudaMemPrefetchAsync(otPool.mem, min_cap, MASTER_GPU, _s));
			PFLDONE(2, 5);
		}
		#endif
		ot = (OT*)otPool.mem;
		LOGERR("Exclusively scanning histogram failed");
		sync(_s); // needed for calling the next constructor on host
		new (ot) OT(inf.nDualVars);
		d_occurs = ot->data();
		assignListPtrs << <nBlocks, BLOCK1D, 0, _s >> > (ot, d_hist, d_segs, inf.nDualVars);
		otPool.cap = min_cap;
	}
	else
		assignListPtrs << <nBlocks, BLOCK1D, 0, _s >> > (ot, d_hist, d_segs, inf.nDualVars);
	if (gopts.sync_always) {
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
	assert(cnfPool.mem);
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

void cuMM::freeVars()
{
	PFLOGN2(2, " Freeing up fixed unified memory..");
	FREE(varsPool);
	d_units = NULL;
	PFLENDING(2, 5, "(remaining: %lld)", cap);
}

void cuMM::freeCNF()
{
	PFLOGN2(2, " Freeing up CNF unified memory..");
	FREE(cnfPool);
	d_cnf_mem = NULL, d_refs_mem = NULL;
	PFLENDING(2, 5, "(remaining: %lld)", cap);
}

void cuMM::freeOT()
{
	PFLOGN2(2, " Freeing up OT unified memory..");
	FREE(otPool);
	PFLENDING(2, 5, "(remaining: %lld)", cap);
}

void cuMM::freeFixed()
{
	PFLOGN2(2, " Freeing up device memory..");
	if (auxPool.mem) {
		DFREE(auxPool);
		d_scatter = NULL, d_stencil = NULL;
		nscatters = 0;
	}
	if (histPool.mem) {
		DFREE(histPool);
		d_segs = NULL, d_hist = NULL;
	}
	if (litsPool.mem) {
		DFREE(litsPool);
		litsPool.size = 0;
	}
	PFLENDING(2, 5, "(remaining: %lld)", dcap);
	dcap = 0;
}

void cuMM::freePinned()
{
	if (pinnedPool.mem) {
		CHECK(cudaFreeHost(pinnedPool.mem)), pinnedPool.mem = NULL;
		pinnedPool.cap = 0;
		pinned_cnf = NULL;
	}
}

void cuMM::breakMirror()
{
	if (hcnfPool.mem) {
		std::free(hcnfPool.mem), hcnfPool.mem = NULL;
		hcnfPool.cap = 0;
	}
}
