/***********************************************************************[memory.cu]
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
#include "memory.cuh"
#include "options.cuh"
#include "definitions.hpp"
#include <cub/device/device_scan.cuh>
#include <cub/device/device_select.cuh>

using namespace cub;
using namespace ParaFROST;

//=============================//
//    CUDA memory management   //
//=============================//
constexpr size_t HC_SREFSIZE = sizeof(S_REF);
constexpr size_t HC_OTSIZE	 = sizeof(OT);
constexpr size_t HC_OLSIZE   = sizeof(OL);
constexpr size_t HC_CNFSIZE  = sizeof(CNF);
constexpr size_t HC_VARSIZE  = sizeof(uint32);
constexpr size_t HC_VECSIZE  = sizeof(cuVecU);

template<class T>
__global__ 
void memset_k(T* mem, T val, size_t size)
{
	for_parallel_x(tid, size) { mem[tid] = val; }
}

__global__ 
void resizeCNF_k(CNF* cnf, const S_REF d_size, const uint32 cs_size)
{
	cnf->resize(d_size, cs_size);
}

__global__ 
void assignListPtrs(OT* __restrict__ ot, const uint32* __restrict__ hist, const S_REF* __restrict__ segs, const uint32 size)
{
	for_parallel_x(tid, size) {
		assert(segs[tid] < UINT32_MAX);
		(*ot)[tid].alloc(ot->data(segs[tid]), hist[tid]);
	}
}

void cuMM::cuMemSetAsync(addr_t mem, const Byte& val, const size_t& size)
{
	grid_t nThreads(BLOCK1D);
	OPTIMIZEBLOCKS(uint32(size), nThreads, 0);
	memset_k<Byte> << <nBlocks, nThreads >> > (mem, val, size);
	if (gopts.sync_always) {
		LASTERR("CUDA memory set failed");
		SYNCALL;
	}
}

void cuMM::resizeCNFAsync(CNF* dcnf, const S_REF& data_size, const uint32& cs_size)
{
	assert(dcnf);
	assert(data_size);
	assert(cs_size);
	resizeCNF_k << <1, 1 >> > (dcnf, data_size, cs_size);
	if (gopts.sync_always) {
		LASTERR("Resizing CNF failed");
		SYNC(0);
	}
}

uint32* cuMM::resizeLits(const size_t& min_lits)
{
	assert(min_lits);
	const size_t min_cap = min_lits * HC_VARSIZE;
	if (litsPool.cap < min_cap) {
		DFREE(litsPool);
		assert(litsPool.mem == NULL);
		if (!hasDeviceMem(min_cap, "Literals")) return NULL;
		CHECK(cudaMalloc((void**)&litsPool.mem, min_cap));
		litsPool.cap = min_cap;
		litsPool.size = min_lits;
	}
	return litsPool.mem;
}

bool cuMM::allocHist(cuHist& cuhist, const bool& proofEnabled)
{
	assert(inf.nDualVars == V2L(inf.maxVar + 1ULL));
	const size_t segBytes = inf.nDualVars * HC_SREFSIZE;
	const size_t histBytes = inf.nDualVars * HC_VARSIZE;
	const size_t varsBytes = (inf.maxVar + 1) * HC_VARSIZE;
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
	const size_t uintVec_sz = inf.maxVar * HC_VARSIZE;
	const size_t varsize = inf.maxVar + 1;
	const size_t scores_sz = varsize * HC_VARSIZE;
	const size_t resolved_sz = resolvedCap * HC_VARSIZE;
	size_t min_cap = HC_VECSIZE * 3;                             // headers: (elected + units + resolved) 
	min_cap += uintVec_sz * 3 + scores_sz + resolved_sz + varsize; // data:    (elected + units + eligible) + scores + resolved + eliminated
	assert(min_cap);
	if (!hasUnifiedMem(min_cap, "Fixed")) return false;
	CHECK(cudaMallocManaged((void**)&varsPool.mem, min_cap));
	CHECK(cudaMemsetAsync(varsPool.mem, 0, min_cap));
	addr_t ea = varsPool.mem, end = ea + min_cap;
	vars->elected = (cuVecU*)ea, ea += HC_VECSIZE;
	vars->units = (cuVecU*)ea, ea += HC_VECSIZE;
	vars->resolved = (cuVecU*)ea, ea += HC_VECSIZE;
	uint32* uintPtr = (uint32*)ea;
	vars->electedData = uintPtr;
	vars->electedSize = (uint32*)((addr_t)vars->elected + sizeof(uint32*));
	SYNCALL; // sync. cudaMemsetAsync
	vars->elected->alloc(uintPtr, inf.maxVar), uintPtr += inf.maxVar, d_units = uintPtr;
	vars->units->alloc(uintPtr, inf.maxVar), uintPtr += inf.maxVar;
	vars->eligible = uintPtr, uintPtr += inf.maxVar;
	vars->scores = uintPtr, uintPtr += varsize;
	vars->resolved->alloc(uintPtr, uint32(resolvedCap)), uintPtr += resolvedCap;
	Byte* bytePtr = (Byte*)uintPtr;
	vars->eliminated = bytePtr, bytePtr += varsize;
	assert(bytePtr == end);
	varsPool.cap = min_cap;
	if (isMemAdviseSafe) {
		LOGN2(2, " Advising GPU driver to favor global over system memory in %s call..", __func__);
		addr_t tmpPtr = ea + uintVec_sz; // skip elected
		CHECK(cudaMemAdvise(tmpPtr, end - tmpPtr, cudaMemAdviseSetPreferredLocation, MASTER_GPU));
		CHECK(cudaMemPrefetchAsync(tmpPtr, end - tmpPtr, MASTER_GPU));
		LOGDONE(2, 5);
	}
	return true;
}

bool cuMM::allocPinned(VARS* vars, cuHist& cuhist)
{
	assert(vars);
	assert(inf.nDualVars == V2L(inf.maxVar + 1ULL));
	const size_t elimBytes = inf.maxVar + 1;
	const size_t unitBytes = inf.maxVar * HC_VARSIZE;
	const size_t histBytes = inf.nDualVars * HC_VARSIZE;
	size_t min_cap = HC_CNFSIZE + elimBytes + unitBytes + histBytes;
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
		LOGWARNING("Pinned memory allocation failure due to %s", cudaGetErrorString(retVal));
		return false;
	}
	addr_t ea = pinnedPool.mem;
	pinned_cnf = (CNF*)ea, ea += HC_CNFSIZE;
	cuhist.h_hist = (uint32*)ea, ea += histBytes;
	vars->cachedUnits = (uint32*)ea, ea += unitBytes;
	vars->cachedEliminated = ea, ea += elimBytes;
    assert(ea == pinnedPool.mem + min_cap);
	pinnedPool.cap = min_cap;
	return true;
}

bool cuMM::allocAux(const size_t& clsCap)
{
	assert(clsCap && clsCap <= UINT32_MAX);
	const size_t scatterBytes = clsCap * HC_SREFSIZE;
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
	const size_t csBytes = clsCap * HC_SREFSIZE;
	const size_t dataBytes = clsCap * SCLAUSESIZE + litsCap * SBUCKETSIZE;
	assert(dataBytes % SBUCKETSIZE == 0);
	const size_t min_cap = HC_CNFSIZE + dataBytes + csBytes;
	assert(min_cap);
	if (cnfPool.cap == 0) {
		assert(cnf == NULL);
		assert(cnfPool.mem == NULL);
		if (!hasUnifiedMem(min_cap, "CNF")) return false;
		CHECK(cudaMallocManaged((void**)&cnfPool.mem, min_cap));
		if (isMemAdviseSafe) {
			LOGN2(2, " Advising GPU driver to favor global over system memory in %s call..", __func__);
			CHECK(cudaMemAdvise(cnfPool.mem, min_cap, cudaMemAdviseSetPreferredLocation, MASTER_GPU));
			LOGDONE(2, 5);
		}
		cnf = (CNF*)cnfPool.mem;
		const S_REF data_cap = S_REF(dataBytes / SBUCKETSIZE);
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
		SYNC(0);
		if (isMemAdviseSafe) {
			LOGN2(2, " Advising GPU driver to favor global over system memory in %s call..", __func__);
			CHECK(cudaMemAdvise(newMem, min_cap, cudaMemAdviseSetPreferredLocation, MASTER_GPU));
			CHECK(cudaMemPrefetchAsync(newMem, min_cap, MASTER_GPU));
			LOGDONE(2, 5);
		}
		CNF* tmp_cnf = (CNF*)newMem;
		const S_REF data_cap = S_REF(dataBytes / SBUCKETSIZE);
		new (tmp_cnf) CNF(data_cap, uint32(clsCap));
		d_cnf_mem = tmp_cnf->data().mem;
		d_refs_mem = tmp_cnf->refsData();
		compactCNF(cnf, tmp_cnf);
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
	grid_t nThreads(BLOCK1D);
	OPTIMIZEBLOCKS(inf.nDualVars, nThreads, 0);
	const size_t min_cap = HC_OTSIZE + inf.nDualVars * HC_OLSIZE + min_lits * HC_SREFSIZE;
	assert(min_cap);
	if (otPool.cap < min_cap) { // realloc
		FREE(otPool);
		assert(otPool.mem == NULL);
		if (!hasUnifiedMem(min_cap, "OT")) return false;
		CHECK(cudaMallocManaged((void**)&otPool.mem, min_cap));
		if (isMemAdviseSafe) {
			LOGN2(2, " Advising GPU driver to favor global over system memory in %s call..", __func__);
			CHECK(cudaMemAdvise(otPool.mem, min_cap, cudaMemAdviseSetPreferredLocation, MASTER_GPU));
			CHECK(cudaMemPrefetchAsync(otPool.mem, min_cap, MASTER_GPU, _s));
			LOGDONE(2, 5);
		}
		ot = (OT*)otPool.mem;
		LASTERR("Exclusively scanning histogram failed");
		SYNC(_s); // needed for calling the next constructor on host
		new (ot) OT(inf.nDualVars);
		d_occurs = ot->data();
		assignListPtrs << <nBlocks, nThreads, 0, _s >> > (ot, d_hist, d_segs, inf.nDualVars);
		otPool.cap = min_cap;
	}
	else
		assignListPtrs << <nBlocks, nThreads, 0, _s >> > (ot, d_hist, d_segs, inf.nDualVars);
	if (gopts.sync_always) {
		LASTERR("Occurrence lists allocation failed");
		SYNC(_s);
	}
	return true;
}

void cuMM::createMirror(CNF*& hcnf, const size_t& clsCap, const size_t& litsCap)
{
	assert(clsCap && clsCap <= UINT32_MAX);
	assert(litsCap && litsCap <= UINT32_MAX);
	assert(litsCap >= clsCap);
	const size_t csBytes = clsCap * HC_SREFSIZE;
	const size_t dataBytes = clsCap * SCLAUSESIZE + litsCap * SBUCKETSIZE;
	assert(dataBytes % SBUCKETSIZE == 0);
	const size_t min_cap = HC_CNFSIZE + dataBytes + csBytes;
	assert(min_cap);
	if (hcnfPool.cap < min_cap) {
		hcnfPool.cap = min_cap;
		pfralloc(hcnfPool.mem, hcnfPool.cap);
		hcnf = (CNF*)hcnfPool.mem;
	}
	const S_REF data_cap = S_REF(dataBytes / SBUCKETSIZE);
	new (hcnf) CNF(data_cap, uint32(clsCap));
}

void cuMM::mirrorCNF(CNF*& hcnf)
{
	assert(cnfPool.cap);
	assert(cnfPool.mem);
	CHECK(cudaMemcpy(hcnf, cnfPool.mem, HC_CNFSIZE, cudaMemcpyDeviceToHost));
	const size_t csBytes = hcnf->size() * HC_SREFSIZE;
	const size_t dataBytes = hcnf->data().size * SBUCKETSIZE;
	const size_t min_cap = HC_CNFSIZE + dataBytes + csBytes;
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
	LOGN2(2, "  freeing up fixed unified memory..");
	FREE(varsPool);
	d_units = NULL;
	LOGENDING(2, 5, "(remaining: %lld)", cap);
}

void cuMM::freeCNF()
{
	LOGN2(2, "  freeing up CNF unified memory..");
	FREE(cnfPool);
	d_cnf_mem = NULL, d_refs_mem = NULL;
	LOGENDING(2, 5, "(remaining: %lld)", cap);
}

void cuMM::freeOT()
{
	LOGN2(2, "  freeing up occurrence table unified memory..");
	FREE(otPool);
	LOGENDING(2, 5, "(remaining: %lld)", cap);
}

void cuMM::freeFixed()
{
	LOGN2(2, "  freeing up histogram and auxiliary memory..");
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
	dcap = 0;
    LOGDONE(2, 5);
}

void cuMM::freePinned()
{
	LOGN2(2, "  freeing up pinned memory..");
	if (pinnedPool.mem) {
		CHECK(cudaFreeHost(pinnedPool.mem));
		pinnedPool.mem = NULL;
		pinnedPool.cap = 0;
		pinned_cnf = NULL;
	}
    LOGDONE(2, 5);
}

void cuMM::breakMirror()
{
	LOGN2(2, "  freeing up CNF host memory..");
	if (hcnfPool.mem) {
		std::free(hcnfPool.mem);
		hcnfPool.mem = NULL;
		hcnfPool.cap = 0;
	}
    LOGDONE(2, 5);
}

bool cuMM::checkMemAdvice()
{
	int concurrentManaged = 0;
	cudaDeviceGetAttribute(&concurrentManaged, cudaDevAttrConcurrentManagedAccess, MASTER_GPU);
	isMemAdviseSafe = concurrentManaged ? true : false;
	return isMemAdviseSafe;
}
