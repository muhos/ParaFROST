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

#if defined(USE_CUARENA) && defined(USE_DEVICE_CNF)
__global__
void initCNF_k(CNF* cnf, const S_REF data_cap, const uint32 cs_cap)
{
	new (cnf) CNF(data_cap, cs_cap);
}

__global__
void initOT_k(OT* ot, const uint32 nlists)
{
	new (ot) OT(nlists);
}

__global__
void initVars_k(
	cuVecU* elected, uint32* electedData,
	cuVecU* units, uint32* unitsData,
	cuVecU* resolved, uint32* resolvedData,
	const uint32 varsize, const uint32 resolvedCap)
{
	elected->alloc(electedData, varsize);
	units->alloc(unitsData, varsize);
	resolved->alloc(resolvedData, resolvedCap);
}
#endif

__global__ 
void assignListPtrs(OT* __restrict__ ot, const uint32* __restrict__ hist, const S_REF* __restrict__ segs, const uint32 size)
{
	for_parallel_x(tid, size) {
		assert(segs[tid] < UINT32_MAX);
		(*ot)[tid].alloc(ot->data(segs[tid]), hist[tid]);
	}
}

cuMM::cuMM()
    : litsPool(),                      // Initialize cuLits
      hcnfPool(), cnfPool(), otPool(), // Initialize cuPool instances
      varsPool(), histPool(), auxPool(),
      pinnedPool(),
      cutimer(),                       // Initialize cuTIMER
      pinned_cnf(nullptr), d_refs_mem(nullptr), d_scatter(nullptr),
      d_segs(nullptr), d_occurs(nullptr), d_hist(nullptr),
      d_cnf_mem(nullptr), d_stencil(nullptr), d_vstate(nullptr),
      nscatters(0), device_cnf(false),
      _compacttime(0.0f), _tot(0), _free(0),
      cap(0), dcap(0), maxcap(0), penalty(0),
      isMemAdviseSafe(false)
{
#ifdef USE_CUARENA
    arena_ready = false;
#endif
}

#ifdef USE_CUARENA
bool cuMM::initDeviceArena(const size_t& numCls, const size_t& numLits, const size_t& resolvedCap, const bool& proofEnabled)
{
	if (arena_ready) return true;
	(void)resolvedCap; // only consumed by the USE_DEVICE_CNF dynamic-sizing path below
	// Stable region: histogram, auxiliary, and variable data (excluding literals)
	const size_t varsize      = inf.maxVar + 1;
	const size_t segBytes     = inf.maxDualVars * HC_SREFSIZE;
	const size_t histBytes    = inf.maxDualVars * HC_VARSIZE;
	const size_t vorgBytes    = varsize * HC_VARSIZE;
	const size_t vstateBytes  = varsize * sizeof(VSTATE);
	size_t hist_cap = segBytes + histBytes + vorgBytes + vstateBytes;
	if (proofEnabled) hist_cap += inf.maxDualVars;
	const size_t scatterBytes = numCls * HC_SREFSIZE;
	const size_t aux_cap = scatterBytes + numCls;
	size_t stable_cap = arena.align_up(hist_cap) + arena.align_up(aux_cap);
#if defined(USE_DEVICE_CNF)
	// Fixed vars pool (allocVars): headers + data + scores + resolved + eliminated.
	const size_t vars_fixed = HC_VECSIZE * 3
		+ varsize * HC_VARSIZE * 4 + resolvedCap * HC_VARSIZE + varsize;
	stable_cap += arena.align_up(vars_fixed); // vars lives in the stable region
#endif
	assert(stable_cap);
	// Estimated dynamic peak (lower bound used for sanity/clamping only).
	const size_t lits_dyn_cap = 2 * numLits * HC_VARSIZE;
#if defined(USE_DEVICE_CNF)
	// CNF: two blocks needed at peak (old + new during realloc)
	const size_t cnf_dyn_cap = HC_CNFSIZE + numCls * (SCLAUSESIZE + HC_SREFSIZE) + numLits * SBUCKETSIZE;
	// OT: one block (header + OL array + literal refs)
	const size_t ot_dyn_cap  = HC_OTSIZE + inf.maxDualVars * HC_OLSIZE + numLits * HC_SREFSIZE;
	const size_t dyn_est = arena.align_up(lits_dyn_cap) + 2 * arena.align_up(cnf_dyn_cap) + arena.align_up(ot_dyn_cap);
#else
	const size_t dyn_est = arena.align_up(lits_dyn_cap);
#endif
	size_t dyn_cap = dyn_est;
#if defined(USE_DEVICE_CNF)
	// In the device-CNF path CNF and OT live inside the arena's dynamic region, so
	// their peak can far exceed the first-round estimate (literals/OT blow up across
	// simplification rounds). The CACHER (thrust/cub/moderngpu scratch) also draws
	// from the dynamic region. The fixed vars pool is already part of stable_cap, so
	// hand the rest of the leftover VRAM to the dynamic region, keeping only a margin.
	size_t gpu_free = 0, gpu_tot = 0;
	CHECK(cudaMemGetInfo(&gpu_free, &gpu_tot));
	const size_t reserve  = size_t(penalty) + (256 * MBYTE);              // driver overhead + slack
	if (gpu_free > stable_cap + reserve) {
		const size_t avail = gpu_free - stable_cap - reserve;
		if (avail > dyn_cap) dyn_cap = avail;                            // grow to use leftover VRAM
	}
	else
		LOGWARNING("cuArena: tight VRAM (free %.3f MB, stable %.3f MB, reserve %.3f MB) -> using estimated dynamic %.3f MB",
			double(gpu_free) / MBYTE, double(stable_cap) / MBYTE, double(reserve) / MBYTE, double(dyn_cap) / MBYTE);
#endif
	const size_t total_cap = stable_cap + dyn_cap;
	LOG2(2, " Creating cuArena device pool (stable %.3f MB, dynamic %.3f MB)..",
		double(stable_cap) / MBYTE, double(dyn_cap) / MBYTE);
	try {
		if (!arena.create_gpu_pool(total_cap, cuArena::GPUMemoryType::Device, 0, stable_cap)) {
			LOGWARNING("cuArena: failed to create device pool (%.3f MB) -> falling back to cudaMalloc",
				double(total_cap) / MBYTE);
			return false;
		}
	} catch (const cuArena::gpu_memory_error& e) {
		LOGWARNING("cuArena: device pool exception (%s) -> falling back to cudaMalloc", e.what());
		return false;
	}
	arena_ready = true;
	LOG2(2, " cuArena device pool ready");
	const size_t elimBytes = varsize * sizeof(Byte);
	const size_t unitBytes = varsize * HC_VARSIZE;
	const size_t pinned_cap = HC_CNFSIZE + elimBytes + unitBytes + histBytes;
	size_t cpu_cap = arena.align_up(pinned_cap) + arena.alignment();
	if (proofEnabled) {
		// Proof stream host mirror: matches the device proof (<= numLits*sizeof(uint32)).
		// Generous fixed reserve; the device proof OOMs first, so this never fails alone.
		const size_t proof_pinned = HC_VECSIZE + 2 * numLits * HC_VARSIZE;
		cpu_cap += arena.align_up(proof_pinned) + arena.alignment();
	}
	if (!arena.create_cpu_pool(cpu_cap, cuArena::CPUMemoryType::Pinned))
		LOGWARNING("cuArena: failed to create pinned pool (%.3f MB) -> falling back to cudaHostAlloc",
			double(cpu_cap) / MBYTE);
	return true;
}
#endif

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
		if (!hasDeviceMem(min_cap, "Literals", cuArena::Region::Dynamic)) return NULL;
	#ifdef USE_CUARENA
		if (arena_ready) {
			try {
				litsPool.mem = arena.allocate<uint32>(min_lits, cuArena::Region::Dynamic);
			}
			catch (const cuArena::gpu_memory_error&) {
				LOGN2(2, "  cuArena: dynamic region fragmented, compacting and retrying..");
				arena.compact_gpu_dynamic(nullptr);
				try {
					litsPool.mem = arena.allocate<uint32>(min_lits, cuArena::Region::Dynamic);
				}
				catch (const cuArena::gpu_memory_error& e2) {
					LOGWARNING("cuArena: Literals alloc failed after compact (%s)", e2.what());
					undoDeviceMem(min_cap);
					return NULL;
				}
				LOGENDING(2, 5, "(%.3f MB reclaimed)", double(arena.gpu_available()) / MBYTE);
			}
		}
		else {
			CHECK(cudaMalloc((void**)&litsPool.mem, min_cap));
		}
	#else
		CHECK(cudaMalloc((void**)&litsPool.mem, min_cap));
	#endif
		litsPool.cap = min_cap;
		litsPool.size = min_lits;
	}
	return litsPool.mem;
}

bool cuMM::allocHist(cuHist& cuhist, const bool& proofEnabled)
{
	assert(inf.maxDualVars == V2L(inf.maxVar + 1ULL));
	const size_t varsize     = inf.maxVar + 1;
	const size_t segBytes    = inf.maxDualVars * HC_SREFSIZE;
	const size_t histBytes   = inf.maxDualVars * HC_VARSIZE;
	const size_t varsBytes   = varsize * HC_VARSIZE;
	const size_t vstateBytes = varsize * sizeof(VSTATE);
	size_t min_cap = segBytes + histBytes + varsBytes + vstateBytes;
	if (proofEnabled)
		min_cap += inf.maxDualVars;
	assert(min_cap);
	if (histPool.cap < min_cap) {
		DFREE(histPool);
		assert(histPool.mem == NULL);
		if (!hasDeviceMem(min_cap, "Histogram")) return false;
	#ifdef USE_CUARENA
		try { histPool.mem = arena.allocate<Byte>(min_cap, cuArena::Region::Stable); }
		catch (const cuArena::gpu_memory_error& e) {
			LOGWARNING("cuArena: Histogram alloc failed (%s)", e.what()); undoDeviceMem(min_cap); return false;
		}
	#else
		CHECK(cudaMalloc((void**)&histPool.mem, min_cap));
	#endif
		// NOTE: d_segs, d_hist used internally by OT allocation and externally
		//       by BVE for calculating resolvents offsets (memory reuse)
		//		 lbyte is used for proof byte counting
		addr_t ea = histPool.mem;
		cuhist.d_segs = d_segs = (S_REF*)ea, ea += segBytes;
		cuhist.d_hist = d_hist = (uint32*)ea, ea += histBytes;
		cuhist.d_vorg = (uint32*)ea, ea += varsBytes;
		d_vstate = (VSTATE*)ea, ea += vstateBytes;
		if (proofEnabled) {
			cuhist.d_lbyte = ea;
			ea += inf.maxDualVars;
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
	const size_t varsize = inf.maxVar + 1;
	const size_t uintVec_sz = varsize * HC_VARSIZE;
	const size_t scores_sz = varsize * HC_VARSIZE;
	const size_t resolved_sz = resolvedCap * HC_VARSIZE;
	size_t min_cap = HC_VECSIZE * 3;                             // headers: (elected + units + resolved) 
	min_cap += uintVec_sz * 3 + scores_sz + resolved_sz + varsize; // data:    (elected + units + eligible) + scores + resolved + eliminated
	assert(min_cap);
#if defined(USE_CUARENA) && defined(USE_DEVICE_CNF)
	if (!hasDeviceMem(min_cap, "Fixed", cuArena::Region::Stable)) return false;
	try { varsPool.mem = arena.allocate<Byte>(min_cap, cuArena::Region::Stable); }
	catch (const cuArena::gpu_memory_error& e) {
		LOGWARNING("cuArena: Fixed alloc failed (%s)", e.what()); undoDeviceMem(min_cap); return false;
	}
#else
	if (!hasUnifiedMem(min_cap, "Fixed")) return false;
	CHECK(cudaMallocManaged((void**)&varsPool.mem, min_cap));
#endif
	CHECK(cudaMemset(varsPool.mem, 0, min_cap));
	addr_t ea = varsPool.mem, end = ea + min_cap;
	vars->elected = (cuVecU*)ea, ea += HC_VECSIZE;
	vars->units = (cuVecU*)ea, ea += HC_VECSIZE;
	vars->resolved = (cuVecU*)ea, ea += HC_VECSIZE;
	uint32* uintPtr = (uint32*)ea;
	vars->electedData = uintPtr;
	vars->electedSize = (uint32*)((addr_t)vars->elected + sizeof(uint32*));
	uintPtr += varsize;
	vars->unitsData = uintPtr;
	vars->unitsSize = (uint32*)((addr_t)vars->units + sizeof(uint32*));
	uintPtr += varsize;
	vars->eligible = uintPtr, uintPtr += varsize;
	vars->scores = uintPtr, uintPtr += varsize;
	uint32* resolvedData = uintPtr; uintPtr += resolvedCap;
	Byte* bytePtr = (Byte*)uintPtr;
	vars->eliminated = bytePtr, bytePtr += varsize;
	assert(bytePtr == end);
	varsPool.cap = min_cap;
#if defined(USE_CUARENA) && defined(USE_DEVICE_CNF)
	initVars_k<<<1, 1>>>(vars->elected, vars->electedData, vars->units, vars->unitsData,
		vars->resolved, resolvedData, uint32(varsize), uint32(resolvedCap));
	LASTERR("Vars device init failed");
	SYNC(0);
#else
	vars->elected->alloc(vars->electedData, varsize);
	vars->units->alloc(vars->unitsData, varsize);
	vars->resolved->alloc(resolvedData, uint32(resolvedCap));
	if (isMemAdviseSafe) {
		LOGN2(2, " Advising GPU driver to favor global over system memory in %s call..", __func__);
		addr_t tmpPtr = ea + uintVec_sz; // skip elected
		cudaMemLocation loc{cudaMemLocationTypeDevice, MASTER_GPU};
		CHECK(cudaMemAdvise(tmpPtr, end - tmpPtr, cudaMemAdviseSetPreferredLocation, loc));
		CHECK(cudaMemPrefetchAsync(tmpPtr, end - tmpPtr, loc, 0, 0));
		LOGDONE(2, 5);
	}
#endif
	return true;
}

bool cuMM::allocPinned(VARS* vars, cuHist& cuhist)
{
	assert(vars);
	assert(inf.maxDualVars == V2L(inf.maxVar + 1ULL));
	const size_t varsize = inf.maxVar + 1;
	const size_t elimBytes = varsize * sizeof(Byte);
	const size_t unitBytes = varsize * HC_VARSIZE;
	const size_t histBytes = inf.maxDualVars * HC_VARSIZE;
	size_t min_cap = HC_CNFSIZE + elimBytes + unitBytes + histBytes;
	assert(min_cap);
#ifdef USE_CUARENA
	const bool arena_pinned = arena_ready && arena.cpu_capacity();
#else
	const bool arena_pinned = false;
#endif
	if (pinnedPool.cap) {
		assert(pinnedPool.mem);
	#ifdef USE_CUARENA
		if (arena_pinned) arena.deallocate_pinned(pinnedPool.mem);
		else CHECK(cudaFreeHost(pinnedPool.mem));
	#else
		CHECK(cudaFreeHost(pinnedPool.mem));
	#endif
		pinnedPool.mem = NULL;
		pinnedPool.cap = 0;
	}
	assert(pinnedPool.mem == NULL);
#ifdef USE_CUARENA
	if (arena_pinned) {
		try { pinnedPool.mem = arena.allocate_pinned<Byte>(min_cap); }
		catch (const cuArena::cpu_memory_error& e) {
			LOGWARNING("cuArena: pinned allocation failure (%s)", e.what());
			return false;
		}
	}
	else {
		cudaError_t retVal = cudaHostAlloc((void**)&pinnedPool.mem, min_cap, cudaHostAllocDefault);
		if (retVal != cudaSuccess || retVal == cudaErrorMemoryAllocation) {
			LOGWARNING("Pinned memory allocation failure due to %s", cudaGetErrorString(retVal));
			return false;
		}
	}
#else
	cudaError_t retVal = cudaHostAlloc((void**)&pinnedPool.mem, min_cap, cudaHostAllocDefault);
	if (retVal != cudaSuccess || retVal == cudaErrorMemoryAllocation) {
		LOGWARNING("Pinned memory allocation failure due to %s", cudaGetErrorString(retVal));
		return false;
	}
#endif
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
	#ifdef USE_CUARENA
		try { auxPool.mem = arena.allocate<Byte>(min_cap, cuArena::Region::Stable); }
		catch (const cuArena::gpu_memory_error& e) {
			LOGWARNING("cuArena: Auxiliary alloc failed (%s)", e.what()); undoDeviceMem(min_cap); return false;
		}
	#else
		CHECK(cudaMalloc((void**)&auxPool.mem, min_cap));
	#endif
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
#if defined(USE_CUARENA) && defined(USE_DEVICE_CNF)
	if (cnfPool.cap == 0) {
		assert(cnf == NULL);
		assert(cnfPool.mem == NULL);
		if (!hasDeviceMem(min_cap, "CNF", cuArena::Region::Dynamic)) return false;
		try { cnfPool.mem = arena.allocate<Byte>(min_cap, cuArena::Region::Dynamic); }
		catch (const cuArena::gpu_memory_error& e) {
			LOGWARNING("cuArena: CNF alloc failed (%s)", e.what()); undoDeviceMem(min_cap); return false;
		}
		cnf = (CNF*)cnfPool.mem;
		const S_REF data_cap = S_REF(dataBytes / SBUCKETSIZE);
		initCNF_k<<<1, 1>>>(cnf, data_cap, uint32(clsCap));
		LASTERR("CNF device init failed");
		SYNC(0);
		cacheCNFPtr(cnf);
		SYNC(0);
		d_cnf_mem = pinned_cnf->data().mem;
		d_refs_mem = pinned_cnf->refsData();
		cnfPool.cap = min_cap;
		device_cnf = true;
	}
	else {
		assert(cnf);
		assert(cnfPool.mem);
		if (!hasDeviceMem(min_cap, "CNF", cuArena::Region::Dynamic)) return false;
		addr_t newMem = NULL;
		try { newMem = arena.allocate<Byte>(min_cap, cuArena::Region::Dynamic); }
		catch (const cuArena::gpu_memory_error&) {
			LOGN2(2, "  cuArena: dynamic region fragmented, compacting and retrying CNF..");
			arena.compact_gpu_dynamic(nullptr);
			try { newMem = arena.allocate<Byte>(min_cap, cuArena::Region::Dynamic); }
			catch (const cuArena::gpu_memory_error& e2) {
				LOGWARNING("cuArena: CNF realloc failed after compact (%s)", e2.what());
				undoDeviceMem(min_cap); return false;
			}
		}
		CNF* tmp_cnf = (CNF*)newMem;
		const S_REF data_cap = S_REF(dataBytes / SBUCKETSIZE);
		initCNF_k<<<1, 1>>>(tmp_cnf, data_cap, uint32(clsCap));
		LASTERR("CNF device reinit failed");
		SYNC(0);
		// Cache OLD cnf into pinned_cnf so compactCNF sees the correct old size.
		// Get new cnf's data/refs pointers via a separate local copy.
		cacheCNFPtr(cnf);
		CNF new_cnf_hdr;
		CHECK(cudaMemcpy(&new_cnf_hdr, tmp_cnf, sizeof(CNF), cudaMemcpyDeviceToHost));
		SYNC(0);
		d_cnf_mem = new_cnf_hdr.data().mem;
		d_refs_mem = new_cnf_hdr.refsData();
		compactCNF(cnf, tmp_cnf);
		DFREE(cnfPool);
		cnfPool.mem = newMem;
		cnfPool.cap = min_cap;
		cnf = tmp_cnf;
	}
#else
	if (cnfPool.cap == 0) {
		assert(cnf == NULL);
		assert(cnfPool.mem == NULL);
		if (!hasUnifiedMem(min_cap, "CNF")) return false;
		CHECK(cudaMallocManaged((void**)&cnfPool.mem, min_cap));
		if (isMemAdviseSafe) {
			LOGN2(2, " Advising GPU driver to favor global over system memory in %s call..", __func__);
			cudaMemLocation loc{cudaMemLocationTypeDevice, MASTER_GPU};
			CHECK(cudaMemAdvise(cnfPool.mem, min_cap, cudaMemAdviseSetPreferredLocation, loc));
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
			cudaMemLocation loc{cudaMemLocationTypeDevice, MASTER_GPU};
			CHECK(cudaMemAdvise(newMem, min_cap, cudaMemAdviseSetPreferredLocation, loc));
			CHECK(cudaMemPrefetchAsync(newMem, min_cap, loc, 0, 0));
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
#endif
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
	DeviceScan::ExclusiveSum(NULL, ebytes, d_hist, d_segs, inf.maxDualVars, _s);
	DeviceScan::ExclusiveSum(tmp, ebytes, d_hist, d_segs, inf.maxDualVars, _s);
	grid_t nThreads(BLOCK1D);
	OPTIMIZEBLOCKS(inf.maxDualVars, nThreads, 0);
	const size_t min_cap = HC_OTSIZE + inf.maxDualVars * HC_OLSIZE + min_lits * HC_SREFSIZE;
	assert(min_cap);
#if defined(USE_CUARENA) && defined(USE_DEVICE_CNF)
	if (otPool.cap < min_cap) { // realloc
		if (otPool.mem) DFREE(otPool);
		assert(otPool.mem == NULL);
		if (!hasDeviceMem(min_cap, "OT", cuArena::Region::Dynamic)) return false;
		try { otPool.mem = arena.allocate<Byte>(min_cap, cuArena::Region::Dynamic); }
		catch (const cuArena::gpu_memory_error&) {
			LOGN2(2, "  cuArena: dynamic region fragmented, compacting and retrying OT..");
			arena.compact_gpu_dynamic(nullptr);
			try { otPool.mem = arena.allocate<Byte>(min_cap, cuArena::Region::Dynamic); }
			catch (const cuArena::gpu_memory_error& e2) {
				LOGWARNING("cuArena: OT alloc failed after compact (%s)", e2.what());
				undoDeviceMem(min_cap);
				return false;
			}
		}
		ot = (OT*)otPool.mem;
		LASTERR("Exclusively scanning histogram failed");
		SYNC(_s);
		initOT_k<<<1, 1>>>(ot, inf.maxDualVars);
		LASTERR("OT device init failed");
		SYNC(0);
		OT ot_hdr;
		CHECK(cudaMemcpy(&ot_hdr, ot, sizeof(OT), cudaMemcpyDeviceToHost));
		d_occurs = ot_hdr.data();
		assignListPtrs<<<nBlocks, nThreads, 0, _s>>>(ot, d_hist, d_segs, inf.maxDualVars);
		otPool.cap = min_cap;
	}
	else
		assignListPtrs<<<nBlocks, nThreads, 0, _s>>>(ot, d_hist, d_segs, inf.maxDualVars);
#else
	if (otPool.cap < min_cap) { // realloc
		FREE(otPool);
		assert(otPool.mem == NULL);
		if (!hasUnifiedMem(min_cap, "OT")) return false;
		CHECK(cudaMallocManaged((void**)&otPool.mem, min_cap));
		if (isMemAdviseSafe) {
			LOGN2(2, " Advising GPU driver to favor global over system memory in %s call..", __func__);
			cudaMemLocation loc{cudaMemLocationTypeDevice, MASTER_GPU};
			CHECK(cudaMemAdvise(otPool.mem, min_cap, cudaMemAdviseSetPreferredLocation, loc));
			CHECK(cudaMemPrefetchAsync(otPool.mem, min_cap, loc, 0, _s));
			LOGDONE(2, 5);
		}
		ot = (OT*)otPool.mem;
		LASTERR("Exclusively scanning histogram failed");
		SYNC(_s); // needed for calling the next constructor on host
		new (ot) OT(inf.maxDualVars);
		d_occurs = ot->data();
		assignListPtrs<<<nBlocks, nThreads, 0, _s>>>(ot, d_hist, d_segs, inf.maxDualVars);
		otPool.cap = min_cap;
	}
	else
		assignListPtrs<<<nBlocks, nThreads, 0, _s>>>(ot, d_hist, d_segs, inf.maxDualVars);
#endif
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
#if defined(USE_CUARENA) && defined(USE_DEVICE_CNF)
	LOGN2(2, " freeing up fixed device memory..");
	DFREE(varsPool);
#else
	LOGN2(2, " freeing up fixed unified memory..");
	FREE(varsPool);
#endif
	LOGDONE(2, 5);;
}

void cuMM::freeCNF()
{
	if (device_cnf) {
		LOGN2(2, " freeing up CNF device memory..");
		DFREE(cnfPool);
		LOGDONE(2, 5);
	}
	else {
		LOGN2(2, " freeing up CNF unified memory..");
		FREE(cnfPool);
		LOGDONE(2, 5);;
	}
	d_cnf_mem = NULL, d_refs_mem = NULL;
}

void cuMM::freeOT()
{
	if (device_cnf) {
		LOGN2(2, " freeing up OT device memory..");
		DFREE(otPool);
		LOGDONE(2, 5);
	}
	else {
		LOGN2(2, " freeing up occurrence table unified memory..");
		FREE(otPool);
		LOGDONE(2, 5);;
	}
}

void cuMM::freeFixed()
{
	LOGN2(2, " freeing up histogram and auxiliary memory..");
	if (auxPool.mem) {
		DFREE(auxPool);
		d_scatter = NULL, d_stencil = NULL;
		nscatters = 0;
	}
	if (histPool.mem) {
		DFREE(histPool);
		d_segs = NULL, d_hist = NULL, d_vstate = NULL;
	}
	if (litsPool.mem) {
		DFREE(litsPool);
		litsPool.size = 0;
	}
    LOGDONE(2, 5);
}

void cuMM::freePinned()
{
	LOGN2(2, " freeing up pinned memory..");
	if (pinnedPool.mem) {
	#ifdef USE_CUARENA
		if (arena_ready && arena.cpu_capacity()) arena.deallocate_pinned(pinnedPool.mem);
		else CHECK(cudaFreeHost(pinnedPool.mem));
	#else
		CHECK(cudaFreeHost(pinnedPool.mem));
	#endif
		pinnedPool.mem = NULL;
		pinnedPool.cap = 0;
		pinned_cnf = NULL;
	}
    LOGDONE(2, 5);
}

void cuMM::breakMirror()
{
	LOGN2(2, " freeing up CNF host memory..");
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
