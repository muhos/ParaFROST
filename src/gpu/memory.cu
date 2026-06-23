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
      nscatters(0),
      _compacttime(0.0f), _tot(0), _free(0),
      cap(0), dcap(0), penalty(0)
{
    arena_ready = false;
}

bool cuMM::initDeviceArena(const size_t& numCls, const size_t& numLits, const size_t& resolvedCap, const bool& proofEnabled)
{
	if (arena_ready) return true;
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
	// Fixed vars pool (allocVars): headers + data + scores + resolved + eliminated.
	const size_t vars_fixed = HC_VECSIZE * 3
		+ varsize * HC_VARSIZE * 4 + resolvedCap * HC_VARSIZE + varsize;
	stable_cap += arena.align_up(vars_fixed); // vars lives in the stable region
	assert(stable_cap);
	// Estimated dynamic peak (lower bound used for sanity/clamping only).
	const size_t lits_dyn_cap = 2 * numLits * HC_VARSIZE;
	// CNF: two blocks needed at peak (old + new during realloc)
	const size_t cnf_dyn_cap = HC_CNFSIZE + numCls * (SCLAUSESIZE + HC_SREFSIZE) + numLits * SBUCKETSIZE;
	// OT: one block (header + OL array + literal refs)
	const size_t ot_dyn_cap  = HC_OTSIZE + inf.maxDualVars * HC_OLSIZE + numLits * HC_SREFSIZE;
	const size_t dyn_est = arena.align_up(lits_dyn_cap) + 2 * arena.align_up(cnf_dyn_cap) + arena.align_up(ot_dyn_cap);
	size_t dyn_cap = dyn_est;
	const size_t gpu_free = getFreeMemory();
	const size_t reserve  = size_t(penalty) + (256 * MBYTE);              // driver overhead
	if (gpu_free > stable_cap + reserve) {
		const size_t avail = gpu_free - stable_cap - reserve;
		if (avail > dyn_cap) dyn_cap = avail;                             // grow to use leftover VRAM
	}
	else
		LOGWARNING("cuArena: low VRAM (free %.3f MB, stable %.3f MB, reserve %.3f MB) -> using estimated dynamic %.3f MB",
			double(gpu_free) / MBYTE, double(stable_cap) / MBYTE, double(reserve) / MBYTE, double(dyn_cap) / MBYTE);
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
		try { histPool.mem = arena.allocate<Byte>(min_cap, cuArena::Region::Stable); }
		catch (const cuArena::gpu_memory_error& e) {
			LOGWARNING("cuArena: Histogram alloc failed (%s)", e.what()); undoDeviceMem(min_cap); return false;
		}
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
	if (!hasDeviceMem(min_cap, "Fixed", cuArena::Region::Stable)) return false;
	try { varsPool.mem = arena.allocate<Byte>(min_cap, cuArena::Region::Stable); }
	catch (const cuArena::gpu_memory_error& e) {
		LOGWARNING("cuArena: Fixed alloc failed (%s)", e.what()); undoDeviceMem(min_cap); return false;
	}
	CHECK(cudaMemset(varsPool.mem, 0, min_cap));
	addr_t ea = varsPool.mem;
	#if !defined(NDEBUG)
	addr_t end = ea + min_cap;
	#endif
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
	#if !defined(NDEBUG)
	assert(bytePtr == end);
	#endif
	varsPool.cap = min_cap;
	initVars_k<<<1, 1>>>(vars->elected, vars->electedData, vars->units, vars->unitsData,
		vars->resolved, resolvedData, uint32(varsize), uint32(resolvedCap));
	LASTERR("Vars device init failed");
	SYNC(0);
	return true;
}

bool cuMM::allocPinned(VARS* vars, cuHist& cuhist)
{
	assert(vars);
	assert(inf.maxDualVars == V2L(inf.maxVar + 1ULL));
	const size_t varsize = inf.maxVar + 1;
	const size_t elimBytes = varsize * sizeof(Byte);
	const size_t histBytes = inf.maxDualVars * HC_VARSIZE;
	size_t min_cap = HC_CNFSIZE + elimBytes + histBytes;
	assert(min_cap);
	const bool arena_pinned = arena_ready && arena.cpu_capacity();
	if (pinnedPool.cap) {
		assert(pinnedPool.mem);
		if (arena_pinned) arena.deallocate_pinned(pinnedPool.mem);
		else CHECK(cudaFreeHost(pinnedPool.mem));
		pinnedPool.mem = NULL;
		pinnedPool.cap = 0;
	}
	assert(pinnedPool.mem == NULL);
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
	addr_t ea = pinnedPool.mem;
	pinned_cnf = (CNF*)ea, ea += HC_CNFSIZE;
	cuhist.h_hist = (uint32*)ea, ea += histBytes;
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
		try { auxPool.mem = arena.allocate<Byte>(min_cap, cuArena::Region::Stable); }
		catch (const cuArena::gpu_memory_error& e) {
			LOGWARNING("cuArena: Auxiliary alloc failed (%s)", e.what()); undoDeviceMem(min_cap); return false;
		}
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
		if (!hasDeviceMem(min_cap, "CNF", cuArena::Region::Dynamic)) return false;
		try { cnfPool.mem = arena.allocate<Byte>(min_cap, cuArena::Region::Dynamic); }
		catch (const cuArena::gpu_memory_error& e) {
			LOGWARNING("cuArena: CNF alloc failed (%s)", e.what()); undoDeviceMem(min_cap); return false;
		}
		cnf = (CNF*)cnfPool.mem;
		const S_REF data_cap = S_REF(dataBytes / SBUCKETSIZE);
		initCNF_k<<<1, 1>>>(cnf, data_cap, uint32(clsCap));
		LASTERR("CNF device init failed");
		cacheCNFPtr(cnf);
		SYNC(0);
		d_cnf_mem = pinned_cnf->data().mem;
		d_refs_mem = pinned_cnf->refsData();
		cnfPool.cap = min_cap;
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
		cacheCNFPtr(cnf);
		CNF new_cnf_hdr;
		CHECK(cudaMemcpy(&new_cnf_hdr, tmp_cnf, sizeof(CNF), cudaMemcpyDeviceToHost));
		d_cnf_mem = new_cnf_hdr.data().mem;
		d_refs_mem = new_cnf_hdr.refsData();
		compactCNF(cnf, tmp_cnf);
		DFREE(cnfPool);
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
	DeviceScan::ExclusiveSum(NULL, ebytes, d_hist, d_segs, inf.maxDualVars, _s);
	DeviceScan::ExclusiveSum(tmp, ebytes, d_hist, d_segs, inf.maxDualVars, _s);
	grid_t nThreads(BLOCK1D);
	OPTIMIZEBLOCKS(inf.maxDualVars, nThreads, 0);
	const size_t min_cap = HC_OTSIZE + inf.maxDualVars * HC_OLSIZE + min_lits * HC_SREFSIZE;
	assert(min_cap);
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
		initOT_k<<<1, 1>>>(ot, inf.maxDualVars);
		LASTERR("OT device init failed");
		OT ot_hdr;
		CHECK(cudaMemcpy(&ot_hdr, ot, sizeof(OT), cudaMemcpyDeviceToHost));
		d_occurs = ot_hdr.data();
		assignListPtrs<<<nBlocks, nThreads, 0, _s>>>(ot, d_hist, d_segs, inf.maxDualVars);
		otPool.cap = min_cap;
	}
	else
		assignListPtrs<<<nBlocks, nThreads, 0, _s>>>(ot, d_hist, d_segs, inf.maxDualVars);
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

void cuMM::freeDevice()
{
	LOGN2(2, " Freeing up device memory..");
	DFREE(varsPool);
	DFREE(cnfPool);
	DFREE(otPool);
	DFREE(auxPool);
	DFREE(histPool);
	DFREE(litsPool);
	d_cnf_mem = NULL, d_refs_mem = NULL;
	d_occurs = NULL, d_scatter = NULL;
	d_stencil = NULL, d_vstate = NULL;
	d_segs = NULL, d_hist = NULL;
	nscatters = 0;
	litsPool.size = 0;
	LOGENDING(2, 5, "(%.3f MB freed)", double(_free) / MBYTE);
}

void cuMM::freePinned()
{
	LOGN2(2, " Freeing up pinned memory..");
	const size_t pinned_cap = pinnedPool.cap;
	if (pinnedPool.mem) {
		if (arena_ready && arena.cpu_capacity()) arena.deallocate_pinned(pinnedPool.mem);
		else CHECK(cudaFreeHost(pinnedPool.mem));
		pinnedPool.mem = NULL;
		pinnedPool.cap = 0;
		pinned_cnf = NULL;
	}
    LOGENDING(2, 5, "(%.3f MB freed)", double(pinned_cap) / MBYTE);
}

void cuMM::breakMirror()
{
	LOGN2(2, " Freeing up CNF host memory..");
	const size_t mirror_cap = hcnfPool.cap;
	if (hcnfPool.mem) {
		std::free(hcnfPool.mem);
		hcnfPool.mem = NULL;
		hcnfPool.cap = 0;
	}
    LOGENDING(2, 5, "(%.3f MB freed)", double(mirror_cap) / MBYTE);
}

