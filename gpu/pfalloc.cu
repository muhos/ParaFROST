
#include "pfsolve.h"
//========================//
//	GPU memory management //
//========================//

__global__ void set_cnf_ptrs(CNF* cnf, uint32 clsCap, uint64 litsCap) { 
	assert(blockDim.x * gridDim.x == 1);
	cnf->allocMem(clsCap, litsCap); 
}
__global__ void copy_cls_k(CNF* dest, CNF* src, uint32 clsCap, uint64 litsCap) {
	assert(blockDim.x * gridDim.x == 1);
	dest->allocMem(clsCap, litsCap);
	uint32 sz = src->size();
	uint64 nLits = 0;
	dest->resize(sz);
	for (uint32 i = 0; i < sz; i++) {
		uint32 c_sz = (*src)[i].size();
		assert((*src)[i].status() == ORIGINAL);
		(*dest)[i].set_ptr(dest->data(nLits));
		(*dest)[i].set_status(ORIGINAL);
		(*dest)[i].set_sig((*src)[i].sig());
		(*dest)[i].resize(c_sz);
		nLits += c_sz;
	}
	assert(nLits == src->numLits());
	dest->resizeData(nLits);
}
__global__ void copy_cnf_k(CNF* dest, CNF* src)
{
	uint64 i = blockDim.x * blockIdx.x + threadIdx.x;
	uint64 stride = blockDim.x * gridDim.x;
	while (i < src->numLits()) { *dest->data(i) = *src->data(i); i += stride; }
}
__global__ void copy_cnf_k(CNF* dest, CNF* src, uint32 clsCap, uint64 litsCap) {
	assert(blockDim.x * gridDim.x == 1);
	dest->allocMem(clsCap, litsCap);
	dest->copyFrom(src); 
}
__global__ void ot_ptrs_k(OT* ot, uint32* hist, int64 nLists, int64 nEntries, uint32 size)
{
	ot->allocMem(nLists, nEntries);
	int64 idx = 0;
	for (uint32 v = 0; v < size; v++) {
		uint32 p = V2D(v + 1), n = NEG(p);
		uint32 ps = hist[p], ns = hist[n];
		(*ot)[p].allocList(ot->data(idx), ps), idx += ps;
		(*ot)[n].allocList(ot->data(idx), ns), idx += ns;
	}
	assert(idx <= ot->capacity());
}

void cuMM::allocPV(PV* pv, const uint32& maxVars) {
	assert(gMemPV == NULL);
	assert(gMemSol == NULL);
	assert(pv->pVars == NULL);
	assert(pv->sol == NULL);
	size_t sol_sz = sizeof(GSOL);
	size_t cuVec_sz = sizeof(cuVec<uint32>);
	size_t uintVec_sz = cuVec_sz + (size_t)maxVars * sizeof(uint32);
	size_t litStVec_sz = maxVars * sizeof(LIT_ST);
	gMemPV_sz = uintVec_sz, gMemSol_sz = sol_sz + uintVec_sz + litStVec_sz;
	assert(gMemPV_sz > 0);
	assert(gMemSol_sz > 0);
	CHECK(cudaMallocManaged((void**)&gMemPV, gMemPV_sz));
	CHECK(cudaMallocManaged((void**)&gMemSol, gMemSol_sz));
	pv->pVars = (cuVec<uint32>*)gMemPV;
	pv->pVars->alloc((uint32*)(gMemPV + cuVec_sz), maxVars);
	pv->sol = (GSOL*)gMemSol;
	pv->sol->assigns = (cuVec<uint32>*)(gMemSol + sol_sz);
	pv->sol->assigns->alloc((uint32*)(gMemSol + sol_sz + cuVec_sz), maxVars);
	pv->sol->value = (LIT_ST*)(gMemSol + sol_sz + uintVec_sz);
	mem_set(pv->sol->value, LIT_ST(UNDEFINED), maxVars);
	cap += gMemPV_sz + gMemSol_sz;
}
void cuMM::allocVO(OCCUR** occurs, SCORE** scores, const uint32& maxVars) {
	assert(gMemVOrd == NULL);
	assert(*occurs == NULL);
	assert(*scores == NULL);
	gMemVOrd_sz = (size_t)maxVars * (sizeof(OCCUR) + sizeof(SCORE));
	assert(gMemVOrd_sz > 0);
	CHECK(cudaMallocManaged((void**)&gMemVOrd, gMemVOrd_sz));
	*occurs = (OCCUR*)gMemVOrd;
	*scores = (SCORE*)(*occurs + maxVars);
	assert(addr_t(*scores + maxVars) == gMemVOrd + gMemVOrd_sz);
	cap += gMemVOrd_sz;
}
void cuMM::allocStats(GSTATS** gstats, const uint32& maxVars) {
	assert(gMemStats == NULL);
	assert(*gstats == NULL);
	gMemStats_sz = sizeof(GSTATS) + maxVars * sizeof(Byte);
	CHECK(cudaMallocManaged((void**)&gMemStats, gMemStats_sz));
	*gstats = (GSTATS*)gMemStats;
	(*gstats)->seen = gMemStats + sizeof(GSTATS);
	cap += gMemStats_sz;
}
CNF* cuMM::resizeCNF(CNF* cnf, const uint32& clsCap, const uint64& litsCap, const bool& _pc) {
	CNF* tmp_cnf = NULL;
	if (gMemCNF_sz == 0) {
		assert(cnf == NULL);
		assert(gMemCNF == NULL);
		gMemCNF_sz = sizeof(CNF) + (size_t)clsCap * sizeof(SCLAUSE) + litsCap * sizeof(uint32);
		assert(gMemCNF_sz > 0);
		CHECK(cudaMallocManaged((void**)&gMemCNF, gMemCNF_sz));
		tmp_cnf = (CNF*)gMemCNF;
		set_cnf_ptrs << <1, 1 >> > (tmp_cnf, clsCap, litsCap);
		LOGERR("CNF pointers assignment failed");
		CHECK(cudaDeviceSynchronize());
	}
	else {
		assert(cnf != NULL);
		assert(gMemCNF != NULL);
		cap -= gMemCNF_sz;
		gMemCNF_sz = sizeof(CNF) + (size_t)clsCap * sizeof(SCLAUSE) + litsCap * sizeof(uint32);
		assert(gMemCNF_sz > 0);
		addr_t gMemCNF_tmp = NULL;
		CHECK(cudaMallocManaged((void**)&gMemCNF_tmp, gMemCNF_sz));
		tmp_cnf = (CNF*)gMemCNF_tmp;
		if (_pc) {
			copy_cls_k << <1, 1 >> > (tmp_cnf, cnf, clsCap, litsCap);
			int nBlocks = maxGPUThreads / BLOCK1D;
			copy_cnf_k << <nBlocks, BLOCK1D >> > (tmp_cnf, cnf);
		}
		else copy_cnf_k << <1, 1 >> > (tmp_cnf, cnf, clsCap, litsCap);
		LOGERR("Moving CNF failed");
		CHECK(cudaDeviceSynchronize());
		CHECK(cudaFree(gMemCNF));
		gMemCNF = gMemCNF_tmp;
	}
	cap += gMemCNF_sz;
	return tmp_cnf;
}
OT* cuMM::resizeOT(uint32* hist, const uint32& maxVars, const int64& maxEntries) {
	assert(hist != NULL);
	OT* tmp_ot = NULL;
	if (gMemOT_sz != 0) {
		assert(gMemOT != NULL);
		cap -= gMemOT_sz;
		CHECK(cudaFree(gMemOT));
		gMemOT = NULL;
	}
	assert(gMemOT == NULL);
	size_t tabSize = V2D(maxVars + 1ULL);
	gMemOT_sz = sizeof(OT) + tabSize * sizeof(OL) + maxEntries * sizeof(uint32);
	CHECK(cudaMallocManaged((void**)&gMemOT, gMemOT_sz));
	tmp_ot = (OT*)gMemOT;
	ot_ptrs_k << <1, 1 >> > (tmp_ot, hist, tabSize, maxEntries, maxVars);
	LOGERR("Assigning OT pointers failed");
	CHECK(cudaDeviceSynchronize());
	cap += gMemOT_sz;
	return tmp_ot;
}