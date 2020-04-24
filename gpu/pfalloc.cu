
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
		(*ot)[p].alloc(ot->data(idx), ps), idx += ps;
		(*ot)[n].alloc(ot->data(idx), ns), idx += ns;
	}
	assert(idx <= ot->capacity());
}
__global__ void clear_cls(CNF* cnf)
{
	uint32 i = blockDim.x * blockIdx.x + threadIdx.x;
	while (i < cnf->size()) { (*cnf)[i].~SCLAUSE(); i += blockDim.x * gridDim.x; }
}
__global__ void clear_lists(OT* ot)
{
	int64 v = blockDim.x * blockIdx.x + threadIdx.x;
	while (v < ot->size()) { (*ot)[v].~cuVec(); v += blockDim.x * gridDim.x; }
}

bool cuMM::allocStats(GSTATS*& gstats, const uint32& maxVars) {
	assert(gMemStats == NULL);
	assert(gstats == NULL);
	gMemStats_sz = sizeof(GSTATS) + maxVars * sizeof(Byte);
	assert(gMemStats_sz > 0);
	cap += gMemStats_sz;
	// memory query
	size_t _free = 0, _tot = 0;
	CHECK(cudaMemGetInfo(&_free, &_tot));
	int64 _used = (_tot - _free) + cap;
	if (verb > 1) printf("c | Used/total GPU memory after %s call = %lld/%zd MB\n", __func__, _used / MBYTE, _tot / MBYTE);
	if (_used >= int64(_tot)) {
		printf("c | WARNING - not enough GPU memory for %s (used/total = %lld/%zd MB), simplifications will be terminated\n", __func__, _used / MBYTE, _tot / MBYTE);
		return false;
	}
	CHECK(cudaMallocManaged((void**)&gMemStats, gMemStats_sz));
	gstats = (GSTATS*)gMemStats;
	gstats->seen = gMemStats + sizeof(GSTATS);
	return true;
}
bool cuMM::allocPV(PV* pv, const uint32& maxVars) {
	assert(gMemPV == NULL);
	assert(gMemSol == NULL);
	assert(pv->pVars == NULL);
	assert(pv->sol == NULL);
	size_t sol_sz = sizeof(GSOL);
	size_t cuVec_sz = sizeof(cuVecU);
	size_t uintVec_sz = cuVec_sz + (size_t)maxVars * sizeof(uint32);
	size_t litStVec_sz = maxVars * sizeof(LIT_ST);
	gMemPV_sz = uintVec_sz, gMemSol_sz = sol_sz + uintVec_sz + litStVec_sz;
	assert(gMemPV_sz > 0);
	assert(gMemSol_sz > 0);
	cap += gMemPV_sz + gMemSol_sz;
	// memory query
	size_t _free = 0, _tot = 0;
	CHECK(cudaMemGetInfo(&_free, &_tot));
	int64 _used = (_tot - _free) + cap;
	if (verb > 1) printf("c | Used/total GPU memory after %s call = %lld/%zd MB\n", __func__, _used / MBYTE, _tot / MBYTE);
	if (_used >= int64(_tot)) {
		printf("c | WARNING - not enough GPU memory for %s (used/total = %lld/%zd MB), simplifications will be terminated\n", __func__, _used / MBYTE, _tot / MBYTE);
		return false;
	}
	CHECK(cudaMallocManaged((void**)&gMemPV, gMemPV_sz));
	CHECK(cudaMallocManaged((void**)&gMemSol, gMemSol_sz));
	pv->pVars = (cuVecU*)gMemPV;
	pv->pVars->alloc((uint32*)(gMemPV + cuVec_sz), maxVars);
	pv->sol = (GSOL*)gMemSol;
	pv->sol->head = 0;
	pv->sol->assigns = (cuVecU*)(gMemSol + sol_sz);
	pv->sol->assigns->alloc((uint32*)(gMemSol + sol_sz + cuVec_sz), maxVars);
	pv->sol->value = (LIT_ST*)(gMemSol + sol_sz + uintVec_sz);
	mem_set(pv->sol->value, LIT_ST(UNDEFINED), maxVars);
	return true;
}
bool cuMM::allocVO(OCCUR*& occurs, SCORE*& scores, const uint32& maxVars) {
	assert(gMemVOrd == NULL);
	assert(occurs == NULL);
	assert(scores == NULL);
	gMemVOrd_sz = (size_t)maxVars * (sizeof(OCCUR) + sizeof(SCORE));
	assert(gMemVOrd_sz > 0);
	cap += gMemVOrd_sz;
	// memory query
	size_t _free = 0, _tot = 0;
	CHECK(cudaMemGetInfo(&_free, &_tot));
	int64 _used = (_tot - _free) + cap;
	if (verb > 1) printf("c | Used/total GPU memory after %s call = %lld/%zd MB\n", __func__, _used / MBYTE, _tot / MBYTE);
	if (_used >= int64(_tot)) {
		printf("c | WARNING - not enough GPU memory for %s (used/total = %lld/%zd MB), simplifications will be terminated\n", __func__, _used / MBYTE, _tot / MBYTE);
		return false;
	}
	CHECK(cudaMallocManaged((void**)&gMemVOrd, gMemVOrd_sz));
	occurs = (OCCUR*)gMemVOrd;
	scores = (SCORE*)(occurs + maxVars);
	assert(addr_t(scores + maxVars) == gMemVOrd + gMemVOrd_sz);
	return true;
}
bool cuMM::resizeCNF(CNF*& cnf, const uint32& clsCap, const uint64& litsCap, const int& phase) {
	CNF* tmp_cnf = NULL;
	if (gMemCNF_sz == 0) {
		assert(cnf == NULL);
		assert(gMemCNF == NULL);
		assert(phase == 0);
		gMemCNF_sz = sizeof(CNF) + (size_t)clsCap * sizeof(SCLAUSE) + litsCap * sizeof(uint32);
		assert(gMemCNF_sz > 0);
		cap += gMemCNF_sz;
		// memory query
		size_t _free = 0, _tot = 0;
		CHECK(cudaMemGetInfo(&_free, &_tot));
		int64 _used = (_tot - _free) + cap;
		if (verb > 1) printf("c | Used/total GPU memory after %s call = %lld/%zd MB\n", __func__, _used / MBYTE, _tot / MBYTE);
		if (_used >= int64(_tot)) {
			printf("c | WARNING - not enough GPU memory for %s (used/total = %lld/%zd MB), simplifications will be terminated\n", __func__, _used / MBYTE, _tot / MBYTE);
			return false;
		}
		CHECK(cudaMallocManaged((void**)&gMemCNF, gMemCNF_sz));
		tmp_cnf = (CNF*)gMemCNF;
		set_cnf_ptrs << <1, 1 >> > (tmp_cnf, clsCap, litsCap);
		LOGERR("CNF pointers assignment failed");
		CHECK(cudaDeviceSynchronize());
	}
	else {
		assert(cnf != NULL);
		assert(gMemCNF != NULL);
		assert(phase >= 0);
		cap -= gMemCNF_sz;
		gMemCNF_sz = sizeof(CNF) + (size_t)clsCap * sizeof(SCLAUSE) + litsCap * sizeof(uint32);
		assert(gMemCNF_sz > 0);
		cap += gMemCNF_sz;
		// memory query
		size_t _free = 0, _tot = 0;
		CHECK(cudaMemGetInfo(&_free, &_tot));
		int64 _used = (_tot - _free) + cap;
		if (verb > 1) printf("c | Used/total GPU memory after %s call = %lld/%zd MB\n", __func__, _used / MBYTE, _tot / MBYTE);
		if (_used >= int64(_tot)) {
			printf("c | WARNING - not enough GPU memory for %s (used/total = %lld/%zd MB), simplifications will be terminated\n", __func__, _used / MBYTE, _tot / MBYTE);
			return false;
		}
		addr_t gMemCNF_tmp = NULL;
		CHECK(cudaMallocManaged((void**)&gMemCNF_tmp, gMemCNF_sz));
		tmp_cnf = (CNF*)gMemCNF_tmp;
		if (phase == 0) {
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
	cnf = tmp_cnf;
	return true;
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
	assert(gMemOT_sz > 0);
	cap += gMemOT_sz;
	// memory query
	size_t _free = 0, _tot = 0;
	CHECK(cudaMemGetInfo(&_free, &_tot));
	int64 _used = (_tot - _free) + cap;
	if (verb > 1) printf("c | Used/total GPU memory after %s call = %lld/%zd MB\n", __func__, _used / MBYTE, _tot / MBYTE);
	if (_used >= int64(_tot)) {
		printf("c | WARNING - not enough GPU memory for %s (used/total = %lld/%zd MB), simplifications will be terminated\n", __func__, _used / MBYTE, _tot / MBYTE);
		return NULL;
	}
	CHECK(cudaMallocManaged((void**)&gMemOT, gMemOT_sz));
	tmp_ot = (OT*)gMemOT;
	ot_ptrs_k << <1, 1 >> > (tmp_ot, hist, tabSize, maxEntries, maxVars);
	LOGERR("Assigning OT pointers failed");
	CHECK(cudaDeviceSynchronize());
	return tmp_ot;
}
cuMM::~cuMM() {
	if (gMemOT != NULL) CHECK(cudaFree(gMemOT)), gMemOT = NULL;
	if (gMemPV != NULL) CHECK(cudaFree(gMemPV)), gMemPV = NULL;
	if (gMemCNF != NULL) CHECK(cudaFree(gMemCNF)), gMemCNF = NULL;
	if (gMemSol != NULL) CHECK(cudaFree(gMemSol)), gMemSol = NULL;
	if (gMemVOrd != NULL) CHECK(cudaFree(gMemVOrd)), gMemVOrd = NULL;
	if (gMemStats != NULL) CHECK(cudaFree(gMemStats)), gMemStats = NULL;
	cap = 0;
}
CNF::~CNF() {
	int nBlocks = MIN((n_cls + BLOCK1D - 1) / BLOCK1D, maxGPUThreads / BLOCK1D);
	clear_cls << <nBlocks, BLOCK1D >> > (this);
	LOGERR("Clearing clauses failed");
	CHECK(cudaDeviceSynchronize());
	cls = NULL, lits = NULL;
}
OT::~OT() {
	int nBlocks = MIN((V2D(nOrgVars() + 1) + BLOCK1D - 1) / BLOCK1D, maxGPUThreads / BLOCK1D);
	clear_lists << <nBlocks, BLOCK1D >> > (this);
	LOGERR("Clearing occurrence lists failed");
	CHECK(cudaDeviceSynchronize());
	lists = NULL, occurs = NULL;
}