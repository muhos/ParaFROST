
#include "pfsolve.h"
#include "pfgrid.cuh"
#include "pfatomics.cuh"

_PFROST_D_ uint32* OT::data(const uint32& offset) { return occurs + atomicAdd(&maxEntries, offset); }
_PFROST_D_ void CNF::newClause(SCLAUSE& src) {
	assert(src.size() > 0);
	S_REF newC = cls + atomicAggInc(&n_cls);
	newC->set_ptr(lits + atomicAdd(&n_lits, src.size()));
	newC->copyFrom(src);
}
__global__ void resetOTCap_k(OT* ot) { ot->resetCap(); }
__global__ void assignOTPtrs_k(OT* ot, uint64 tabSize) { ot->assignPtrs(tabSize); }
__global__ void assignListPtrs_k(OT* ot, uint32* hist, uint32 size)
{
	uint32 v = global_tx();
	while (v < size) {
		uint32 p = V2D(v + 1), n = NEG(p);
		uint32 ps = hist[p], ns = hist[n];
		uint32* head = ot->data(ps + ns);
		(*ot)[p].assignPtr(head, ps);
		(*ot)[n].assignPtr(head + ps, ns);
		v += stride_x();
	}
}
__global__ void shrink_k(CNF* dest, CNF* src)
{
	uint32 i = global_tx();
	while (i < src->size()) {
		SCLAUSE& s = (*src)[i];
		if (s.status() == LEARNT || s.status() == ORIGINAL) dest->newClause(s);
		i += stride_x();
	}
}
__global__ void assignCNFPtrs(CNF* cnf, uint32 clsCap, uint64 litsCap) { cnf->assignPtrs(clsCap, litsCap); }
//========================//
//	GPU memory management //
//========================//
bool cuMM::allocPV(PV* pv) {
	assert(gMemPVars == NULL);
	assert(gMemVSt == NULL);
	assert(pv->pVars == NULL);
	assert(pv->units == NULL);
	assert(pv->scores == NULL);
	assert(pv->gsts == NULL);
	size_t gsts_sz = sizeof(GSTATS);
	size_t cuVec_sz = sizeof(cuVecU);
	size_t uintVec_sz = nOrgVars() * sizeof(uint32);
	size_t scVec_sz = nOrgVars() * sizeof(SCORE);
	gMemPVars_sz = ((cuVec_sz + uintVec_sz) << 1ULL) + scVec_sz;
	gMemVSt_sz = gsts_sz + nOrgVars() + 1;
	assert(gMemPVars_sz > 0);
	assert(gMemVSt_sz > 0);
	cap += gMemPVars_sz + gMemVSt_sz;
	if (!hasFreeMem(__func__)) return false;
	CHECK(cudaMallocManaged((void**)&gMemVSt, gMemVSt_sz));
	CHECK(cudaMallocManaged((void**)&gMemPVars, gMemPVars_sz));
	pv->gsts = (GSTATS*)gMemVSt, numDelVars = &pv->gsts->numDelVars;
	seen = gMemVSt + gsts_sz, pv->gsts->seen = seen;
	addr_t ea = gMemPVars, end = ea + gMemPVars_sz;
	pv->pVars = (cuVecU*)ea, ea += cuVec_sz;
	pv->units = (cuVecU*)ea, ea += cuVec_sz;
	pv->pVars->assignPtr((uint32*)ea, nOrgVars()), ea += uintVec_sz, rawUnits = (uint32*)ea;
	pv->units->assignPtr(rawUnits, nOrgVars()), ea += uintVec_sz;
	pv->scores = (SCORE*)ea, ea += scVec_sz;
	assert(ea == end);
	if (devProp.major > 5) {
		CHECK(cudaMemAdvise(gMemVSt, gMemVSt_sz, cudaMemAdviseSetPreferredLocation, MASTER_GPU));
		CHECK(cudaMemAdvise(rawUnits, uintVec_sz + scVec_sz, cudaMemAdviseSetPreferredLocation, MASTER_GPU));
		CHECK(cudaMemPrefetchAsync(rawUnits, uintVec_sz + scVec_sz, MASTER_GPU));
	}
	return true;
}
bool cuMM::resizeCNF(CNF*& cnf, const uint32& clsCap, const uint64& litsCap) {
	assert(clsCap > 0);
	assert(litsCap > 0);
	if (gMemCNF_sz == 0) {
		assert(cnf == NULL);
		assert(gMemCNF == NULL);
		gMemCNF_sz = sizeof(CNF) + (size_t)clsCap * sizeof(SCLAUSE) + litsCap * sizeof(uint32);
		assert(gMemCNF_sz > 0);
		cap += gMemCNF_sz;
		if (!hasFreeMem(__func__)) return false;
		CHECK(cudaMallocManaged((void**)&gMemCNF, gMemCNF_sz));
		cnf = (CNF*)gMemCNF, cnf->assignPtrs(clsCap, litsCap);
		rawCls = *cnf, rawLits = cnf->data();
	}
	else {
		assert(cnf != NULL);
		assert(gMemCNF != NULL);
		size_t new_gMemCNF_sz = sizeof(CNF) + (size_t)clsCap * sizeof(SCLAUSE) + litsCap * sizeof(uint32);
		assert(new_gMemCNF_sz > 0);
		if (gMemCNF_sz >= new_gMemCNF_sz) cnf->shrink();
		else {
			cap -= gMemCNF_sz, cap += new_gMemCNF_sz;
			if (!hasFreeMem(__func__)) return false;
			addr_t gMemCNF_tmp = NULL;
			CHECK(cudaMallocManaged((void**)&gMemCNF_tmp, new_gMemCNF_sz));
			CNF* tmp_cnf = (CNF*)gMemCNF_tmp;
			tmp_cnf->assignPtrs(clsCap, litsCap), rawCls = *tmp_cnf, rawLits = tmp_cnf->data();
			tmp_cnf->copyFrom(cnf);
			CHECK(cudaFree(gMemCNF));
			gMemCNF = gMemCNF_tmp, gMemCNF_sz = new_gMemCNF_sz, cnf = tmp_cnf;
			prefetchCNF();
		}
	}
	return true;
}
bool cuMM::resizeOTAsync(OT*& ot, uint32* d_hist, const int64& entriesCap, const cudaStream_t& _s) {
	assert(d_hist != NULL);
	assert(entriesCap > 0);
	size_t tabSize = nDualVars();
	size_t new_gMemOT_sz = sizeof(OT) + tabSize * sizeof(OL) + entriesCap * sizeof(uint32);
	assert(new_gMemOT_sz > 0);
	if (gMemOT_sz < new_gMemOT_sz) { // realloc
		if (gMemOT_sz != 0) {
			assert(gMemOT != NULL);
			assert(ot != NULL);
			CHECK(cudaFree(gMemOT));
			gMemOT = NULL;
		}
		assert(gMemOT == NULL);
		cap -= gMemOT_sz, cap += new_gMemOT_sz;
		if (!hasFreeMem(__func__)) return false;
		CHECK(cudaMallocManaged((void**)&gMemOT, new_gMemOT_sz));
		if (devProp.major > 5) {
			CHECK(cudaMemAdvise(gMemOT, new_gMemOT_sz, cudaMemAdviseSetPreferredLocation, MASTER_GPU));
			CHECK(cudaMemPrefetchAsync(gMemOT, new_gMemOT_sz, MASTER_GPU));
		}
		ot = (OT*)gMemOT;
		assignOTPtrs_k << <1, 1 >> > (ot, tabSize);
		assignListPtrs_k << <otBlocks, BLOCK1D >> > (ot, d_hist, nOrgVars());
		gMemOT_sz = new_gMemOT_sz;
	}
	else assignListPtrs_k << <otBlocks, BLOCK1D, 0, _s >> > (ot, d_hist, nOrgVars());
	return true;
}
void cuMM::resetOTCapAsync(OT* ot, const cudaStream_t& _s) {
	resetOTCap_k << <1, 1, 0, _s >> > (ot);
}
void cuMM::freeHostCNF() { if (hMemCNF != NULL) free(hMemCNF), hMemCNF = NULL; }
CNF* cuMM::allocHostCNF() {
	assert(gMemCNF != NULL);
	assert(hMemCNF == NULL);
	assert(hMemCNF_sz == 0);
	hMemCNF_sz = gMemCNF_sz;
	assert(hMemCNF_sz > 0);
	hMemCNF = (addr_t)malloc(hMemCNF_sz);
	CHECK(cudaMemcpy(hMemCNF, gMemCNF, sizeof(CNF), cudaMemcpyDeviceToHost));
	CNF* cnf = (CNF*)hMemCNF;
	cnf->assignPtrs();
	return cnf;
}
cuMM::~cuMM() {
	rawLits = NULL, rawUnits = NULL, rawCls = NULL;
	freeHostCNF();
	if (gMemOT != NULL) CHECK(cudaFree(gMemOT)), gMemOT = NULL;
	if (gMemCNF != NULL) CHECK(cudaFree(gMemCNF)), gMemCNF = NULL;
	if (gMemVSt != NULL) CHECK(cudaFree(gMemVSt)), gMemVSt = NULL;
	if (gMemPVars != NULL) CHECK(cudaFree(gMemPVars)), gMemPVars = NULL;
	cap = 0;
}