
#include "pfdevice.cuh"

template<class T>
__global__ void memset_k(T* mem, T val, size_t size)
{
	size_t tid = global_tx();
	while (tid < size) { mem[tid] = val; tid += stride_x(); }
}

__global__ void sts_reset(GSTATS* gsts, uint32 size)
{
	size_t tid = global_tx();
	while (tid < size) { 
		if (tid == 0) gsts->numDelVars = 0, gsts->numClauses = 0, gsts->numLits = 0;
		gsts->seen[tid] = 0;
		tid += stride_x();
	}
}

__global__ void reset_ot_k(OT* ot)
{
	uint64 tid = global_tx();
	while (tid < ot->size()) { (*ot)[tid].clear(); tid += stride_x(); }
}

__global__ void reduce_ot(CNF* cnf, OT* ot)
{
	uint64 tid = global_tx();
	while (tid < ot->size()) { reduceOL(*cnf, (*ot)[tid]); tid += stride_x(); }
}

__global__ void reduce_ot_p(CNF* cnf, OT* ot, cuVecU* pVars)
{
	uint32 tid = global_tx();
	while (tid < pVars->size()) { 
		assert(pVars->at(tid)); 
		reduceOL(*cnf, (*ot)[V2D(pVars->at(tid))]); 
		tid += stride_x(); 
	}
}

__global__ void reduce_ot_n(CNF* cnf, OT* ot, cuVecU* pVars)
{
	uint32 tid = global_tx();
	while (tid < pVars->size()) { 
		assert(pVars->at(tid));
		reduceOL(*cnf, (*ot)[NEG(V2D(pVars->at(tid)))]);
		tid += stride_x(); 
	}
}

__global__ void create_ot_k(CNF* cnf, OT* ot)
{
	uint32 tid = global_tx();
	while (tid < cnf->size()) {
		SCLAUSE& c = (*cnf)[tid];
		if (c.status() == ORIGINAL || c.status() == LEARNT) {
#pragma unroll
			for (LIT_POS k = 0; k < c.size(); k++) (*ot)[c[k]].insert(tid);
		}
		tid += stride_x();
	}
}

__global__ void assign_scores(SCORE* scores, uint32* hist, uint32 size)
{
	uint32 tid = global_tx();
	while (tid < size) {
		uint32 p = V2D(tid + 1), ps = hist[p], ns = hist[NEG(p)];
		scores[tid].v = tid + 1;
		scores[tid].sc = (!ps || !ns) ? (ps | ns) : (ps * ns);
		tid += stride_x();
	}
}

__global__ void assign_scores(SCORE* scores, uint32* hist, OT* ot, uint32 size)
{
	uint32 tid = global_tx();
	while (tid < size) {
		uint32 p = V2D(tid + 1), n = NEG(p), ps = (*ot)[p].size(), ns = (*ot)[n].size();
		hist[p] = ps, hist[n] = ns;
		scores[tid].v = tid + 1;
		scores[tid].sc = (!ps || !ns) ? (ps | ns) : (ps * ns);
		tid += stride_x();
	}
}

__global__ void calc_sig_k(CNF* cnf, uint32 offset, uint32 size)
{
	uint32 tid = global_tx() + offset;
	while (tid < size) { calcSig((*cnf)[tid]); tid += stride_x(); }
}

__global__ void copy_k(uint32* dest, CNF* src, int64 size)
{
	int64 tid = global_tx();
	while (tid < size) { dest[tid] = *src->data(tid); tid += stride_x(); }
}

__global__ void copy_k(CNF* dest, CNF* src, int64 size)
{
	int64 tid = global_tx();
	while (tid < size) { *dest->data(tid) = *src->data(tid); tid += stride_x(); }
}

__global__ void copy_if_k(uint32* dest, CNF* src, GSTATS* gsts)
{
	uint32 tid = global_tx();
	while (tid < src->size()) {
		SCLAUSE& c = (*src)[tid];
		if (c.status() == ORIGINAL || c.status() == LEARNT) {
			uint64 lits_idx = atomicAdd(&gsts->numLits, c.size());
#pragma unroll
			for (LIT_POS k = 0; k < c.size(); k++) dest[lits_idx++] = c[k];
		}
		tid += stride_x();
	}
}

__global__ void calc_added_cls_k(CNF* cnf, OT* ot, cuVecU* pVars, GSTATS* gsts)
{
	uint32* sh_rCls = SharedMemory<uint32>();
	uint32 tid = global_tx_off();
	uint32 x, p, n, nCls = 0;
	while (tid < pVars->size()) {
		x = (*pVars)[tid];
		assert(x);
		p = V2D(x), n = NEG(p);
		calcResolvents(x, *cnf, (*ot)[p], (*ot)[n], nCls);
		if (tid + blockDim.x < pVars->size()) {
			x = (*pVars)[tid + blockDim.x];
			assert(x);
			p = V2D(x), n = NEG(p);
			calcResolvents(x, *cnf, (*ot)[p], (*ot)[n], nCls);
		}
		tid += stride_x_off();
	}
	loadShared(sh_rCls, nCls, pVars->size());
	sharedReduce(sh_rCls, nCls);
	warpReduce(sh_rCls, nCls);
	if (threadIdx.x == 0) atomicAdd(&gsts->numClauses, nCls);
}

__global__ void calc_added_all_k(CNF* cnf, OT* ot, cuVecU* pVars, GSTATS* gsts)
{
	uint32* sh_rCls = SharedMemory<uint32>();
	uint64* sh_rLits = (uint64*)(sh_rCls + blockDim.x);
	uint32 tid = global_tx_off();
	uint32 x, p, n, nCls = 0;
	uint64 nLits = 0;
	while (tid < pVars->size()) {
		x = (*pVars)[tid];
		assert(x);
		p = V2D(x), n = NEG(p);
		calcResolvents(x, *cnf, (*ot)[p], (*ot)[n], nCls, nLits);
		if (tid + blockDim.x < pVars->size()) {
			x = (*pVars)[tid + blockDim.x];
			assert(x);
			p = V2D(x), n = NEG(p);
			calcResolvents(x, *cnf, (*ot)[p], (*ot)[n], nCls, nLits);
		}
		tid += stride_x_off();
	}
	loadShared(sh_rCls, nCls, sh_rLits, nLits, pVars->size());
	sharedReduce(sh_rCls, nCls, sh_rLits, nLits);
	warpReduce(sh_rCls, nCls, sh_rLits, nLits);
	if (threadIdx.x == 0) {
		atomicAdd(&gsts->numClauses, nCls);
		atomicAdd(&gsts->numLits, nLits);
	}
}

__global__ void cnt_del_vars(GSTATS* gsts, uint32 size)
{
	uint32* sh_delVars = SharedMemory<uint32>();
	uint32 tid = global_tx_off();
	uint32 nVarsDeleted = 0;
	while (tid < size) {
		if (!gsts->seen[tid + 1]) nVarsDeleted++;
		uint32 off = tid + blockDim.x + 1;
		if (off < size && !gsts->seen[off]) nVarsDeleted++;
		tid += stride_x_off();
	}
	loadShared(sh_delVars, nVarsDeleted, size);
	sharedReduce(sh_delVars, nVarsDeleted);
	warpReduce(sh_delVars, nVarsDeleted);
	if (threadIdx.x == 0) atomicAdd(&gsts->numDelVars, nVarsDeleted);
}

__global__ void mark_vars(CNF* cnf, GSTATS* gsts)
{
	uint32 tid = global_tx();
	while (tid < cnf->size()) {
		SCLAUSE& c = (*cnf)[tid];
		if (c.status() == ORIGINAL || c.status() == LEARNT) {
#pragma unroll
			for (LIT_POS k = 0; k < c.size(); k++) { assert(c[k]); gsts->seen[ABS(c[k])] = 1; }
		}
		tid += stride_x();
	}
}

__global__ void cnt_reds(CNF* cnf, GSTATS* gsts)
{
	uint32* sh_rCls = SharedMemory<uint32>();
	uint64* sh_rLits = (uint64*)(sh_rCls + blockDim.x);
	uint32 tid = global_tx_off();
	uint32 nCls = 0;
	uint64 nLits = 0;
	while (tid < cnf->size()) {
		SCLAUSE& c1 = (*cnf)[tid];
		if (c1.status() == LEARNT || c1.status() == ORIGINAL) {
			CL_LEN cl_size = c1.size();
			nCls++, nLits += cl_size;
#pragma unroll
			for (LIT_POS k = 0; k < cl_size; k++) { assert(c1[k]); gsts->seen[ABS(c1[k])] = 1; }
		}
		if (tid + blockDim.x < cnf->size()) {
			SCLAUSE& c2 = (*cnf)[tid + blockDim.x];
			if (c2.status() == LEARNT || c2.status() == ORIGINAL) {
				CL_LEN cl_size = c2.size();
				nCls++, nLits += cl_size;
#pragma unroll
				for (LIT_POS k = 0; k < cl_size; k++) { assert(c2[k]); gsts->seen[ABS(c2[k])] = 1; }
			}
		}
		tid += stride_x_off();
	}
	loadShared(sh_rCls, nCls, sh_rLits, nLits, cnf->size());
	sharedReduce(sh_rCls, nCls, sh_rLits, nLits);
	warpReduce(sh_rCls, nCls, sh_rLits, nLits);
	if (threadIdx.x == 0) {
		atomicAdd(&gsts->numClauses, nCls);
		atomicAdd(&gsts->numLits, nLits);
	}
}

__global__ void cnt_lits(CNF* cnf, GSTATS* gsts)
{
	uint64* sh_rLits = SharedMemory<uint64>();
	uint32 tid = global_tx_off();
	uint64 nLits = 0;
	while (tid < cnf->size()) {
		SCLAUSE& c1 = (*cnf)[tid];
		if (c1.status() == LEARNT || c1.status() == ORIGINAL) nLits += c1.size();
		if (tid + blockDim.x < cnf->size()) {
			SCLAUSE& c2 = (*cnf)[tid + blockDim.x];
			if (c2.status() == LEARNT || c2.status() == ORIGINAL) nLits += c2.size();
		}
		tid += stride_x_off();
	}
	loadShared(sh_rLits, nLits, cnf->size());
	sharedReduce(sh_rLits, nLits);
	warpReduce(sh_rLits, nLits);
	if (threadIdx.x == 0) atomicAdd(&gsts->numLits, nLits);
}

__global__ void cnt_cls_lits(CNF* cnf, GSTATS* gsts)
{
	uint32* sh_rCls = SharedMemory<uint32>();
	uint64* sh_rLits = (uint64*)(sh_rCls + blockDim.x);
	uint32 tid = global_tx_off();
	uint32 nCls = 0;
	uint64 nLits = 0;
	while (tid < cnf->size()) {
		SCLAUSE& c1 = (*cnf)[tid];
		if (c1.status() == LEARNT || c1.status() == ORIGINAL) nCls++, nLits += c1.size();
		if (tid + blockDim.x < cnf->size()) {
			SCLAUSE& c2 = (*cnf)[tid + blockDim.x];
			if (c2.status() == LEARNT || c2.status() == ORIGINAL) nCls++, nLits += c2.size();
		}
		tid += stride_x_off();
	}
	loadShared(sh_rCls, nCls, sh_rLits, nLits, cnf->size());
	sharedReduce(sh_rCls, nCls, sh_rLits, nLits);
	warpReduce(sh_rCls, nCls, sh_rLits, nLits);
	if (threadIdx.x == 0) {
		atomicAdd(&gsts->numClauses, nCls);
		atomicAdd(&gsts->numLits, nLits);
	}
}

__global__ void ve_k(CNF *cnf, OT* ot, cuVecU* pVars, cuVecU* units)
{
	uint32 tx = threadIdx.x;
	uint32 tid = global_tx();
	__shared__ uint32 defs[BLVE * FAN_LMT];
	__shared__ uint32 outs[BLVE * SH_MAX_BVE_OUT];
	while (tid < pVars->size()) {
		assert((*pVars)[tid]);
		uint32 x = (*pVars)[tid], p = V2D(x), n = NEG(p);
		if ((*ot)[p].size() == 0 || (*ot)[n].size() == 0) { // pure
			deleteClauses(*cnf, (*ot)[p], (*ot)[n]);
			(*pVars)[tid] = 0, (*ot)[p].clear(true), (*ot)[n].clear(true);
		}
		else if ((*ot)[p].size() == 1 || (*ot)[n].size() == 1) {
			if (resolve_x(x, *cnf, (*ot)[p], (*ot)[n], units, &outs[tx * SH_MAX_BVE_OUT]))
				(*pVars)[tid] = 0, (*ot)[p].clear(true), (*ot)[n].clear(true);
		}
		else if (gateReasoning_x(p, *cnf, (*ot)[p], (*ot)[n], units, &defs[tx * FAN_LMT], &outs[tx * SH_MAX_BVE_OUT])
			  || resolve_x(x, *cnf, (*ot)[p], (*ot)[n], units, &outs[tx * SH_MAX_BVE_OUT])) {
			(*pVars)[tid] = 0, (*ot)[p].clear(true), (*ot)[n].clear(true);
		}
		tid += stride_x();
	}
}

__global__ void hse_k(CNF *cnf, OT* ot, cuVecU* pVars, cuVecU* units)
{
	uint32 tid = global_tx();
	__shared__ uint32 sh_cls[BLHSE * SH_MAX_HSE_IN];
	while (tid < pVars->size()) {
		assert(pVars->at(tid));
		uint32 p = V2D(pVars->at(tid)), n = NEG(p);
		self_sub_x(p, *cnf, (*ot)[p], (*ot)[n], units, &sh_cls[threadIdx.x * SH_MAX_HSE_IN]);
		tid += stride_x();
	}
}

__global__ void bce_k(CNF *cnf, OT* ot, cuVecU* pVars)
{
	uint32 tid = global_tx();
	__shared__ uint32 sh_cls[BLBCE * SH_MAX_BCE_IN];
	while (tid < pVars->size()) {
		assert((*pVars)[tid]);
		uint32 x = (*pVars)[tid], p = V2D(x), n = NEG(p);
		blocked_x(x, *cnf, (*ot)[p], (*ot)[n], &sh_cls[threadIdx.x * SH_MAX_BCE_IN]);
		tid += stride_x();
	}
}

__global__ void hre_k(CNF *cnf, OT* ot, cuVecU* pVars)
{
	uint32 gid = global_ty();
	uint32* smem = SharedMemory<uint32>();
	uint32* m_c = smem + warpSize * SH_MAX_HRE_IN + threadIdx.y * SH_MAX_HRE_OUT; // shared memory for resolvent
	while (gid < pVars->size()) {
		assert(pVars->at(gid));
		uint32 p = V2D(pVars->at(gid));
		OL& poss = (*ot)[p];
		// do merging and apply forward equality check (on-the-fly) over resolvents
#pragma unroll
		for (uint32 i = 0; i < poss.size(); i++) {
			if ((*cnf)[poss[i]].status() == DELETED) continue;
			uint32 pos_size = 0;
			if (threadIdx.x == 0) pos_size = (*cnf)[poss[i]].size();
			pos_size = __shfl_sync(FULLWARP, pos_size, 0);
			assert(pos_size == (*cnf)[poss[i]].size());
			if (pos_size <= SH_MAX_HRE_IN) { // use shared memory for positives
				uint32* sh_pos = smem + threadIdx.y * SH_MAX_HRE_IN;
				if (threadIdx.x == 0) (*cnf)[poss[i]].shareTo(sh_pos);
				OL& negs = (*ot)[NEG(p)];
#pragma unroll
				for (uint32 j = 0; j < negs.size(); j++) {
					SCLAUSE& neg = (*cnf)[negs[j]];
					if (neg.status() == DELETED || (pos_size + neg.size() - 2) > SH_MAX_HRE_OUT) continue;
					CL_LEN m_c_size = 0;
					uint32 m_c_sig = 0;
					if (threadIdx.x == 0) {
						assert(warpSize == blockDim.x);
						m_c_size = merge(pVars->at(gid), sh_pos, pos_size, neg, m_c);
						calcSig(m_c, m_c_size, m_c_sig);
					}
					m_c_size = __shfl_sync(FULLWARP, m_c_size, 0);
					m_c_sig = __shfl_sync(FULLWARP, m_c_sig, 0);
					forward_equ(*cnf, *ot, m_c, m_c_sig, m_c_size);
				}
			}
			else { // use global memory
				OL& negs = (*ot)[NEG(p)];
#pragma unroll
				for (uint32 j = 0; j < negs.size(); j++) {
					SCLAUSE& neg = (*cnf)[negs[j]];
					if (neg.status() == DELETED || (pos_size + neg.size() - 2) > SH_MAX_HRE_OUT) continue;
					CL_LEN m_c_size = 0;
					uint32 m_c_sig = 0;
					if (threadIdx.x == 0) {
						assert(warpSize == blockDim.x);
						m_c_size = merge(pVars->at(gid), (*cnf)[poss[i]], neg, m_c);
						calcSig(m_c, m_c_size, m_c_sig);
					}
					m_c_size = __shfl_sync(FULLWARP, m_c_size, 0);
					m_c_sig = __shfl_sync(FULLWARP, m_c_sig, 0);
					forward_equ(*cnf, *ot, m_c, m_c_sig, m_c_size);
				} 
			} 
		} 
		gid += stride_y();
	}
}
//==============================================//
//          ParaFROST Wrappers/helpers          //
//==============================================//
void cuMemSetAsync(addr_t mem, const Byte& val, const size_t& size)
{
	uint32 nBlocks = MIN(uint32((size + BLOCK1D - 1) / BLOCK1D), maxGPUThreads / BLOCK1D);
	memset_k<Byte> << <nBlocks, BLOCK1D >> > (mem, val, size);
}
void stsResetAsync(GSTATS* gsts)
{
	uint32 nBlocks = MIN((nOrgVars() + BLOCK1D) / BLOCK1D, maxGPUThreads / BLOCK1D);
	sts_reset << <nBlocks, BLOCK1D >> > (gsts, nOrgVars() + 1);
}
void copy(CNF* dest, CNF* src, const int64& size)
{
	uint32 nBlocks = MIN(uint32((size + BLOCK1D - 1) / BLOCK1D), maxGPUThreads / BLOCK1D);
	copy_k << <nBlocks, BLOCK1D >> > (dest, src, size);
	LOGERR("Copying failed");
	CHECK(cudaDeviceSynchronize());
}
void copy(uint32* dest, CNF* src, const int64& size)
{
	uint32 nBlocks = MIN(uint32((size + BLOCK1D - 1) / BLOCK1D), maxGPUThreads / BLOCK1D);
	copy_k << <nBlocks, BLOCK1D >> > (dest, src, size);
	LOGERR("Copying failed");
	CHECK(cudaDeviceSynchronize());
}
void copyIf(uint32* dest, CNF* src, GSTATS* gsts)
{
	gsts->numLits = 0;
	uint32 nBlocks = MIN((nClauses() + BLOCK1D - 1) / BLOCK1D, maxGPUThreads / BLOCK1D);
	copy_if_k << <nBlocks, BLOCK1D >> > (dest, src, gsts);
	LOGERR("Copying failed");
	CHECK(cudaDeviceSynchronize());
}
void calcVarScores(SCORE* scores, uint32* hist)
{
	uint32 nBlocks = MIN((nOrgVars() + BLOCK1D - 1) / BLOCK1D, maxGPUThreads / BLOCK1D);
	assign_scores << <nBlocks, BLOCK1D >> > (scores, hist, nOrgVars());
	LOGERR("Assigning scores failed");
	CHECK(cudaDeviceSynchronize());
}
void calcVarScores(SCORE* scores, uint32* hist, OT* ot)
{
	uint32 nBlocks = MIN((nOrgVars() + BLOCK1D - 1) / BLOCK1D, maxGPUThreads / BLOCK1D);
	assign_scores << <nBlocks, BLOCK1D >> > (scores, hist, ot, nOrgVars());
	LOGERR("Assigning scores failed");
	CHECK(cudaDeviceSynchronize());
}
void calcAdded(CNF* cnf, OT* ot, PV* pv)
{
	assert(pv->numPVs > 0);
	pv->gsts->numClauses = 0;
	uint32 nBlocks = MIN((pv->numPVs + (BLUB << 1) - 1) / (BLUB << 1), maxGPUThreads / (BLUB << 1));
	uint32 smemSize = BLUB * sizeof(uint32);
	calc_added_cls_k << <nBlocks, BLUB, smemSize >> > (cnf, ot, pv->pVars, pv->gsts);
	LOGERR("Added clauses calculation failed");
	CHECK(cudaDeviceSynchronize());
	cnf_stats.max_added_cls = pv->gsts->numClauses;
}
void calcSigAsync(CNF* cnf, const uint32& offset, const uint32& size, const cudaStream_t& _s)
{
	uint32 nBlocks = MIN((size + BLOCK1D - 1) / BLOCK1D, maxGPUThreads / BLOCK1D);
	calc_sig_k << <nBlocks, BLOCK1D, 0, _s >> > (cnf, offset, size);
}
void reduceOTAsync(CNF* cnf, OT* ot, const bool& p)
{
	assert(cnf != NULL);
	assert(ot != NULL);
	uint32 nBlocks = MIN(uint32((nDualVars() + BLOCK1D - 1) / BLOCK1D), maxGPUThreads / BLOCK1D);
	reduce_ot << <nBlocks, BLOCK1D >> > (cnf, ot);
	if (p) {
		LOGERR("Occurrence table reduction failed");
		CHECK(cudaDeviceSynchronize());
		cout << "c |==========================|" << endl;
		cout << "c |==== occurrence table ====|" << endl;
		ot->print();
		cout << "c |==========================|" << endl;
	}
}
void reduceOT(CNF* cnf, OT* ot, PV* pv, cudaStream_t* streams, const bool& p)
{
	assert(pv->numPVs);
	assert(cnf != NULL);
	assert(ot != NULL);
	uint32 nBlocks = MIN((pv->numPVs + BLOCK1D - 1) / BLOCK1D, maxGPUThreads / BLOCK1D);
	reduce_ot_p << <nBlocks, BLOCK1D, 0, streams[0] >> > (cnf, ot, pv->pVars);
	reduce_ot_n << <nBlocks, BLOCK1D, 0, streams[1] >> > (cnf, ot, pv->pVars);
	LOGERR("Occurrence table reduction failed");
	CHECK(cudaDeviceSynchronize());
	assert(ot->accViolation());
	if (p) {
		cout << "c |==========================|" << endl;
		cout << "c |==== occurrence table ====|" << endl;
		ot->print();
		cout << "c |==========================|" << endl;
	}
}
void createOT(CNF* cnf, OT* ot, const bool& p)
{
	assert(cnf != NULL);
	assert(ot != NULL);
	uint32 cnf_size = nClauses() + (maxAddedCls() >> 1);
	uint32 rstGridSize = MIN(uint32((nDualVars() + BLOCK1D - 1) / BLOCK1D), maxGPUThreads / BLOCK1D);
	uint32 otGridSize = MIN((cnf_size + BLOCK1D - 1) / BLOCK1D, maxGPUThreads / BLOCK1D);
	reset_ot_k << <rstGridSize, BLOCK1D >> > (ot);
	create_ot_k << <otGridSize, BLOCK1D >> > (cnf, ot);
	LOGERR("Occurrence table creation failed");
	CHECK(cudaDeviceSynchronize());
	assert(ot->accViolation());
	if (p) {
		cout << "c |==========================|" << endl;
		cout << "c |==== occurrence table ====|" << endl;
		ot->print();
		cout << "c |==========================|" << endl;
	}
}
void createOTAsync(CNF* cnf, OT* ot, const bool& p)
{
	assert(cnf != NULL);
	assert(ot != NULL);
	uint32 cnf_size = nClauses() + (maxAddedCls() >> 1);
	uint32 rstGridSize = MIN(uint32((nDualVars() + BLOCK1D - 1) / BLOCK1D), maxGPUThreads / BLOCK1D);
	uint32 otGridSize = MIN((cnf_size + BLOCK1D - 1) / BLOCK1D, maxGPUThreads / BLOCK1D);
	reset_ot_k << <rstGridSize, BLOCK1D >> > (ot);
	create_ot_k << <otGridSize, BLOCK1D >> > (cnf, ot);
	if (p) {
		LOGERR("Occurrence table creation failed");
		CHECK(cudaDeviceSynchronize());
		assert(ot->accViolation());
		cout << "c |==========================|" << endl;
		cout << "c |==== occurrence table ====|" << endl;
		ot->print();
		cout << "c |==========================|" << endl;
	}
}
void ve(CNF *cnf, OT* ot, PV *pv)
{   
	assert(pv->numPVs > 0);
	uint32 nBlocks = MIN((pv->numPVs + BLVE - 1) / BLVE, maxGPUThreads / BLVE);
#if VE_DBG
	ve_k << <1, 1 >> > (cnf, ot, pv->pVars, pv->units);
#else
	ve_k << <nBlocks, BLVE >> > (cnf, ot, pv->pVars, pv->units);
#endif
	LOGERR("Parallel BVE failed");
	CHECK(cudaDeviceSynchronize());
}
void hse(CNF *cnf, OT* ot, PV *pv)
{
	assert(pv->numPVs > 0);
	uint32 nBlocks = MIN((pv->numPVs + BLHSE - 1) / BLHSE, maxGPUThreads / BLHSE);
#if SS_DBG
	hse_k << <1, 1 >> > (cnf, ot, pv->pVars, pv->units);
#else
	hse_k << <nBlocks, BLHSE >> > (cnf, ot, pv->pVars, pv->units);
#endif
	LOGERR("Parallel HSE failed");
	CHECK(cudaDeviceSynchronize());
}
void bce(CNF *cnf, OT* ot, PV *pv)
{
	assert(pv->numPVs > 0);
	uint32 nBlocks = MIN((pv->numPVs + BLBCE - 1) / BLBCE, maxGPUThreads / BLBCE);
	bce_k << <nBlocks, BLBCE >> > (cnf, ot, pv->pVars);
	LOGERR("Parallel BCE failed");
	CHECK(cudaDeviceSynchronize());
}
void hre(CNF *cnf, OT* ot, PV *pv)
{
	assert(pv->numPVs > 0);
	dim3 block2D(devProp.warpSize, devProp.warpSize), grid2D(1, 1, 1);
	grid2D.y = MIN((pv->numPVs + block2D.y - 1) / block2D.y, maxGPUThreads / block2D.y);
	uint32 smemSize = devProp.warpSize * (SH_MAX_HRE_IN + SH_MAX_HRE_OUT) * sizeof(uint32);
	hre_k << <grid2D, block2D, smemSize >> > (cnf, ot, pv->pVars);
	LOGERR("HRE Elimination failed");
	CHECK(cudaDeviceSynchronize());
}
void evalReds(CNF* cnf, GSTATS* gsts)
{
	stsResetAsync(gsts);
	uint32 cnf_sz = nClauses() + (maxAddedCls() >> 1);
	uint32 nBlocks1 = MIN((cnf_sz + (BLOCK1D << 1) - 1) / (BLOCK1D << 1), maxGPUThreads / (BLOCK1D << 1));
	uint32 smemSize1 = BLOCK1D * (sizeof(uint64) + sizeof(uint32));
	cnt_reds << <nBlocks1, BLOCK1D, smemSize1 >> > (cnf, gsts);
	uint32 nBlocks2 = MIN((nOrgVars() + (BLOCK1D << 1) - 1) / (BLOCK1D << 1), maxGPUThreads / (BLOCK1D << 1));
	uint32 smemSize2 = BLOCK1D * sizeof(uint32);
	cnt_del_vars << <nBlocks2, BLOCK1D, smemSize2 >> > (gsts, nOrgVars());
	LOGERR("Counting reductions failed");
	CHECK(cudaDeviceSynchronize());
	cnf_stats.n_del_vars = gsts->numDelVars;
	cnf_stats.n_cls_after = gsts->numClauses;
	cnf_stats.n_lits_after = gsts->numLits;
}
void countLits(CNF* cnf, GSTATS* gsts)
{
	gsts->numLits = 0;
	uint32 cnf_sz = nClauses() + (maxAddedCls() >> 1);
	uint32 nBlocks = MIN((cnf_sz + (BLOCK1D << 1) - 1) / (BLOCK1D << 1), maxGPUThreads / (BLOCK1D << 1));
	uint32 smemSize = BLOCK1D * sizeof(uint64);
	cnt_lits << <nBlocks, BLOCK1D, smemSize >> > (cnf, gsts);
	LOGERR("Counting literals failed");
	CHECK(cudaDeviceSynchronize());
	cnf_stats.n_lits_after = gsts->numLits;
}
void countCls(CNF* cnf, GSTATS* gsts)
{
	gsts->numClauses = 0;
	gsts->numLits = 0;
	uint32 cnf_sz = nClauses() + (maxAddedCls() >> 1);
	uint32 nBlocks = MIN((cnf_sz + (BLOCK1D << 1) - 1) / (BLOCK1D << 1), maxGPUThreads / (BLOCK1D << 1));
	uint32 smemSize = BLOCK1D * (sizeof(uint64) + sizeof(uint32));
	cnt_cls_lits << <nBlocks, BLOCK1D, smemSize >> > (cnf, gsts);
	LOGERR("Counting clauses-literals failed");
	CHECK(cudaDeviceSynchronize());
	cnf_stats.n_cls_after = gsts->numClauses;
	cnf_stats.n_lits_after = gsts->numLits;
}
void countDelVars(CNF* cnf, GSTATS* gsts, uint32* dvarsEA)
{
	stsResetAsync(gsts);
	uint32 cnf_sz = nClauses() + (maxAddedCls() >> 1);
	uint32 nBlocks1 = MIN((cnf_sz + BLOCK1D - 1) / BLOCK1D, maxGPUThreads / BLOCK1D);
	mark_vars << <nBlocks1, BLOCK1D >> > (cnf, gsts);
	uint32 nBlocks2 = MIN((nOrgVars() + (BLOCK1D << 1) - 1) / (BLOCK1D << 1), maxGPUThreads / (BLOCK1D << 1));
	uint32 smemSize = BLOCK1D * sizeof(uint32);
	cnt_del_vars << <nBlocks2, BLOCK1D, smemSize >> > (gsts, nOrgVars());
	CHECK(cudaMemcpyAsync(&cnf_stats.n_del_vars, dvarsEA, sizeof(uint32), cudaMemcpyDeviceToHost));
	CHECK(cudaStreamSynchronize(0));
}