
#include "pfdevice.cuh"

template<class T>
__global__ void memset_k(T* mem, T val, size_t sz)
{
	size_t i = blockDim.x * blockIdx.x + threadIdx.x;
	while (i < sz) { mem[i] = val; i += blockDim.x * gridDim.x; }
}

__global__ void reset_ot_k(OT* ot)
{
	int64 v = blockDim.x * blockIdx.x + threadIdx.x;
	while (v < ot->size()) { (*ot)[v].clear(); v += blockDim.x * gridDim.x; }
}

__global__ void create_ot_k(CNF* cnf, OT* ot, uint32 size)
{
	uint32 i = blockDim.x * blockIdx.x + threadIdx.x;
	while (i < size) {
		SCLAUSE& c = (*cnf)[i];
		if (c.status() == ORIGINAL || c.status() == LEARNT) {
#pragma unroll
			for (LIT_POS l = 0; l < c.size(); l++) (*ot)[c[l]].push(i);
		}
		i += blockDim.x * gridDim.x;
	}
}

__global__ void create_ot_k(CNF* cnf, OT* ot)
{
	uint32 i = blockDim.x * blockIdx.x + threadIdx.x;
	while (i < cnf->size()) {
		SCLAUSE& c = (*cnf)[i];
		if (c.status() == ORIGINAL || c.status() == LEARNT) {
#pragma unroll
			for (LIT_POS l = 0; l < c.size(); l++) (*ot)[c[l]].push(i);
		}
		i += blockDim.x * gridDim.x;
	}
}

__global__ void assign_scores(OCCUR* occurs, SCORE* scores, uint32* hist, uint32 size)
{
	uint32 v = blockDim.x * blockIdx.x + threadIdx.x;
	while (v < size) {
		uint32 p = V2D(v + 1), n = NEG(p);
		uint32 ps = hist[p], ns = hist[n];
		occurs[v].ps = ps, occurs[v].ns = ns;
		scores[v].v = v;
		scores[v].sc = (!ps || !ns) ? (ps | ns) : (ps * ns);
		v += blockDim.x * gridDim.x;
	}
}

__global__ void calc_sig_k(CNF* cnf, uint32 offset, uint32 size)
{
	uint32 i = blockDim.x * blockIdx.x + threadIdx.x + offset;
	while (i < size) { (*cnf)[i].calcSig(); i += blockDim.x * gridDim.x; }
}

__global__ void copy_k(uint32* dest, CNF* src, int64 size)
{
	int64 i = blockDim.x * blockIdx.x + threadIdx.x;
	int64 stride = blockDim.x * gridDim.x;
	while (i < size) { dest[i] = *src->data(i); i += stride; }
}

__global__ void copy_if_k(uint32* dest, CNF* src, GSTATS* gstats)
{
	uint32 i = blockDim.x * blockIdx.x + threadIdx.x;
	uint32 stride = blockDim.x * gridDim.x;
	while (i < src->size()) {
		SCLAUSE& c = (*src)[i];
		if (c.status() == ORIGINAL || c.status() == LEARNT) {
			uint64 lits_idx = atomicAdd(&gstats->numLits, c.size());
#pragma unroll
			for (LIT_POS l = 0; l < c.size(); l++) dest[lits_idx++] = c[l];
		}
		i += stride;
	}
}

__global__ void calc_added_cls_k(CNF* cnf, OT* ot, cuVec<uint32>* pVars, GSTATS* gstats)
{
	uint32* sh_rCls = SharedMemory<uint32>();
	uint32 tx = threadIdx.x;
	uint32 v = blockIdx.x * (blockDim.x << 1) + tx;
	uint32 stride = (blockDim.x << 1) * gridDim.x;
	uint32 x, p, n, nCls = 0;
	while (v < pVars->size()) {
		x = (*pVars)[v], p = V2D(x + 1), n = NEG(p);
		calcResolvents(x + 1, *cnf, (*ot)[p], (*ot)[n], nCls);
		if (v + blockDim.x < pVars->size()) {
			x = (*pVars)[v + blockDim.x], p = V2D(x + 1), n = NEG(p);
			calcResolvents(x + 1, *cnf, (*ot)[p], (*ot)[n], nCls);
		}
		v += stride;
	}
	// write local sum in shared memory
	sh_rCls[tx] = (tx < pVars->size()) ? nCls : 0;
	__syncthreads();
	// do reduction in shared memory
	if ((blockDim.x >= 512) && (tx < 256)) sh_rCls[tx] = nCls = nCls + sh_rCls[tx + 256];
	__syncthreads();
	if ((blockDim.x >= 256) && (tx < 128)) sh_rCls[tx] = nCls = nCls + sh_rCls[tx + 128];
	__syncthreads();
	if ((blockDim.x >= 128) && (tx < 64)) sh_rCls[tx] = nCls = nCls + sh_rCls[tx + 64];
	__syncthreads();
	if (tx < warpSize) {
		if (blockDim.x >= 64) nCls += sh_rCls[tx + warpSize];
		for (int i = (warpSize >> 1); i > 0; i >>= 1)
			nCls += __shfl_xor_sync(0xffffffff, nCls, i, warpSize);
	}
	// write result for this block to global mem
	if (tx == 0) atomicAdd(&gstats->numClauses, nCls);
}

__global__ void calc_added_all_k(CNF* cnf, OT* ot, cuVec<uint32>* pVars, GSTATS* gstats)
{
	uint32* sh_rCls = SharedMemory<uint32>();
	uint64* sh_rLits = (uint64*)(sh_rCls + blockDim.x);
	uint32 tx = threadIdx.x;
	uint32 v = blockIdx.x * (blockDim.x << 1) + tx;
	uint32 stride = (blockDim.x << 1) * gridDim.x;
	uint32 x, p, n, nCls = 0;
	uint64 nLits = 0;
	while (v < pVars->size()) {
		x = (*pVars)[v], p = V2D(x + 1), n = NEG(p);
		calcResolvents(x + 1, *cnf, (*ot)[p], (*ot)[n], nCls, nLits);
		if (v + blockDim.x < pVars->size()) {
			x = (*pVars)[v + blockDim.x], p = V2D(x + 1), n = NEG(p);
			calcResolvents(x + 1, *cnf, (*ot)[p], (*ot)[n], nCls, nLits);
		}
		v += stride;
	}
	// write local sum in shared memory
	sh_rCls[tx] = (tx < pVars->size()) ? nCls : 0;
	sh_rLits[tx] = (tx < pVars->size()) ? nLits : 0;
	__syncthreads();
	// do reduction in shared memory
	if ((blockDim.x >= 512) && (tx < 256)) {
		sh_rCls[tx] = nCls = nCls + sh_rCls[tx + 256];
		sh_rLits[tx] = nLits = nLits + sh_rLits[tx + 256];
	}
	__syncthreads();
	if ((blockDim.x >= 256) && (tx < 128)) {
		sh_rCls[tx] = nCls = nCls + sh_rCls[tx + 128];
		sh_rLits[tx] = nLits = nLits + sh_rLits[tx + 128];
	}
	__syncthreads();
	if ((blockDim.x >= 128) && (tx < 64)) {
		sh_rCls[tx] = nCls = nCls + sh_rCls[tx + 64];
		sh_rLits[tx] = nLits = nLits + sh_rLits[tx + 64];
	}
	__syncthreads();
	if (tx < warpSize) {
		if (blockDim.x >= 64) {
			nCls += sh_rCls[tx + warpSize];
			nLits += sh_rLits[tx + warpSize];
		}
		for (int i = (warpSize >> 1); i > 0; i >>= 1) {
			nCls += __shfl_xor_sync(0xffffffff, nCls, i, warpSize);
			nLits += __shfl_xor_sync(0xffffffff, nLits, i, warpSize);
		}
	}
	// write result for this block to global mem
	if (tx == 0) {
		atomicAdd(&gstats->numClauses, nCls);
		atomicAdd(&gstats->numLits, nLits);
	}
}

__global__ void cnt_del_vars(GSTATS* gstats, uint32 size)
{
	uint32* sh_delVars = SharedMemory<uint32>();
	uint32 tx = threadIdx.x;
	uint32 i = blockIdx.x * (blockDim.x << 1) + tx;
	uint32 stride = (blockDim.x << 1) * gridDim.x;
	uint32 nDelVars = 0;
	while (i < size) {
		if (!gstats->seen[i]) nDelVars++;
		if (i + blockDim.x < size && !gstats->seen[i + blockDim.x]) nDelVars++;
		i += stride;
	}
	sh_delVars[tx] = (tx < size) ? nDelVars : 0;
	__syncthreads();
	if ((blockDim.x >= 512) && (tx < 256)) sh_delVars[tx] = nDelVars = nDelVars + sh_delVars[tx + 256];
	__syncthreads();
	if ((blockDim.x >= 256) && (tx < 128)) sh_delVars[tx] = nDelVars = nDelVars + sh_delVars[tx + 128];
	__syncthreads();
	if ((blockDim.x >= 128) && (tx < 64)) sh_delVars[tx] = nDelVars = nDelVars + sh_delVars[tx + 64];
	__syncthreads();
	if (tx < warpSize) {
		if (blockDim.x >= 64) nDelVars += sh_delVars[tx + warpSize];
		for (int i = (warpSize >> 1); i > 0; i >>= 1)
			nDelVars += __shfl_xor_sync(0xffffffff, nDelVars, i, warpSize);
	}
	if (tx == 0) atomicAdd(&gstats->numDelVars, nDelVars);
}

__global__ void cnt_reds(CNF* cnf, GSTATS* gstats)
{
	uint32* sh_rCls = SharedMemory<uint32>();
	uint64* sh_rLits = (uint64*)(sh_rCls + blockDim.x);
	uint32 tx = threadIdx.x;
	uint32 i = blockIdx.x * (blockDim.x << 1) + tx;
	uint32 stride = (blockDim.x << 1) * gridDim.x;
	uint32 nCls = 0;
	uint64 nLits = 0;
	while (i < cnf->size()) {
		SCLAUSE& c1 = (*cnf)[i], &c2 = (*cnf)[i + blockDim.x];
		if (c1.status() == LEARNT || c1.status() == ORIGINAL) {
			CL_LEN cl_size = c1.size();
			nCls++, nLits += cl_size;
			for (LIT_POS k = 0; k < cl_size; k++) { assert(c1[k]); gstats->seen[V2X(c1[k])] = 1; }
		}
		if (i + blockDim.x < cnf->size() && (c2.status() == LEARNT || c2.status() == ORIGINAL)) {
			CL_LEN cl_size = c2.size();
			nCls++, nLits += cl_size;
			for (LIT_POS k = 0; k < cl_size; k++) { assert(c2[k]); gstats->seen[V2X(c2[k])] = 1; }
		}
		i += stride;
	}
	// write local sum in shared memory
	sh_rCls[tx] = (tx < cnf->size()) ? nCls : 0;
	sh_rLits[tx] = (tx < cnf->size()) ? nLits : 0;
	__syncthreads();
	// do reduction in shared memory
	if ((blockDim.x >= 512) && (tx < 256)) {
		sh_rCls[tx] = nCls = nCls + sh_rCls[tx + 256];
		sh_rLits[tx] = nLits = nLits + sh_rLits[tx + 256];
	}
	__syncthreads();
	if ((blockDim.x >= 256) && (tx < 128)) {
		sh_rCls[tx] = nCls = nCls + sh_rCls[tx + 128];
		sh_rLits[tx] = nLits = nLits + sh_rLits[tx + 128];
	}
	__syncthreads();
	if ((blockDim.x >= 128) && (tx < 64)) {
		sh_rCls[tx] = nCls = nCls + sh_rCls[tx + 64];
		sh_rLits[tx] = nLits = nLits + sh_rLits[tx + 64];
	}
	__syncthreads();
	if (tx < warpSize) {
		if (blockDim.x >= 64) {
			nCls += sh_rCls[tx + warpSize];
			nLits += sh_rLits[tx + warpSize];
		}
		for (int i = (warpSize >> 1); i > 0; i >>= 1) {
			nCls += __shfl_xor_sync(0xffffffff, nCls, i, warpSize);
			nLits += __shfl_xor_sync(0xffffffff, nLits, i, warpSize);
		}
	}
	// write result for this block to global mem
	if (tx == 0) {
		atomicAdd(&gstats->numClauses, nCls);
		atomicAdd(&gstats->numLits, nLits);
	}
}

__global__ void cnt_lits(CNF* cnf, GSTATS* gstats)
{
	uint64* sh_rLits = SharedMemory<uint64>();
	uint32 tx = threadIdx.x;
	uint32 i = blockIdx.x * (blockDim.x << 1) + tx;
	uint32 stride = (blockDim.x << 1) * gridDim.x;
	uint64 nLits = 0;
	while (i < cnf->size()) {
		SCLAUSE& c1 = (*cnf)[i], &c2 = (*cnf)[i + blockDim.x];
		if (c1.status() == LEARNT || c1.status() == ORIGINAL) nLits += c1.size();
		if (i + blockDim.x < cnf->size() && (c2.status() == LEARNT || c2.status() == ORIGINAL)) nLits += c2.size();
		i += stride;
	}
	sh_rLits[tx] = (tx < cnf->size()) ? nLits : 0;
	__syncthreads();
	if ((blockDim.x >= 512) && (tx < 256)) sh_rLits[tx] = nLits = nLits + sh_rLits[tx + 256];
	__syncthreads();
	if ((blockDim.x >= 256) && (tx < 128)) sh_rLits[tx] = nLits = nLits + sh_rLits[tx + 128];
	__syncthreads();
	if ((blockDim.x >= 128) && (tx < 64)) sh_rLits[tx] = nLits = nLits + sh_rLits[tx + 64];
	__syncthreads();
	if (tx < warpSize) {
		if (blockDim.x >= 64) nLits += sh_rLits[tx + warpSize];
		for (int i = (warpSize >> 1); i > 0; i >>= 1)
			nLits += __shfl_xor_sync(0xffffffff, nLits, i, warpSize);
	}
	if (tx == 0) atomicAdd(&gstats->numLits, nLits);
}

__global__ void cnt_cls_lits(CNF* cnf, GSTATS* gstats)
{
	uint32* sh_rCls = SharedMemory<uint32>();
	uint64* sh_rLits = (uint64*)(sh_rCls + blockDim.x);
	uint32 tx = threadIdx.x;
	uint32 i = blockIdx.x * (blockDim.x << 1) + tx;
	uint32 stride = (blockDim.x << 1) * gridDim.x;
	uint32 nCls = 0;
	uint64 nLits = 0;
	while (i < cnf->size()) {
		SCLAUSE& c1 = (*cnf)[i], &c2 = (*cnf)[i + blockDim.x];
		if (c1.status() == LEARNT || c1.status() == ORIGINAL) { nCls++, nLits += c1.size(); }
		if (i + blockDim.x < cnf->size() && 
			(c2.status() == LEARNT || c2.status() == ORIGINAL)) { nCls++, nLits += c2.size(); }
		i += stride;
	}
	sh_rCls[tx] = (tx < cnf->size()) ? nCls : 0;
	sh_rLits[tx] = (tx < cnf->size()) ? nLits : 0;
	__syncthreads();
	if ((blockDim.x >= 512) && (tx < 256)) {
		sh_rCls[tx] = nCls = nCls + sh_rCls[tx + 256];
		sh_rLits[tx] = nLits = nLits + sh_rLits[tx + 256];
	}
	__syncthreads();
	if ((blockDim.x >= 256) && (tx < 128)) {
		sh_rCls[tx] = nCls = nCls + sh_rCls[tx + 128];
		sh_rLits[tx] = nLits = nLits + sh_rLits[tx + 128];
	}
	__syncthreads();
	if ((blockDim.x >= 128) && (tx < 64)) {
		sh_rCls[tx] = nCls = nCls + sh_rCls[tx + 64];
		sh_rLits[tx] = nLits = nLits + sh_rLits[tx + 64];
	}
	__syncthreads();
	if (tx < warpSize) {
		if (blockDim.x >= 64) {
			nCls += sh_rCls[tx + warpSize];
			nLits += sh_rLits[tx + warpSize];
		}
		for (int i = (warpSize >> 1); i > 0; i >>= 1) {
			nCls += __shfl_xor_sync(0xffffffff, nCls, i, warpSize);
			nLits += __shfl_xor_sync(0xffffffff, nLits, i, warpSize);
		}
	}
	// write result for this block to global mem
	if (tx == 0) {
		atomicAdd(&gstats->numClauses, nCls);
		atomicAdd(&gstats->numLits, nLits);
	}
}

__global__ void ve_k(CNF *cnf, OT* ot, cuVec<uint32>* pVars, GSOL* sol)
{
	uint32 tx = threadIdx.x;
	uint32 v = blockDim.x * blockIdx.x + threadIdx.x;
	__shared__ uint32 defs[BLSIMP * FAN_LMT];
	__shared__ uint32 outs[BLSIMP * MAX_SH_RES_LEN];
	while (v < pVars->size()) {
		uint32 x = (*pVars)[v], p = V2D(x + 1), n = NEG(p);
		assert(sol->value[x] == UNDEFINED);
		if ((*ot)[p].size() == 0 || (*ot)[n].size() == 0) { // pure
			deleteClauses(*cnf, (*ot)[p], (*ot)[n]);
			sol->value[x] = (*ot)[p].size() ? 1 : 0;
			(*ot)[p].clear(), (*ot)[n].clear();
		}
		else if ((*ot)[p].size() == 1 || (*ot)[n].size() == 1) {
			resolve_x(x + 1, *cnf, (*ot)[p], (*ot)[n], sol, &outs[tx * MAX_SH_RES_LEN]);
			sol->value[x] = 1;
			(*ot)[p].clear(), (*ot)[n].clear();
		}
		else if (gateReasoning_x(p, *cnf, (*ot)[p], (*ot)[n], sol, &defs[tx * FAN_LMT], &outs[tx * MAX_SH_RES_LEN])
			  || resolve_x(x + 1, *cnf, (*ot)[p], (*ot)[n], sol, &outs[tx * MAX_SH_RES_LEN], true)) {
			sol->value[x] = 1;
			(*ot)[p].clear(), (*ot)[n].clear();
		}
		v += blockDim.x * gridDim.x;
	}
}

__global__ void hse_k(CNF *cnf, OT* ot, cuVec<uint32>* pVars, GSOL* sol)
{
	uint32 v = blockDim.x * blockIdx.x + threadIdx.x;
	__shared__ uint32 sh_cls[BLSIMP * SHARED_CL_LEN];
	while (v < pVars->size()) {
		assert(sol->value[pVars->at(v)] == UNDEFINED);
		uint32 p = V2D(pVars->at(v) + 1), n = NEG(p);
		self_sub_x(p, *cnf, (*ot)[p], (*ot)[n], sol, &sh_cls[threadIdx.x * SHARED_CL_LEN]);
		v += blockDim.x * gridDim.x;
	}
}

__global__ void bce_k(CNF *cnf, OT* ot, cuVec<uint32>* pVars, GSOL* sol)
{
	uint32 v = blockDim.x * blockIdx.x + threadIdx.x;
	__shared__ uint32 sh_cls[BLSIMP * CL_MAX_LEN_BCE];
	while (v < pVars->size()) {
		uint32 x = (*pVars)[v], p = V2D(x + 1), n = NEG(p);
		assert(sol->value[x] == UNDEFINED);
		blocked_x(x + 1, *cnf, (*ot)[p], (*ot)[n], &sh_cls[threadIdx.x * CL_MAX_LEN_BCE]);
		v += blockDim.x * gridDim.x;
	}
}

__global__ void hre_k(CNF *cnf, OT* ot, cuVec<uint32>* pVars, GSOL* sol)
{
	uint32 v = blockDim.y * blockIdx.y + threadIdx.y;
	uint32* smem = SharedMemory<uint32>();
	uint32* m_c = smem + warpSize * SH_MAX_HRE_IN + threadIdx.y * SH_MAX_HRE_OUT; // shared memory for resolvent
	while (v < pVars->size()) {
		assert(sol->value[pVars->at(v)] == UNDEFINED);
		uint32 p = V2D(pVars->at(v) + 1);
		OL& poss = (*ot)[p], &negs = (*ot)[NEG(p)];
		// do merging and apply forward equality check (on-the-fly) over resolvents
#pragma unroll
		for (uint32 i = 0; i < poss.size(); i++) {
			if ((*cnf)[poss[i]].status() == DELETED) continue;
			uint32 pos_size = (*cnf)[poss[i]].size();
			if (pos_size <= SH_MAX_HRE_IN) { // use shared memory for positives
				uint32* sh_pos = smem + threadIdx.y * SH_MAX_HRE_IN;
				if (threadIdx.x == 0) (*cnf)[poss[i]].shareTo(sh_pos);
#pragma unroll
				for (uint32 j = 0; j < negs.size(); j++) {
					SCLAUSE& neg = (*cnf)[negs[j]];
					if (neg.status() == DELETED || (pos_size + neg.size() - 2) > SH_MAX_HRE_OUT) continue;
					CL_LEN m_c_size = 0;
					uint32 m_c_sig = 0;
					if (threadIdx.x == 0) {
						assert(warpSize == blockDim.x);
						m_c_size = merge(pVars->at(v) + 1, sh_pos, pos_size, neg, m_c);
						calcSig(m_c, m_c_size, m_c_sig);
					}
					m_c_size = __shfl_sync(0xffffffff, m_c_size, 0);
					m_c_sig = __shfl_sync(0xffffffff, m_c_sig, 0);
					forward_equ(*cnf, *ot, m_c, m_c_sig, m_c_size);
				}
			}
			else { // use global memory
#pragma unroll
				for (uint32 j = 0; j < negs.size(); j++) {
					SCLAUSE& neg = (*cnf)[negs[j]];
					if (neg.status() == DELETED || (pos_size + neg.size() - 2) > SH_MAX_HRE_OUT) continue;
					CL_LEN m_c_size = 0;
					uint32 m_c_sig = 0;
					if (threadIdx.x == 0) {
						assert(warpSize == blockDim.x);
						m_c_size = merge(pVars->at(v) + 1, (*cnf)[poss[i]], neg, m_c);
						calcSig(m_c, m_c_size, m_c_sig);
					}
					m_c_size = __shfl_sync(0xffffffff, m_c_size, 0);
					m_c_sig = __shfl_sync(0xffffffff, m_c_sig, 0);
					forward_equ(*cnf, *ot, m_c, m_c_sig, m_c_size);
				} 
			} 
		} 
		v += blockDim.y * gridDim.y;
	}
}

__global__ void prop_k() 
{

}
//==============================================//
//          ParaFROST Wrappers/helpers          //
//==============================================//
void mem_set(addr_t mem, const Byte& val, const size_t& size)
{
	int nBlocks = MIN((size + BLOCK1D - 1) / BLOCK1D, maxGPUThreads / BLOCK1D);
	memset_k<Byte> << <nBlocks, BLOCK1D >> > (mem, val, size);
	LOGERR("Memory set failed");
	CHECK(cudaDeviceSynchronize());
}
void mem_set(LIT_ST* mem, const LIT_ST& val, const size_t& size)
{
	int nBlocks = MIN((size + BLOCK1D - 1) / BLOCK1D, maxGPUThreads / BLOCK1D);
	memset_k<LIT_ST> << <nBlocks, BLOCK1D >> > (mem, val, size);
	LOGERR("Memory set failed");
	CHECK(cudaDeviceSynchronize());
}
void copy(uint32* dest, CNF* src, const int64& size)
{
	int nBlocks = MIN((size + BLOCK1D - 1) / BLOCK1D, maxGPUThreads / BLOCK1D);
	copy_k << <nBlocks, BLOCK1D >> > (dest, src, size);
	LOGERR("Copying failed");
	CHECK(cudaDeviceSynchronize());
}
void copyIf(uint32* dest, CNF* src, GSTATS* gstats)
{
	gstats->numLits = 0;
	int nBlocks = MIN((nClauses() + maxAddedCls() + BLOCK1D - 1) / BLOCK1D, maxGPUThreads / BLOCK1D);
	copy_if_k << <nBlocks, BLOCK1D >> > (dest, src, gstats);
	LOGERR("Copying failed");
	CHECK(cudaDeviceSynchronize());
}
void calc_vscores(OCCUR* occurs, SCORE* scores, uint32* histogram)
{
	int nBlocks = MIN((nOrgVars() + BLOCK1D - 1) / BLOCK1D, maxGPUThreads / BLOCK1D);
	assign_scores << <nBlocks, BLOCK1D >> > (occurs, scores, histogram, nOrgVars());
	LOGERR("Assigning scores failed");
	CHECK(cudaDeviceSynchronize());
}
void calc_added(CNF* cnf, OT* ot, PV* pv, GSTATS* gstats)
{
	assert(pv->numPVs > 0);
	gstats->numClauses = 0;
	int nBlocks = MIN((pv->numPVs + (BLUB << 1) - 1) / (BLUB << 1), maxGPUThreads / (BLUB << 1));
	int smemSize = BLUB * sizeof(uint32);
	calc_added_cls_k << <nBlocks, BLUB, smemSize >> > (cnf, ot, pv->pVars, gstats);
	LOGERR("Added clauses calculation failed");
	CHECK(cudaDeviceSynchronize());
	cnf_stats.max_added_cls = gstats->numClauses;
	cnf_stats.max_added_lits = cnf_stats.max_added_cls * MAX_GL_RES_LEN;
	printf("c | added cls = %d\n", maxAddedCls()), printf("c | added lits = %zd\n", maxAddedLits());
}
void calc_sig(CNF* cnf, const uint32& offset, const uint32& size)
{
	int nBlocks = MIN((size + BLOCK1D - 1) / BLOCK1D, maxGPUThreads / BLOCK1D);
	calc_sig_k << <nBlocks, BLOCK1D >> > (cnf, offset, size);
	LOGERR("Calculating signatures failed");
	CHECK(cudaDeviceSynchronize());
}
void create_ot(CNF* cnf, OT* ot, const bool& p)
{
	assert(cnf != NULL);
	assert(ot != NULL);
	int rstGridSize = MIN((V2D(nOrgVars() + 1) + BLOCK1D - 1) / BLOCK1D, maxGPUThreads / BLOCK1D);
	int otGridSize = MIN((nClauses() + maxAddedCls() + BLOCK1D - 1) / BLOCK1D, maxGPUThreads / BLOCK1D);
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
void ve(CNF *cnf, OT* ot, PV *pv)
{   
	assert(pv->numPVs > 0);
	int nBlocks = MIN((pv->numPVs + BLSIMP - 1) / BLSIMP, maxGPUThreads / BLSIMP);
#if VE_DBG
	ve_k << <1, 1 >> > (cnf, ot, pv->pVars, pv->sol);
#else
	ve_k << <nBlocks, BLSIMP >> > (cnf, ot, pv->pVars, pv->sol);
#endif
	//if (prop() == UNSAT) killSolver(UNSAT);
	LOGERR("Parallel BVE failed");
	CHECK(cudaDeviceSynchronize());
}
void hse(CNF *cnf, OT* ot, PV *pv)
{
	assert(pv->numPVs > 0);
	int nBlocks = MIN((pv->numPVs + BLSIMP - 1) / BLSIMP, maxGPUThreads / BLSIMP);
#if SS_DBG
	hse_k << <1, 1 >> > (cnf, ot, pv->pVars, pv->sol);
#else
	hse_k << <nBlocks, BLSIMP >> > (cnf, ot, pv->pVars, pv->sol);
#endif
	//if (prop() == UNSAT) killSolver(UNSAT);
	LOGERR("Parallel HSE failed");
	CHECK(cudaDeviceSynchronize());
}
void bce(CNF *cnf, OT* ot, PV *pv)
{
	assert(pv->numPVs > 0);
	int nBlocks = MIN((pv->numPVs + BLSIMP - 1) / BLSIMP, maxGPUThreads / BLSIMP);
	bce_k << <nBlocks, BLSIMP >> > (cnf, ot, pv->pVars, pv->sol);
	LOGERR("Parallel BCE failed");
	CHECK(cudaDeviceSynchronize());
}
void hre(CNF *cnf, OT* ot, PV *pv)
{
	assert(pv->numPVs > 0);
	dim3 block2D(devProp.warpSize, devProp.warpSize), grid2D(1, 1, 1);
	grid2D.y = MIN((pv->numPVs + block2D.y - 1) / block2D.y, maxGPUThreads / block2D.y);
	int smemSize = devProp.warpSize * (SH_MAX_HRE_IN + SH_MAX_HRE_OUT) * sizeof(uint32);
	hre_k << <grid2D, block2D, smemSize >> > (cnf, ot, pv->pVars, pv->sol);
	LOGERR("HRE Elimination failed");
	CHECK(cudaDeviceSynchronize());
}
void evalReds(CNF* cnf, GSTATS* gstats)
{
	gstats->numDelVars = 0;
	gstats->numClauses = 0;
	gstats->numLits = 0;
	mem_set(gstats->seen, 0, nOrgVars());
	uint32 cnf_sz = nClauses() + maxAddedCls();
	int nBlocks1 = MIN((cnf_sz + (BLOCK1D << 1) - 1) / (BLOCK1D << 1), maxGPUThreads / (BLOCK1D << 1));
	int nBlocks2 = MIN((nOrgVars() + (BLOCK1D << 1) - 1) / (BLOCK1D << 1), maxGPUThreads / (BLOCK1D << 1));
	int smemSize1 = BLOCK1D * (sizeof(uint64) + sizeof(uint32));
	int smemSize2 = BLOCK1D * sizeof(uint32);
	cnt_reds << <nBlocks1, BLOCK1D, smemSize1 >> > (cnf, gstats);
	cnt_del_vars << <nBlocks2, BLOCK1D, smemSize2 >> > (gstats, nOrgVars());
	LOGERR("Counting reductions failed");
	CHECK(cudaDeviceSynchronize());
	cnf_stats.n_del_vars = gstats->numDelVars;
	cnf_stats.n_cls_after = gstats->numClauses;
	cnf_stats.n_lits_after = gstats->numLits;
}
void countLits(CNF* cnf, GSTATS* gstats)
{
	gstats->numLits = 0;
	uint32 cnf_sz = nClauses() + maxAddedCls();
	int nBlocks = MIN((cnf_sz + (BLOCK1D << 1) - 1) / (BLOCK1D << 1), maxGPUThreads / (BLOCK1D << 1));
	int smemSize = BLOCK1D * sizeof(uint64);
	cnt_lits << <nBlocks, BLOCK1D, smemSize >> > (cnf, gstats);
	LOGERR("Counting literals failed");
	CHECK(cudaDeviceSynchronize());
	cnf_stats.n_lits_after = gstats->numLits;
}
void countCls(CNF* cnf, GSTATS* gstats)
{
	gstats->numClauses = 0;
	gstats->numLits = 0;
	uint32 cnf_sz = nClauses() + maxAddedCls();
	int nBlocks = MIN((cnf_sz + (BLOCK1D << 1) - 1) / (BLOCK1D << 1), maxGPUThreads / (BLOCK1D << 1));
	int smemSize = BLOCK1D * (sizeof(uint64) + sizeof(uint32));
	cnt_cls_lits << <nBlocks, BLOCK1D, smemSize >> > (cnf, gstats);
	LOGERR("Counting clauses-literals failed");
	CHECK(cudaDeviceSynchronize());
	cnf_stats.n_cls_after = gstats->numClauses;
	cnf_stats.n_lits_after = gstats->numLits;
}