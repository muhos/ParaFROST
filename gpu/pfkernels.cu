/***********************************************************************[pfkernels.cu]
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

#include "pfdevice.cuh"
#include "pfbve.cuh"
#include "pfhse.cuh"

namespace pFROST {

	namespace SIGmA {

		template<class T>
		__global__ void memset_k(T* mem, T val, size_t size)
		{
			size_t tid = global_tx();
			while (tid < size) { mem[tid] = val; tid += stride_x(); }
		}

		__global__ void reset_stats(GSTATS* gstats) { gstats->numDelVars = 0, gstats->numClauses = 0, gstats->numLits = 0; }

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
				S_REF r = cnf->ref(tid);
				SCLAUSE& c = (*cnf)[r];
				if (c.original() || c.learnt()) {
					uint32* lit = c, * cend = c.end();
#pragma unroll
					while (lit != cend) (*ot)[*lit++].insert(r);
				}
				tid += stride_x();
			}
		}

		__global__ void assign_scores(uint32* eligible, uint32* scores, uint32* hist, uint32 size)
		{
			uint32 tid = global_tx();
			while (tid < size) {
				uint32 v = tid + 1;
				uint32 p = V2D(v), ps = hist[p], ns = hist[NEG(p)];
				eligible[tid] = v;
				scores[v] = rscore(ps, ns);
				tid += stride_x();
			}
		}

		__global__ void assign_scores(uint32* eligible, uint32* scores, uint32* hist, OT* ot, uint32 size)
		{
			uint32 tid = global_tx();
			while (tid < size) {
				uint32 v = tid + 1;
				uint32 p = V2D(v), n = NEG(p), ps = (*ot)[p].size(), ns = (*ot)[n].size();
				hist[p] = ps, hist[n] = ns;
				eligible[tid] = v;
				scores[v] = rscore(ps, ns);
				tid += stride_x();
			}
		}

		__global__ void calc_sig_k(CNF* cnf, uint32 offset, uint32 size)
		{
			uint32 tid = global_tx() + offset;
			while (tid < size) { calcSig(cnf->clause(tid)); tid += stride_x(); }
		}

		__global__ void copy_if_k(uint32* dest, CNF* src, GSTATS* gstats)
		{
			uint32 tid = global_tx();
			while (tid < src->size()) {
				SCLAUSE& c = src->clause(tid);
				if (c.original() || c.learnt()) {
					uint32* d = dest + atomicAdd(&gstats->numLits, c.size());
					uint32* s = c, *cend = c.end();
#pragma unroll
					while (s != cend) *d++ = *s++;
				}
				tid += stride_x();
			}
		}

		__global__ void copy_if_k(CNF* dest, CNF* src)
		{
			uint32 i = global_tx();
			while (i < src->size()) {
				SCLAUSE& s = src->clause(i);
				if (s.original() || s.learnt())
					dest->insert(s);
				i += stride_x();
			}
		}

		__global__ void cnt_reds(CNF* cnf, GSTATS* gstats)
		{
			uint32* sh_rCls = SharedMemory<uint32>();
			uint32* sh_rLits = sh_rCls + blockDim.x;
			uint32 tid = global_tx_off();
			uint32 nCls = 0;
			uint32 nLits = 0;
			while (tid < cnf->size()) {
				SCLAUSE& c1 = cnf->clause(tid);
				if (c1.original() || c1.learnt())
					nCls++, nLits += c1.size();
				if (tid + blockDim.x < cnf->size()) {
					SCLAUSE& c2 = cnf->clause(tid + blockDim.x);
					if (c2.original() || c2.learnt())
						nCls++, nLits += c2.size();
				}
				tid += stride_x_off();
			}
			loadShared(sh_rCls, nCls, sh_rLits, nLits, cnf->size());
			sharedReduce(sh_rCls, nCls, sh_rLits, nLits);
			warpReduce(sh_rCls, nCls, sh_rLits, nLits);
			if (threadIdx.x == 0) {
				atomicAdd(&gstats->numClauses, nCls);
				atomicAdd(&gstats->numLits, nLits);
			}
		}

		__global__ void cnt_cls(CNF* cnf, GSTATS* gstats)
		{
			uint32* sh_rCls = SharedMemory<uint32>();
			uint32 tid = global_tx_off();
			uint32 nCls = 0;
			while (tid < cnf->size()) {
				SCLAUSE& c1 = cnf->clause(tid);
				if (c1.original() || c1.learnt()) nCls++;
				uint32 off = tid + blockDim.x;
				if (off < cnf->size()) {
					SCLAUSE& c2 = cnf->clause(off);
					if (c2.original() || c2.learnt()) nCls++;
				}
				tid += stride_x_off();
			}
			loadShared(sh_rCls, nCls, cnf->size());
			sharedReduce(sh_rCls, nCls);
			warpReduce(sh_rCls, nCls);
			if (threadIdx.x == 0) atomicAdd(&gstats->numClauses, nCls);
		}

		__global__ void cnt_lits(CNF* cnf, GSTATS* gstats)
		{
			uint32* sh_rLits = SharedMemory<uint32>();
			uint32 tid = global_tx_off();
			uint32 nLits = 0;
			while (tid < cnf->size()) {
				SCLAUSE& c1 = cnf->clause(tid);
				if (c1.original() || c1.learnt()) nLits += c1.size();
				uint32 off = tid + blockDim.x;
				if (off < cnf->size()) {
					SCLAUSE& c2 = cnf->clause(off);
					if (c2.original() || c2.learnt()) nLits += c2.size();
				}
				tid += stride_x_off();
			}
			loadShared(sh_rLits, nLits, cnf->size());
			sharedReduce(sh_rLits, nLits);
			warpReduce(sh_rLits, nLits);
			if (threadIdx.x == 0) atomicAdd(&gstats->numLits, nLits);
		}

		__global__ void cnt_cls_lits(CNF* cnf, GSTATS* gstats)
		{
			uint32* sh_rCls = SharedMemory<uint32>();
			uint32* sh_rLits = sh_rCls + blockDim.x;
			uint32 tid = global_tx_off();
			uint32 nCls = 0;
			uint32 nLits = 0;
			while (tid < cnf->size()) {
				SCLAUSE& c1 = cnf->clause(tid);
				if (c1.original() || c1.learnt()) nCls++, nLits += c1.size();
				uint32 off = tid + blockDim.x;
				if (off < cnf->size()) {
					SCLAUSE& c2 = cnf->clause(off);
					if (c2.original() || c2.learnt()) nCls++, nLits += c2.size();
				}
				tid += stride_x_off();
			}
			loadShared(sh_rCls, nCls, sh_rLits, nLits, cnf->size());
			sharedReduce(sh_rCls, nCls, sh_rLits, nLits);
			warpReduce(sh_rCls, nCls, sh_rLits, nLits);
			if (threadIdx.x == 0) {
				atomicAdd(&gstats->numClauses, nCls);
				atomicAdd(&gstats->numLits, nLits);
			}
		}

		__global__ void ve_k(CNF* cnf, OT* ot, cuVecU* pVars, cuVecU* units, cuVecU* resolved)
		{
			uint32 tx = threadIdx.x;
			uint32 tid = global_tx();
			__shared__ uint32 defs[BLVE * FAN_LMT];
			__shared__ uint32 outs[BLVE * SH_MAX_BVE_OUT];
			while (tid < pVars->size()) {
				uint32& x = (*pVars)[tid];
				assert(x);
				assert(!ELIMINATED(x));
				uint32 p = V2D(x), n = NEG(p);
				uint32 pOrgs = (*ot)[p].size(), nOrgs = (*ot)[n].size();
				// try pure-literal elimination
				if (!pOrgs || !nOrgs) {
					uint32 psLits = 0, nsLits = 0;
					countLitsBefore(*cnf, (*ot)[p], psLits);
					countLitsBefore(*cnf, (*ot)[n], nsLits);
					toblivion(*cnf, (*ot)[p], (*ot)[n], resolved, pOrgs, nOrgs, psLits, nsLits, p);
					x |= MELTING_MASK;
				}
				// try simple resolution
				else if ((pOrgs == 1 || nOrgs == 1) &&
					resolve_x(x, pOrgs, nOrgs, *cnf, (*ot)[p], (*ot)[n], units, resolved, &outs[tx * SH_MAX_BVE_OUT])) x |= MELTING_MASK;
				// try gate-equivalence reasoning, otherwise resolution
				else if (gateReasoning_x(p, pOrgs, nOrgs, *cnf, (*ot)[p], (*ot)[n], units, resolved, &defs[tx * FAN_LMT], &outs[tx * SH_MAX_BVE_OUT]) ||
					resolve_x(x, pOrgs, nOrgs, *cnf, (*ot)[p], (*ot)[n], units, resolved, &outs[tx * SH_MAX_BVE_OUT])) x |= MELTING_MASK;
				tid += stride_x();
			}
		}

		__global__ void in_ve_k(CNF* cnf, OT* ot, cuVecU* pVars, cuVecU* units, cuVecU* resolved)
		{
			uint32 tx = threadIdx.x;
			uint32 tid = global_tx();
			__shared__ uint32 defs[BLVE * FAN_LMT];
			__shared__ uint32 outs[BLVE * SH_MAX_BVE_OUT];
			while (tid < pVars->size()) {
				uint32& x = (*pVars)[tid];
				assert(x);
				assert(!ELIMINATED(x));
				uint32 p = V2D(x), n = NEG(p);
				uint32 pOrgs = 0, nOrgs = 0;
				countOrgs(*cnf, (*ot)[p], pOrgs), countOrgs(*cnf, (*ot)[n], nOrgs);
				// try pure-literal elimination
				if (!pOrgs || !nOrgs) {
					uint32 psLits = 0, nsLits = 0;
					countLitsBefore(*cnf, (*ot)[p], psLits);
					countLitsBefore(*cnf, (*ot)[n], nsLits);
					toblivion(*cnf, (*ot)[p], (*ot)[n], resolved, pOrgs, nOrgs, psLits, nsLits, p);
					x |= MELTING_MASK;
				}
				// try simple resolution
				else if ((pOrgs == 1 || nOrgs == 1) &&
					resolve_x(x, pOrgs, nOrgs, *cnf, (*ot)[p], (*ot)[n], units, resolved, &outs[tx * SH_MAX_BVE_OUT])) x |= MELTING_MASK;
				// try gate-equivalence reasoning, otherwise resolution
				else if (gateReasoning_x(p, pOrgs, nOrgs, *cnf, (*ot)[p], (*ot)[n], units, resolved, &defs[tx * FAN_LMT], &outs[tx * SH_MAX_BVE_OUT]) ||
					resolve_x(x, pOrgs, nOrgs , *cnf, (*ot)[p], (*ot)[n], units, resolved, &outs[tx * SH_MAX_BVE_OUT])) x |= MELTING_MASK;
				tid += stride_x();
			}
		}

		__global__ void hse_k(CNF* cnf, OT* ot, cuVecU* pVars, cuVecU* units, uint32 limit)
		{
			uint32 tid = global_tx();
			__shared__ uint32 sh_cls[BLHSE * SH_MAX_HSE_IN];
			while (tid < pVars->size()) {
				uint32 x = (*pVars)[tid];
				assert(x);
				assert(!ELIMINATED(x));
				uint32 p = V2D(x), n = NEG(p);
				if ((*ot)[p].size() <= limit && (*ot)[n].size() <= limit)
					self_sub_x(p, *cnf, (*ot)[p], (*ot)[n], units, &sh_cls[threadIdx.x * SH_MAX_HSE_IN]);
				tid += stride_x();
			}
		}

		__global__ void bce_k(CNF* cnf, OT* ot, cuVecU* pVars, uint32 limit)
		{
			uint32 tid = global_tx();
			__shared__ uint32 sh_cls[BLBCE * SH_MAX_BCE_IN];
			while (tid < pVars->size()) {
				uint32 x = (*pVars)[tid];
				assert(x);
				assert(!ELIMINATED(x));
				uint32 p = V2D(x), n = NEG(p);
				if ((*ot)[p].size() <= limit && (*ot)[n].size() <= limit)
					blocked_x(x, *cnf, (*ot)[p], (*ot)[n], &sh_cls[threadIdx.x * SH_MAX_BCE_IN]);
				tid += stride_x();
			}
		}

		__global__ void hre_k(CNF* cnf, OT* ot, cuVecU* pVars, uint32 limit)
		{
			uint32 gid = global_ty();
			uint32* smem = SharedMemory<uint32>();
			uint32* m_c = smem + warpSize * SH_MAX_HRE_IN + threadIdx.y * SH_MAX_HRE_OUT; // shared memory for resolvent
			while (gid < pVars->size()) {
				assert(pVars->at(gid));
				assert(!ELIMINATED(pVars->at(gid)));
				uint32 p = V2D(pVars->at(gid)), n = NEG(p);
				// do merging and apply forward equality check (on-the-fly) over resolvents
				if ((*ot)[p].size() <= limit && (*ot)[n].size() <= limit) {
					for (uint32* i = (*ot)[p]; i != (*ot)[p].end(); i++) {
						SCLAUSE& pos = (*cnf)[*i];
						if (pos.deleted() || pos.learnt()) continue;
						if (pos.size() <= SH_MAX_HRE_IN) { // use shared memory for positives
							uint32* sh_pos = smem + threadIdx.y * SH_MAX_HRE_IN;
							if (threadIdx.x == 0) pos.shareTo(sh_pos);
							for (uint32* j = (*ot)[n]; j != (*ot)[n].end(); j++) {
								SCLAUSE& neg = (*cnf)[*j];
								if (neg.deleted() || neg.learnt() || (pos.size() + neg.size() - 2) > SH_MAX_HRE_OUT) continue;
								int m_len = 0;
								if (threadIdx.x == 0) m_len = merge(pVars->at(gid), sh_pos, pos.size(), neg, m_c);
								m_len = __shfl_sync(FULLWARP, m_len, 0);
								if (m_len) forward_equ(*cnf, *ot, m_c, m_len);
							}
						}
						else { // use global memory
							for (uint32* j = (*ot)[n]; j != (*ot)[n].end(); j++) {
								SCLAUSE& neg = (*cnf)[*j];
								if (neg.deleted() || neg.learnt() || (pos.size() + neg.size() - 2) > SH_MAX_HRE_OUT) continue;
								int m_len = 0;
								if (threadIdx.x == 0) m_len = merge(pVars->at(gid), pos, neg, m_c);
								m_len = __shfl_sync(FULLWARP, m_len, 0);
								if (m_len) forward_equ(*cnf, *ot, m_c, m_len);
							}
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
		void copyIf(uint32* dest, CNF* src, GSTATS* gstats)
		{
			reset_stats << <1, 1 >> > (gstats);
			uint32 nBlocks = MIN((inf.nClauses + BLOCK1D - 1) / BLOCK1D, maxGPUThreads / BLOCK1D);
			copy_if_k << <nBlocks, BLOCK1D >> > (dest, src, gstats);
			LOGERR("Copying literals failed");
			CHECK(cudaDeviceSynchronize());
		}
		void shrinkSimp(CNF* dest, CNF* src)
		{
			uint32 cnf_size = inf.nClauses + (inf.maxAddedCls >> 1);
			uint32 nBlocks = MIN((cnf_size + BLOCK1D - 1) / BLOCK1D, maxGPUThreads / BLOCK1D);
			copy_if_k << <nBlocks, BLOCK1D >> > (dest, src);
			LOGERR("Copying CNF failed");
			CHECK(cudaDeviceSynchronize());
		}
		void calcScores(VARS* vars, uint32* hist)
		{
			uint32 nBlocks = MIN((inf.maxVar + BLOCK1D - 1) / BLOCK1D, maxGPUThreads / BLOCK1D);
			assign_scores << <nBlocks, BLOCK1D >> > (vars->eligible, vars->scores, hist, inf.maxVar);
			LOGERR("Assigning scores failed");
			CHECK(cudaDeviceSynchronize());
		}
		void calcScores(VARS* vars, uint32* hist, OT* ot)
		{
			uint32 nBlocks = MIN((inf.maxVar + BLOCK1D - 1) / BLOCK1D, maxGPUThreads / BLOCK1D);
			assign_scores << <nBlocks, BLOCK1D >> > (vars->eligible, vars->scores, hist, ot, inf.maxVar);
			LOGERR("Assigning scores failed");
			CHECK(cudaDeviceSynchronize());
		}
		void calcSigCNFAsync(CNF* cnf, const uint32& offset, const uint32& size, const cudaStream_t& _s)
		{
			assert(size);
			uint32 nBlocks = MIN((size + BLOCK1D - 1) / BLOCK1D, maxGPUThreads / BLOCK1D);
			calc_sig_k << <nBlocks, BLOCK1D, 0, _s >> > (cnf, offset, size);
		}
		void calcSigCNF(CNF* cnf, const uint32& size)
		{
			assert(size);
			uint32 nBlocks = MIN((size + BLOCK1D - 1) / BLOCK1D, maxGPUThreads / BLOCK1D);
			calc_sig_k << <nBlocks, BLOCK1D >> > (cnf, 0, size);
			LOGERR("Signature calculation failed");
			CHECK(cudaDeviceSynchronize());
		}
		void reduceOTAsync(CNF* cnf, OT* ot, const bool& p)
		{
			assert(cnf != NULL);
			assert(ot != NULL);
			uint32 nBlocks = MIN(uint32((inf.nDualVars + BLOCK1D - 1) / BLOCK1D), maxGPUThreads / BLOCK1D);
			reduce_ot << <nBlocks, BLOCK1D >> > (cnf, ot);
			if (p) {
				LOGERR("Occurrence table reduction failed");
				CHECK(cudaDeviceSynchronize());
				PFLOGR('=', 30);
				PFLOG0("\toccurrence table");
				ot->print();
				PFLOGR('=', 30);
			}
		}
		void reduceOT(CNF* cnf, OT* ot, VARS* vars, cudaStream_t* streams, const bool& p)
		{
			assert(vars->numPVs);
			assert(cnf != NULL);
			assert(ot != NULL);
			uint32 nBlocks = MIN((vars->numPVs + BLOCK1D - 1) / BLOCK1D, maxGPUThreads / BLOCK1D);
			reduce_ot_p << <nBlocks, BLOCK1D, 0, streams[0] >> > (cnf, ot, vars->pVars);
			reduce_ot_n << <nBlocks, BLOCK1D, 0, streams[1] >> > (cnf, ot, vars->pVars);
			LOGERR("Occurrence table reduction failed");
			CHECK(cudaDeviceSynchronize());
			assert(ot->accViolation());
			if (p) {	
				PFLOGR('=', 30);
				PFLOG0("\toccurrence table");
				ot->print();
				PFLOGR('=', 30);
			}
		}
		void createOT(CNF* cnf, OT* ot, const bool& p)
		{
			assert(cnf != NULL);
			assert(ot != NULL);
			uint32 cnf_size = inf.nClauses + (inf.maxAddedCls >> 1);
			uint32 rstGridSize = MIN(uint32((inf.nDualVars + BLOCK1D - 1) / BLOCK1D), maxGPUThreads / BLOCK1D);
			uint32 otGridSize = MIN((cnf_size + BLOCK1D - 1) / BLOCK1D, maxGPUThreads / BLOCK1D);
			reset_ot_k << <rstGridSize, BLOCK1D >> > (ot);
			create_ot_k << <otGridSize, BLOCK1D >> > (cnf, ot);
			LOGERR("Occurrence table creation failed");
			CHECK(cudaDeviceSynchronize());
			assert(ot->accViolation());
			if (p) {
				PFLOGR('=', 30);
				PFLOG0("\toccurrence table");
				ot->print();
				PFLOGR('=', 30);
			}
		}
		void createOTAsync(CNF* cnf, OT* ot, const bool& p)
		{
			assert(cnf != NULL);
			assert(ot != NULL);
			uint32 cnf_size = inf.nClauses + (inf.maxAddedCls >> 1);
			uint32 rstGridSize = MIN(uint32((inf.nDualVars + BLOCK1D - 1) / BLOCK1D), maxGPUThreads / BLOCK1D);
			uint32 otGridSize = MIN((cnf_size + BLOCK1D - 1) / BLOCK1D, maxGPUThreads / BLOCK1D);
			reset_ot_k << <rstGridSize, BLOCK1D >> > (ot);
			create_ot_k << <otGridSize, BLOCK1D >> > (cnf, ot);
			if (p) {
				LOGERR("Occurrence table creation failed");
				CHECK(cudaDeviceSynchronize());
				assert(ot->accViolation());
				PFLOGR('=', 30);
				PFLOG0("\toccurrence table");
				ot->print();
				PFLOGR('=', 30);
			}
		}
		void ve(CNF* cnf, OT* ot, VARS* vars, const bool& in)
		{
			assert(vars->numPVs);
#if VE_DBG
			putchar('\n');
			if (in) in_ve_k << <1, 1 >> > (cnf, ot, vars->pVars, vars->units, vars->resolved);
			else	   ve_k << <1, 1 >> > (cnf, ot, vars->pVars, vars->units, vars->resolved);
#else
			uint32 nBlocks = MIN((vars->numPVs + BLVE - 1) / BLVE, maxGPUThreads / BLVE);
			if (in)	in_ve_k << <nBlocks, BLVE >> > (cnf, ot, vars->pVars, vars->units, vars->resolved);
			else	   ve_k << <nBlocks, BLVE >> > (cnf, ot, vars->pVars, vars->units, vars->resolved);
#endif
			LOGERR("Parallel BVE failed");
			CHECK(cudaDeviceSynchronize());
		}
		void hse(CNF* cnf, OT* ot, VARS* vars, const uint32& limit)
		{
			assert(vars->numPVs);
			assert(limit);
#if SS_DBG
			putchar('\n');
			hse_k << <1, 1 >> > (cnf, ot, vars->pVars, vars->units, limit);
#else
			uint32 nBlocks = MIN((vars->numPVs + BLHSE - 1) / BLHSE, maxGPUThreads / BLHSE);
			hse_k << <nBlocks, BLHSE >> > (cnf, ot, vars->pVars, vars->units, limit);
#endif
			LOGERR("Parallel HSE failed");
			CHECK(cudaDeviceSynchronize());
		}
		void bce(CNF* cnf, OT* ot, VARS* vars, const uint32& limit)
		{
			assert(vars->numPVs);
			assert(limit);
			uint32 nBlocks = MIN((vars->numPVs + BLBCE - 1) / BLBCE, maxGPUThreads / BLBCE);
			bce_k << <nBlocks, BLBCE >> > (cnf, ot, vars->pVars, limit);
			LOGERR("Parallel BCE failed");
			CHECK(cudaDeviceSynchronize());
		}
		void hre(CNF* cnf, OT* ot, VARS* vars, const uint32& limit)
		{
			assert(vars->numPVs);
			assert(limit);
			dim3 block2D(devProp.warpSize, devProp.warpSize), grid2D(1, 1, 1);
			grid2D.y = MIN((vars->numPVs + block2D.y - 1) / block2D.y, maxGPUThreads / block2D.y);
			uint32 smemSize = devProp.warpSize * (SH_MAX_HRE_IN + SH_MAX_HRE_OUT) * sizeof(uint32);
			hre_k << <grid2D, block2D, smemSize >> > (cnf, ot, vars->pVars, limit);
			LOGERR("HRE Elimination failed");
			CHECK(cudaDeviceSynchronize());
		}
		void countMelted(LIT_ST* vstate)
		{
			inf.n_del_vars_after = 0;
			for (uint32 v = 1; v <= inf.maxVar; v++) 
				if (vstate[v] == MELTED)
					inf.n_del_vars_after++;
			assert(inf.n_del_vars_after >= inf.maxMelted);
			inf.n_del_vars_after -= inf.maxMelted;
			inf.maxMelted += inf.n_del_vars_after;
		}
		void evalReds(CNF* cnf, GSTATS* gstats, LIT_ST* vstate)
		{
			reset_stats << <1, 1 >> > (gstats);
			uint32 cnf_sz = inf.nClauses + (inf.maxAddedCls >> 1); // approximate cnf size 
			uint32 nBlocks1 = MIN((cnf_sz + (BLOCK1D << 1) - 1) / (BLOCK1D << 1), maxGPUThreads / (BLOCK1D << 1));
			uint32 smemSize1 = BLOCK1D * sizeof(uint32) * 2;
			cnt_reds << <nBlocks1, BLOCK1D, smemSize1 >> > (cnf, gstats);
			countMelted(vstate);
			LOGERR("Counting reductions failed");
			CHECK(cudaDeviceSynchronize());
			GSTATS hstats;
			CHECK(cudaMemcpy(&hstats, gstats, sizeof(GSTATS), cudaMemcpyDeviceToHost)); // avoids unified memory migration on large scale
			inf.n_cls_after = hstats.numClauses;
			inf.n_lits_after = hstats.numLits;
		}
		void countFinal(CNF* cnf, GSTATS* gstats, LIT_ST* vstate)
		{
			reset_stats << <1, 1 >> > (gstats);
			uint32 cnf_sz = inf.nClauses + (inf.maxAddedCls >> 1);
			uint32 nBlocks = MIN((cnf_sz + (BLOCK1D << 1) - 1) / (BLOCK1D << 1), maxGPUThreads / (BLOCK1D << 1));
			uint32 smemSize = BLOCK1D * sizeof(uint32);
			cnt_cls << <nBlocks, BLOCK1D, smemSize >> > (cnf, gstats);
			countMelted(vstate);
			LOGERR("Counting clauses failed");
			CHECK(cudaDeviceSynchronize());
			GSTATS hstats;
			CHECK(cudaMemcpy(&hstats, gstats, sizeof(GSTATS), cudaMemcpyDeviceToHost));
			inf.n_cls_after = hstats.numClauses;
		}
		void countCls(CNF* cnf, GSTATS* gstats)
		{
			reset_stats << <1, 1>> > (gstats);
			uint32 cnf_sz = inf.nClauses + (inf.maxAddedCls >> 1);
			uint32 nBlocks = MIN((cnf_sz + (BLOCK1D << 1) - 1) / (BLOCK1D << 1), maxGPUThreads / (BLOCK1D << 1));
			uint32 smemSize = BLOCK1D * sizeof(uint32);
			cnt_cls << <nBlocks, BLOCK1D, smemSize>> > (cnf, gstats);
			LOGERR("Counting clauses failed");
			CHECK(cudaDeviceSynchronize());
			GSTATS hstats;
			CHECK(cudaMemcpy(&hstats, gstats, sizeof(GSTATS), cudaMemcpyDeviceToHost));
			inf.n_cls_after = hstats.numClauses;
		}
		void countLits(CNF* cnf, GSTATS* gstats)
		{
			reset_stats << <1, 1 >> > (gstats);
			uint32 cnf_sz = inf.nClauses + (inf.maxAddedCls >> 1);
			uint32 nBlocks = MIN((cnf_sz + (BLOCK1D << 1) - 1) / (BLOCK1D << 1), maxGPUThreads / (BLOCK1D << 1));
			uint32 smemSize = BLOCK1D * sizeof(uint32);
			cnt_lits << <nBlocks, BLOCK1D, smemSize >> > (cnf, gstats);
			LOGERR("Counting literals failed");
			CHECK(cudaDeviceSynchronize());
			GSTATS hstats;
			CHECK(cudaMemcpy(&hstats, gstats, sizeof(GSTATS), cudaMemcpyDeviceToHost));
			inf.n_lits_after = hstats.numLits;
		}
		void countAll(CNF* cnf, GSTATS* gstats)
		{
			reset_stats << <1, 1 >> > (gstats);
			uint32 cnf_sz = inf.nClauses + (inf.maxAddedCls >> 1);
			uint32 nBlocks = MIN((cnf_sz + (BLOCK1D << 1) - 1) / (BLOCK1D << 1), maxGPUThreads / (BLOCK1D << 1));
			uint32 smemSize = BLOCK1D * (sizeof(uint32) + sizeof(uint32));
			cnt_cls_lits << <nBlocks, BLOCK1D, smemSize >> > (cnf, gstats);
			LOGERR("Counting clauses-literals failed");
			CHECK(cudaDeviceSynchronize());
			GSTATS hstats;
			CHECK(cudaMemcpy(&hstats, gstats, sizeof(GSTATS), cudaMemcpyDeviceToHost));
			inf.n_cls_after = hstats.numClauses;
			inf.n_lits_after = hstats.numLits;
		}

	}
}