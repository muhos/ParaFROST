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

#include "pfsimplify.cuh"
#include "pfmemory.cuh"
#include "pfdevice.cuh"
#include "pfbve.cuh"
#include "pfhse.cuh"
#include "pfere.cuh"
#include <cub/device/device_scan.cuh>

using namespace cub;

namespace pFROST {

	namespace SIGmA {

		template<class T>
		__global__ void memset_k(T* mem, T val, size_t size)
		{
			size_t tid = global_tx();
			while (tid < size) { mem[tid] = val; tid += stride_x(); }
		}

		__global__ void reset_stats(GSTATS* gstats) { gstats->numDelVars = 0, gstats->numClauses = 0, gstats->numLits = 0; }

		__global__ void prep_cnf_k(CNF* cnf)
		{
			uint32 tid = global_tx();
			while (tid < cnf->size()) { 
				SCLAUSE& c = cnf->clause(tid);
				devSort(c.data(), c.size());
				calcSig(c); 
				tid += stride_x();
			}
		}

		__global__ void reset_ot_k(OT* ot)
		{
			uint32 tid = global_tx();
			while (tid < ot->size()) { (*ot)[tid].clear(); tid += stride_x(); }
		}

		__global__ void reduce_ot(const CNF* __restrict__ cnf, OT* __restrict__ ot)
		{
			uint32 tid = global_tx();
			while (tid < ot->size()) { reduceOL(*cnf, (*ot)[tid]); tid += stride_x(); }
		}

		__global__ void sort_ot_p(const CNF* __restrict__ cnf, OT* __restrict__ ot, const cuVecU* __restrict__ pVars)
		{
			uint32 tid = global_tx();
			while (tid < pVars->size()) {
				const uint32 x = pVars->at(tid), p = V2L(x);
				assert(x);
				OL& ol = (*ot)[p];
				devSort(ol.data(), ol.size(), CNF_CMP_KEY(cnf));
				tid += stride_x();
			}
		}

		__global__ void sort_ot_n(const CNF* __restrict__ cnf, OT* __restrict__ ot, const cuVecU* __restrict__ pVars)
		{
			uint32 tid = global_tx();
			while (tid < pVars->size()) {
				const uint32 x = pVars->at(tid), n = NEG(V2L(x));
				assert(x);
				OL& ol = (*ot)[n];
				devSort(ol.data(), ol.size(), CNF_CMP_KEY(cnf));
				tid += stride_x();
			}
		}

		__global__ void create_ot_k(CNF* __restrict__ cnf, OT* __restrict__ ot)
		{
			uint32 tid = global_tx();
			while (tid < cnf->size()) {
				const S_REF r = cnf->ref(tid);
				SCLAUSE& c = (*cnf)[r];
				if (c.original() || c.learnt()) {
#pragma unroll
					forall_clause(c, lit) (*ot)[*lit].insert(r);
				}
				tid += stride_x();
			}
		}

		__global__ void assign_scores(uint32* __restrict__ eligible, uint32* __restrict__ scores, const uint32* __restrict__ hist, uint32 size)
		{
			uint32 tid = global_tx();
			while (tid < size) {
				const uint32 v = tid + 1;
				const uint32 p = V2L(v), ps = hist[p], ns = hist[NEG(p)];
				eligible[tid] = v;
				scores[v] = ps * ns;
				tid += stride_x();
			}
		}

		__global__ void assign_scores(uint32* __restrict__ eligible, uint32* __restrict__ scores, uint32* __restrict__ hist, const OT* __restrict__ ot, uint32 size)
		{
			uint32 tid = global_tx();
			while (tid < size) {
				const uint32 v = tid + 1;
				const uint32 p = V2L(v), n = NEG(p), ps = (*ot)[p].size(), ns = (*ot)[n].size();
				hist[p] = ps, hist[n] = ns;
				eligible[tid] = v;
				scores[v] = ps * ns;
				tid += stride_x();
			}
		}

		__global__ void copy_if_k(uint32* __restrict__ dest, CNF* __restrict__ src, GSTATS* __restrict__ gstats)
		{
			uint32 tid = global_tx();
			while (tid < src->size()) {
				SCLAUSE& c = src->clause(tid);
				if (c.original() || c.learnt()) {
					uint32* d = dest + atomicAdd(&gstats->numLits, c.size());
#pragma unroll
					forall_clause(c, s) { *d++ = *s; }
				}
				tid += stride_x();
			}
		}

		__global__ void cnt_reds(const CNF* __restrict__ cnf, GSTATS* __restrict__ gstats)
		{
			uint32* sh_rCls = SharedMemory<uint32>();
			uint32* sh_rLits = sh_rCls + blockDim.x;
			uint32 tid = global_tx_off();
			uint32 nCls = 0;
			uint32 nLits = 0;
			while (tid < cnf->size()) {
				const SCLAUSE& c1 = cnf->clause(tid);
				if (c1.original() || c1.learnt())
					nCls++, nLits += c1.size();
				if (tid + blockDim.x < cnf->size()) {
					const SCLAUSE& c2 = cnf->clause(tid + blockDim.x);
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

		__global__ void cnt_cls(const CNF* __restrict__ cnf, GSTATS* __restrict__ gstats)
		{
			uint32* sh_rCls = SharedMemory<uint32>();
			uint32 tid = global_tx_off();
			uint32 nCls = 0;
			while (tid < cnf->size()) {
				const SCLAUSE& c1 = cnf->clause(tid);
				if (c1.original() || c1.learnt()) nCls++;
				uint32 off = tid + blockDim.x;
				if (off < cnf->size()) {
					const SCLAUSE& c2 = cnf->clause(off);
					if (c2.original() || c2.learnt()) nCls++;
				}
				tid += stride_x_off();
			}
			loadShared(sh_rCls, nCls, cnf->size());
			sharedReduce(sh_rCls, nCls);
			warpReduce(sh_rCls, nCls);
			if (threadIdx.x == 0) atomicAdd(&gstats->numClauses, nCls);
		}

		__global__ void cnt_lits(const CNF* __restrict__ cnf, GSTATS* __restrict__ gstats)
		{
			uint32* sh_rLits = SharedMemory<uint32>();
			uint32 tid = global_tx_off();
			uint32 nLits = 0;
			while (tid < cnf->size()) {
				const SCLAUSE& c1 = cnf->clause(tid);
				if (c1.original() || c1.learnt()) nLits += c1.size();
				uint32 off = tid + blockDim.x;
				if (off < cnf->size()) {
					const SCLAUSE& c2 = cnf->clause(off);
					if (c2.original() || c2.learnt()) nLits += c2.size();
				}
				tid += stride_x_off();
			}
			loadShared(sh_rLits, nLits, cnf->size());
			sharedReduce(sh_rLits, nLits);
			warpReduce(sh_rLits, nLits);
			if (threadIdx.x == 0) atomicAdd(&gstats->numLits, nLits);
		}

		__global__ void cnt_cls_lits(const CNF* __restrict__ cnf, GSTATS* __restrict__ gstats)
		{
			uint32* sh_rCls = SharedMemory<uint32>();
			uint32* sh_rLits = sh_rCls + blockDim.x;
			uint32 tid = global_tx_off();
			uint32 nCls = 0;
			uint32 nLits = 0;
			while (tid < cnf->size()) {
				const SCLAUSE& c1 = cnf->clause(tid);
				if (c1.original() || c1.learnt()) nCls++, nLits += c1.size();
				uint32 off = tid + blockDim.x;
				if (off < cnf->size()) {
					const SCLAUSE& c2 = cnf->clause(off);
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

		__global__ void ve_k(CNF* __restrict__ cnfptr, OT* __restrict__ otptr, cuVecU* __restrict__ pVars, cuVecU* __restrict__ units, cuVecU* __restrict__ resolved, const uint32* __restrict__ vorg)
		{
			uint32 tid = global_tx();
			__shared__ uint32 outs[BLVE * SH_MAX_BVE_OUT];
			while (tid < pVars->size()) {
				uint32& x = (*pVars)[tid];
				assert(x);
				assert(!ELIMINATED(x));
				const uint32 p = V2L(x), n = NEG(p);
				CNF& cnf = *cnfptr;
				OT& ot = *otptr;
				OL& poss = ot[p], &negs = ot[n];
				uint32 pOrgs = 0, nOrgs = 0;
				countOrgs(cnf, poss, pOrgs), countOrgs(cnf, negs, nOrgs);
				bool elim = false;
				// pure-literal elimination
				if (!pOrgs || !nOrgs) {
					toblivion(p, vorg, pOrgs, nOrgs, cnf, poss, negs, resolved), elim = true;
				}
				// Equiv/NOT-gate Reasoning
				else if (uint32 def = find_equ_gate(p, cnf, poss, negs)) {
					saveResolved(p, vorg, pOrgs, nOrgs, cnf, poss, negs, resolved); // must be called before substitution
					substitute_single(p, def, cnf, poss, negs, units), elim = true;
				}
				else {
					assert(pOrgs && nOrgs);
					const uint32 nClsBefore = pOrgs + nOrgs;
					uint32 *shared_outs = outs + threadIdx.x * SH_MAX_BVE_OUT;
					uint32 nAddedCls, nAddedLits;
					// AND-gate Reasoning
					if (nOrgs < SH_MAX_BVE_OUT && find_ao_gate(n, nClsBefore, cnf, ot, shared_outs, nAddedCls, nAddedLits)) {
						if (nAddedCls) substitute_x(x, nAddedCls, nAddedLits, cnf, poss, negs, units, shared_outs);
						toblivion(p, vorg, pOrgs, nOrgs, cnf, poss, negs, resolved), elim = true;
					}
					// OR-gate Reasoning
					else if (pOrgs < SH_MAX_BVE_OUT && find_ao_gate(p, nClsBefore, cnf, ot, shared_outs, nAddedCls, nAddedLits)) {
						if (nAddedCls) substitute_x(x, nAddedCls, nAddedLits, cnf, poss, negs, units, shared_outs);
						toblivion(p, vorg, pOrgs, nOrgs, cnf, poss, negs, resolved), elim = true;
					}
					// ITE-gate Reasoning
					else if (find_ite_gate(p, nClsBefore, cnf, ot, nAddedCls, nAddedLits)
						|| find_ite_gate(n, nClsBefore, cnf, ot, nAddedCls, nAddedLits)) {
						if (nAddedCls) substitute_x(x, nAddedCls, nAddedLits, cnf, poss, negs, units, shared_outs);
						toblivion(p, vorg, pOrgs, nOrgs, cnf, poss, negs, resolved), elim = true;
					}
					// XOR-gate Reasoning
					else if (find_xor_gate(p, nClsBefore, cnf, ot, shared_outs, nAddedCls, nAddedLits)) {
						if (nAddedCls) substitute_x(x, nAddedCls, nAddedLits, cnf, poss, negs, units, shared_outs);
						toblivion(p, vorg, pOrgs, nOrgs, cnf, poss, negs, resolved), elim = true;
					}
					// n-by-m resolution
					else if (resolve(x, nClsBefore, cnf, poss, negs, nAddedCls, nAddedLits)) {
						if (nAddedCls) resolve_x(x, nAddedCls, nAddedLits, cnf, poss, negs, units, shared_outs);
						toblivion(p, vorg, pOrgs, nOrgs, cnf, poss, negs, resolved);
						elim = true;
					}
				}
				if (elim) x |= MELTING_MASK;
				tid += stride_x();
			}
		}

		__global__ void in_ve_k_1(CNF* __restrict__ cnfptr, OT* __restrict__ otptr, cuVecU* __restrict__ pVars, cuVecU* __restrict__ units, cuVecU* __restrict__ resolved, const uint32* __restrict__ vorg, uint32* __restrict__ type, uint32* __restrict__ rpos, S_REF* __restrict__ rref)
		{
			uint32 tid = global_tx();
			uint32* outs = SharedMemory<uint32>();
			while (tid < pVars->size()) {
				uint32& x = (*pVars)[tid];
				assert(x);
				assert(!ELIMINATED(x));
				const uint32 p = V2L(x), n = NEG(p);
				CNF& cnf = *cnfptr;
				OT& ot = *otptr;
				OL& poss = ot[p], & negs = ot[n];
				uint32 pOrgs = 0, nOrgs = 0;
				countOrgs(cnf, poss, pOrgs), countOrgs(cnf, negs, nOrgs);
				// pure-literal elimination
				if (!pOrgs || !nOrgs) {
					toblivion(p, vorg, pOrgs, nOrgs, cnf, poss, negs, resolved);
					type[tid] = 0, rref[tid] = 0, rpos[tid] = 0, x |= MELTING_MASK;
				}
				// Equiv/NOT-gate Reasoning
				else if (uint32 def = find_equ_gate(p, cnf, poss, negs)) {
					saveResolved(p, vorg, pOrgs, nOrgs, cnf, poss, negs, resolved); // must be called before substitution
					substitute_single(p, def, cnf, poss, negs, units);
					type[tid] = 0, rref[tid] = 0, rpos[tid] = 0, x |= MELTING_MASK;
				}
				else {
					assert(pOrgs && nOrgs);
					const uint32 nClsBefore = pOrgs + nOrgs;
					uint32* shared_outs = outs + threadIdx.x * SH_MAX_BVE_OUT1;
					uint32 elimType = 0, nAddedCls = 0, nAddedLits = 0;
					//=====================
					// check resolvability 
					//=====================
					// AND/OR-gate Reasoning
					if ((nOrgs < SH_MAX_BVE_OUT1 && find_ao_gate(n, nClsBefore, cnf, ot, shared_outs, nAddedCls, nAddedLits))
					||  (pOrgs < SH_MAX_BVE_OUT1 && find_ao_gate(p, nClsBefore, cnf, ot, shared_outs, nAddedCls, nAddedLits))) elimType = AOIX_MASK;
					// ITE-gate Reasoning
					else if (find_ite_gate(p, nClsBefore, cnf, ot, nAddedCls, nAddedLits)
						||	 find_ite_gate(n, nClsBefore, cnf, ot, nAddedCls, nAddedLits)) elimType = AOIX_MASK;
					// XOR-gate Reasoning
					else if (find_xor_gate(p, nClsBefore, cnf, ot, shared_outs, nAddedCls, nAddedLits)) elimType = AOIX_MASK;
					// n-by-m resolution
					else if (!nAddedCls && resolve(x, nClsBefore, cnf, poss, negs, nAddedCls, nAddedLits)) elimType = RES_MASK;
					//=====================
					// check addibility 
					//=====================
					if (!nAddedCls) { // eliminated without resolvents
						toblivion(p, vorg, pOrgs, nOrgs, cnf, poss, negs, resolved);
						type[tid] = 0, rref[tid] = 0, rpos[tid] = 0, x |= MELTING_MASK;
					}
					else if (elimType) { // can be eliminated with resolvents in next phase
						assert(nAddedLits >= nAddedCls);
						assert(elimType < TYPE_MASK);
						assert(nAddedCls <= ADDEDCLS_MAX);
						assert(nAddedLits <= ADDEDLITS_MAX);
						// save elimination info.
						type[tid] = ENCODEVARINFO(elimType, nAddedCls, nAddedLits);
						rpos[tid] = nAddedCls, rref[tid] = nAddedLits + dc_nbuckets * nAddedCls;			
					}
					else  // cannot be eliminated
						type[tid] = 0, rref[tid] = 0, rpos[tid] = 0;
				}
				tid += stride_x();
			}
		}

		__global__ void ve_k_1(CNF* __restrict__ cnfptr, OT* __restrict__ otptr, cuVecU* __restrict__ pVars, cuVecU* __restrict__ units, cuVecU* __restrict__ resolved, const uint32* __restrict__ vorg, uint32* __restrict__ type, uint32* __restrict__ rpos, S_REF* __restrict__ rref)
		{
			uint32 tid = global_tx();
			uint32* outs = SharedMemory<uint32>();
			while (tid < pVars->size()) {
				uint32& x = (*pVars)[tid];
				assert(x);
				assert(!ELIMINATED(x));
				const uint32 p = V2L(x), n = NEG(p);
				CNF& cnf = *cnfptr;
				OT& ot = *otptr;
				OL& poss = ot[p], & negs = ot[n];
				const uint32 pOrgs = poss.size(), nOrgs = negs.size();
				// pure-literal elimination
				if (!pOrgs || !nOrgs) {
					toblivion(p, vorg, pOrgs, nOrgs, cnf, poss, negs, resolved);
					type[tid] = 0, rref[tid] = 0, rpos[tid] = 0, x |= MELTING_MASK;
				}
				// Equiv/NOT-gate Reasoning
				else if (uint32 def = find_equ_gate(p, cnf, poss, negs)) {
					saveResolved(p, vorg, pOrgs, nOrgs, cnf, poss, negs, resolved); // must be called before substitution
					substitute_single(p, def, cnf, poss, negs, units);
					type[tid] = 0, rref[tid] = 0, rpos[tid] = 0, x |= MELTING_MASK;
				}
				else {
					assert(pOrgs && nOrgs);
					const uint32 nClsBefore = pOrgs + nOrgs;
					uint32* shared_outs = outs + threadIdx.x * SH_MAX_BVE_OUT1;
					uint32 elimType = 0, nAddedCls = 0, nAddedLits = 0;
					//=====================
					// check resolvability 
					//=====================
					// AND/OR-gate Reasoning
					if ((nOrgs < SH_MAX_BVE_OUT1 && find_ao_gate(n, nClsBefore, cnf, ot, shared_outs, nAddedCls, nAddedLits))
					||  (pOrgs < SH_MAX_BVE_OUT1 && find_ao_gate(p, nClsBefore, cnf, ot, shared_outs, nAddedCls, nAddedLits))) elimType = AOIX_MASK;
					// ITE-gate Reasoning
					else if (find_ite_gate(p, nClsBefore, cnf, ot, nAddedCls, nAddedLits) 
						  || find_ite_gate(n, nClsBefore, cnf, ot, nAddedCls, nAddedLits)) elimType = AOIX_MASK;
					// XOR-gate Reasoning
					else if (find_xor_gate(p, nClsBefore, cnf, ot, shared_outs, nAddedCls, nAddedLits)) elimType = AOIX_MASK;
					// n-by-m resolution
					else if (!nAddedCls && resolve(x, nClsBefore, cnf, poss, negs, nAddedCls, nAddedLits)) elimType = RES_MASK;
					//=====================
					// check addibility 
					//=====================
					if (!nAddedCls) { // eliminated without resolvents
						toblivion(p, vorg, pOrgs, nOrgs, cnf, poss, negs, resolved);
						type[tid] = 0, rref[tid] = 0, rpos[tid] = 0, x |= MELTING_MASK;
					}
					else if (elimType) { // can be eliminated with resolvents in next phase
						assert(nAddedLits >= nAddedCls);
						assert(elimType < TYPE_MASK);
						assert(nAddedCls <= ADDEDCLS_MAX);
						assert(nAddedLits <= ADDEDLITS_MAX);
						// save elimination info.
						type[tid] = ENCODEVARINFO(elimType, nAddedCls, nAddedLits);
						rpos[tid] = nAddedCls, rref[tid] = nAddedLits + dc_nbuckets * nAddedCls;
					}
					else  // cannot be eliminated
						type[tid] = 0, rref[tid] = 0, rpos[tid] = 0;
				}
				tid += stride_x();
			}
		}

		__global__ void ve_k_2(CNF* __restrict__ cnf, OT* __restrict__ ot, cuVecU* __restrict__ pVars, cuVecU* __restrict__ units, cuVecU* __restrict__ resolved, const uint32* __restrict__ vorg, const uint32* __restrict__ type, const uint32* __restrict__ rpos, const S_REF* __restrict__ rref)
		{
			uint32 tid = global_tx();
			uint32* outs = SharedMemory<uint32>();
			while (tid < pVars->size()) {
				uint32& x = (*pVars)[tid];
				assert(x);
				const uint32 xinfo = type[tid];
				const uint32 elimType = RECOVERTYPE(xinfo);
				assert(elimType < TYPE_MASK);
				if (elimType) {
					assert(!ELIMINATED(x));
					const uint32 p = V2L(x);
					const uint32 nAddedCls = RECOVERADDEDCLS(xinfo);
					const uint32 nAddedLits = RECOVERADDEDLITS(xinfo);
					assert(nAddedCls && nAddedCls <= ADDEDCLS_MAX);
					assert(nAddedLits && nAddedLits <= ADDEDLITS_MAX);
					const uint32 added_pos = rpos[tid];
					const S_REF added_ref = rref[tid];
					OL& poss = (*ot)[p], &negs = (*ot)[NEG(p)];
					if (IS_RES(elimType)) {
						if (memorySafe(tid, x, nAddedCls, nAddedLits, added_pos, added_ref, cnf)) {
							if (nAddedCls) saveResolved(p, vorg, *cnf, poss, negs, resolved);
							resolve_x(x, nAddedCls, nAddedLits, added_pos, added_ref, *cnf, poss, negs, units, outs + threadIdx.x * SH_MAX_BVE_OUT2);
							x |= MELTING_MASK;
						}
					}
					else {
						assert(IS_AOIX(elimType));
						if (memorySafe(tid, x, nAddedCls, nAddedLits, added_pos, added_ref, cnf)) {
							if (nAddedCls) saveResolved(p, vorg, *cnf, poss, negs, resolved);
							substitute_x(x, nAddedCls, nAddedLits, added_pos, added_ref, *cnf, poss, negs, units, outs + threadIdx.x * SH_MAX_BVE_OUT2);
							x |= MELTING_MASK;
						}
						else freezeClauses(*cnf, poss, negs);
					}
				}
				tid += stride_x();
			}
		}

		__global__ void roll_back_gates(CNF* __restrict__ cnf, OT* __restrict__ ot, cuVecU* __restrict__ pVars, uint32* __restrict__ type)
		{
			uint32 tid = global_tx();
			while (tid < pVars->size()) {
				const uint32 x = (*pVars)[tid];
				assert(x);
				const uint32 xinfo = type[tid];
				const uint32 elimType = RECOVERTYPE(xinfo);
				assert(elimType < TYPE_MASK);
				if (elimType && IS_AOIX(elimType)) {
					assert(!ELIMINATED(x));
					const uint32 p = V2L(x);
					freezeClauses(*cnf, (*ot)[p], (*ot)[NEG(p)]);
				}
				tid += stride_x();
			}
		}

		__global__ void hse_k(CNF* __restrict__ cnf, OT* __restrict__ ot, const cuVecU* __restrict__ pVars, cuVecU* __restrict__ units)
		{
			uint32 tid = global_tx();
			__shared__ uint32 sh_cls[BLHSE * SH_MAX_HSE_IN];
			while (tid < pVars->size()) {
				const uint32 x = (*pVars)[tid];
				assert(x);
				assert(!ELIMINATED(x));
				const uint32 p = V2L(x), n = NEG(p);
				if ((*ot)[p].size() <= dc_limits[0] && (*ot)[n].size() <= dc_limits[0])
					subsume_x(p, *cnf, (*ot)[p], (*ot)[n], units, sh_cls + threadIdx.x * SH_MAX_HSE_IN);
				tid += stride_x();
			}
		}

		__global__ void bce_k(CNF* __restrict__ cnf, OT* __restrict__ ot, const cuVecU* __restrict__ pVars, cuVecU* __restrict__ resolved, const uint32* __restrict__ vorg)
		{
			uint32 tid = global_tx();
			__shared__ uint32 sh_cls[BLBCE * SH_MAX_BCE_IN];
			while (tid < pVars->size()) {
				const uint32 x = (*pVars)[tid];
				assert(x);
				assert(!ELIMINATED(x));
				const uint32 p = V2L(x), n = NEG(p);
				if ((*ot)[p].size() <= dc_limits[1] && (*ot)[n].size() <= dc_limits[1])
					blocked_x(x, vorg, *cnf, (*ot)[p], (*ot)[n], resolved, sh_cls + threadIdx.x * SH_MAX_BCE_IN);
				tid += stride_x();
			}
		}

		__global__ void ere_k(CNF* __restrict__ cnf, OT* __restrict__ ot, const cuVecU* __restrict__ pVars)
		{
			uint32 gid = global_ty();
			uint32* smem = SharedMemory<uint32>();
			while (gid < pVars->size()) {
				const uint32 v = pVars->at(gid);
				assert(v);
				assert(!ELIMINATED(v));
				const uint32 p = V2L(v), n = NEG(p);
				OL& poss = (*ot)[p], & negs = (*ot)[n];
				// do merging and apply forward equality check (on-the-fly) over resolvents
				if (poss.size() <= dc_limits[2] && negs.size() <= dc_limits[2]) {
					forall_occurs(poss, i) {
						SCLAUSE& pos = (*cnf)[*i];
						if (pos.deleted()) continue;
						forall_occurs(negs, j) {
							SCLAUSE& neg = (*cnf)[*j];
							if (neg.deleted() || (pos.size() + neg.size() - 2) > SH_MAX_ERE_OUT) continue;
							uint32* m_c = smem + threadIdx.y * SH_MAX_ERE_OUT; // shared memory for resolvent
							int m_len = 0;
							if ((m_len = merge_ere(v, pos, neg, m_c)) > 1) {
								CL_ST type;
								if (pos.learnt() || neg.learnt()) type = LEARNT;
								else type = ORIGINAL;
								forward_equ(*cnf, *ot, m_c, m_len, type);
							}
						}
					}
				}
				gid += stride_y();
			}
		}
		//======================================================//
		//                GPU Wrappers Definitions              //
		//======================================================//
		void initConstants(cuLimit culimit)
		{
			CHECK(cudaMemcpyToSymbol(dc_limits, &culimit.limits, sizeof(uint32) * NLIMITS, 0, cudaMemcpyHostToDevice));
		}
		void copyIf(uint32* dest, CNF* src, GSTATS* gstats)
		{
			if (profile_gpu) cutimer->start();
			reset_stats << <1, 1 >> > (gstats);
			uint32 nBlocks = std::min((inf.nClauses + BLOCK1D - 1) / BLOCK1D, maxGPUTereads / BLOCK1D);
			copy_if_k << <nBlocks, BLOCK1D >> > (dest, src, gstats);
			if (profile_gpu) cutimer->stop(), cutimer->vo += cutimer->gpuTime();
			LOGERR("Copying literals failed");
			syncAll();
		}
		void calcScores(VARS* vars, uint32* hist)
		{
			if (profile_gpu) cutimer->start();
			uint32 nBlocks = std::min((inf.maxVar + BLOCK1D - 1) / BLOCK1D, maxGPUTereads / BLOCK1D);
			assign_scores << <nBlocks, BLOCK1D >> > (vars->eligible, vars->scores, hist, inf.maxVar);
			if (profile_gpu) cutimer->stop(), cutimer->vo += cutimer->gpuTime();
			LOGERR("Assigning scores failed");
			syncAll();
		}
		void calcScores(VARS* vars, uint32* hist, OT* ot)
		{
			if (profile_gpu) cutimer->start();
			uint32 nBlocks = std::min((inf.maxVar + BLOCK1D - 1) / BLOCK1D, maxGPUTereads / BLOCK1D);
			assign_scores << <nBlocks, BLOCK1D >> > (vars->eligible, vars->scores, hist, ot, inf.maxVar);
			if (profile_gpu) cutimer->stop(), cutimer->vo += cutimer->gpuTime();
			LOGERR("Assigning scores failed");
			syncAll();
		}
		void countMelted(VSTATE* vstate)
		{
			inf.n_del_vars_after = 0;
			forall_variables(v) {
				if (MELTED(vstate[v].state))
					inf.n_del_vars_after++;
			}
			assert(inf.n_del_vars_after >= inf.maxMelted);
			inf.n_del_vars_after -= inf.maxMelted;
			inf.maxMelted += inf.n_del_vars_after;
		}
		void countFinal(CNF* cnf, GSTATS* gstats, VSTATE* vstate)
		{
			reset_stats << <1, 1 >> > (gstats);
			const uint32 cnf_sz = inf.nClauses + (inf.nClauses >> 1);
			uint32 nBlocks = std::min((cnf_sz + (BLOCK1D << 1) - 1) / (BLOCK1D << 1), maxGPUTereads / (BLOCK1D << 1));
			uint32 smemSize = BLOCK1D * (sizeof(uint32) + sizeof(uint32));
			cnt_cls_lits << <nBlocks, BLOCK1D, smemSize >> > (cnf, gstats);
			countMelted(vstate);
			if (unified_access || sync_always) {
				LOGERR("Final CNF counting failed");
				syncAll();
			}
			GSTATS hstats;
			CHECK(cudaMemcpy(&hstats, gstats, sizeof(GSTATS), cudaMemcpyDeviceToHost));
			inf.n_cls_after = hstats.numClauses;
			inf.n_lits_after = hstats.numLits;
		}
		void countCls(CNF* cnf, GSTATS* gstats)
		{
			reset_stats << <1, 1 >> > (gstats);
			const uint32 cnf_sz = inf.nClauses;
			uint32 nBlocks = std::min((cnf_sz + (BLOCK1D << 1) - 1) / (BLOCK1D << 1), maxGPUTereads / (BLOCK1D << 1));
			uint32 smemSize = BLOCK1D * sizeof(uint32);
			cnt_cls << <nBlocks, BLOCK1D, smemSize >> > (cnf, gstats);
			GSTATS hstats;
			CHECK(cudaMemcpy(&hstats, gstats, sizeof(GSTATS), cudaMemcpyDeviceToHost));
			inf.n_cls_after = hstats.numClauses;
		}
		void countLits(CNF* cnf, GSTATS* gstats)
		{
			reset_stats << <1, 1 >> > (gstats);
			const uint32 cnf_sz = inf.nClauses;
			uint32 nBlocks = std::min((cnf_sz + (BLOCK1D << 1) - 1) / (BLOCK1D << 1), maxGPUTereads / (BLOCK1D << 1));
			uint32 smemSize = BLOCK1D * sizeof(uint32);
			cnt_lits << <nBlocks, BLOCK1D, smemSize >> > (cnf, gstats);
			GSTATS hstats;
			CHECK(cudaMemcpy(&hstats, gstats, sizeof(GSTATS), cudaMemcpyDeviceToHost));
			inf.n_lits_after = hstats.numLits;
		}
		void countAll(CNF* cnf, GSTATS* gstats)
		{
			reset_stats << <1, 1 >> > (gstats);
			const uint32 cnf_sz = inf.nClauses + (inf.nClauses >> 1);
			uint32 nBlocks = std::min((cnf_sz + (BLOCK1D << 1) - 1) / (BLOCK1D << 1), maxGPUTereads / (BLOCK1D << 1));
			uint32 smemSize = BLOCK1D * (sizeof(uint32) + sizeof(uint32));
			cnt_cls_lits << <nBlocks, BLOCK1D, smemSize >> > (cnf, gstats);
			GSTATS hstats;
			CHECK(cudaMemcpy(&hstats, gstats, sizeof(GSTATS), cudaMemcpyDeviceToHost));
			inf.n_cls_after = hstats.numClauses;
			inf.n_lits_after = hstats.numLits;
		}
		void evalReds(CNF* cnf, GSTATS* gstats, VSTATE* vstate)
		{
			reset_stats << <1, 1 >> > (gstats);
			const uint32 cnf_sz = inf.nClauses + (inf.nClauses >> 1);
			uint32 nBlocks1 = std::min((cnf_sz + (BLOCK1D << 1) - 1) / (BLOCK1D << 1), maxGPUTereads / (BLOCK1D << 1));
			uint32 smemSize1 = BLOCK1D * sizeof(uint32) * 2;
			cnt_reds << <nBlocks1, BLOCK1D, smemSize1 >> > (cnf, gstats);
			countMelted(vstate);
			GSTATS hstats;
			CHECK(cudaMemcpy(&hstats, gstats, sizeof(GSTATS), cudaMemcpyDeviceToHost)); // avoids unified memory migration on large scale
			inf.n_cls_after = hstats.numClauses;
			inf.n_lits_after = hstats.numLits;
		}
		void cuMemSetAsync(addr_t mem, const Byte& val, const size_t& size)
		{
			uint32 nBlocks = std::min(uint32((size + BLOCK1D - 1) / BLOCK1D), maxGPUTereads / BLOCK1D);
			memset_k<Byte> << <nBlocks, BLOCK1D >> > (mem, val, size);
			if (sync_always) {
				LOGERR("CUDA memory set failed");
				syncAll();
			}
		}
		void prepareCNFAsync(CNF* cnf, const cudaStream_t& _s)
		{
			assert(inf.nClauses);
			if (profile_gpu) cutimer->start(_s);
			uint32 nBlocks = std::min((inf.nClauses + BLOCK1D - 1) / BLOCK1D, maxGPUTereads / BLOCK1D);
			prep_cnf_k << <nBlocks, BLOCK1D, 0, _s >> > (cnf);
			if (profile_gpu) cutimer->stop(_s), cutimer->sig += cutimer->gpuTime();
			if (sync_always) {
				LOGERR("Signature calculation failed");
				syncAll();
			}
		}
		void reduceOTAsync(CNF* cnf, OT* ot, const bool& p)
		{
			assert(cnf != NULL);
			assert(ot != NULL);
			if (profile_gpu) cutimer->start();
			uint32 nBlocks = std::min(uint32((inf.nDualVars + BLOCK1D - 1) / BLOCK1D), maxGPUTereads / BLOCK1D);
			reduce_ot << <nBlocks, BLOCK1D >> > (cnf, ot);
			if (p || sync_always) {
				LOGERR("Occurrence table reduction failed");
				syncAll();
				if (p) {
					PFLRULER('=', 30);
					PFLOG0("\toccurrence table");
					ot->print();
					PFLRULER('=', 30);
				}
			}
			if (profile_gpu) cutimer->stop(), cutimer->rot += cutimer->gpuTime();
		}
		void createOTAsync(CNF* cnf, OT* ot, const bool& p)
		{
			assert(cnf != NULL);
			assert(ot != NULL);
			if (profile_gpu) cutimer->start();
			uint32 rstGridSize = std::min(uint32((inf.nDualVars + BLOCK1D - 1) / BLOCK1D), maxGPUTereads / BLOCK1D);
			reset_ot_k << <rstGridSize, BLOCK1D >> > (ot);
			uint32 otGridSize = std::min((inf.nClauses + BLOCK1D - 1) / BLOCK1D, maxGPUTereads / BLOCK1D);
			create_ot_k << <otGridSize, BLOCK1D >> > (cnf, ot);
			if (p || sync_always) {
				LOGERR("Occurrence table creation failed");
				syncAll();
				assert(ot->accViolation());
				if (p) {
					PFLRULER('=', 30);
					PFLOG0("\toccurrence table");
					ot->print();
					PFLRULER('=', 30);
				}
			}
			if (profile_gpu) cutimer->stop(), cutimer->cot += cutimer->gpuTime();
		}
		void sortOTAsync(CNF* cnf, OT* ot, VARS* vars, cudaStream_t* streams)
		{
			assert(cnf != NULL);
			assert(ot != NULL);
			assert(vars->numPVs);
			uint32 nBlocks = std::min((vars->numPVs + BLSORT - 1) / BLSORT, maxGPUTereads / BLSORT);
			cudaStream_t s1, s2;
			if (profile_gpu) s1 = s2 = 0, cutimer->start();
			else s1 = streams[0], s2 = streams[1];
			sort_ot_p << <nBlocks, BLSORT, 0, s1 >> > (cnf, ot, vars->pVars);
			sort_ot_n << <nBlocks, BLSORT, 0, s2 >> > (cnf, ot, vars->pVars);
			if (sync_always) {
				LOGERR("Sorting OT failed");
				syncAll();
			}
			if (profile_gpu) cutimer->stop(), cutimer->sot += cutimer->gpuTime();
		}
		void veRollBack(CNF* cnf, OT* ot, VARS* vars, uint32* type)
		{
			// undo marked gate clauses if phase-2 in 've' ran out of memory which is an extreme case
			uint32 nBlocks = std::min((vars->numPVs + BLVE2 - 1) / BLVE2, maxGPUTereads / BLVE2);
			roll_back_gates << <nBlocks, BLVE2 >> > (cnf, ot, vars->pVars, type);
			LOGERR("BVE Roll-back failed");
			sync();
		}
		void veAsync(CNF* cnf, OT* ot, VARS* vars, cudaStream_t* streams, cuMM& cumm, const cuHist& cuhist, const bool& in)
		{
			assert(vars->numPVs);
			if (profile_gpu) cutimer->start();
			if (!atomic_ve) {
				const uint32 cs_offset = cumm.pinnedCNF()->size();
				const S_REF data_offset = cumm.pinnedCNF()->data().size;
				S_REF* rref = cuhist.d_segs;
				uint32* type = cuhist.d_hist, *rpos = type + inf.maxVar;
				uint32* vorg = cuhist.d_vorg;
				// Phase-1
				uint32 smSize1 = (BLVE1 * SH_MAX_BVE_OUT1) * sizeof(uint32);
				#if VE_DBG
				if (in) in_ve_k_1 << <1, 1, smSize1 >> > (cnf, ot, vars->pVars, vars->units, vars->resolved, vorg, type, rpos, rref);
				else	   ve_k_1 << <1, 1, smSize1 >> > (cnf, ot, vars->pVars, vars->units, vars->resolved, vorg, type, rpos, rref);
				#else
				uint32 nBlocks1 = std::min((vars->numPVs + BLVE1 - 1) / BLVE1, maxGPUTereads / BLVE1);
				if (in) in_ve_k_1 << <nBlocks1, BLVE1, smSize1 >> > (cnf, ot, vars->pVars, vars->units, vars->resolved, vorg, type, rpos, rref);
				else	   ve_k_1 << <nBlocks1, BLVE1, smSize1 >> > (cnf, ot, vars->pVars, vars->units, vars->resolved, vorg, type, rpos, rref);
				#endif
				// Phase-2
				size_t tb1 = 0, tb2 = 0;
				void* ts1 = NULL, *ts2 = NULL;
				DeviceScan::ExclusiveScan(NULL, tb1, rpos, rpos, Sum(), cs_offset, vars->numPVs);
				DeviceScan::ExclusiveScan(NULL, tb2, rref, rref, Sum(), data_offset, vars->numPVs);
				LOGERR("BVE Phase-1 failed"); 
				sync(); //sync phase-1
				if (tb1 > cumm.literalsCap() || tb2 > cumm.scatterCap()) {
					addr_t tmpmem = cumm.allocTemp(tb1 + tb2);
					if (tmpmem == NULL) {
						veRollBack(cnf, ot, vars, type);
						return; 
					}
					ts1 = tmpmem, ts2 = tmpmem + tb1;
				}
				else ts1 = cuhist.d_lits, ts2 = cumm.scatter();
				DeviceScan::ExclusiveScan(ts1, tb1, rpos, rpos, Sum(), cs_offset, vars->numPVs, streams[0]);
				DeviceScan::ExclusiveScan(ts2, tb2, rref, rref, Sum(), data_offset, vars->numPVs, streams[1]);
				uint32 nBlocks2 = std::min((vars->numPVs + BLVE2 - 1) / BLVE2, maxGPUTereads / BLVE2);
				uint32 smSize2 = (BLVE2 * SH_MAX_BVE_OUT2) * sizeof(uint32);
				LOGERR("BVE Phase-2 failed"); 
				sync(streams[0]), sync(streams[1]); // sync phase-2
				#if VE_DBG
				ve_k_2 << <1, 1, smSize2 >> > (cnf, ot, vars->pVars, vars->units, vars->resolved, vorg, type, rpos, rref);
				#else
				// Phase-3
				ve_k_2 << <nBlocks2, BLVE2, smSize2 >> > (cnf, ot, vars->pVars, vars->units, vars->resolved, vorg, type, rpos, rref);
				#endif
			}
			else {
				#if VE_DBG
				ve_k << <1, 1 >> > (cnf, ot, vars->pVars, vars->units, vars->resolved, cuhist.d_vorg);
				#else
				uint32 nBlocks1 = std::min((vars->numPVs + BLVE - 1) / BLVE, maxGPUTereads / BLVE);
				ve_k << <nBlocks1, BLVE >> > (cnf, ot, vars->pVars, vars->units, vars->resolved, cuhist.d_vorg);
				#endif
			}
			if (profile_gpu) cutimer->stop(), cutimer->ve += cutimer->gpuTime();
			LOGERR("BVE Elimination failed"); 
			syncAll();
		}
		void hseAsync(CNF* cnf, OT* ot, VARS* vars)
		{
			assert(vars->numPVs);
			if (profile_gpu) cutimer->start();
#if SS_DBG
			putchar('\n');
			hse_k << <1, 1 >> > (cnf, ot, vars->pVars, vars->units);
#else
			uint32 nBlocks = std::min((vars->numPVs + BLHSE - 1) / BLHSE, maxGPUTereads / BLHSE);
			hse_k << <nBlocks, BLHSE >> > (cnf, ot, vars->pVars, vars->units);
#endif
			if (profile_gpu) cutimer->stop(), cutimer->hse += cutimer->gpuTime();
			if (sync_always) { 
				LOGERR("HSE Elimination failed");
				syncAll();
			}
		}
		void bceAsync(CNF* cnf, OT* ot, VARS* vars, const uint32* vorg)
		{
			assert(vars->numPVs);
			if (profile_gpu) cutimer->start();
			uint32 nBlocks = std::min((vars->numPVs + BLBCE - 1) / BLBCE, maxGPUTereads / BLBCE);
			bce_k << <nBlocks, BLBCE >> > (cnf, ot, vars->pVars, vars->resolved, vorg);
			if (profile_gpu) cutimer->stop(), cutimer->bce += cutimer->gpuTime();
			if (sync_always) {
				LOGERR("BCE Elimination failed");
				syncAll();
			}
		}
		void ereAsync(CNF* cnf, OT* ot, VARS* vars)
		{
			assert(vars->numPVs);
			if (profile_gpu) cutimer->start();
			dim3 block2D(devProp.warpSize, devProp.warpSize), grid2D(1, 1, 1);
			grid2D.y = std::min((vars->numPVs + block2D.y - 1) / block2D.y, maxGPUTereads / block2D.y);
			uint32 smemSize = devProp.warpSize * SH_MAX_ERE_OUT * sizeof(uint32);
			ere_k << <grid2D, block2D, smemSize >> > (cnf, ot, vars->pVars);
			if (profile_gpu) cutimer->stop(), cutimer->ere += cutimer->gpuTime();
			if (sync_always) {
				LOGERR("ERE Elimination failed");
				syncAll();
			}
		}

	}
}