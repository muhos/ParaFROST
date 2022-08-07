/***********************************************************************[kernels.cu]
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

#include "bounded.cuh"
#include "subsume.cuh"
#include "blocked.cuh"
#include "redundancy.cuh"
#include "simplify.cuh"
#include <cub/device/device_scan.cuh>

using namespace cub;

namespace pFROST {

	// ===== for CNF counting =====
	#define MAXREDUCEBLOCKS 256
	uint32 hostCBlocks[MAXREDUCEBLOCKS], hostLBlocks[MAXREDUCEBLOCKS];
	__device__ uint32 devCBlocks[MAXREDUCEBLOCKS], devLBlocks[MAXREDUCEBLOCKS];
	__device__ uint32 gcounter;
	__device__ int    lastEliminatedID;
	__global__ void reset_id() { lastEliminatedID = -1; }
	__global__ void print_id() { printf("c lastEliminatedID = %d\n", lastEliminatedID); }
	__global__ void reset_counter() { gcounter = 0; }
	__global__ void print_counter() { printf("c gcounter = %d\n", gcounter); }
	__global__ void check_counter(const uint32 checksum) { assert(checksum == gcounter); }
	//=============================

	__global__ void printCMem()
	{
		printf("c   Number of buckets per SCLAUSE = %d\n", dc_nbuckets);
		printf("c   Variable mapping array %p\n", dc_ptrs->d_vorg);
		printf("c   Literal  bytes array %p\n", dc_ptrs->d_lbyte);
		printf("c   Limits[%d] = {", NLIMITS);
		for (int i = 0; i < NLIMITS; i++)
			printf("%d, ", dc_options[i]);
		printf("}\n");
	}

	template<class T>
	__global__ void memset_k(T* mem, T val, size_t size)
	{
		size_t tid = global_tx;
		while (tid < size) { mem[tid] = val; tid += stride_x; }
	}

	__global__ void prep_cnf_k(CNF* cnf)
	{
		gridtype tid = global_tx;
		while (tid < cnf->size()) {
			SCLAUSE& c = cnf->clause(tid);
			devSort(c.data(), c.size());
			calcSig(c);
			tid += stride_x;
		}
	}

	__global__ void reset_ot_k(OT* ot)
	{
		gridtype tid = global_tx;
		while (tid < ot->size()) {
			(*ot)[tid].clear();
			tid += stride_x;
		}
	}

	__global__ void create_ot_k(CNF* __restrict__ cnf, OT* __restrict__ ot_ptr)
	{
		gridtype tid = global_tx;
		while (tid < cnf->size()) {
			const S_REF r = cnf->ref(tid);
			SCLAUSE& c = (*cnf)[r];
			if (c.original() || c.learnt()) {
				OT& ot = *ot_ptr;
#pragma unroll
				forall_clause(c, lit) {
					ot[*lit].insert(r);
				}
			}
			tid += stride_x;
		}
	}

	__global__ void reduce_ot(const CNF* __restrict__ cnf, OT* __restrict__ ot)
	{
		gridtype tid = global_tx;
		while (tid < ot->size()) {
			reduceOL(*cnf, (*ot)[tid]);
			tid += stride_x;
		}
	}

	__global__ void sort_ot_p(const CNF* __restrict__ cnf, OT* __restrict__ ot, const cuVecU* __restrict__ pVars)
	{
		gridtype tid = global_tx;
		while (tid < pVars->size()) {
			const uint32 x = pVars->at(tid), p = V2L(x);
			assert(x);
			OL& ol = (*ot)[p];
			devSort(ol.data(), ol.size(), OLIST_CMP(cnf));
			tid += stride_x;
		}
	}

	__global__ void sort_ot_n(const CNF* __restrict__ cnf, OT* __restrict__ ot, const cuVecU* __restrict__ pVars)
	{
		gridtype tid = global_tx;
		while (tid < pVars->size()) {
			const uint32 x = pVars->at(tid), n = NEG(V2L(x));
			assert(x);
			OL& ol = (*ot)[n];
			devSort(ol.data(), ol.size(), OLIST_CMP(cnf));
			tid += stride_x;
		}
	}

	__global__ void copy_if_k(uint32* __restrict__ dest, CNF* __restrict__ src)
	{
		gridtype tid = global_tx;
		while (tid < src->size()) {
			SCLAUSE& c = src->clause(tid);
			if (c.original() || c.learnt()) {
				uint32* d = dest + atomicAdd(&gcounter, c.size());
#pragma unroll
				forall_clause(c, s) { *d++ = *s; }
			}
			tid += stride_x;
		}
	}

	__global__ void cnt_proof_verify(const uint32* __restrict__ literals, const uint32 numLits)
	{
		gridtype tid = 0;
		while (tid < numLits) {
			addr_t lbyte = dc_ptrs->d_lbyte;
			const uint32 lit = literals[tid];
			if (lit & 0xF0000000) lbyte[lit] = 5;
			else if (lit & 0x0FE00000) lbyte[lit] = 4;
			else if (lit & 0x001FC000) lbyte[lit] = 3;
			else if (lit & 0x00003F80) lbyte[lit] = 2;
			else lbyte[lit] = 1;
			printf(" literal(%d) has %d bytes of its original\n", SIGN(lit) ? -int(ABS(lit)) : ABS(lit), lbyte[lit]);
			gcounter += lbyte[lit];
			tid++;
		}
		printf(" total = %d\n", gcounter);
	}

	__global__ void cnt_proof(const uint32* __restrict__ literals, const uint32 numLits)
	{
		uint32* sh_bytes = SharedMemory<uint32>();
		gridtype tid = global_tx_off;
		uint32 nbytes = 0;
		while (tid < numLits) {
			addr_t lbyte = dc_ptrs->d_lbyte;
			uint32 lit = literals[tid];
			countBytes(lit, lbyte[lit], nbytes);
			gridtype off = tid + blockDim.x;
			if (off < numLits) {
				lit = literals[off];
				countBytes(lit, lbyte[lit], nbytes);
			}
			tid += stride_x_off;
		}
		loadShared(sh_bytes, nbytes, numLits);
		sharedReduce(sh_bytes, nbytes);
		warpReduce(sh_bytes, nbytes);
		if (!threadIdx.x) devLBlocks[blockIdx.x] = nbytes;
	}

	__global__ void cnt_cls(const CNF* __restrict__ cnf)
	{
		uint32* sh_rCls = SharedMemory<uint32>();
		gridtype tid = global_tx_off;
		uint32 nCls = 0;
		while (tid < cnf->size()) {
			const SCLAUSE& c1 = cnf->clause(tid);
			if (c1.original() || c1.learnt()) nCls++;
			uint32 off = tid + blockDim.x;
			if (off < cnf->size()) {
				const SCLAUSE& c2 = cnf->clause(off);
				if (c2.original() || c2.learnt()) nCls++;
			}
			tid += stride_x_off;
		}
		loadShared(sh_rCls, nCls, cnf->size());
		sharedReduce(sh_rCls, nCls);
		warpReduce(sh_rCls, nCls);
		if (!threadIdx.x) devCBlocks[blockIdx.x] = nCls;
	}

	__global__ void cnt_lits(const CNF* __restrict__ cnf)
	{
		uint32* sh_rLits = SharedMemory<uint32>();
		gridtype tid = global_tx_off;
		uint32 nLits = 0;
		while (tid < cnf->size()) {
			const SCLAUSE& c1 = cnf->clause(tid);
			if (c1.original() || c1.learnt()) nLits += c1.size();
			uint32 off = tid + blockDim.x;
			if (off < cnf->size()) {
				const SCLAUSE& c2 = cnf->clause(off);
				if (c2.original() || c2.learnt()) nLits += c2.size();
			}
			tid += stride_x_off;
		}
		loadShared(sh_rLits, nLits, cnf->size());
		sharedReduce(sh_rLits, nLits);
		warpReduce(sh_rLits, nLits);
		if (!threadIdx.x) devLBlocks[blockIdx.x] = nLits;
	}

	__global__ void cnt_cls_lits(const CNF* __restrict__ cnf)
	{
		uint32* sh_rCls = SharedMemory<uint32>();
		uint32* sh_rLits = sh_rCls + blockDim.x;
		gridtype tid = global_tx_off;
		uint32 nCls = 0;
		uint32 nLits = 0;
		while (tid < cnf->size()) {
			const SCLAUSE& c1 = cnf->clause(tid);
			if (c1.original() || c1.learnt()) nCls++, nLits += c1.size();
			gridtype off = tid + blockDim.x;
			if (off < cnf->size()) {
				const SCLAUSE& c2 = cnf->clause(off);
				if (c2.original() || c2.learnt()) nCls++, nLits += c2.size();
			}
			tid += stride_x_off;
		}
		loadShared(sh_rCls, nCls, sh_rLits, nLits, cnf->size());
		sharedReduce(sh_rCls, nCls, sh_rLits, nLits);
		warpReduce(sh_rCls, nCls, sh_rLits, nLits);
		if (!threadIdx.x) {
			gridtype bx = blockIdx.x;
			devCBlocks[bx] = nCls;
			devLBlocks[bx] = nLits;
		}
	}

	__global__ void mapfrozen_k(const uint32* __restrict__ frozen, uint32* __restrict__ varcore, const uint32 size)
	{
		gridtype tid = global_tx;
		while (tid < size) {
			assert(frozen[tid] && frozen[tid] < NOVAR);
			varcore[frozen[tid]] = tid;
			tid += stride_x;
		}
	}

	__global__ void assign_scores(
		uint32* __restrict__ eligible,
		uint32* __restrict__ scores,
		const uint32* __restrict__ hist,
		uint32 size)
	{
		gridtype tid = global_tx;
		while (tid < size) {
			const uint32 v = tid + 1;
			const uint32 p = V2L(v), ps = hist[p], ns = hist[NEG(p)];
			eligible[tid] = v;
			scores[v] = ps * ns;
			tid += stride_x;
		}
	}

	__global__ void assign_scores(
		uint32* __restrict__ eligible,
		uint32* __restrict__ scores,
		uint32* __restrict__ hist,
		const OT* __restrict__ ot,
		uint32 size)
	{
		gridtype tid = global_tx;
		while (tid < size) {
			const uint32 v = tid + 1;
			const uint32 p = V2L(v), n = NEG(p), ps = (*ot)[p].size(), ns = (*ot)[n].size();
			hist[p] = ps, hist[n] = ns;
			eligible[tid] = v;
			scores[v] = ps * ns;
			tid += stride_x;
		}
	}

	__global__ void ve_k(
		CNF* __restrict__ cnfptr,
		OT* __restrict__ otptr,
		const cuVecU* __restrict__ pVars,
		Byte* __restrict__ eliminated,
		cuVecU* __restrict__ units,
		cuVecU* __restrict__ resolved,
		cuVecB* __restrict__ proof)
	{
		gridtype tid = global_tx;
		__shared__ uint32 outs[BLVE * SH_MAX_BVE_OUT];
		while (tid < pVars->size()) {
			const uint32 x = pVars->at(tid);
			assert(x);
			assert(!ELIMINATED(eliminated[x]));
			assert(!IS_ADDING(eliminated[x]));
			const uint32 p = V2L(x), n = NEG(p);
			CNF& cnf = *cnfptr;
			OT& ot = *otptr;
			OL& poss = ot[p], & negs = ot[n];
			uint32 pOrgs = 0, nOrgs = 0;
			countOrgs(cnf, poss, pOrgs), countOrgs(cnf, negs, nOrgs);
			bool elim = false;
			// pure-literal elimination
			if (!pOrgs || !nOrgs) {
				toblivion(p, pOrgs, nOrgs, cnf, poss, negs, resolved), elim = true;
			}
			// Equiv/NOT-gate Reasoning
			else if (uint32 def = find_equ_gate(p, cnf, poss, negs)) {
				saveResolved(p, pOrgs, nOrgs, cnf, poss, negs, resolved); // must be called before substitution
				substitute_single(p, def, cnf, poss, negs, units, proof), elim = true;
			}
			else {
				assert(pOrgs && nOrgs);
				const uint32 nClsBefore = pOrgs + nOrgs;
				uint32* shared_outs = outs + threadIdx.x * SH_MAX_BVE_OUT;
				uint32 nElements = 0, nAddedCls = 0, nAddedLits = 0;
				// AND-gate Reasoning
				if (nOrgs < SH_MAX_BVE_OUT && find_ao_gate(n, nClsBefore, cnf, ot, shared_outs, nElements, nAddedCls, nAddedLits)) {
					if (nAddedCls) substitute_x(x, nElements, nAddedCls, nAddedLits, cnf, poss, negs, units, proof, shared_outs);
					toblivion(p, pOrgs, nOrgs, cnf, poss, negs, resolved), elim = true;
				}
				// OR-gate Reasoning
				else if (!nAddedCls && pOrgs < SH_MAX_BVE_OUT && find_ao_gate(p, nClsBefore, cnf, ot, shared_outs, nElements, nAddedCls, nAddedLits)) {
					if (nAddedCls) substitute_x(x, nElements, nAddedCls, nAddedLits, cnf, poss, negs, units, proof, shared_outs);
					toblivion(p, pOrgs, nOrgs, cnf, poss, negs, resolved), elim = true;
				}
				// ITE-gate Reasoning
				else if (find_ite_gate(p, nClsBefore, cnf, ot, nElements, nAddedCls, nAddedLits)) {
					if (nAddedCls) substitute_x(x, nElements, nAddedCls, nAddedLits, cnf, poss, negs, units, proof, shared_outs);
					toblivion(p, pOrgs, nOrgs, cnf, poss, negs, resolved), elim = true;
				}
				else if (!nAddedCls && find_ite_gate(n, nClsBefore, cnf, ot, nElements, nAddedCls, nAddedLits)) {
					if (nAddedCls) substitute_x(x, nElements, nAddedCls, nAddedLits, cnf, poss, negs, units, proof, shared_outs);
					toblivion(p, pOrgs, nOrgs, cnf, poss, negs, resolved), elim = true;
				}
				// XOR-gate Reasoning
				else if (find_xor_gate(p, nClsBefore, cnf, ot, shared_outs, nElements, nAddedCls, nAddedLits)) {
					if (nAddedCls) substitute_x(x, nElements, nAddedCls, nAddedLits, cnf, poss, negs, units, proof, shared_outs);
					toblivion(p, pOrgs, nOrgs, cnf, poss, negs, resolved), elim = true;
				}
				else if (!nAddedCls && find_xor_gate(n, nClsBefore, cnf, ot, shared_outs, nElements, nAddedCls, nAddedLits)) {
					if (nAddedCls) substitute_x(x, nElements, nAddedCls, nAddedLits, cnf, poss, negs, units, proof, shared_outs);
					toblivion(p, pOrgs, nOrgs, cnf, poss, negs, resolved), elim = true;
				}
				// n-by-m resolution
				else if (!nAddedCls && resolve(x, nClsBefore, cnf, poss, negs, nElements, nAddedCls, nAddedLits)) {
					if (nAddedCls) resolve_x(x, nElements, nAddedCls, nAddedLits, cnf, poss, negs, units, proof, shared_outs);
					toblivion(p, pOrgs, nOrgs, cnf, poss, negs, resolved);
					elim = true;
				}
			}
			if (elim) ELIMINATE(eliminated[x]);
			tid += stride_x;
		}
	}

	__global__ void in_ve_k_1(
		CNF* __restrict__ cnfptr,
		OT* __restrict__ otptr,
		const cuVecU* __restrict__ pVars,
		Byte* __restrict__ eliminated,
		cuVecU* __restrict__ units,
		cuVecU* __restrict__ resolved,
		cuVecB* __restrict__ proof,
		const uint32* __restrict__ varcore,
		uint32* __restrict__ ucnt,
		uint32* __restrict__ type,
		uint32* __restrict__ rpos,
		S_REF* __restrict__ rref)
	{
		gridtype tid = global_tx;
		uint32* outs = SharedMemory<uint32>();
		while (tid < pVars->size()) {
			const uint32 x = pVars->at(tid);
			assert(x);
			assert(!ELIMINATED(eliminated[x]));
			assert(!IS_ADDING(eliminated[x]));
			const uint32 p = V2L(x), n = NEG(p);
			CNF& cnf = *cnfptr;
			OT& ot = *otptr;
			OL& poss = ot[p], & negs = ot[n];
			uint32 pOrgs = 0, nOrgs = 0;
			countOrgs(cnf, poss, pOrgs), countOrgs(cnf, negs, nOrgs);
			// pure-literal elimination
			if (!pOrgs || !nOrgs) {
				toblivion(p, pOrgs, nOrgs, cnf, poss, negs, resolved);
				ucnt[tid] = 0, type[tid] = 0, rref[tid] = 0, rpos[tid] = 0, ELIMINATE(eliminated[x]);
			}
			// Equiv/NOT-gate Reasoning
			else if (uint32 def = find_equ_gate(p, cnf, poss, negs)) {
				saveResolved(p, pOrgs, nOrgs, cnf, poss, negs, resolved); // must be called before substitution
				substitute_single(p, def, cnf, poss, negs, units, proof);
				ucnt[tid] = 0, type[tid] = 0, rref[tid] = 0, rpos[tid] = 0, ELIMINATE(eliminated[x]);
			}
			else {
				assert(pOrgs && nOrgs);
				const uint32 nClsBefore = pOrgs + nOrgs;
				uint32 elimType = 0, nElements = 0, nAddedCls = 0, nAddedLits = 0; // nElements: compresses '#units' and '#proof bytes'  
				//=====================
				// check resolvability 
				//=====================
				if (nClsBefore > 2) {
					// AND/OR-gate Reasoning
					if (nOrgs < SH_MAX_BVE_OUT1 &&
						find_ao_gate(n, nClsBefore, cnf, ot, &outs[threadIdx.x * SH_MAX_BVE_OUT1], nElements, nAddedCls, nAddedLits))
						elimType = AOIX_MASK;
					else if (!nAddedCls && pOrgs < SH_MAX_BVE_OUT1 &&
						find_ao_gate(p, nClsBefore, cnf, ot, &outs[threadIdx.x * SH_MAX_BVE_OUT1], nElements, nAddedCls, nAddedLits))
						elimType = AOIX_MASK;
				}
				if (!elimType && nClsBefore > 3) {
					// ITE-gate Reasoning
					if (find_ite_gate(p, nClsBefore, cnf, ot, nElements, nAddedCls, nAddedLits))
						elimType = AOIX_MASK;
					else if (!nAddedCls &&
						find_ite_gate(n, nClsBefore, cnf, ot, nElements, nAddedCls, nAddedLits))
						elimType = AOIX_MASK;
					// XOR-gate Reasoning
					else if (find_xor_gate(p, nClsBefore, cnf, ot, &outs[threadIdx.x * SH_MAX_BVE_OUT1], nElements, nAddedCls, nAddedLits))
						elimType = AOIX_MASK;
					else if (!nAddedCls &&
						find_xor_gate(n, nClsBefore, cnf, ot, &outs[threadIdx.x * SH_MAX_BVE_OUT1], nElements, nAddedCls, nAddedLits))
						elimType = AOIX_MASK;
				}
				// fun-tab reasoning
				if (varcore && !elimType && nClsBefore > 2 &&
					find_fun_gate(p, nClsBefore, varcore, cnf, ot, &outs[threadIdx.x * SH_MAX_BVE_OUT1], nElements, nAddedCls, nAddedLits))
					elimType = CORE_MASK;
				// n-by-m resolution
				else if (!elimType && !nAddedCls &&
					resolve(x, nClsBefore, cnf, poss, negs, nElements, nAddedCls, nAddedLits))
					elimType = RES_MASK;
				//=====================
				// check addibility 
				//=====================
				if (!nAddedCls) { // eliminated without resolvents
					toblivion(p, pOrgs, nOrgs, cnf, poss, negs, resolved);
					ucnt[tid] = 0, type[tid] = 0, rref[tid] = 0, rpos[tid] = 0, ELIMINATE(eliminated[x]);
				}
				else if (elimType) { // can be eliminated with resolvents in next phase
					assert(nAddedLits >= nAddedCls);
					assert(elimType <= TYPE_MASK);
					assert(nAddedCls <= ADDEDCLS_MAX);
					assert(nAddedLits <= ADDEDLITS_MAX);
					// save elimination info.
					type[tid] = ENCODEVARINFO(elimType, nAddedCls, nAddedLits);
					ucnt[tid] = nElements, rpos[tid] = nAddedCls, rref[tid] = nAddedLits + dc_nbuckets * nAddedCls;
				}
				else  // cannot be eliminated
					ucnt[tid] = 0, type[tid] = 0, rref[tid] = 0, rpos[tid] = 0;
			}
			tid += stride_x;
		}
	}

	__global__ void ve_k_1(
		CNF* __restrict__ cnfptr,
		OT* __restrict__ otptr,
		const cuVecU* __restrict__ pVars,
		Byte* __restrict__ eliminated,
		cuVecU* __restrict__ units,
		cuVecU* __restrict__ resolved,
		cuVecB* __restrict__ proof,
		const uint32* __restrict__ varcore,
		uint32* __restrict__ ucnt,
		uint32* __restrict__ type,
		uint32* __restrict__ rpos,
		S_REF* __restrict__ rref)
	{
		gridtype tid = global_tx;
		uint32* outs = SharedMemory<uint32>();
		while (tid < pVars->size()) {
			const uint32 x = pVars->at(tid);
			assert(x);
			assert(!ELIMINATED(eliminated[x]));
			assert(!IS_ADDING(eliminated[x]));
			const uint32 p = V2L(x), n = NEG(p);
			CNF& cnf = *cnfptr;
			OT& ot = *otptr;
			OL& poss = ot[p], & negs = ot[n];
			const uint32 pOrgs = poss.size();
			const uint32 nOrgs = negs.size();
			// pure-literal elimination
			if (!pOrgs || !nOrgs) {
				toblivion(p, pOrgs, nOrgs, cnf, poss, negs, resolved);
				ucnt[tid] = 0, type[tid] = 0, rref[tid] = 0, rpos[tid] = 0, ELIMINATE(eliminated[x]);
			}
			// Equiv/NOT-gate Reasoning
			else if (uint32 def = find_equ_gate(p, cnf, poss, negs)) {
				saveResolved(p, pOrgs, nOrgs, cnf, poss, negs, resolved); // must be called before substitution
				substitute_single(p, def, cnf, poss, negs, units, proof);
				ucnt[tid] = 0, type[tid] = 0, rref[tid] = 0, rpos[tid] = 0, ELIMINATE(eliminated[x]);
			}
			else {
				assert(pOrgs && nOrgs);
				const uint32 nClsBefore = pOrgs + nOrgs;
				uint32 elimType = 0, nElements = 0, nAddedCls = 0, nAddedLits = 0;
				//=====================
				// check resolvability 
				//=====================
				if (nClsBefore > 2) {
					// AND/OR-gate Reasoning
					if (nOrgs < SH_MAX_BVE_OUT1 &&
						find_ao_gate(n, nClsBefore, cnf, ot, &outs[threadIdx.x * SH_MAX_BVE_OUT1], nElements, nAddedCls, nAddedLits))
						elimType = AOIX_MASK;
					else if (!nAddedCls && pOrgs < SH_MAX_BVE_OUT1 &&
						find_ao_gate(p, nClsBefore, cnf, ot, &outs[threadIdx.x * SH_MAX_BVE_OUT1], nElements, nAddedCls, nAddedLits))
						elimType = AOIX_MASK;
				}
				if (!elimType && nClsBefore > 3) {
					// ITE-gate Reasoning
					if (find_ite_gate(p, nClsBefore, cnf, ot, nElements, nAddedCls, nAddedLits))
						elimType = AOIX_MASK;
					else if (!nAddedCls &&
						find_ite_gate(n, nClsBefore, cnf, ot, nElements, nAddedCls, nAddedLits))
						elimType = AOIX_MASK;
					// XOR-gate Reasoning
					else if (find_xor_gate(p, nClsBefore, cnf, ot, &outs[threadIdx.x * SH_MAX_BVE_OUT1], nElements, nAddedCls, nAddedLits))
						elimType = AOIX_MASK;
					else if (!nAddedCls &&
						find_xor_gate(n, nClsBefore, cnf, ot, &outs[threadIdx.x * SH_MAX_BVE_OUT1], nElements, nAddedCls, nAddedLits))
						elimType = AOIX_MASK;
				}
				// fun-tab reasoning
				if (varcore && !elimType && nClsBefore > 2 &&
					find_fun_gate(p, nClsBefore, varcore, cnf, ot, &outs[threadIdx.x * SH_MAX_BVE_OUT1], nElements, nAddedCls, nAddedLits))
					elimType = CORE_MASK;
				// n-by-m resolution
				else if (!elimType && !nAddedCls &&
					resolve(x, nClsBefore, cnf, poss, negs, nElements, nAddedCls, nAddedLits))
					elimType = RES_MASK;
				//=====================
				// check addibility 
				//=====================
				if (!nAddedCls) { // eliminated without resolvents
					toblivion(p, pOrgs, nOrgs, cnf, poss, negs, resolved);
					ucnt[tid] = 0, type[tid] = 0, rref[tid] = 0, rpos[tid] = 0, ELIMINATE(eliminated[x]);
				}
				else if (elimType) { // can be eliminated with resolvents in next phase
					assert(nAddedLits >= nAddedCls);
					assert(elimType <= TYPE_MASK);
					assert(nAddedCls <= ADDEDCLS_MAX);
					assert(nAddedLits <= ADDEDLITS_MAX);
					// save elimination info.
					type[tid] = ENCODEVARINFO(elimType, nAddedCls, nAddedLits);
					ucnt[tid] = nElements, rpos[tid] = nAddedCls, rref[tid] = nAddedLits + dc_nbuckets * nAddedCls;
				}
				else  // cannot be eliminated
					ucnt[tid] = 0, type[tid] = 0, rref[tid] = 0, rpos[tid] = 0;
			}
			tid += stride_x;
		}
	}

	__global__ void ve_k_2(
		CNF* __restrict__ cnf,
		OT* __restrict__ ot,
		const cuVecU* __restrict__ pVars,
		Byte* __restrict__ eliminated,
		uint32* __restrict__ eligible,
		cuVecU* __restrict__ units,
		cuVecU* __restrict__ resolved,
		cuVecB* __restrict__ proof,
		const uint32* __restrict__ ucnt,
		const uint32* __restrict__ type,
		const uint32* __restrict__ rpos,
		const S_REF* __restrict__ rref)
	{
		gridtype tid = global_tx;
		uint32* outs = SharedMemory<uint32>();
		while (tid < pVars->size()) {
			const uint32 x = pVars->at(tid);
			assert(x);
			const uint32 xinfo = type[tid];
			const uint32 elimType = RECOVERTYPE(xinfo);
			assert(elimType <= TYPE_MASK);
			if (elimType) {
				assert(!ELIMINATED(eliminated[x]));
				assert(!IS_ADDING(eliminated[x]));
				const uint32 p = V2L(x);
				const uint32 nAddedCls = RECOVERADDEDCLS(xinfo);
				const uint32 nAddedLits = RECOVERADDEDLITS(xinfo);
				assert(nAddedCls && nAddedCls <= ADDEDCLS_MAX);
				assert(nAddedLits && nAddedLits <= ADDEDLITS_MAX);
				const uint32 added_pos = rpos[tid];
				const S_REF added_ref = rref[tid];
				OL& poss = (*ot)[p], & negs = (*ot)[NEG(p)];
				if (memorySafe(tid, x, nAddedCls, nAddedLits, added_pos, added_ref, cnf)) {
					saveResolved(p, *cnf, poss, negs, resolved);
					if (IS_RES(elimType))
						resolve_x(x, ucnt[tid], nAddedCls, nAddedLits, added_pos, added_ref, *cnf, poss, negs, units, proof, outs + threadIdx.x * SH_MAX_BVE_OUT2);
					else if (IS_CORE(elimType))
						coresubstitute_x(x, ucnt[tid], nAddedCls, nAddedLits, added_pos, added_ref, *cnf, poss, negs, units, proof, outs + threadIdx.x * SH_MAX_BVE_OUT2);
					else {
						assert(IS_AOIX(elimType));
						substitute_x(x, ucnt[tid], nAddedCls, nAddedLits, added_pos, added_ref, *cnf, poss, negs, units, proof, outs + threadIdx.x * SH_MAX_BVE_OUT2);
					}
					ELIMINATE(eliminated[x]);
					MARKADDING(eliminated[x]);
					atomicAggMax(&lastEliminatedID, tid);
				}
				else if (!IS_RES(elimType)) freezeClauses(*cnf, poss, negs);
			}
			eligible[tid] = eliminated[x] ? 0 : x;
			tid += stride_x;
		}
	}

	__global__ void resizeCNF_k(CNF* cnf,
		const uint32* __restrict__ type,
		const uint32* __restrict__ rpos,
		const S_REF* __restrict__ rref,
		const int verbose)
	{
		if (lastEliminatedID >= 0) {
			const uint32 lastAdded = type[lastEliminatedID];
			const uint32 lastAddedPos = rpos[lastEliminatedID];
			const S_REF  lastAddedRef = rref[lastEliminatedID];
			assert(lastAdded < NOVAR);
			assert(lastAddedPos < NOVAR);
			assert(lastAddedRef < GNOREF);
			assert(RECOVERTYPE(lastAdded) < TYPE_MASK);
			const uint32 lastAddedCls = RECOVERADDEDCLS(lastAdded);
			const uint32 lastAddedLits = RECOVERADDEDLITS(lastAdded);
			assert(lastAddedCls && lastAddedCls <= ADDEDCLS_MAX);
			assert(lastAddedLits && lastAddedLits <= ADDEDLITS_MAX);
			const S_REF lastAddedBuckets = lastAddedLits + dc_nbuckets * lastAddedCls;
			const S_REF data_size = lastAddedBuckets + lastAddedRef;
			const uint32 cs_size = lastAddedCls + lastAddedPos;
			cnf->resize(data_size, cs_size);
			if (verbose > 1) printf("c   resized CNF to %d clauses and %lld data for a last ID %d\n", cs_size, data_size, lastEliminatedID);
		}
	}

	__global__ void sub_k(
		CNF* __restrict__ cnf,
		OT* __restrict__ otptr,
		cuVecB* __restrict__ proof,
		cuVecU* __restrict__ units,
		const cuVecU* __restrict__ pVars,
		const Byte* __restrict__ eliminated)
	{
		gridtype tid = global_tx;
		uint32* sh_cls = SharedMemory<uint32>();
		while (tid < pVars->size()) {
			const uint32 x = (*pVars)[tid];
			assert(x);
			assert(!ELIMINATED(eliminated[x]));
			const uint32 p = V2L(x), n = NEG(p);
			OT& ot = *otptr;
			if (ot[p].size() <= dc_options[0] && ot[n].size() <= dc_options[0])
				subsume_x(p, *cnf, ot[p], ot[n], units, proof, sh_cls + threadIdx.x * SH_MAX_SUB_IN);
			tid += stride_x;
		}
	}

	__global__ void bce_k(
		CNF* __restrict__ cnf,
		OT* __restrict__ ot,
		cuVecB* __restrict__ proof,
		cuVecU* __restrict__ resolved,
		const cuVecU* __restrict__ pVars,
		const Byte* __restrict__ eliminated)
	{
		gridtype tid = global_tx;
		__shared__ uint32 sh_cls[BLBCE * SH_MAX_BCE_IN];
		while (tid < pVars->size()) {
			const uint32 x = (*pVars)[tid];
			assert(x);
			assert(!ELIMINATED(eliminated[x]));
			const uint32 p = V2L(x), n = NEG(p);
			if ((*ot)[p].size() <= dc_options[1] && (*ot)[n].size() <= dc_options[1])
				blocked_x(x, *cnf, (*ot)[p], (*ot)[n], resolved, proof, sh_cls + threadIdx.x * SH_MAX_BCE_IN);
			tid += stride_x;
		}
	}

	__global__ void ere_k(
		CNF* __restrict__ cnf,
		OT* __restrict__ ot,
		cuVecB* __restrict__ proof,
		const cuVecU* __restrict__ pVars,
		const Byte* __restrict__ eliminated)
	{
		gridtype gid = global_ty;
		uint32* smem = SharedMemory<uint32>();
		while (gid < pVars->size()) {
			const uint32 v = pVars->at(gid);
			assert(v);
			assert(!ELIMINATED(eliminated[v]));
			const uint32 p = V2L(v), n = NEG(p);
			OL& poss = (*ot)[p], & negs = (*ot)[n];
			// do merging and apply forward equality check (on-the-fly) over resolvents
			if (poss.size() <= dc_options[2] && negs.size() <= dc_options[2]) {
				forall_occurs(poss, i) {
					SCLAUSE& pos = (*cnf)[*i];
					if (pos.deleted()) continue;
					forall_occurs(negs, j) {
						SCLAUSE& neg = (*cnf)[*j];
						if (neg.deleted()) continue;
						uint32* m_c = smem + threadIdx.y * SH_MAX_ERE_OUT; // shared memory for resolvent
						int m_len = 0;
						if ((m_len = merge_ere(v, pos, neg, m_c)) > 1) {
							CL_ST type = (pos.learnt() || neg.learnt()) ? LEARNT : ORIGINAL;
							forward_equ(*cnf, *ot, proof, m_c, m_len, type);
						}
					}
				}
			}
			gid += stride_y;
		}
	}
	//======================================================//
	//                GPU Wrappers Definitions              //
	//======================================================//
	void initDevOpts(const cuOptions& cuopt)
	{
		CHECK(cudaMemcpyToSymbol(dc_options, &cuopt.options, sizeof(uint32) * NLIMITS, 0, cudaMemcpyHostToDevice));
	}

	void initDevVorg(const cuHist& cuhist)
	{
		DCPTR ptrs = { cuhist.d_vorg, cuhist.d_lbyte };
		CHECK(cudaMemcpyToSymbol(dc_ptrs, &ptrs, sizeof(DCPTR), 0, cudaMemcpyHostToDevice));
	}

	void initSharedMem()
	{
#if defined(EXTSHMEM)
		if (devProp.major < 7) PFLOGE("extending shared memory size is not supported");
		PFLOGN2(2, " Setting maximum shared memory to 64KB..");
		int maxbytes = 65536; // 64 KB
		cudaFuncSetAttribute(in_ve_k_1, cudaFuncAttributeMaxDynamicSharedMemorySize, maxbytes);
		cudaFuncSetAttribute(ve_k_1, cudaFuncAttributeMaxDynamicSharedMemorySize, maxbytes);
		cudaFuncSetAttribute(ve_k_2, cudaFuncAttributeMaxDynamicSharedMemorySize, maxbytes);
		cudaFuncSetAttribute(sub_k, cudaFuncAttributeMaxDynamicSharedMemorySize, maxbytes);
		cudaFuncSetAttribute(ere_k, cudaFuncAttributeMaxDynamicSharedMemorySize, maxbytes);
		maxGPUSharedMem = maxbytes;
		PFLDONE(2, 5);
#endif
	}

	void printConstants()
	{
		printCMem << <1, 1 >> > ();
		LOGERR("Printing constant memory failed");
		syncAll();
	}

	void cuMemSetAsync(addr_t mem, const Byte& val, const size_t& size)
	{
		OPTIMIZEBLOCKS(uint32(size), BLOCK1D);
		memset_k<Byte> << <nBlocks, BLOCK1D >> > (mem, val, size);
		if (gopts.sync_always) {
			LOGERR("CUDA memory set failed");
			syncAll();
		}
	}

	void prepareCNFAsync(CNF* cnf, const cudaStream_t& _s)
	{
		assert(inf.nClauses);
		if (gopts.profile_gpu) cutimer->start(_s);
		OPTIMIZEBLOCKS(inf.nClauses, BLOCK1D);
		prep_cnf_k << <nBlocks, BLOCK1D, 0, _s >> > (cnf);
		if (gopts.sync_always) {
			LOGERR("Preparing CNF failed");
			syncAll();
		}
		if (gopts.profile_gpu) cutimer->stop(_s), cutimer->sig += cutimer->gpuTime();
	}

	void mapFrozenAsync(VARS* vars, const uint32& size)
	{
		assert(vars->varcore == vars->eligible); // an alies of eligible
		if (gopts.profile_gpu) cutimer->start();
		OPTIMIZEBLOCKS(size, BLOCK1D);
		// 'vars->scores' is an alies for frozen vars on the GPU side
		mapfrozen_k << <nBlocks, BLOCK1D >> > (vars->scores, vars->varcore, size);
		if (gopts.sync_always) {
			LOGERR("Mapping frozen failed");
			syncAll();
		}
		if (gopts.profile_gpu) cutimer->stop(), cutimer->ve += cutimer->gpuTime();
	}
	//=======================
	// histogram related
	//=======================
	void copyIf(uint32* dest, CNF* src)
	{
		if (gopts.profile_gpu) cutimer->start();
		reset_counter << <1, 1 >> > ();
		OPTIMIZEBLOCKS(inf.nClauses, BLOCK1D);
		copy_if_k << <nBlocks, BLOCK1D >> > (dest, src);
		check_counter << <1, 1 >> > (inf.nLiterals);
		if (gopts.profile_gpu) cutimer->stop(), cutimer->vo += cutimer->gpuTime();
		LOGERR("Copying literals failed");
		syncAll();
	}

	void copyIfAsync(uint32* dest, CNF* src)
	{
		if (gopts.profile_gpu) cutimer->start();
		reset_counter << <1, 1 >> > ();
		OPTIMIZEBLOCKS(inf.nClauses, BLOCK1D);
		copy_if_k << <nBlocks, BLOCK1D >> > (dest, src);
		if (gopts.sync_always) {
			check_counter << <1, 1 >> > (inf.nLiterals);
			LOGERR("Copying literals failed");
			syncAll();
		}
		if (gopts.profile_gpu) cutimer->stop(), cutimer->vo += cutimer->gpuTime();
	}

	void calcScores(VARS* vars, uint32* hist)
	{
		if (gopts.profile_gpu) cutimer->start();
		OPTIMIZEBLOCKS(inf.maxVar, BLOCK1D);
		assign_scores << <nBlocks, BLOCK1D >> > (vars->eligible, vars->scores, hist, inf.maxVar);
		if (gopts.profile_gpu) cutimer->stop(), cutimer->vo += cutimer->gpuTime();
		LOGERR("Assigning scores failed");
		syncAll();
	}

	void calcScores(VARS* vars, uint32* hist, OT* ot)
	{
		if (gopts.profile_gpu) cutimer->start();
		OPTIMIZEBLOCKS(inf.maxVar, BLOCK1D);
		assign_scores << <nBlocks, BLOCK1D >> > (vars->eligible, vars->scores, hist, ot, inf.maxVar);
		if (gopts.profile_gpu) cutimer->stop(), cutimer->vo += cutimer->gpuTime();
		LOGERR("Assigning scores failed");
		syncAll();
	}
	//=======================
	// CNF measurements
	//=======================
	template <typename D, typename I>
	inline D seqreduceBlocks(const D* blocks, const I& n)
	{
		D finalcount = 0;
		for (I i = 0; i < n; i++)
			finalcount += blocks[i];
		return finalcount;
	}

	template <typename D, typename I>
	inline void seqreduceBlocks(const D* CBlocks, const D* LBlocks, const I& n)
	{
		inf.n_cls_after = 0;
		inf.n_lits_after = 0;
		for (I i = 0; i < n; i++) {
			inf.n_cls_after += CBlocks[i];
			inf.n_lits_after += LBlocks[i];
		}
	}

	uint32 cuPROOF::count(const uint32* literals, const uint32& numLits)
	{
		if (!proof.checkFile()) PFLOGE("host proof system is not activated");
		if (!literals) return 0;
		if (!numLits) return 0;
		enabled = true;
		PFLOGN2(2, " Counting proof bytes..");
		OPTIMIZEBLOCKS2(numLits, BLOCK1D);
		OPTIMIZESHARED(BLOCK1D, sizeof(uint32));
		syncAll(); // sync any pending kernels or transfers
		if (gopts.profile_gpu) cutimer->start();
		cnt_proof << <nBlocks, BLOCK1D, smemSize >> > (literals, numLits);
		LOGERR("Proof counting failed");
		CHECK(cudaMemcpyFromSymbol(hostLBlocks, devLBlocks, nBlocks * sizeof(uint32)));
		if (gopts.profile_gpu) cutimer->stop(), cutimer->ve += cutimer->gpuTime();
		const uint32 maxcap = seqreduceBlocks(hostLBlocks, nBlocks);
		assert(maxcap && maxcap < (numLits * sizeof(uint32)));
		PFLENDING(2, 5, "(%d bytes)", maxcap);
		return maxcap;
	}

	void parcountCls(CNF* cnf)
	{
		const uint32 cnf_sz = inf.nClauses;
		OPTIMIZEBLOCKS2(cnf_sz, BLOCK1D);
		OPTIMIZESHARED(BLOCK1D, sizeof(uint32));
		cnt_cls << <nBlocks, BLOCK1D, smemSize >> > (cnf);
		CHECK(cudaMemcpyFromSymbol(hostCBlocks, devCBlocks, nBlocks * sizeof(uint32)));
		inf.n_cls_after = seqreduceBlocks(hostCBlocks, nBlocks);
	}

	void parcountLits(CNF* cnf)
	{
		const uint32 cnf_sz = inf.nClauses;
		OPTIMIZEBLOCKS2(cnf_sz, BLOCK1D);
		OPTIMIZESHARED(BLOCK1D, sizeof(uint32));
		cnt_lits << <nBlocks, BLOCK1D, smemSize >> > (cnf);
		CHECK(cudaMemcpyFromSymbol(hostLBlocks, devLBlocks, nBlocks * sizeof(uint32)));
		inf.n_lits_after = seqreduceBlocks(hostLBlocks, nBlocks);
	}

	void parcountAll(CNF* cnf)
	{
		const uint32 cnf_sz = inf.nClauses + (inf.nClauses >> 1);
		OPTIMIZEBLOCKS2(cnf_sz, BLOCK1D);
		OPTIMIZESHARED(BLOCK1D, sizeof(uint32) * 2);
		cnt_cls_lits << <nBlocks, BLOCK1D, smemSize >> > (cnf);
		CHECK(cudaMemcpyFromSymbol(hostCBlocks, devCBlocks, nBlocks * sizeof(uint32)));
		CHECK(cudaMemcpyFromSymbol(hostLBlocks, devLBlocks, nBlocks * sizeof(uint32)));
		seqreduceBlocks(hostCBlocks, hostLBlocks, nBlocks);
	}
	//=======================
	// occurrence table
	//=======================
	void reduceOTAsync(CNF* cnf, OT* ot, const bool& p)
	{
		assert(cnf);
		assert(ot);
		if (gopts.profile_gpu) cutimer->start();
		OPTIMIZEBLOCKS(inf.nDualVars, BLOCK1D);
		reduce_ot << <nBlocks, BLOCK1D >> > (cnf, ot);
		if (p || gopts.sync_always) {
			LOGERR("Occurrence table reduction failed");
			syncAll();
			if (p) {
				PFLRULER('=', 30);
				PFLOG0("\toccurrence table");
				ot->print();
				PFLRULER('=', 30);
			}
		}
		if (gopts.profile_gpu) cutimer->stop(), cutimer->rot += cutimer->gpuTime();
	}

	inline void resetOTAsync(CNF* cnf, OT* ot)
	{
		assert(cnf);
		assert(ot);
		OPTIMIZEBLOCKS(inf.nDualVars, BLOCK1D);
		reset_ot_k << <nBlocks, BLOCK1D >> > (ot);
		if (gopts.sync_always) {
			LOGERR("Occurrence table reset failed");
			syncAll();
			assert(ot->accViolation());
		}
	}

	void createOTAsync(CNF* cnf, OT* ot, const bool& p)
	{
		assert(cnf);
		assert(ot);
		if (gopts.profile_gpu) cutimer->start();
		resetOTAsync(cnf, ot);
		OPTIMIZEBLOCKS(inf.nClauses, BLOCK1D);
		create_ot_k << <nBlocks, BLOCK1D >> > (cnf, ot);
		if (p || gopts.sync_always) {
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
		if (gopts.profile_gpu) cutimer->stop(), cutimer->cot += cutimer->gpuTime();
	}

	void sortOTAsync(CNF* cnf, OT* ot, VARS* vars)
	{
		assert(cnf);
		assert(ot);
		assert(vars->numPVs);
		OPTIMIZEBLOCKS(vars->numPVs, BLSORT);
		if (gopts.profile_gpu) cutimer->start();
		sort_ot_p << <nBlocks, BLSORT >> > (cnf, ot, vars->pVars);
		sort_ot_n << <nBlocks, BLSORT >> > (cnf, ot, vars->pVars);
		if (gopts.sync_always) {
			LOGERR("Sorting OT failed");
			syncAll();
		}
		if (gopts.profile_gpu) cutimer->stop(), cutimer->sot += cutimer->gpuTime();
	}
	//=======================
	// variable elimination
	//=======================
	inline void vePhase1(
		CNF* cnf,
		OT* ot,
		VARS* vars,
		cuVecB* proof,
		uint32* ucnt,
		uint32* type,
		uint32* rpos,
		S_REF* rref,
		const bool& in)
	{
		PFLOGN2(2, "  resetting last eliminated ID to -1.. ");
		reset_id << <1, 1 >> > ();
		PFLDONE(2, 5);
		PFLOGN2(2, "  configuring phase 1 with ");
		OPTIMIZEBLOCKSELIM(vars->numPVs, BLVE1, gopts.ve);
		OPTIMIZESHARED(nThreads, SH_MAX_BVE_OUT1 * sizeof(uint32));
		PFLENDING(2, 5, "(%d/%d ths, %d/%d bls) and %zd KB shared memory",
			nThreads, BLVE1, nBlocks, MAXBLOCKS, smemSize / KBYTE);
#if VE_DBG
		if (in) in_ve_k_1 << <1, 1, smemSize >> > (cnf, ot, vars->pVars, vars->eliminated, vars->units, vars->resolved, proof, vars->varcore, ucnt, type, rpos, rref);
		else	   ve_k_1 << <1, 1, smemSize >> > (cnf, ot, vars->pVars, vars->eliminated, vars->units, vars->resolved, proof, vars->varcore, ucnt, type, rpos, rref);
#else
		if (in) in_ve_k_1 << <nBlocks, nThreads, smemSize >> > (cnf, ot, vars->pVars, vars->eliminated, vars->units, vars->resolved, proof, vars->varcore, ucnt, type, rpos, rref);
		else	   ve_k_1 << <nBlocks, nThreads, smemSize >> > (cnf, ot, vars->pVars, vars->eliminated, vars->units, vars->resolved, proof, vars->varcore, ucnt, type, rpos, rref);
#endif
		LOGERR("BVE Phase-1 failed");
		sync();
	}

	inline void vePhase2(
		VARS* vars,
		uint32* rpos,
		S_REF* rref,
		cudaStream_t* streams,
		cuMM& cumm)
	{
		const uint32 cs_offset = cumm.pinnedCNF()->size();
		const S_REF data_offset = cumm.pinnedCNF()->data().size;
		size_t tb1 = 0, tb2 = 0;
		DeviceScan::ExclusiveScan(NULL, tb1, rpos, rpos, Sum(), cs_offset, vars->numPVs);
		DeviceScan::ExclusiveScan(NULL, tb2, rref, rref, Sum(), data_offset, vars->numPVs);
		size_t tmpcap = tb1 + tb2;
		addr_t ts1 = NULL, ts2 = NULL;
		addr_t tmpmem = (addr_t) ((tmpcap > cumm.scatterCap()) ? cacher.allocate(tmpcap) : cumm.scatter());
		ts1 = tmpmem, ts2 = ts1 + tb1;
		DeviceScan::ExclusiveScan(ts1, tb1, rpos, rpos, Sum(), cs_offset, vars->numPVs, streams[0]);
		DeviceScan::ExclusiveScan(ts2, tb2, rref, rref, Sum(), data_offset, vars->numPVs, streams[1]);
		LOGERR("BVE Phase-2 failed");
		sync(streams[0]);
		sync(streams[1]);
		if (tmpcap > cumm.scatterCap()) {
			assert(tmpmem != (addr_t)cumm.scatter());
			cacher.deallocate(tmpmem);
		}
	}

	inline void vePhase3(
		CNF* cnf,
		OT* ot,
		VARS* vars,
		cuVecB* proof,
		uint32* ucnt,
		uint32* type,
		uint32* rpos,
		S_REF* rref)
	{
		PFLOGN2(2, "  configuring phase 3 with ");
		OPTIMIZEBLOCKSELIM(vars->numPVs, BLVE2, gopts.ve);
		OPTIMIZESHARED(nThreads, SH_MAX_BVE_OUT2 * sizeof(uint32));
		PFLENDING(2, 5, "(%d/%d ths, %d/%d bls) and %zd KB shared memory",
			nThreads, BLVE2, nBlocks, MAXBLOCKS, smemSize / KBYTE);
#if VE_DBG
		ve_k_2 << <1, 1, smemSize >> > (cnf, ot, vars->pVars, vars->eliminated, vars->eligible, vars->units, vars->resolved, proof, ucnt, type, rpos, rref);
#else
		ve_k_2 << <nBlocks, nThreads, smemSize >> > (cnf, ot, vars->pVars, vars->eliminated, vars->eligible, vars->units, vars->resolved, proof, ucnt, type, rpos, rref);
#endif
		if (gopts.sync_always) {
			LOGERR("BVE Phase-3 failed");
			sync();
		}
	}

	void veAsync(
		CNF* cnf,
		OT* ot,
		VARS* vars,
		cudaStream_t* streams,
		cuVecB* proof,
		cuMM& cumm,
		const cuHist& cuhist,
		const bool& in)
	{
		assert(cnf);
		assert(ot);
		assert(vars->numPVs);
		if (gopts.profile_gpu) cutimer->start();
		if (!gopts.ve_atomic) {
			S_REF* rref = cuhist.d_segs;
			uint32* type = cuhist.d_hist;
			uint32* rpos = type + inf.maxVar;
			uint32* ucnt = cumm.resizeLits(inf.maxVar);
			assert(ucnt); // inf.maxVar cannot be larger than nr. of literals
			vePhase1(cnf, ot, vars, proof, ucnt, type, rpos, rref, in);
			vePhase2(vars, rpos, rref, streams, cumm);
			vePhase3(cnf, ot, vars, proof, ucnt, type, rpos, rref);
		}
		else {
#if VE_DBG
			ve_k << <1, 1 >> > (cnf, ot, vars->pVars, vars->units, vars->resolved, proof);
#else
			OPTIMIZEBLOCKS(vars->numPVs, BLVE);
			ve_k << <nBlocks, BLVE >> > (cnf, ot, vars->pVars, vars->eliminated, vars->units, vars->resolved, proof);
#endif
			LOGERR("Atomic BVE failed");
			syncAll();
		}
		if (gopts.profile_gpu) cutimer->stop(), cutimer->ve += cutimer->gpuTime();
	}

	void veResizeCNFAsync(CNF* cnf, const cuHist& cuhist)
	{
		S_REF* rref = cuhist.d_segs;
		uint32* type = cuhist.d_hist, * rpos = type + inf.maxVar;
		resizeCNF_k << <1, 1 >> > (cnf, type, rpos, rref, verbose);
		if (gopts.sync_always) {
			LOGERR("Resizing CNF after BVE failed");
			sync();
		}
	}

	void subAsync(CNF* cnf, OT* ot, VARS* vars, cuVecB* proof)
	{
		assert(cnf);
		assert(ot);
		assert(vars->numPVs);
		if (gopts.profile_gpu) cutimer->start();
		PFLOGN2(2, "  configuring SUB kernel with ");
		OPTIMIZEBLOCKSELIM(vars->numPVs, BLSUB, gopts.sub);
		OPTIMIZESHARED(nThreads, SH_MAX_SUB_IN * sizeof(uint32));
		PFLENDING(2, 5, "(%d/%d ths, %d/%d bls) and %zd KB shared memory",
			nThreads, BLSUB, nBlocks, MAXBLOCKS, smemSize / KBYTE);
#if SS_DBG
		sub_k << <1, 1, smemSize >> > (cnf, ot, proof, vars->units, vars->pVars, vars->eliminated);
#else
		sub_k << <nBlocks, nThreads, smemSize >> > (cnf, ot, proof, vars->units, vars->pVars, vars->eliminated);
#endif
		if (gopts.profile_gpu) cutimer->stop(), cutimer->sub += cutimer->gpuTime();
		if (gopts.sync_always) {
			LOGERR("SUB Elimination failed");
			syncAll();
		}
	}

	void bceAsync(CNF* cnf, OT* ot, VARS* vars, cuVecB* proof)
	{
		assert(cnf);
		assert(ot);
		assert(vars->numPVs);
		if (gopts.profile_gpu) cutimer->start();
		OPTIMIZEBLOCKS(vars->numPVs, BLBCE);
		bce_k << <nBlocks, BLBCE >> > (cnf, ot, proof, vars->resolved, vars->pVars, vars->eliminated);
		if (gopts.profile_gpu) cutimer->stop(), cutimer->bce += cutimer->gpuTime();
		if (gopts.sync_always) {
			LOGERR("BCE Elimination failed");
			syncAll();
		}
	}

	void ereAsync(CNF* cnf, OT* ot, VARS* vars, cuVecB* proof)
	{
		assert(cnf);
		assert(ot);
		assert(vars->numPVs);
		if (gopts.profile_gpu) cutimer->start();
#if	defined(_DEBUG) || defined(DEBUG) || !defined(NDEBUG)
		dim3 block2D(16, devProp.warpSize);
#else 
		dim3 block2D(devProp.warpSize, devProp.warpSize);
#endif
		dim3 grid2D(1, 1, 1);
		PFLOGN2(2, "  configuring ERE kernel with ");
		OPTIMIZEBLOCKSERE(vars->numPVs, block2D, gopts.ere);
		OPTIMIZESHARED(block2D.y, SH_MAX_ERE_OUT * sizeof(uint32));
		grid2D.y = nBlocks;
		PFLENDING(2, 5, "(%d/%d ths, %d/%d bls) and %zd KB shared memory",
			block2D.y, devProp.warpSize, grid2D.y, MAXBLOCKS, smemSize / KBYTE);
		ere_k << <grid2D, block2D, smemSize >> > (cnf, ot, proof, vars->pVars, vars->eliminated);
		if (gopts.profile_gpu) cutimer->stop(), cutimer->ere += cutimer->gpuTime();
		if (gopts.sync_always) {
			LOGERR("ERE Elimination failed");
			syncAll();
		}
	}

}