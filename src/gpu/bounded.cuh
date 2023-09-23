/***********************************************************************[bounded.cuh]
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

#ifndef __BVE_
#define __BVE_

#include "and.cuh"
#include "xor.cuh"
#include "resolve.cuh"
#include "function.cuh"
#include "ifthenelse.cuh"
#include "equivalence.cuh"

namespace ParaFROST {

	//=========== Debugging macros =============//
	// each header above has its own flag

	#define MEMORYSAFE_DBG 0

	#if AND_DBG || XOR_DBG || ITE_DBG || EQU_DBG
	#define SUBST_DBG 1
	#else
	#define SUBST_DBG 0
	#endif

	#if FUN_DBG
	#define CORE_DBG 1
	#else
	#define CORE_DBG 0
	#endif

	#if SUBST_DBG || CORE_DBG || RES_DBG
	#define VE_DBG 1
	#else 
	#define VE_DBG 0
	#endif

	#if VE_DBG
		_PFROST_D_ void printResolvents(const uint32& addedPos, const S_REF& newref, const SCLAUSE& added) 
		{
			#if VE_DBG
				printf("c  C(%d, r: %lld)->", addedPos - 1, newref - (added->size() + DC_NBUCKETS));
				added->print();
			#endif
		}
	#else
		#define printResolvents(DUMMY, ...)
	#endif
	
	//==========================================//

	__device__ int   lastEliminatedID;

	#define ADD_RESOLVENT \
	{ \
		int rsize; \
		if ((ci.size() + cj.size() - 2) > SH_MAX_BVE_OUT2) { \
			uint32 unit; \
			/* must use a dummy "merge" first to avoid data racing on global memory */ \
			if (rsize = merge(x, ci, cj, unit)) {  \
				if (rsize == 1) { \
					assert(unit > 1); \
					assert(udata); \
					*udata++ = unit; \
					if (proof) \
						saveProofUnit(pdata, unit); \
				} \
				else { \
					assert(rsize > 1); \
					assert(!unit); \
					SCLAUSE* added = new (cnf.cref(newref)) SCLAUSE(); \
					merge(x, ci, cj, added); \
					assert(added->size() == rsize); \
					added->markAdded(); \
					refs[addedPos++] = newref, newref += rsize + nbuckets; \
					if (proof) \
						saveProofClause(pdata, *added, PROOF_ADDED); \
					printResolvents(addedPos, newref, added); \
				} \
			} \
		} \
		/* use shared memory */ \
		else if (rsize = merge(x, ci, cj, out_c)) {  \
			assert(rsize <= SH_MAX_BVE_OUT2); \
			if (rsize == 1) { \
				assert(*out_c > 1); \
				assert(udata); \
				*udata++ = *out_c; \
				if (proof) \
					saveProofUnit(pdata, *out_c); \
			} \
			else { \
				uint32 sig = 0; \
				calcSig(out_c, rsize, sig); \
				SCLAUSE* added = cnf.cref(newref); \
				added = new (added) SCLAUSE(out_c, rsize); \
				added->set_sig(sig); \
				added->markAdded(); \
				refs[addedPos++] = newref, newref += rsize + nbuckets; \
				assert(added->isSorted()); \
				assert(added->hasZero() < 0); \
				if (proof) \
					saveProofClause(pdata, out_c, rsize, PROOF_ADDED); \
				printResolvents(addedPos, newref, added); \
			} \
		} \
	}

	#define NEW_RESOLVENTS(REG_NBUCKETS,CHECKSUM,NEWREF,REFS,ENDME,ENDOTHER,USTART,UDATA,PSTART,PDATA,PROOFBYTES) \
		const int REG_NBUCKETS = DC_NBUCKETS; \
		const uint32 CHECKSUM = addedPos + nAddedCls; \
		S_REF NEWREF = addedRef; \
		S_REF* REFS = cnf.refsData(); \
		S_REF* ENDME = me.end(); \
		S_REF* ENDOTHER = other.end(); \
		uint32* USTART = RECOVERADDEDUNITS(nElements) ? units->jump(RECOVERADDEDUNITS(nElements)) : NULL; \
		uint32* UDATA = USTART; \
		const uint32 PROOFBYTES = RECOVERADDEDPROOF(nElements); \
		addr_t PSTART = (proof && PROOFBYTES) ? proof->jump(PROOFBYTES) : NULL; \
		addr_t PDATA = PSTART; \

	#define DELETE_RESOLVED \
	{ \
		assert(checksum == addedPos); \
		assert(udata - RECOVERADDEDUNITS(nElements) == ustart); \
		assert(pdata - proofbytes == pstart); \
		assert((addedRef + S_REF(nAddedLits) + S_REF(nbuckets * nAddedCls)) == newref); \
		toblivion(cnf, me, other); \
	}

	_PFROST_D_ void resolve_x(
		const uint32& x, 
		const uint32& nElements, 
		const uint32& nAddedCls, 
		const uint32& nAddedLits,		
		uint32 addedPos,
		S_REF addedRef,
		CNF& cnf, 
		OL& me, 
		OL& other, 
		cuVecU* units,
		cuVecB* proof,
		uint32* out_c)
	{
		assert(x);
		assert(nAddedCls);
		assert(nAddedLits);
		
		NEW_RESOLVENTS(nbuckets, checksum, newref, refs, endme, endother, ustart, udata, pstart, pdata, proofbytes);

		#if RES_DBG
		printf("c  Resolving %d, proof bytes = %d, units = %d, added = %d, addedPos = %d, addedRef = %lld:\n", 
			x, proofbytes, RECOVERADDEDUNITS(nElements), nAddedCls, addedPos, addedRef);
		#endif

		for (S_REF* i = me; i != endme && addedPos < checksum; i++) {
			SCLAUSE& ci = cnf[*i];
			if (ci.learnt()) continue;
			for (S_REF* j = other; j != endother && addedPos < checksum; j++) {
				SCLAUSE& cj = cnf[*j];
				if (cj.learnt()) continue;
				ADD_RESOLVENT;
			}
		}

		DELETE_RESOLVED;
	}

	_PFROST_D_ void substitute_x(
		const uint32& x, 
		const uint32& nElements, 
		const uint32& nAddedCls, 
		const uint32& nAddedLits,
		uint32 addedPos,
		S_REF addedRef,
		CNF& cnf,
		OL& me, 
		OL& other,
		cuVecU* units, 
		cuVecB* proof,
		uint32* out_c)
	{
		assert(x);
		assert(nAddedCls);
		assert(nAddedLits);

		NEW_RESOLVENTS(nbuckets, checksum, newref, refs, endme, endother, ustart, udata, pstart, pdata, proofbytes);

		#if SUBST_DBG
		printf("c  Substituting %d, proof bytes = %d, units = %d, added = %d, addedPos = %d, addedRef = %lld:\n", 
			x, proofbytes, RECOVERADDEDUNITS(nElements), nAddedCls, addedPos, addedRef);
		#endif

		for (S_REF* i = me; i != endme && addedPos < checksum; i++) {
			SCLAUSE& ci = cnf[*i];
			if (ci.learnt()) continue;
			for (S_REF* j = other; j != endother && addedPos < checksum; j++) {
				SCLAUSE& cj = cnf[*j];
				if (cj.learnt()) continue;
				if (ci.molten() == cj.molten()) continue;
				ADD_RESOLVENT;
			}
		}

		DELETE_RESOLVED;
	}

	_PFROST_D_ void coresubstitute_x(
		const uint32& x, 
		const uint32& nElements, 
		const uint32& nAddedCls, 
		const uint32& nAddedLits,
		uint32 addedPos,
		S_REF addedRef,
		CNF& cnf,
		OL& me, 
		OL& other,
		cuVecU* units, 
		cuVecB* proof,
		uint32* out_c)
	{
		assert(x);
		assert(nAddedCls);
		assert(nAddedLits);
		
		NEW_RESOLVENTS(nbuckets, checksum, newref, refs, endme, endother, ustart, udata, pstart, pdata, proofbytes);

		#if CORE_DBG
		printf("c  Core substituting %d, proof bytes = %d, units = %d, added = %d, addedPos = %d, addedRef = %lld:\n", 
			x, proofbytes, RECOVERADDEDUNITS(nElements), nAddedCls, addedPos, addedRef);
		#endif

		for (S_REF* i = me; i != endme && addedPos < checksum; i++) {
			SCLAUSE& ci = cnf[*i];
			if (ci.learnt()) continue;
			for (S_REF* j = other; j != endother && addedPos < checksum; j++) {
				SCLAUSE& cj = cnf[*j];
				if (cj.learnt() || (ci.molten() && cj.molten())) continue;
				ADD_RESOLVENT;
			}
		}
		
		DELETE_RESOLVED;
	}

	//=========================================================//
	// kernels
	//=========================================================//

	__global__ void reset_id() { lastEliminatedID = -1; }

	__global__ void print_id() { printf("c lastEliminatedID = %d\n", lastEliminatedID); }
	
	__global__ void mapfrozen_k(const uint32* __restrict__ frozen, uint32* __restrict__ varcore, const uint32 size)
	{
		grid_t tid = global_tx;
		while (tid < size) {
			assert(frozen[tid] && frozen[tid] < NOVAR);
			varcore[frozen[tid]] = tid;
			tid += stride_x;
		}
	}

	// Macros for checking applicability of variable elimination
	#define MEMORY_SAFE_DBG \
		if ((addedPos + nAddedCls) > cnf->refs().capacity()) { \
				printf("c Memory violation to clauses capacity(%d) @(tid = %d, v = %d): added clauses: %d, added position = %d\n", cnf->refs().capacity(), tid, x, nAddedCls, addedPos); \
				assert(0); \
		} \
		const S_REF data_size = addedRef + nAddedLits + DC_NBUCKETS * nAddedCls; \
		if (data_size > cnf->data().cap) { \
			const S_REF nAddedBuckets = nAddedLits + DC_NBUCKETS * nAddedCls; \
			printf("c Memory violation to data capacity(%lld) @(tid = %d, v = %d): added buckets: %lld, added reference = %lld\n", cnf->data().cap, tid, x, nAddedBuckets, addedRef); \
			assert(0); \
		}


	#define MEMORY_SAFE ((addedPos + nAddedCls) <= cnf->refs().capacity()) && \
						((addedRef + nAddedLits + DC_NBUCKETS * nAddedCls) <= cnf->data().cap)

	_PFROST_D_ void variable_elimination(
					  const grid_t& tid,
					  const uint32& x,
					  const uint32& p,
					  const uint32& n,
					  const uint32& pOrgs,
					  const uint32& nOrgs,
							CNF&    cnf,
							OT&		ot,
							OL&     poss,
							OL&     negs,
							cuVecU* units,
							cuVecU* resolved,
							cuVecB* proof,
					  const uint32* varcore,
							uint32* outs,
							uint32* ucnt,
							uint32* type,
							uint32* rpos,
							S_REF*  rref,
							Byte*   eliminated)
	{
		uint32 elimType = 0, nElements = 0, nAddedCls = 0, nAddedLits = 0;

		/* pure-literal elimination */
		if (!pOrgs || !nOrgs) {
			toblivion(p, n, pOrgs, nOrgs, cnf, poss, negs, resolved);
			ucnt[tid] = 0, type[tid] = 0, rref[tid] = 0, rpos[tid] = 0, ELIMINATE(eliminated[x]);
		}
		/* Equiv/NOT-gate Reasoning */
		else if (uint32 def = find_equ_gate(p, n, cnf, poss, negs)) {
			saveResolved(p, n, pOrgs, nOrgs, cnf, poss, negs, resolved);
			substitute_single(p, n, def, cnf, poss, negs, units, proof);
			ucnt[tid] = 0, type[tid] = 0, rref[tid] = 0, rpos[tid] = 0, ELIMINATE(eliminated[x]);
		}
		/* simple resolution case */
		else if ((pOrgs == 1 || nOrgs == 1) && countResolvents(x, cnf, poss, negs, nElements, nAddedCls, nAddedLits)) {
			/* can be eliminated with resolvents in next phase */
			if (nAddedCls) {
				assert(nAddedLits >= nAddedCls);
				assert(elimType <= TYPE_MASK);
				assert(nAddedCls <= ADDEDCLS_MAX);
				assert(nAddedLits <= ADDEDLITS_MAX);
				/* save elimination info. */
				type[tid] = ENCODEVARINFO(elimType, nAddedCls, nAddedLits);
				ucnt[tid] = nElements, rpos[tid] = nAddedCls, rref[tid] = nAddedLits + DC_NBUCKETS * nAddedCls;
			}
			else {
				toblivion(p, n, pOrgs, nOrgs, cnf, poss, negs, resolved);
				ucnt[tid] = 0, type[tid] = 0, rref[tid] = 0, rpos[tid] = 0, ELIMINATE(eliminated[x]);
			}
		}
		else {
			assert(pOrgs && nOrgs);
			const uint32 nClsBefore = pOrgs + nOrgs;
			elimType = 0, nElements = 0, nAddedCls = 0, nAddedLits = 0;
			/*=====================*/
			/* check resolvability */
			/*=====================*/
			if (nClsBefore > 2) {
				/* AND/OR-gate Reasoning */
				if (nOrgs < SH_MAX_BVE_OUT1 &&
					find_ao_gate(n, negs, p, poss, nClsBefore, cnf, &outs[threadIdx.x * SH_MAX_BVE_OUT1], nElements, nAddedCls, nAddedLits))
					elimType = AOIX_MASK;
				else if (!nAddedCls && pOrgs < SH_MAX_BVE_OUT1 &&
					find_ao_gate(p, poss, n, negs, nClsBefore, cnf, &outs[threadIdx.x * SH_MAX_BVE_OUT1], nElements, nAddedCls, nAddedLits))
					elimType = AOIX_MASK;
			}
			if (!elimType && nClsBefore > 3) {
				/* ITE-gate Reasoning */
				if (find_ite_gate(p, poss, n, negs, nClsBefore, cnf, ot, nElements, nAddedCls, nAddedLits))
					elimType = AOIX_MASK;
				else if (!nAddedCls &&
					find_ite_gate(n, negs, p, poss, nClsBefore, cnf, ot, nElements, nAddedCls, nAddedLits))
					elimType = AOIX_MASK;
				/* XOR-gate Reasoning */
				else if (find_xor_gate(p, poss, n, negs, nClsBefore, cnf, ot, &outs[threadIdx.x * SH_MAX_BVE_OUT1], nElements, nAddedCls, nAddedLits))
					elimType = AOIX_MASK;
				else if (!nAddedCls &&
					find_xor_gate(n, negs, p, poss, nClsBefore, cnf, ot, &outs[threadIdx.x * SH_MAX_BVE_OUT1], nElements, nAddedCls, nAddedLits))
					elimType = AOIX_MASK;
			}
			/* fun-tab reasoning */
			if (varcore && !elimType && nClsBefore > 2 &&
				find_fun_gate(p, n, nClsBefore, varcore, cnf, ot, &outs[threadIdx.x * SH_MAX_BVE_OUT1], nElements, nAddedCls, nAddedLits))
				elimType = CORE_MASK;
			/* n-by-m resolution */
			else if (!elimType && !nAddedCls &&
				countResolvents(x, nClsBefore, cnf, poss, negs, nElements, nAddedCls, nAddedLits))
				elimType = RES_MASK;
			/*=====================*/
			/* check addibility    */
			/*=====================*/
			/* eliminated without resolvents */
			if (!nAddedCls) {
				toblivion(p, n, pOrgs, nOrgs, cnf, poss, negs, resolved);
				ucnt[tid] = 0, type[tid] = 0, rref[tid] = 0, rpos[tid] = 0, ELIMINATE(eliminated[x]);
			}
			/* can be eliminated with resolvents in next phase */
			else if (elimType) {
				assert(nAddedLits >= nAddedCls);
				assert(elimType <= TYPE_MASK);
				assert(nAddedCls <= ADDEDCLS_MAX);
				assert(nAddedLits <= ADDEDLITS_MAX);
				/* save elimination info. */
				type[tid] = ENCODEVARINFO(elimType, nAddedCls, nAddedLits);
				ucnt[tid] = nElements, rpos[tid] = nAddedCls, rref[tid] = nAddedLits + DC_NBUCKETS * nAddedCls;
			}
			/* cannot be eliminated */
			else
				ucnt[tid] = 0, type[tid] = 0, rref[tid] = 0, rpos[tid] = 0;
		}
	}

	// 3-phase BVE
	__global__ void in_ve_k_1(
		CNF* __restrict__ cnfptr,
		OT* __restrict__ otptr,
		const cuVecU* __restrict__ elected,
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
		grid_t tid = global_tx;
		uint32* outs = SharedMemory<uint32>();
		while (tid < elected->size()) {
			const uint32 x = elected->at(tid);
			assert(x);
			assert(!ELIMINATED(eliminated[x]));
			assert(!IS_ADDING(eliminated[x]));
			const uint32 p = V2L(x), n = NEG(p);
			CNF& cnf = *cnfptr;
			OT& ot = *otptr;
			OL& poss = ot[p], & negs = ot[n];
			assert(devIsSorted(poss.data(), poss.size(), OLIST_CMP(cnfptr)));
			assert(devIsSorted(negs.data(), negs.size(), OLIST_CMP(cnfptr)));
			// nElements: compresses '#units' and '#proof bytes'
			uint32 pOrgs = 0, nOrgs = 0;
			countOrgs(cnf, poss, pOrgs);
			countOrgs(cnf, negs, nOrgs);
			variable_elimination(tid, x, p, n, pOrgs, nOrgs, cnf, ot, poss, negs, 
								 units, resolved, proof, varcore, outs, ucnt, type, rpos, rref, eliminated);
			tid += stride_x;
		}
	}

	__global__ void ve_k_1(
		CNF* __restrict__ cnfptr,
		OT* __restrict__ otptr,
		const cuVecU* __restrict__ elected,
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
		grid_t tid = global_tx;
		uint32* outs = SharedMemory<uint32>();
		while (tid < elected->size()) {
			const uint32 x = elected->at(tid);
			assert(x);
			assert(!ELIMINATED(eliminated[x]));
			assert(!IS_ADDING(eliminated[x]));
			const uint32 p = V2L(x), n = NEG(p);
			CNF& cnf = *cnfptr;
			OT& ot = *otptr;
			OL& poss = ot[p], & negs = ot[n];
			const uint32 pOrgs = poss.size();
			const uint32 nOrgs = negs.size();
			variable_elimination(tid, x, p, n, pOrgs, nOrgs, cnf, ot, poss, negs, 
								 units, resolved, proof, varcore, outs, ucnt, type, rpos, rref, eliminated);
			tid += stride_x;
		}
	}

	__global__ void ve_k_2(
		CNF* __restrict__ cnf,
		OT* __restrict__ ot,
		const cuVecU* __restrict__ elected,
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
		grid_t tid = global_tx;
		uint32* outs = SharedMemory<uint32>();
		while (tid < elected->size()) {
			const uint32 x = elected->at(tid);
			assert(x);
			const uint32 xinfo = type[tid];
			const uint32 elimType = RECOVERTYPE(xinfo);
			assert(elimType <= TYPE_MASK);
			if (elimType) {
				assert(!ELIMINATED(eliminated[x]));
				assert(!IS_ADDING(eliminated[x]));
				const uint32 p = V2L(x), n = NEG(p);
				const uint32 nAddedCls = RECOVERADDEDCLS(xinfo);
				const uint32 nAddedLits = RECOVERADDEDLITS(xinfo);
				assert(nAddedCls && nAddedCls <= ADDEDCLS_MAX);
				assert(nAddedLits && nAddedLits <= ADDEDLITS_MAX);
				const uint32 addedPos = rpos[tid];
				const S_REF addedRef = rref[tid];
				OL& poss = (*ot)[p], & negs = (*ot)[n];
				if (MEMORY_SAFE) {
					saveResolved(p, n, *cnf, poss, negs, resolved);
					if (IS_RES(elimType))
						resolve_x(x, ucnt[tid], nAddedCls, nAddedLits, addedPos, addedRef, *cnf, poss, negs, units, proof, outs + threadIdx.x * SH_MAX_BVE_OUT2);
					else if (IS_CORE(elimType))
						coresubstitute_x(x, ucnt[tid], nAddedCls, nAddedLits, addedPos, addedRef, *cnf, poss, negs, units, proof, outs + threadIdx.x * SH_MAX_BVE_OUT2);
					else {
						assert(IS_AOIX(elimType));
						substitute_x(x, ucnt[tid], nAddedCls, nAddedLits, addedPos, addedRef, *cnf, poss, negs, units, proof, outs + threadIdx.x * SH_MAX_BVE_OUT2);
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
			const S_REF lastAddedBuckets = lastAddedLits + DC_NBUCKETS * lastAddedCls;
			const S_REF data_size = lastAddedBuckets + lastAddedRef;
			const uint32 cs_size = lastAddedCls + lastAddedPos;
			cnf->resize(data_size, cs_size);
			if (verbose > 1) printf("c   resized CNF to %d clauses and %lld data for a last ID %d\n", cs_size, data_size, lastEliminatedID);
		}
	}

} 


#endif