/***********************************************************************[bounded.cuh]
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

#ifndef __GPU_BVE_
#define __GPU_BVE_

#include "and.cuh"
#include "xor.cuh"
#include "function.cuh"
#include "resolve.cuh"
#include "ifthenelse.cuh"
#include "equivalence.cuh"

namespace pFROST {

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
	//===========================================//

	// atomic approach
	_PFROST_D_ void resolve_x(
		const uint32& x, 
		const uint32& nElements, 
		const uint32& nAddedCls, 
		const uint32& nAddedLits, 	
		CNF& cnf, 
		OL& me,
		OL& other,
		cuVecU* units, 
		cuVecB* proof, 
		uint32* out_c)
	{
		assert(x);
		const int nbuckets = dc_nbuckets;
		uint32 checksum = 0;
		S_REF ref;
		S_REF* refs = cnf.jump(ref, nAddedCls, nAddedLits);
		uint32* ustart = RECOVERADDEDUNITS(nElements) ? units->jump(RECOVERADDEDUNITS(nElements)) : NULL;
		uint32* udata = ustart;
		const uint32 proofBytes = RECOVERADDEDPROOF(nElements);
		addr_t pdata = (proof && proofBytes) ? proof->jump(proofBytes) : NULL;
		#pragma unroll
		forall_occurs(me, i) {
			SCLAUSE& ci = cnf[*i];
			if (ci.learnt()) continue;
			S_REF* otherend = other.end();
			#pragma unroll
			for (S_REF* j = other; j != otherend && checksum < nAddedCls; j++) {
				SCLAUSE& cj = cnf[*j];
				if (cj.learnt()) continue;
				if ((ci.size() + cj.size() - 2) > SH_MAX_BVE_OUT) { // use global memory
					uint32 unit;
					int rsize;
					if (rsize = merge(x, ci, cj, unit)) {
						if (rsize == 1) {
							assert(unit > 1);
							assert(udata);
							*udata++ = unit;
							if (proof) saveProofUnit(pdata, unit);
						}
						else {
							assert(rsize > 1);
							assert(!unit);
							SCLAUSE* added = new (cnf.cref(ref)) SCLAUSE();
							merge(x, ci, cj, added);
							assert(added->size() == rsize);
							added->markAdded();
							refs[checksum++] = ref, ref += rsize + nbuckets;
							if (proof) saveProofClause(pdata, *added, PROOF_ADDED);
						}
					}
				}
				else { // use shared memory
					int rsize = 0;
					if (rsize = merge(x, ci, cj, out_c)) {
						assert(rsize <= SH_MAX_BVE_OUT);
						if (rsize == 1) {
							assert(*out_c > 1);
							assert(udata);
							*udata++ = *out_c;
							if (proof) saveProofUnit(pdata, *out_c);
						}
						else {
							assert(rsize > 1);
							uint32 sig = 0;
							calcSig(out_c, rsize, sig);
							SCLAUSE* added = cnf.cref(ref);
							added = new (added) SCLAUSE(out_c, rsize);
							added->set_sig(sig);
							added->markAdded();
							refs[checksum++] = ref, ref += rsize + nbuckets;
							assert(added->isSorted());
							assert(added->hasZero() < 0);
							if (proof) saveProofClause(pdata, out_c, rsize, PROOF_ADDED);
						}
					}
				}
			}
		}
		assert(checksum == nAddedCls);
		assert(udata - RECOVERADDEDUNITS(nElements) == ustart);
	}

	_PFROST_D_ void substitute_x(
		const uint32& x,
		const uint32& nElements, 
		const uint32& nAddedCls, 
		const uint32& nAddedLits, 
		CNF& cnf,
		OL& me, 
		OL& other, 
		cuVecU* units, 
		cuVecB* proof,
		uint32* out_c)
	{
		assert(x);
		const int nbuckets = dc_nbuckets;
		uint32 checksum = 0;
		S_REF ref;
		S_REF* refs = cnf.jump(ref, nAddedCls, nAddedLits);
		uint32* ustart = RECOVERADDEDUNITS(nElements) ? units->jump(RECOVERADDEDUNITS(nElements)) : NULL;
		uint32* udata = ustart;
		const uint32 proofBytes = RECOVERADDEDPROOF(nElements);
		addr_t pdata = (proof && proofBytes) ? proof->jump(proofBytes) : NULL;
		#pragma unroll
		forall_occurs(me, i) {
			SCLAUSE& ci = cnf[*i];
			if (ci.learnt()) continue;
			S_REF* otherend = other.end();
			#pragma unroll
			for (S_REF* j = other; j != otherend && checksum < nAddedCls; j++) {
				SCLAUSE& cj = cnf[*j];
				if (cj.learnt()) continue;
				if (ci.molten() == cj.molten()) continue;
				if ((ci.size() + cj.size() - 2) > SH_MAX_BVE_OUT) { // use global memory
					uint32 unit;
					int rsize;
					if (rsize = merge(x, ci, cj, unit)) {
						if (rsize == 1) {
							assert(unit > 1);
							assert(udata);
							*udata++ = unit;
							if (proof) saveProofUnit(pdata, unit);
						}
						else {
							assert(rsize > 1);
							assert(!unit);
							SCLAUSE* added = new (cnf.cref(ref)) SCLAUSE();
							merge(x, ci, cj, added);
							assert(added->size() == rsize);
							added->markAdded();
							refs[checksum++] = ref, ref += rsize + nbuckets;
							if (proof) saveProofClause(pdata, *added, PROOF_ADDED);
						}
					}
				}
				else { // use shared memory
					int rsize = 0;
					if (rsize = merge(x, ci, cj, out_c)) {
						assert(rsize <= SH_MAX_BVE_OUT);
						if (rsize == 1) {
							assert(*out_c > 1);
							assert(udata);
							*udata++ = *out_c;
							if (proof) saveProofUnit(pdata, *out_c);
						}
						else {
							assert(rsize > 1);
							uint32 sig = 0;
							calcSig(out_c, rsize, sig);
							SCLAUSE* added = cnf.cref(ref);
							added = new (added) SCLAUSE(out_c, rsize);
							added->set_sig(sig);
							added->markAdded();
							refs[checksum++] = ref, ref += rsize + nbuckets;
							assert(added->isSorted());
							assert(added->hasZero() < 0);
							if (proof) saveProofClause(pdata, out_c, rsize, PROOF_ADDED);
						}
					}
				}
			}
		}
		assert(checksum == nAddedCls);
		assert(udata - RECOVERADDEDUNITS(nElements) == ustart);
	}


	// 3-phase approach
	_PFROST_D_ bool memorySafe(
		const uint32& tid,
		const uint32& x,
		const uint32& nAddedCls,
		const uint32& nAddedLits,
		const uint32& addedPos, 
		const S_REF& addedRef,
		CNF* cnf)
	{
		#if MEMORYSAFE_DBG

		if ((addedPos + nAddedCls) > cnf->refs().capacity()) {
			printf("c Memory violation to clauses capacity(%d) @(tid = %d, v = %d): added clauses: %d, added position = %d\n", cnf->refs().capacity(), tid, x, nAddedCls, addedPos);
			return false;
		}
		const S_REF data_size = addedRef + nAddedLits + dc_nbuckets * nAddedCls;
		if (data_size > cnf->data().cap) {
			const S_REF nAddedBuckets = nAddedLits + dc_nbuckets * nAddedCls;
			printf("c Memory violation to data capacity(%lld) @(tid = %d, v = %d): added buckets: %lld, added reference = %lld\n", cnf->data().cap, tid, x, nAddedBuckets, addedRef);
			return false;
		}
		return true;

		#else

		return ((addedPos + nAddedCls) <= cnf->refs().capacity())
			&& ((addedRef + nAddedLits + dc_nbuckets * nAddedCls) <= cnf->data().cap);

		#endif
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
		const int nbuckets = dc_nbuckets;
		const uint32 checksum = addedPos + nAddedCls;
		S_REF newref = addedRef;
		S_REF* refs = cnf.refsData();
		S_REF* endme = me.end();
		S_REF* endother = other.end();
		uint32* ustart = RECOVERADDEDUNITS(nElements) ? units->jump(RECOVERADDEDUNITS(nElements)) : NULL;
		uint32* udata = ustart;
		const uint32 proofBytes = RECOVERADDEDPROOF(nElements);
		addr_t pstart = (proof && proofBytes) ? proof->jump(proofBytes) : NULL;
		addr_t pdata = pstart;

		#if RES_DBG
		printf("c  Resolving %d, proof bytes = %d, units = %d, added = %d, addedPos = %d, addedRef = %lld:\n", 
			x, proofBytes, RECOVERADDEDUNITS(nElements), nAddedCls, addedPos, addedRef);
		#endif

		for (S_REF* i = me; i != endme && addedPos < checksum; i++) {
			SCLAUSE& ci = cnf[*i];
			if (ci.learnt()) continue;
			for (S_REF* j = other; j != endother && addedPos < checksum; j++) {
				SCLAUSE& cj = cnf[*j];
				if (cj.learnt()) continue;
				int rsize;
				if ((ci.size() + cj.size() - 2) > SH_MAX_BVE_OUT2) {
					uint32 unit;
					if (rsize = merge(x, ci, cj, unit)) {  // must use a dummy "merge" first to avoid data racing on global memory
						if (rsize == 1) {
							assert(unit > 1);
							assert(udata);
							*udata++ = unit;
							if (proof) saveProofUnit(pdata, unit);
						}
						else {
							assert(rsize > 1);
							assert(!unit);
							SCLAUSE* added = new (cnf.cref(newref)) SCLAUSE();
							merge(x, ci, cj, added);
							assert(added->size() == rsize);
							added->markAdded();
							refs[addedPos++] = newref, newref += rsize + nbuckets;
							if (proof) saveProofClause(pdata, *added, PROOF_ADDED);
							#if RES_DBG
							printf("c  C(%d, r: %lld)->", addedPos - 1, newref - added->blockSize()), added->print();
							#endif
						}
					}
				}
				else if (rsize = merge(x, ci, cj, out_c)) {  // use shared memory
					assert(rsize <= SH_MAX_BVE_OUT2);
					if (rsize == 1) {
						assert(*out_c > 1);
						assert(udata);
						*udata++ = *out_c;
						if (proof) saveProofUnit(pdata, *out_c);
					}
					else {
						uint32 sig = 0;
						calcSig(out_c, rsize, sig);
						SCLAUSE* added = cnf.cref(newref);
						added = new (added) SCLAUSE(out_c, rsize);
						added->set_sig(sig);
						added->markAdded();
						refs[addedPos++] = newref, newref += rsize + nbuckets;
						assert(added->isSorted());
						assert(added->hasZero() < 0);
						if (proof) saveProofClause(pdata, out_c, rsize, PROOF_ADDED);
						#if RES_DBG
						printf("c  C(%d, r: %lld)->", addedPos - 1, newref - added->blockSize()), added->print();
						#endif	
					}
				}
			}
		}
		assert(checksum == addedPos);
		assert((addedRef + S_REF(nAddedLits) + S_REF(nbuckets * nAddedCls)) == newref);
		assert(udata - RECOVERADDEDUNITS(nElements) == ustart);
		assert(pdata - proofBytes == pstart);
		// delete resolved
		toblivion(cnf, me, other);
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
		const int nbuckets = dc_nbuckets;
		const uint32 checksum = addedPos + nAddedCls;
		S_REF newref = addedRef;
		S_REF* refs = cnf.refsData();
		S_REF* endme = me.end();
		S_REF* endother = other.end();
		uint32* ustart = RECOVERADDEDUNITS(nElements) ? units->jump(RECOVERADDEDUNITS(nElements)) : NULL;
		uint32* udata = ustart;
		const uint32 proofBytes = RECOVERADDEDPROOF(nElements);
		addr_t pstart = (proof && proofBytes) ? proof->jump(proofBytes) : NULL;
		addr_t pdata = pstart;

		#if SUBST_DBG
		printf("c  Substituting %d, proof bytes = %d, units = %d, added = %d, addedPos = %d, addedRef = %lld:\n", 
			x, proofBytes, RECOVERADDEDUNITS(nElements), nAddedCls, addedPos, addedRef);
		#endif

		for (S_REF* i = me; i != endme && addedPos < checksum; i++) {
			SCLAUSE& ci = cnf[*i];
			if (ci.learnt()) continue;
			for (S_REF* j = other; j != endother && addedPos < checksum; j++) {
				SCLAUSE& cj = cnf[*j];
				if (cj.learnt()) continue;
				if (ci.molten() == cj.molten()) continue;
				int rsize;
				if ((ci.size() + cj.size() - 2) > SH_MAX_BVE_OUT2) {
					uint32 unit;
					if (rsize = merge(x, ci, cj, unit)) {  // must use "merge" check first to avoid data racing on global memory
						if (rsize == 1) {
							assert(unit > 1);
							assert(udata);
							*udata++ = unit;
							if (proof) saveProofUnit(pdata, unit);
						}
						else {
							assert(rsize > 1);
							assert(!unit);
							SCLAUSE* added = new (cnf.cref(newref)) SCLAUSE();
							merge(x, ci, cj, added);
							assert(added->size() == rsize);
							added->markAdded();
							refs[addedPos++] = newref, newref += rsize + nbuckets;
							if (proof) saveProofClause(pdata, *added, PROOF_ADDED);
							#if SUBST_DBG
							printf("c  C(%d, r: %lld)->", addedPos - 1, newref - added->blockSize()), added->print();
							#endif
						}
					}
				}
				else if (rsize = merge(x, ci, cj, out_c)) {  // use shared memory
					assert(rsize <= SH_MAX_BVE_OUT2);
					if (rsize == 1) {
						assert(*out_c > 1);
						assert(udata);
						*udata++ = *out_c;
						if (proof) saveProofUnit(pdata, *out_c);
					}
					else {
						assert(rsize > 1);
						uint32 sig = 0;
						calcSig(out_c, rsize, sig);
						SCLAUSE* added = cnf.cref(newref);
						added = new (added) SCLAUSE(out_c, rsize);
						added->set_sig(sig);
						added->markAdded();
						refs[addedPos++] = newref, newref += rsize + nbuckets;
						assert(added->isSorted());
						assert(added->hasZero() < 0);
						if (proof) saveProofClause(pdata, out_c, rsize, PROOF_ADDED);
						#if SUBST_DBG
						printf("c  C(%d, r: %lld)->", addedPos - 1, newref - added->blockSize()), added->print();
						#endif
					}
				}
			}
		}
		assert(checksum == addedPos);
		assert(udata - RECOVERADDEDUNITS(nElements) == ustart);
		assert(pdata - proofBytes == pstart);
		assert((addedRef + S_REF(nAddedLits) + S_REF(nbuckets * nAddedCls)) == newref);
		// delete resolved
		toblivion(cnf, me, other);
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
		const int nbuckets = dc_nbuckets;
		const uint32 checksum = addedPos + nAddedCls;
		S_REF newref = addedRef;
		S_REF* refs = cnf.refsData();
		S_REF* endme = me.end();
		S_REF* endother = other.end();
		uint32* ustart = RECOVERADDEDUNITS(nElements) ? units->jump(RECOVERADDEDUNITS(nElements)) : NULL;
		uint32* udata = ustart;
		const uint32 proofBytes = RECOVERADDEDPROOF(nElements);
		addr_t pstart = (proof && proofBytes) ? proof->jump(proofBytes) : NULL;
		addr_t pdata = pstart;

		#if CORE_DBG
		printf("c  Core substituting %d, proof bytes = %d, units = %d, added = %d, addedPos = %d, addedRef = %lld:\n", 
			x, proofBytes, RECOVERADDEDUNITS(nElements), nAddedCls, addedPos, addedRef);
		#endif

		for (S_REF* i = me; i != endme && addedPos < checksum; i++) {
			SCLAUSE& ci = cnf[*i];
			if (ci.learnt()) continue;
			for (S_REF* j = other; j != endother && addedPos < checksum; j++) {
				SCLAUSE& cj = cnf[*j];
				if (cj.learnt() || (ci.molten() && cj.molten())) continue;
				int rsize;
				if ((ci.size() + cj.size() - 2) > SH_MAX_BVE_OUT2) {
					uint32 unit;
					if (rsize = merge(x, ci, cj, unit)) {
						if (rsize == 1) {
							assert(unit > 1);
							assert(udata);
							*udata++ = unit;
							if (proof) saveProofUnit(pdata, unit);
						}
						else {
							assert(rsize > 1);
							assert(!unit);
							SCLAUSE* added = new (cnf.cref(newref)) SCLAUSE();
							merge(x, ci, cj, added);
							assert(added->size() == rsize);
							added->markAdded();
							refs[addedPos++] = newref, newref += rsize + nbuckets;
							if (proof) saveProofClause(pdata, *added, PROOF_ADDED);
							#if CORE_DBG
							printf("c  C(%d, r: %lld)->", addedPos - 1, newref - added->blockSize()), added->print();
							#endif
						}
					}
				}
				else if (rsize = merge(x, ci, cj, out_c)) { // use shared memory
					assert(rsize <= SH_MAX_BVE_OUT2);
					if (rsize == 1) {
						assert(*out_c > 1);
						assert(udata);
						*udata++ = *out_c;
						if (proof) saveProofUnit(pdata, *out_c);
					}
					else {
						assert(rsize > 1);
						uint32 sig = 0;
						calcSig(out_c, rsize, sig);
						SCLAUSE* added = cnf.cref(newref);
						added = new (added) SCLAUSE(out_c, rsize);
						added->set_sig(sig);
						added->markAdded();
						refs[addedPos++] = newref, newref += rsize + nbuckets;
						assert(added->isSorted());
						assert(added->hasZero() < 0);
						if (proof) saveProofClause(pdata, out_c, rsize, PROOF_ADDED);
						#if CORE_DBG
						printf("c  C(%d, r: %lld)->", addedPos - 1, newref - added->blockSize()), added->print();
						#endif
					}
				}
			}
		}
		assert(checksum == addedPos);
		assert(udata - RECOVERADDEDUNITS(nElements) == ustart);
		assert(pdata - proofBytes == pstart);
		assert((addedRef + S_REF(nAddedLits) + S_REF(nbuckets * nAddedCls)) == newref);
		// delete resolved
		toblivion(cnf, me, other);
	}

} // parafrost namespace


#endif