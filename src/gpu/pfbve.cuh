/***********************************************************************[pfbve.cuh]
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

#ifndef __SIGMA_BVE_
#define __SIGMA_BVE_

#include "pfdevice.cuh"
#include "pfresolve.cuh"
#include "pfand.cuh"
#include "pfequ.cuh"
#include "pfite.cuh"
#include "pfxor.cuh"

namespace pFROST {

	namespace SIGmA {

		// atomic approach
		_PFROST_D_ void resolve_x(const uint32& x, const uint32& nAddedCls, const uint32& nAddedLits, CNF& cnf, OL& me, OL& other, cuVecU* units, uint32* out_c)
		{
			assert(x);
			const int nbuckets = dc_nbuckets;
			uint32 checksum = 0;
			S_REF ref;
			S_REF* cs = cnf.jump(ref, nAddedCls, nAddedLits);
#pragma unroll
			for (S_REF* i = me; i != me.end(); i++) {
				SCLAUSE& ci = cnf[*i];
				if (ci.learnt()) continue;
#pragma unroll
				for (S_REF* j = other; j != other.end() && checksum < nAddedCls; j++) {
					SCLAUSE& cj = cnf[*j];
					if (cj.learnt()) continue;
					if ((ci.size() + cj.size() - 2) > SH_MAX_BVE_OUT) { // use global memory
						uint32 unit;
						int rsize;
						if (rsize = merge(x, ci, cj, unit)) {
							if (rsize == 1) {
								assert(unit);
								units->push(unit);
							}
							else {
								assert(rsize > 1);
								assert(!unit);
								SCLAUSE* added = new (cnf.cref(ref)) SCLAUSE();
								merge(x, ci, cj, added);
								assert(added->size() == rsize);
								added->markAdded();
								cs[checksum++] = ref, ref += rsize + nbuckets;
							}
						}
					}
					else { // use shared memory
						int merged_sz = 0;
						if (merged_sz = merge(x, ci, cj, out_c)) {
							assert(merged_sz <= SH_MAX_BVE_OUT);
							if (merged_sz == 1) units->push(*out_c);
							else {
								assert(merged_sz > 1);
								uint32 sig = 0;
								calcSig(out_c, merged_sz, sig);
								SCLAUSE* added = cnf.cref(ref);
								added = new (added) SCLAUSE(out_c, merged_sz);
								added->set_sig(sig);
								added->markAdded();
								cs[checksum++] = ref, ref += merged_sz + nbuckets;
								assert(added->isSorted());
								assert(added->hasZero() < 0);
							}
						}
					}
				}
			}
			assert(checksum == nAddedCls);
		}

		_PFROST_D_ void substitute_x(const uint32& x, const uint32& nAddedCls, const uint32& nAddedLits, CNF& cnf, OL& me, OL& other, cuVecU* units, uint32* out_c)
		{
			assert(x);
			const int nbuckets = dc_nbuckets;
			uint32 checksum = 0;
			S_REF ref;
			S_REF* cs = cnf.jump(ref, nAddedCls, nAddedLits);
#pragma unroll
			for (S_REF* i = me; i != me.end(); i++) {
				SCLAUSE& ci = cnf[*i];
				if (ci.learnt()) continue;
#pragma unroll
				for (S_REF* j = other; j != other.end() && checksum < nAddedCls; j++) {
					SCLAUSE& cj = cnf[*j];
					if (cj.learnt()) continue;
					if (ci.molten() == cj.molten()) continue;
					if ((ci.size() + cj.size() - 2) > SH_MAX_BVE_OUT) { // use global memory
						uint32 unit;
						int rsize;
						if (rsize = merge(x, ci, cj, unit)) {
							if (rsize == 1) {
								assert(unit);
								units->push(unit);
							}
							else {
								assert(rsize > 1);
								assert(!unit);
								SCLAUSE* added = new (cnf.cref(ref)) SCLAUSE();
								merge(x, ci, cj, added);
								assert(added->size() == rsize);
								added->markAdded();
								cs[checksum++] = ref, ref += rsize + nbuckets;
							}
						}
					}
					else { // use shared memory
						int merged_sz = 0;
						if (merged_sz = merge(x, ci, cj, out_c)) {
							assert(merged_sz <= SH_MAX_BVE_OUT);
							if (merged_sz == 1) units->push(*out_c);
							else {
								assert(merged_sz > 1);
								uint32 sig = 0;
								calcSig(out_c, merged_sz, sig);
								SCLAUSE* added = cnf.cref(ref);
								added = new (added) SCLAUSE(out_c, merged_sz);
								added->set_sig(sig);
								added->markAdded();
								cs[checksum++] = ref, ref += merged_sz + nbuckets;
								assert(added->isSorted());
								assert(added->hasZero() < 0);
							}
						}
					}
				}
			}
			assert(checksum == nAddedCls);
		}


		// 3-phase approach
		_PFROST_D_ bool memorySafe(const uint32& tid, const uint32& x, const uint32& nAddedCls, const uint32& nAddedLits, const uint32& addedPos, const S_REF& addedRef, CNF* cnf)
		{
			const uint32 cs_size = addedPos + nAddedCls;
			if (cs_size > cnf->refs().capacity()) {
				#if VE_DBG
				printf("Memory violation to clauses capacity(%d) @(tid = %d, v = %d): added clauses: %d, added position = %d\n", cnf->refs().capacity(), tid, x, nAddedCls, addedPos);
				#endif
				return false;
			}
			const S_REF data_size = addedRef + nAddedLits + dc_nbuckets * nAddedCls;
			if (data_size > cnf->data().cap) {
				#if VE_DBG
				const S_REF nAddedBuckets = nAddedLits + dc_nbuckets * nAddedCls;
				printf("Memory violation to data capacity(%lld) @(tid = %d, v = %d): added buckets: %lld, added reference = %lld\n", cnf->data().cap, tid, x, nAddedBuckets, addedRef);
				#endif
				return false;
			}
			return true;
		}

		_PFROST_D_ void freezeClauses(CNF& cnf, OL& poss, OL& negs)
		{
#pragma unroll
			forall_occurs(poss, i) {
				SCLAUSE& c = cnf[*i];
				if (c.original() && c.molten())
					c.freeze();
			}
#pragma unroll
			forall_occurs(negs, i) {
				SCLAUSE& c = cnf[*i];
				if (c.original() && c.molten())
					c.freeze();
			}
		}

		_PFROST_D_ void resolve_x(const uint32& x, const uint32& nAddedCls, const uint32& nAddedLits, uint32 addedPos, S_REF addedRef,
			CNF& cnf, OL& me, OL& other, cuVecU* units, uint32* out_c)
		{
			assert(x);
			assert(nAddedCls);
			assert(nAddedLits);
#if VE_DBG
			printf("c  Resolving %d, added = %d, addedPos = %d, addedRef = %lld:\n", x, nAddedCls, addedPos, addedRef);
#endif
			const int nbuckets = dc_nbuckets;
			const uint32 checksum = addedPos + nAddedCls;
			S_REF newref = addedRef;
			S_REF* cs = cnf.refsData();
			for (S_REF* i = me; i != me.end() && addedPos < checksum; i++) {
				SCLAUSE& ci = cnf[*i];
				if (ci.learnt()) continue;
				for (S_REF* j = other; j != other.end() && addedPos < checksum; j++) {
					SCLAUSE& cj = cnf[*j];
					if (cj.learnt()) continue;
					if ((ci.size() + cj.size() - 2) > SH_MAX_BVE_OUT2) {
						// must use "merge" check first to avoid data racing on global memory
						uint32 unit;
						int rsize;
						if (rsize = merge(x, ci, cj, unit)) {
							if (rsize == 1) {
								assert(unit);
								units->push(unit);
							}
							else {
								assert(rsize > 1);
								assert(!unit);
								SCLAUSE* added = new (cnf.cref(newref)) SCLAUSE();
								merge(x, ci, cj, added);
								assert(added->size() == rsize);
								added->markAdded();
								cs[addedPos++] = newref, newref += rsize + nbuckets;
#if VE_DBG
								printf("c  C(%d, r: %lld)->", addedPos - 1, newref - added->blockSize()), added->print();

#endif
							}
						}
					}
					else { // use shared memory
						int merged_sz;
						if (merged_sz = merge(x, ci, cj, out_c)) {
							assert(merged_sz <= SH_MAX_BVE_OUT2);
							if (merged_sz == 1) units->push(*out_c);
							else {
								assert(merged_sz > 1);
								uint32 sig = 0;
								calcSig(out_c, merged_sz, sig);
								SCLAUSE* added = cnf.cref(newref);
								added = new (added) SCLAUSE(out_c, merged_sz);
								added->set_sig(sig);
								added->markAdded();
								cs[addedPos++] = newref, newref += merged_sz + nbuckets;
								assert(added->isSorted());
								assert(added->hasZero() < 0);
#if VE_DBG
								printf("c  C(%d, r: %lld)->", addedPos - 1, newref - added->blockSize()), added->print();
#endif
							}
						}
					}
				}
			}
			assert(checksum == addedPos);
			assert((addedRef + S_REF(nAddedLits) + S_REF(dc_nbuckets * nAddedCls)) == newref);
			// delete resolved
			toblivion(cnf, me, other);
		}

		_PFROST_D_ void substitute_x(const uint32& x, const uint32& nAddedCls, const uint32& nAddedLits, uint32 addedPos, S_REF addedRef,
			CNF& cnf, OL& me, OL& other, cuVecU* units, uint32* out_c)
		{
			assert(x);
			assert(nAddedCls);
			assert(nAddedLits);
#if VE_DBG
			printf("c  Substituting %d, added = %d, addedPos = %d, addedRef = %lld:\n", x, nAddedCls, addedPos, addedRef);
#endif
			const int nbuckets = dc_nbuckets;
			const uint32 checksum = addedPos + nAddedCls;
			S_REF newref = addedRef;
			S_REF* cs = cnf.refsData();
			for (S_REF* i = me; i != me.end() && addedPos < checksum; i++) {
				SCLAUSE& ci = cnf[*i];
				if (ci.learnt()) continue;
				for (S_REF* j = other; j != other.end() && addedPos < checksum; j++) {
					SCLAUSE& cj = cnf[*j];
					if (cj.learnt()) continue;
					if (ci.molten() == cj.molten()) continue;
					if ((ci.size() + cj.size() - 2) > SH_MAX_BVE_OUT2) {
						// must use "merge" check first to avoid data racing on global memory
						uint32 unit;
						int rsize;
						if (rsize = merge(x, ci, cj, unit)) {
							if (rsize == 1) {
								assert(unit);
								units->push(unit);
							}
							else {
								assert(rsize > 1);
								assert(!unit);
								SCLAUSE* added = new (cnf.cref(newref)) SCLAUSE();
								merge(x, ci, cj, added);
								assert(added->size() == rsize);
								added->markAdded();
								cs[addedPos++] = newref, newref += rsize + nbuckets;
#if VE_DBG
								printf("c  C(%d, r: %lld)->", addedPos - 1, newref - added->blockSize()), added->print();
#endif
							}
						}
					}
					else { // use shared memory
						int merged_sz;
						if (merged_sz = merge(x, ci, cj, out_c)) {
							assert(merged_sz <= SH_MAX_BVE_OUT2);
							if (merged_sz == 1) units->push(*out_c);
							else {
								assert(merged_sz > 1);
								uint32 sig = 0;
								calcSig(out_c, merged_sz, sig);
								SCLAUSE* added = cnf.cref(newref);
								added = new (added) SCLAUSE(out_c, merged_sz);
								added->set_sig(sig);
								added->markAdded();
								cs[addedPos++] = newref, newref += merged_sz + nbuckets;
								assert(added->isSorted());
								assert(added->hasZero() < 0);
#if VE_DBG
								printf("c  C(%d, r: %lld)->", addedPos - 1, newref - added->blockSize()), added->print();
#endif
							}
						}
					}
				}
			}
			assert(checksum == addedPos);
			assert((addedRef + S_REF(nAddedLits) + S_REF(dc_nbuckets * nAddedCls)) == newref);
			// delete resolved
			toblivion(cnf, me, other);
		}


	} // sigma namespace
} // parafrost namespace


#endif