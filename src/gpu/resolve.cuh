/***********************************************************************[resolve.cuh]
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

#ifndef __RESOLVE_
#define __RESOLVE_

#include "elimination.cuh"

namespace ParaFROST {

	#define RES_DBG 0

_PFROST_D_ bool countResolvents(
		const uint32& x,
		CNF& cnf,
		OL& me,
		OL& other,
		uint32& nElements,
		uint32& nAddedCls,
		uint32& nAddedLits)
	{
		assert(x);
		assert(!nElements);
		assert(!nAddedCls);
		assert(!nAddedLits);
		const int rlimit = kOpts->ve_clause_max;
		// check if proof bytes has to be calculated
		if (kOpts->proof_en) {
			uint32 proofBytes = 0;
			forall_occurs(me, i) {
				SCLAUSE& ci = cnf[*i];
				if (ci.learnt()) continue;
				forall_occurs(other, j) {
					SCLAUSE& cj = cnf[*j];
					if (cj.learnt()) continue;
					const int rsize = mergeProof(x, ci, cj, proofBytes);
					if (rsize == 1)
						nElements++;
					else if (rsize) {
						#if RES_UNBOUNDED
						++nAddedCls;
						#else
						if (rlimit && rsize > rlimit) return false;
						++nAddedCls;
						#endif
						nAddedLits += rsize;
					}
				}
			}

			// GUARD for compressed proof size and #units
			if (nElements > ADDEDCLS_MAX || proofBytes > ADDEDPROOF_MAX)
				return false;

			nElements = ENCODEPROOFINFO(nElements, proofBytes);

			#if PROOF_DBG
			printf("c  Variable %d: counted %d units and %d proof bytes\n", x, nElements, proofBytes);
			#endif
		}
		else { // no proof
			forall_occurs(me, i) {
				SCLAUSE& ci = cnf[*i];
				if (ci.learnt()) continue;
				forall_occurs(other, j) {
					SCLAUSE& cj = cnf[*j];
					if (cj.learnt()) continue;
					const int rsize = merge(x, ci, cj);
					if (rsize == 1)
						nElements++;
					else if (rsize) {
						#if RES_UNBOUNDED
						++nAddedCls;
						#else
						if (rlimit && rsize > rlimit) return false;
						++nAddedCls;
						#endif
						nAddedLits += rsize;
					}
				}
			}
		}

		// GUARD for compressed variable limits
		if (nAddedCls > ADDEDCLS_MAX || nAddedLits > ADDEDLITS_MAX)
			return false;

		#if RES_DBG
		printf("c  Resolving(%d) ==> added = %d, deleted = %d\n", x, nAddedCls, me.size() + other.size());
		print_clause_set(cnf, me, other);
		#endif	

		return true;
	}

	_PFROST_D_ bool countResolvents(
		const uint32& x,
		const uint32& nClsBefore,
		CNF& cnf,
		OL& me,
		OL& other,
		uint32& nElements,
		uint32& nAddedCls,
		uint32& nAddedLits)
	{
		assert(x);
		assert(nClsBefore);
		assert(checkMolten(cnf, me, other));

		nElements = 0, nAddedCls = 0, nAddedLits = 0;

		const int rlimit = kOpts->ve_clause_max;

		// check if proof bytes has to be calculated
		if (kOpts->proof_en) {
			uint32 proofBytes = 0;
			forall_occurs(me, i) {
				SCLAUSE& ci = cnf[*i];
				if (ci.learnt()) continue;
				forall_occurs(other, j) {
					SCLAUSE& cj = cnf[*j];
					if (cj.learnt()) continue;
					const int rsize = mergeProof(x, ci, cj, proofBytes);
					if (rsize == 1)
						nElements++;
					else if (rsize) {
						#if RES_UNBOUNDED
						++nAddedCls;
						#else
						if (++nAddedCls > nClsBefore || (rlimit && rsize > rlimit)) return false;
						#endif
						nAddedLits += rsize;
					}
				}
			}

			// GUARD for compressed proof size and #units
			if (nElements > ADDEDCLS_MAX || proofBytes > ADDEDPROOF_MAX) return false;
			nElements = ENCODEPROOFINFO(nElements, proofBytes);

			#if PROOF_DBG
			printf("c  Variable %d: counted %d units and %d proof bytes\n", x, nElements, proofBytes);
			#endif
		}
		else { // no proof
			forall_occurs(me, i) {
				SCLAUSE& ci = cnf[*i];
				if (ci.learnt()) continue;
				forall_occurs(other, j) {
					SCLAUSE& cj = cnf[*j];
					if (cj.learnt()) continue;
					const int rsize = merge(x, ci, cj);
					if (rsize == 1)
						nElements++;
					else if (rsize) {
						#if RES_UNBOUNDED
						++nAddedCls;
						#else
						if (++nAddedCls > nClsBefore || (rlimit && rsize > rlimit)) return false;
						#endif
						nAddedLits += rsize;
					}
				}
			}
		}

		// GUARD for compressed variable limits
		if (nAddedCls > ADDEDCLS_MAX || nAddedLits > ADDEDLITS_MAX) return false;

		// check bound on literals
		if (kOpts->ve_lbound_en) {
			uint32 nLitsBefore = 0;
			countLitsBefore(cnf, me, nLitsBefore);
			countLitsBefore(cnf, other, nLitsBefore);
			if (nAddedLits > nLitsBefore) return false;
		}

		#if RES_DBG
		printf("c  Resolving(%d) ==> added = %d, deleted = %d\n", x, nAddedCls, me.size() + other.size());
		print_clause_set(cnf, me, other);
		#endif	

		return true;
	}

}

#endif