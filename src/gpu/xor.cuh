/***********************************************************************[xor.cuh]
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

#ifndef __SIGMA_XOR_
#define __SIGMA_XOR_

#include "device.cuh"

namespace pFROST {

	namespace SIGmA {

#define COUNTFLIPS(X) \
	{ \
		while (__popc(++X) & 1); \
	}

		_PFROST_D_ void freeze_arities(CNF& cnf, OL& list)
		{
#pragma unroll
			forall_occurs(list, i) {
				SCLAUSE& c = cnf[*i];
				if (c.size() > 2 && c.molten()) c.freeze();
			}
		}

		_PFROST_D_ void freeze_arities(CNF& cnf, OL& me, OL& other)
		{
			freeze_arities(cnf, me);
			freeze_arities(cnf, other);
		}

		_PFROST_D_ void shareXORClause(SCLAUSE& c, uint32* shared)
		{
#pragma unroll
			forall_clause(c, k)
				* shared++ = *k;
		}

		_PFROST_D_ bool checkArity(SCLAUSE& c, uint32* literals, const int& size)
		{
			assert(literals);
			const uint32* end = literals + size;
			forall_clause(c, i) {
				const uint32 lit = *i;
				uint32* j;
				for (j = literals; j != end; j++) {
					if (lit == *j) {
						break;
					}
				}
				if (j == end) return false;
			}
			return true;
		}

		_PFROST_D_ bool makeArity(CNF& cnf, OT& ot, uint32& parity, uint32* literals, const int& size)
		{
			const uint32 oldparity = parity;
			COUNTFLIPS(parity);
			for (int k = 0; k < size; k++) {
				const uint32 bit = (1UL << k);
				if (NEQUAL(parity & bit, oldparity & bit))
					literals[k] = FLIP(literals[k]);
			}
			// search for an arity clause
			assert(size > 2);
			uint32 best = *literals;
			assert(best > 1);
			int minsize = ot[best].size();
			for (int k = 1; k < size; k++) {
				const uint32 lit = literals[k];
				assert(lit > 1);
				int lsize = ot[lit].size();
				if (lsize < minsize) {
					minsize = lsize;
					best = literals[k];
				}
			}
			forall_occurs(ot[best], i) {
				SCLAUSE& c = cnf[*i];
				if (c.original() && c.size() == size && checkArity(c, literals, size)) {
					assert(c.original());
					c.melt();
					return true;
				}
			}
			return false;
		}

		_PFROST_D_ bool find_xor_gate(const uint32& dx, const uint32& nOrgCls, CNF& cnf, OT& ot, uint32* out_c, uint32& nAddedCls, uint32& nAddedLits)
		{
			assert(dx > 1);
			const int xor_max_arity = dc_limits[3];
			const uint32 fx = FLIP(dx), v = ABS(dx);
			assert(checkMolten(cnf, ot[dx], ot[fx]));
			OL& itarget = ot[dx];
			OL& otarget = ot[fx];
			forall_occurs(itarget, i) {
				SCLAUSE& ci = cnf[*i];
				if (ci.original()) {
					const int size = ci.size();
					const int arity = size - 1; // XOR arity
					if (size < 3 || arity > xor_max_arity) continue;
					assert(arity <= SH_MAX_BVE_OUT1);
					// share to out_c
					shareXORClause(ci, out_c);
					// find arity clauses
					uint32 parity = 0;
					int itargets = (1 << arity);
					while (--itargets && makeArity(cnf, ot, parity, out_c, size));
					assert(parity < (1UL << size)); // overflow check
					assert(itargets >= 0);
					if (itargets)
						freeze_arities(cnf, itarget, otarget);
					else {
						ci.melt();
						// check resolvability
						nAddedCls = 0, nAddedLits = 0;
						if (countSubstituted(v, nOrgCls, cnf, itarget, otarget, nAddedCls, nAddedLits)) {
							freeze_arities(cnf, itarget, otarget);
							break;
						}
						// can be substituted
#if VE_DBG
						printf("c  Gate %d = XOR(", ABS(dx));
						for (int k = 0; k < size; k++) {
							printf(" %d", ABS(out_c[k]));
							if (k < arity - 1) printf(",");
						}
						printf(" ) found ==> added = %d, deleted = %d\n", nAddedCls, itarget.size() + otarget.size());
						printGate(cnf, itarget, otarget);
#endif
						return true;
					}
				} // ci-original block
			} // itarget-for block
			return false;
		}

	} // sigma namespace
} // parafrost namespace


#endif