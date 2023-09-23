/***********************************************************************[xor.cuh]
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

#ifndef __GPU_XOR_
#define __GPU_XOR_

#include "elimination.cuh"

namespace ParaFROST {

	#define XOR_DBG 0

	#define COUNTFLIPS(X) \
	{ \
		while (__popc(++X) & 1); \
	}

	_PFROST_D_ void freeze_arities(CNF& cnf, OL& list)
	{
		forall_occurs(list, i) {
			SCLAUSE& c = cnf[*i];
			if (c.size() > 2 && c.molten()) 
				c.freeze();
		}
	}

	_PFROST_D_ void freeze_arities(CNF& cnf, OL& me, OL& other)
	{
		forall_occurs(me, i) {
			SCLAUSE& c = cnf[*i];
			if (c.size() > 2 && c.molten()) 
				c.freeze();
		}
		forall_occurs(other, i) {
			SCLAUSE& c = cnf[*i];
			if (c.size() > 2 && c.molten()) 
				c.freeze();
		}
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

	_PFROST_D_ bool find_xor_gate(
		const uint32& dx, OL& dx_list, 
		const uint32& fx, OL& fx_list,
		const uint32& nOrgCls,
		CNF& cnf,
		OT& ot, 
		uint32* out_c, 
		uint32& nElements,
		uint32& nAddedCls,
		uint32& nAddedLits)
	{
		assert(dx > 1);
		assert(dx_list.size());
		assert(fx_list.size());

		if (cnf[dx_list.back()].size() == 2 ||
			cnf[fx_list.back()].size() == 2)
			return false;

		const int maxarity = kOpts->xor_max_arity;

		int arity = cnf[*dx_list].size() - 1;
		if (arity > maxarity) return false;

		assert(checkMolten(cnf, ot[dx], ot[fx]));

		const uint32 v = ABS(dx);

		forall_occurs(dx_list, i) {
			SCLAUSE& ci = cnf[*i];
			if (ci.original()) {
				const int size = ci.size();
				const int arity = size - 1; // XOR arity
				if (size < 3 || arity > maxarity) continue;

				// share to out_c
				uint32* shared = out_c;
				forall_clause(ci, k)
					*shared++ = *k;

				// find arity clauses
				uint32 parity = 0;
				int itargets = (1 << arity);
				while (--itargets && makeArity(cnf, ot, parity, out_c, size));

				assert(parity < (1UL << size)); // overflow check
				assert(itargets >= 0);

				if (itargets)
					freeze_arities(cnf, dx_list, fx_list);
				else {

					ci.melt();

					// check resolvability
					nElements = 0, nAddedCls = 0, nAddedLits = 0;
					if (countSubstituted(v, nOrgCls, cnf, dx_list, fx_list, nElements, nAddedCls, nAddedLits)) {
						freeze_arities(cnf, dx_list, fx_list);
						break;
					}

					// can be substituted
					#if XOR_DBG
					printf("c  Gate %d = XOR(", ABS(dx));
					for (int k = 0; k < arity; k++) {
						printf(" %d", ABS(out_c[k]));
						if (k < arity - 1) printf(",");
					}
					printf(" ) found ==> added = %d, deleted = %d\n", nAddedCls, dx_list.size() + fx_list.size());
					printGate(cnf, dx_list, fx_list);
					#endif

					return true;
				}
			} // ci-original block
		} // dx_list-for block

		return false;
	}

} 


#endif