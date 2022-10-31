/***********************************************************************[equ.cuh]
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

#ifndef __EQU_
#define __EQU_

#include "elimination.cuh"

namespace ParaFROST {

	#define EQU_DBG 0

	_PFROST_D_ void substitute_single(
		const uint32& dx,
		const uint32& def,
		SCLAUSE& org, 
		uint32& nUnits)
	{
		assert(dx > 1);
		assert(def != dx);
		assert(org.original());
		int n = 0;
		#pragma unroll
		for (int i = 0; i < org.size(); i++) {
			const uint32 lit = org[i];
			if (lit == dx) org[n++] = def;
			else if (NEQUAL(lit, def)) org[n++] = lit;
		}
		org.resize(n);
		if (n == 1) nUnits++;
		else {
			devSort(org.data(), n);
			calcSig(org);
		}
	}

	_PFROST_D_ void substitute_single(
		const uint32& p, 
		const uint32& def,
		CNF& cnf, 
		OL& poss,
		OL& negs, 
		cuVecU* units,
		cuVecB* proof)
	{
		assert(def > 1);
		assert(!SIGN(p));
		uint32 n = NEG(p), def_f = FLIP(def);
		// substitute negatives 
		uint32 nNegUnits = 0;
		#pragma unroll
		forall_occurs(negs, i) {
			SCLAUSE& c = cnf[*i];
			if (c.learnt() || c.molten() || c.has(def)) c.markDeleted(); // learnt or tautology
			else substitute_single(n, def_f, c, nNegUnits);
		}
		// substitute positives
		uint32 nPosUnits = 0;
		#pragma unroll
		forall_occurs(poss, i) {
			SCLAUSE& c = cnf[*i];
			if (c.learnt() || c.molten() || c.has(def_f)) c.markDeleted(); // learnt or tautology
			else substitute_single(p, def, c, nPosUnits);
		}
		// add units
		if (nPosUnits || nNegUnits) {
			uint32* ustart = units->jump(nPosUnits + nNegUnits);
			uint32* udata = ustart;
			if (nNegUnits) appendUnits(cnf, negs, udata);
			if (nPosUnits) appendUnits(cnf, poss, udata);
			assert(udata == ustart + nPosUnits + nNegUnits);
		}
		// add proof
		if (proof) {
			uint32 bytes = 0;
			countProofOrg(cnf, poss, bytes);
			countProofOrg(cnf, negs, bytes);
			addr_t pstart = proof->jump(bytes);
			addr_t pdata = pstart;
			assert(pdata);
			addProof(pdata, cnf, negs);
			addProof(pdata, cnf, poss);
			assert(pdata == pstart + bytes);
		}
	}

	_PFROST_D_ uint32 find_sfanin(const uint32& gate_out, CNF& cnf, OL& list)
	{
		assert(gate_out > 1);
		uint32 imp = 0;
		int nImps = 0;
		forall_occurs(list, i) {
			SCLAUSE& c = cnf[*i];
			if (c.original() && c.size() == 2) {
				imp = FLIP(c[0] ^ c[1] ^ gate_out);
				c.melt(); // mark as gate clause
				nImps++;
			}
			if (nImps > 1) return 0; // cannot be a single-input gate
		}
		return imp;
	}

	_PFROST_D_ uint32 find_equ_gate(const uint32& p, CNF& cnf, OL& poss, OL& negs)
	{
		assert(p > 1);
		assert(!SIGN(p));
		assert(checkMolten(cnf, poss, negs));
		uint32 first;
		if (first = find_sfanin(p, cnf, poss)) {
			uint32 second = NEG(p), def = first;
			if (second < def) first = second, second = def;
			assert(first < second);
			forall_occurs(negs, i) {
				SCLAUSE& c = cnf[*i];
				if (c.original() && c.size() == 2 && c[0] == first && c[1] == second) {
					assert(def == first || def == second);
					c.melt(); // mark as fanout clause
					#if EQU_DBG
					printf("c  Gate %d = -/+%d found\n", ABS(p), ABS(def));
					pClauseSet(cnf, poss, negs);
					#endif
					return def;
				}
			}
		}
		freezeBinaries(cnf, poss);
		return 0;
	}

} // parafrost namespace


#endif